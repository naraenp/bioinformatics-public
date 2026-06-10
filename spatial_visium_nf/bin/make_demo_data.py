#!/usr/bin/env python3
"""Generate a tiny self-contained spatial dataset that exercises the whole
pipeline — with a *known, planted* ground truth the deconvolution must recover.

The real pipeline runs on a 10x Visium breast-cancer section plus the Wu et al.
2021 scRNA-seq atlas as the reference (see fetch_real_data.sh), which is ~0.6 GB
and needs the full scanpy/squidpy stack. For CI, local iteration, and the smoke
test we instead synthesize:

  * a scRNA-seq REFERENCE of K cell types, each with its own marker genes, as a
    10x-style sparse MTX bundle + a metadata.csv carrying the cell-type labels;
  * a Visium SPATIAL sample laid out on a grid, where each spot is a *known*
    mixture of the K cell types — the proportions are planted with spatial
    structure (each cell type dominates a region), so (a) the NNLS deconvolution
    must recover them within tolerance and (b) each type's marker genes are
    spatially autocorrelated and must surface as spatially variable (Moran's I).

The demo writes the SAME on-disk formats as the real inputs, so load_spatial.py
and load_reference.py run an identical code path on toy and real data. The
planted proportions are the spatial analogue of the planted DE genes in
plant_rnaseq_nf / aml_rnaseq_nf.

Writes into --outdir:
    reference/matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz, metadata.csv
    spatial/filtered_feature_bc_matrix.h5
    spatial/spatial/tissue_positions_list.csv, scalefactors_json.json
    truth/proportions.csv     planted per-spot cell-type proportions
    truth/markers.csv         planted marker genes per cell type
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np

# Breast-cancer-flavoured cell types, so the demo rehearses the real biology.
TYPE_NAMES = ["Tcell", "Bcell", "Myeloid", "Epithelial", "Stroma"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--genes", type=int, default=200)
    p.add_argument("--markers", type=int, default=12,
                   help="marker genes per cell type")
    p.add_argument("--marker-fold", type=float, default=8.0,
                   help="expression fold-up of a type's markers in its signature")
    p.add_argument("--cells-per-type", type=int, default=300)
    p.add_argument("--ref-depth", type=int, default=3000,
                   help="mean counts per reference cell")
    p.add_argument("--grid", type=int, default=14,
                   help="spatial grid side; spots = grid*grid")
    p.add_argument("--spot-depth", type=int, default=8000,
                   help="mean counts per spatial spot (spots mix many cells)")
    p.add_argument("--sigma", type=float, default=3.0,
                   help="spatial spread of each cell type's domain (in spots)")
    p.add_argument("--mix-noise", type=float, default=12.0,
                   help="Dirichlet concentration; higher = tighter to the planted mean")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


# ---------------------------------------------------------------------------
# Pure-numpy core (unit-tested without the heavy env)
# ---------------------------------------------------------------------------

def build_signatures(rng, n_genes, n_types, markers_per_type, marker_fold):
    """Per-cell-type expression-rate signatures with distinct marker blocks.

    Returns (S, marker_idx) where S is (n_genes, n_types), each column an
    expression-rate vector normalised to sum 1, and marker_idx maps each type
    index to its block of marker-gene indices.
    """
    base = rng.uniform(0.5, 2.0, size=n_genes)            # shared baseline rate
    S = np.tile(base[:, None], (1, n_types)).astype(float)
    marker_idx: dict[int, np.ndarray] = {}
    used = 0
    for k in range(n_types):
        idx = np.arange(used, used + markers_per_type)
        marker_idx[k] = idx
        S[idx, k] *= marker_fold                          # up in this type only
        used += markers_per_type
    # mild per-type, per-gene jitter so non-marker genes aren't identical
    S *= rng.uniform(0.9, 1.1, size=S.shape)
    S /= S.sum(axis=0, keepdims=True)                     # column = rate profile
    return S, marker_idx


def planted_proportions(rng, coords, n_types, sigma, concentration):
    """Per-spot cell-type proportions with spatial structure.

    Each type is assigned a domain centre; a spot's mean proportion of a type
    falls off as a Gaussian of distance to that centre (softmax-normalised).
    Dirichlet noise around the mean keeps it realistic but recoverable.
    """
    coords = np.asarray(coords, float)
    lo, hi = coords.min(0), coords.max(0)
    # Spread K centres across the tissue (corners + centre for K=5).
    centres = lo + rng.uniform(0.15, 0.85, size=(n_types, 2)) * (hi - lo)
    d2 = ((coords[:, None, :] - centres[None, :, :]) ** 2).sum(-1)   # spots x K
    logits = -d2 / (2.0 * sigma ** 2)
    mean = np.exp(logits - logits.max(1, keepdims=True))
    mean /= mean.sum(1, keepdims=True)
    P = np.array([rng.dirichlet(concentration * m + 1e-3) for m in mean])
    return P


def mix_spot_expression(P, S):
    """Expected spot expression rates = proportions @ signatures^T.

    P is (spots, K), S is (genes, K); returns (spots, genes) rate matrix with
    each row normalised to sum 1.
    """
    rate = P @ S.T                       # spots x genes
    rate /= rate.sum(1, keepdims=True)
    return rate


def poisson_counts(rng, rate, depth):
    """Sample integer counts ~ Poisson(depth * rate) row-wise."""
    return rng.poisson(depth * rate).astype(np.int64)


# ---------------------------------------------------------------------------
# I/O — emit the same formats as the real inputs
# ---------------------------------------------------------------------------

def _gzip_text(path: Path, text: str) -> None:
    with gzip.open(path, "wt") as fh:
        fh.write(text)


def write_reference_mtx(outdir: Path, counts, barcodes, gene_ids, gene_names,
                        cell_types) -> None:
    """Write a CellRanger-style sparse bundle (features x cells) + metadata.csv,
    matching the Wu et al. GSE176078 reference layout the real loader reads."""
    from scipy import sparse
    from scipy.io import mmwrite

    refdir = outdir / "reference"
    refdir.mkdir(parents=True, exist_ok=True)
    mat = sparse.csc_matrix(counts.T)                    # features x cells
    # mmwrite needs a real path; write then gzip.
    tmp = refdir / "matrix.mtx"
    mmwrite(tmp, mat, field="integer")
    with open(tmp, "rb") as fh, gzip.open(refdir / "matrix.mtx.gz", "wb") as gz:
        gz.write(fh.read())
    tmp.unlink()
    _gzip_text(refdir / "barcodes.tsv.gz", "\n".join(barcodes) + "\n")
    _gzip_text(refdir / "features.tsv.gz",
               "\n".join(f"{i}\t{n}\tGene Expression" for i, n in
                         zip(gene_ids, gene_names)) + "\n")
    lines = ["barcode,cell_type"]
    lines += [f"{b},{t}" for b, t in zip(barcodes, cell_types)]
    (refdir / "metadata.csv").write_text("\n".join(lines) + "\n")


def write_spatial_10x_h5(path: Path, counts, barcodes, gene_ids, gene_names) -> None:
    """Write a 10x v3 filtered_feature_bc_matrix.h5 (CSC, features x barcodes)."""
    import h5py
    from scipy import sparse

    mat = sparse.csc_matrix(counts.T.astype(np.int32))   # features x barcodes
    n_features, n_bc = mat.shape
    with h5py.File(path, "w") as f:
        g = f.create_group("matrix")
        g.create_dataset("barcodes", data=np.array(barcodes, dtype="S"))
        g.create_dataset("data", data=mat.data.astype(np.int32))
        g.create_dataset("indices", data=mat.indices.astype(np.int64))
        g.create_dataset("indptr", data=mat.indptr.astype(np.int64))
        g.create_dataset("shape", data=np.array([n_features, n_bc], dtype=np.int32))
        fg = g.create_group("features")
        fg.create_dataset("id", data=np.array(gene_ids, dtype="S"))
        fg.create_dataset("name", data=np.array(gene_names, dtype="S"))
        fg.create_dataset("feature_type",
                          data=np.array(["Gene Expression"] * n_features, dtype="S"))
        fg.create_dataset("genome", data=np.array(["demo"] * n_features, dtype="S"))


def write_spatial_positions(spatial_dir: Path, barcodes, coords) -> None:
    """tissue_positions_list.csv (Space Ranger 1.x: no header) + scalefactors."""
    spatial_dir.mkdir(parents=True, exist_ok=True)
    # barcode, in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres
    px = 100.0                                            # pixels per array unit
    rows = []
    for bc, (r, c) in zip(barcodes, coords):
        rows.append(f"{bc},1,{int(r)},{int(c)},{r * px:.1f},{c * px:.1f}")
    (spatial_dir / "tissue_positions_list.csv").write_text("\n".join(rows) + "\n")
    scale = {"tissue_hires_scalef": 1.0, "tissue_lowres_scalef": 1.0,
             "fiducial_diameter_fullres": px, "spot_diameter_fullres": px * 0.8}
    (spatial_dir / "scalefactors_json.json").write_text(json.dumps(scale))


# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)

    n_types = len(TYPE_NAMES)
    n_genes = max(args.genes, n_types * args.markers + 20)

    # Gene names: marker blocks first (MK_<type>_<j>), then housekeeping (HK_i).
    S, marker_idx = build_signatures(rng, n_genes, n_types, args.markers,
                                     args.marker_fold)
    gene_names = [""] * n_genes
    for k, idx in marker_idx.items():
        for j, gi in enumerate(idx):
            gene_names[gi] = f"MK_{TYPE_NAMES[k]}_{j:02d}"
    for gi in range(n_genes):
        if not gene_names[gi]:
            gene_names[gi] = f"HK_{gi:04d}"
    gene_ids = [f"GENE{gi:05d}" for gi in range(n_genes)]

    # ---- reference: cells_per_type cells of each type ----------------------
    ref_counts, ref_barcodes, ref_types = [], [], []
    for k, tname in enumerate(TYPE_NAMES):
        rate = S[:, k] / S[:, k].sum()
        cnt = poisson_counts(rng, np.tile(rate, (args.cells_per_type, 1)),
                             args.ref_depth)
        ref_counts.append(cnt)
        ref_barcodes += [f"{tname}_cell{i:04d}" for i in range(args.cells_per_type)]
        ref_types += [tname] * args.cells_per_type
    ref_counts = np.vstack(ref_counts)
    write_reference_mtx(out, ref_counts, ref_barcodes, gene_ids, gene_names,
                        ref_types)

    # ---- spatial: grid of spots, planted mixtures --------------------------
    coords = [(r, c) for r in range(args.grid) for c in range(args.grid)]
    spot_barcodes = [f"spot_{r:02d}x{c:02d}-1" for (r, c) in coords]
    P = planted_proportions(rng, coords, n_types, args.sigma, args.mix_noise)
    rate = mix_spot_expression(P, S)
    spot_counts = poisson_counts(rng, rate, args.spot_depth)

    spatial_root = out / "spatial"
    spatial_root.mkdir(parents=True, exist_ok=True)
    write_spatial_10x_h5(spatial_root / "filtered_feature_bc_matrix.h5",
                         spot_counts, spot_barcodes, gene_ids, gene_names)
    write_spatial_positions(spatial_root / "spatial", spot_barcodes, coords)

    # ---- ground truth for the self-check -----------------------------------
    truth = out / "truth"
    truth.mkdir(parents=True, exist_ok=True)
    header = "barcode," + ",".join(TYPE_NAMES)
    prop_lines = [header]
    for bc, p in zip(spot_barcodes, P):
        prop_lines.append(bc + "," + ",".join(f"{v:.6f}" for v in p))
    (truth / "proportions.csv").write_text("\n".join(prop_lines) + "\n")
    mk_lines = ["cell_type,marker_genes"]
    for k, tname in enumerate(TYPE_NAMES):
        names = ";".join(gene_names[gi] for gi in marker_idx[k])
        mk_lines.append(f"{tname},{names}")
    (truth / "markers.csv").write_text("\n".join(mk_lines) + "\n")

    print(f"demo data -> {out}")
    print(f"  reference: {ref_counts.shape[0]} cells x {n_genes} genes, "
          f"{n_types} types ({args.cells_per_type}/type)")
    print(f"  spatial:   {len(spot_barcodes)} spots ({args.grid}x{args.grid}) "
          f"x {n_genes} genes")
    print(f"  planted proportions -> {truth/'proportions.csv'}; "
          f"{args.markers} markers/type spatially structured (sigma={args.sigma})")


if __name__ == "__main__":
    main()
