#!/usr/bin/env python3
"""NORMALIZE — library-size normalisation, log1p, and a shared HVG space.

Puts the spatial spots and the reference cells into the SAME feature space and
the same normalisation so the downstream signature + NNLS deconvolution are
comparable:

  * restrict both to the genes they share;
  * pick highly variable genes on the reference (where the cell-type structure
    lives), intersected with the shared set;
  * store, per object, X = library-size-normalised LINEAR expression over the
    HVGs (target_sum=1e4) — the space NNLS deconvolution runs in — plus a
    layers['lognorm'] = log1p(X) used for HVG selection, Moran's I and plots.

Writes spatial_norm.h5ad, reference_norm.h5ad and hvgs.txt.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import scanpy as sc


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--spatial", required=True, type=Path, help="spatial_qc.h5ad")
    p.add_argument("--reference", required=True, type=Path, help="reference.h5ad")
    p.add_argument("--n-hvg", type=int, default=2000)
    p.add_argument("--target-sum", type=float, default=1e4)
    p.add_argument("--out-spatial", required=True, type=Path)
    p.add_argument("--out-reference", required=True, type=Path)
    p.add_argument("--hvg-out", required=True, type=Path)
    return p.parse_args()


def lin_and_log(adata, genes, target_sum):
    """Return a copy restricted to `genes` with X = linear-normalised counts and
    layers['lognorm'] = log1p of it."""
    sub = adata[:, genes].copy()
    sub.layers["counts"] = sub.X.copy()
    sc.pp.normalize_total(sub, target_sum=target_sum)   # linear, library-size
    sub.layers["lognorm"] = sub.X.copy()
    sc.pp.log1p(sub.layers["lognorm"])                  # log1p of the linear layer
    return sub


def main() -> None:
    args = parse_args()
    spatial = sc.read_h5ad(args.spatial)
    reference = sc.read_h5ad(args.reference)

    shared = spatial.var_names.intersection(reference.var_names)
    if len(shared) < 10:
        raise SystemExit(f"only {len(shared)} shared genes between spatial and "
                         "reference — check the two inputs use the same gene IDs")

    # HVGs from the reference (carries the cell-type signal), within shared genes.
    ref_shared = reference[:, shared].copy()
    sc.pp.normalize_total(ref_shared, target_sum=args.target_sum)
    sc.pp.log1p(ref_shared)
    n_top = min(args.n_hvg, ref_shared.n_vars)
    sc.pp.highly_variable_genes(ref_shared, n_top_genes=n_top)
    hvgs = ref_shared.var_names[ref_shared.var["highly_variable"]].tolist()

    spatial_n = lin_and_log(spatial, hvgs, args.target_sum)
    reference_n = lin_and_log(reference, hvgs, args.target_sum)

    args.out_spatial.parent.mkdir(parents=True, exist_ok=True)
    spatial_n.write_h5ad(args.out_spatial)
    reference_n.write_h5ad(args.out_reference)
    args.hvg_out.write_text("\n".join(hvgs) + "\n")
    print(f"{len(shared)} shared genes -> {len(hvgs)} HVGs")
    print(f"wrote {args.out_spatial} ({spatial_n.n_obs} spots) and "
          f"{args.out_reference} ({reference_n.n_obs} cells)")


if __name__ == "__main__":
    main()
