#!/usr/bin/env python3
"""LOAD_REFERENCE — read an annotated scRNA-seq reference into an AnnData.

Reads a CellRanger-style sparse bundle (matrix.mtx.gz, barcodes.tsv.gz,
features.tsv.gz) plus a metadata CSV that carries the per-cell cell-type label,
applies a light per-cell count filter, and writes reference.h5ad with
obs['cell_type'].

The real reference is the Wu et al. 2021 breast-cancer atlas (GSE176078);
fetch_real_data.sh arranges its files into the standard 10x names this reads.
The demo writes the same layout (see make_demo_data.py).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dir", required=True, type=Path,
                   help="10x mtx directory (matrix/barcodes/features .tsv.gz)")
    p.add_argument("--metadata", required=True, type=Path,
                   help="CSV indexed by barcode with a cell-type column")
    p.add_argument("--celltype-col", default="cell_type")
    p.add_argument("--min-counts", type=int, default=200,
                   help="drop reference cells below this total count")
    p.add_argument("--out", required=True, type=Path, help="reference.h5ad")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    adata = sc.read_10x_mtx(args.dir)
    adata.var_names_make_unique()

    meta = pd.read_csv(args.metadata, index_col=0)
    if args.celltype_col not in meta.columns:
        raise SystemExit(
            f"--celltype-col '{args.celltype_col}' not in metadata columns "
            f"{list(meta.columns)}")
    common = adata.obs_names.intersection(meta.index)
    if len(common) == 0:
        raise SystemExit("no barcodes shared between matrix and metadata")
    adata = adata[common].copy()
    adata.obs["cell_type"] = meta.loc[adata.obs_names, args.celltype_col].astype(str).to_numpy()

    sc.pp.filter_cells(adata, min_counts=args.min_counts)
    adata = adata[~adata.obs["cell_type"].isin(["nan", "NA", ""])].copy()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out)
    n_types = adata.obs["cell_type"].nunique()
    print(f"wrote {args.out}: {adata.n_obs} cells x {adata.n_vars} genes, "
          f"{n_types} cell types")


if __name__ == "__main__":
    main()
