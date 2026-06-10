#!/usr/bin/env python3
"""LOAD_SPATIAL — read a 10x Visium sample into an AnnData (.h5ad).

Reads the Space Ranger filtered feature-barcode matrix (filtered_feature_bc_matrix.h5)
plus the spot positions (spatial/tissue_positions[_list].csv), keeps the in-tissue
spots, and attaches the array (row, col) grid coordinates as obsm['spatial'].

Works identically on the real breast-cancer download and on the synthetic demo,
which writes the same files (see make_demo_data.py).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5", required=True, type=Path,
                   help="filtered_feature_bc_matrix.h5")
    p.add_argument("--positions", required=True, type=Path,
                   help="spatial/tissue_positions[_list].csv")
    p.add_argument("--out", required=True, type=Path, help="spatial.h5ad")
    return p.parse_args()


POS_COLS = ["barcode", "in_tissue", "array_row", "array_col",
            "pxl_row_in_fullres", "pxl_col_in_fullres"]


def read_positions(path: Path) -> pd.DataFrame:
    """Read tissue positions, tolerating the headered (SpaceRanger 2.x) and
    header-less list (1.x) variants."""
    head = path.read_text().splitlines()[0]
    has_header = head.lower().startswith("barcode")
    df = pd.read_csv(path, header=0 if has_header else None)
    df.columns = POS_COLS[:df.shape[1]]
    return df.set_index("barcode")


def main() -> None:
    args = parse_args()
    adata = sc.read_10x_h5(args.h5)
    adata.var_names_make_unique()

    pos = read_positions(args.positions)
    common = adata.obs_names.intersection(pos.index)
    adata = adata[common].copy()
    pos = pos.loc[adata.obs_names]
    for c in ["in_tissue", "array_row", "array_col"]:
        adata.obs[c] = pos[c].to_numpy()
    # keep in-tissue spots only
    adata = adata[adata.obs["in_tissue"] == 1].copy()
    # obsm['spatial'] as (x, y) = (array_col, array_row) for plotting
    adata.obsm["spatial"] = adata.obs[["array_col", "array_row"]].to_numpy()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out)
    print(f"wrote {args.out}: {adata.n_obs} in-tissue spots x {adata.n_vars} genes")


if __name__ == "__main__":
    main()
