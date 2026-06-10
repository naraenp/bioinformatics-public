#!/usr/bin/env python3
"""BUILD_SIGNATURE — per-cell-type reference signature matrix.

Averages the library-size-normalised (linear) reference expression within each
cell type to form a genes x cell-types signature matrix. This is the basis the
NNLS deconvolution regresses each spot onto, and it is deliberately transparent:
a cell type's signature is just the mean expression of its cells.

Input:  reference_norm.h5ad (X = linear normalised over HVGs, obs['cell_type'])
Output: signature.tsv (rows = genes, cols = cell types)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def mean_by_label(X: np.ndarray, labels) -> tuple[np.ndarray, list[str]]:
    """Column-mean of X (cells x genes) within each label.

    Returns (M, types) where M is (genes x n_types) and types is the sorted
    unique label list. Pure numpy so it is unit-testable without scanpy.
    """
    labels = np.asarray(labels)
    types = sorted(set(labels.tolist()))
    M = np.zeros((X.shape[1], len(types)), dtype=float)
    for j, t in enumerate(types):
        M[:, j] = X[labels == t].mean(axis=0)
    return M, types


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", required=True, type=Path, help="reference_norm.h5ad")
    # load_reference.py canonicalises whatever source column was chosen onto
    # obs['cell_type'], so this stays 'cell_type' regardless of the real dataset's
    # own label name (e.g. Wu et al.'s 'celltype_major').
    p.add_argument("--celltype-col", default="cell_type")
    p.add_argument("--out", required=True, type=Path, help="signature.tsv")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    import scanpy as sc

    adata = sc.read_h5ad(args.reference)
    X = adata.X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    M, types = mean_by_label(X, adata.obs[args.celltype_col].to_numpy())
    sig = pd.DataFrame(M, index=adata.var_names, columns=types)
    sig.index.name = "gene"

    args.out.parent.mkdir(parents=True, exist_ok=True)
    sig.to_csv(args.out, sep="\t")
    print(f"wrote {args.out}: {sig.shape[0]} genes x {sig.shape[1]} cell types "
          f"({', '.join(types)})")


if __name__ == "__main__":
    main()
