#!/usr/bin/env python3
"""SVG — spatially variable genes by Moran's I on a kNN spot graph.

Builds a symmetric k-nearest-neighbour graph over the spot coordinates and
scores each gene's spatial autocorrelation with Moran's I, computed
transparently from first principles (no black box) and vectorised across genes.
Significance is assessed by a seeded permutation test on the spot ordering.

Squidpy's sq.gr.spatial_autocorr is the production equivalent; this hand-rolled
version keeps the core CI light and the maths legible — in line with the
pipeline's "transparent over black-box" stance.

Input:  spatial_norm.h5ad (uses layers['lognorm']; obsm['spatial'])
Output: svg.tsv (gene, morans_I, pval, rank), sorted by Moran's I.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def knn_weights(coords: np.ndarray, k: int) -> np.ndarray:
    """Symmetric binary kNN adjacency (n x n, zero diagonal). Pure scipy."""
    from scipy.spatial import cKDTree

    coords = np.asarray(coords, float)
    n = coords.shape[0]
    k = min(k, n - 1)
    _, idx = cKDTree(coords).query(coords, k=k + 1)      # includes self
    W = np.zeros((n, n), dtype=float)
    rows = np.repeat(np.arange(n), k)
    cols = idx[:, 1:].ravel()
    W[rows, cols] = 1.0
    W = np.maximum(W, W.T)                                # symmetrise
    np.fill_diagonal(W, 0.0)
    return W


def morans_i(V: np.ndarray, W: np.ndarray) -> np.ndarray:
    """Moran's I for each column of V (spots x genes) under weights W.

    I = (n / W_sum) * (z' W z) / (z' z), with z the mean-centred values.
    Vectorised across genes. Pure numpy.
    """
    n = V.shape[0]
    Z = V - V.mean(axis=0, keepdims=True)
    w_sum = W.sum()
    num = np.einsum("ij,jg,ig->g", W, Z, Z)              # z' W z per gene
    den = (Z * Z).sum(axis=0)                            # z' z per gene
    den = np.where(den == 0, np.nan, den)
    return (n / w_sum) * (num / den)


def permutation_pvals(V, W, observed, n_perm, seed):
    """One-sided permutation p-values: fraction of permuted I >= observed."""
    rng = np.random.default_rng(seed)
    n = V.shape[0]
    ge = np.zeros(V.shape[1], dtype=int)
    for _ in range(n_perm):
        perm = rng.permutation(n)
        ge += (morans_i(V[perm], W) >= observed).astype(int)
    return (ge + 1) / (n_perm + 1)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--spatial", required=True, type=Path, help="spatial_norm.h5ad")
    p.add_argument("--k", type=int, default=6, help="neighbours per spot")
    p.add_argument("--n-perm", type=int, default=100, help="permutations for p-values")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", required=True, type=Path, help="svg.tsv")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    import pandas as pd
    import scanpy as sc

    adata = sc.read_h5ad(args.spatial)
    coords = np.asarray(adata.obsm["spatial"], float)
    X = adata.layers["lognorm"] if "lognorm" in adata.layers else adata.X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)

    W = knn_weights(coords, args.k)
    I = morans_i(X, W)
    pvals = permutation_pvals(X, W, I, args.n_perm, args.seed)

    df = pd.DataFrame({"gene": adata.var_names, "morans_I": I, "pval": pvals})
    df = df.dropna(subset=["morans_I"]).sort_values(
        "morans_I", ascending=False, kind="mergesort").reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    top = df.head(5)
    print(f"wrote {args.out}: {len(df)} genes scored (k={args.k}, "
          f"{args.n_perm} perms)")
    print("  top SVGs: " + ", ".join(f"{g}({i:.2f})"
          for g, i in zip(top["gene"], top["morans_I"])))


if __name__ == "__main__":
    main()
