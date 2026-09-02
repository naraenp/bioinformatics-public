#!/usr/bin/env python3
"""DECONVOLVE — per-spot cell-type deconvolution by non-negative least squares.

A Visium spot holds several cells, so its expression is a mixture. For each spot
we solve `minimise ||S w - x||_2` subject to `w >= 0` (scipy.optimize.nnls) and
normalise the weights to proportions.

Each gene is first divided by its mean signature level. Without this, the fit is
dominated by a few very high-expression genes (immunoglobulins here) and collapses
onto one or two cell types on real data. It is weighted least squares with one
deterministic weight per gene; RCTD and cell2location are the principled fix, since
they model counts and return uncertainty.

Input:  spatial_norm.h5ad (linear normalised over HVGs) + signature.tsv
Output: proportions.tsv (spot x cell type, + dominant_type, coords)
Self-check: --truth reports MAE / correlation vs planted proportions; --check
exits non-zero if the MAE exceeds --tol.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def gene_scale_factors(S: np.ndarray, floor: float = 1e-6) -> np.ndarray:
    """Per-gene weighting factor: the mean signature level across cell types.

    Dividing both S and the spot expression by this puts every gene on a
    comparable scale before the least-squares fit, so high-expression genes do
    not dominate the residual. Floored to avoid dividing by ~zero. Pure numpy."""
    g = np.asarray(S, dtype=float).mean(axis=1)
    return np.maximum(g, floor)


def nnls_deconvolve(X: np.ndarray, S: np.ndarray) -> np.ndarray:
    """Proportions for each row of X (spots x genes) against signatures S
    (genes x types). Returns (spots x types) rows summing to 1. Pure scipy."""
    from scipy.optimize import nnls

    n, k = X.shape[0], S.shape[1]
    P = np.zeros((n, k), dtype=float)
    for i in range(n):
        w, _ = nnls(S, X[i])
        tot = w.sum()
        P[i] = w / tot if tot > 0 else np.full(k, 1.0 / k)
    return P


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--spatial", required=True, type=Path, help="spatial_norm.h5ad")
    p.add_argument("--signature", required=True, type=Path, help="signature.tsv")
    p.add_argument("--out", required=True, type=Path, help="proportions.tsv")
    p.add_argument("--truth", type=Path, default=None,
                   help="planted proportions CSV for the self-check")
    p.add_argument("--check", action="store_true",
                   help="exit non-zero if MAE vs --truth exceeds --tol")
    p.add_argument("--tol", type=float, default=0.12,
                   help="max tolerated mean absolute error in the self-check")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    import scanpy as sc

    adata = sc.read_h5ad(args.spatial)
    sig = pd.read_csv(args.signature, sep="\t", index_col=0)

    # Align spots and signatures onto their shared genes, in a common order.
    genes = adata.var_names.intersection(sig.index)
    if len(genes) < 5:
        raise SystemExit("too few shared genes between spatial data and signature")
    adata = adata[:, genes]
    sig = sig.loc[genes]
    X = adata.X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    S = sig.to_numpy()
    types = list(sig.columns)

    # Inverse-mean per-gene weighting (see module docstring) so the fit is not
    # dominated by a few very high-expression genes on real cross-platform data.
    g = gene_scale_factors(S)
    P = nnls_deconvolve(X / g[None, :], S / g[:, None])
    prop = pd.DataFrame(P, index=adata.obs_names, columns=types)
    prop.index.name = "barcode"
    prop["dominant_type"] = prop[types].idxmax(axis=1)
    for c in ["array_row", "array_col"]:
        if c in adata.obs:
            prop[c] = adata.obs[c].to_numpy()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    prop.to_csv(args.out, sep="\t")
    print(f"wrote {args.out}: {prop.shape[0]} spots x {len(types)} cell types")
    counts = prop["dominant_type"].value_counts()
    print("  dominant-type spots: " + ", ".join(f"{t}={n}" for t, n in counts.items()))

    if args.truth is not None:
        truth = pd.read_csv(args.truth, index_col=0)
        common_bc = prop.index.intersection(truth.index)
        common_t = [t for t in types if t in truth.columns]
        if len(common_bc) == 0 or len(common_t) == 0:
            raise SystemExit("truth file shares no barcodes/types with the result")
        A = prop.loc[common_bc, common_t].to_numpy()
        B = truth.loc[common_bc, common_t].to_numpy()
        mae = float(np.abs(A - B).mean())
        corr = float(np.corrcoef(A.ravel(), B.ravel())[0, 1])
        per_type = {t: round(float(np.abs(A[:, j] - B[:, j]).mean()), 4)
                    for j, t in enumerate(common_t)}
        print(f"  [self-check] vs planted: MAE={mae:.4f}  corr={corr:.4f}")
        print(f"  [self-check] per-type MAE: {per_type}")
        if args.check and mae > args.tol:
            raise SystemExit(
                f"SELF-CHECK FAILED: MAE {mae:.4f} > tol {args.tol}")
        if args.check:
            print(f"  [self-check] PASS (MAE {mae:.4f} <= tol {args.tol})")


if __name__ == "__main__":
    main()
