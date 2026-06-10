#!/usr/bin/env python3
"""DECONVOLVE — per-spot cell-type deconvolution by non-negative least squares.

For each Visium spot we model its (library-size-normalised, linear) expression
as a non-negative mixture of the reference cell-type signatures and solve

    minimise || S w - x ||_2   subject to  w >= 0

with scipy.optimize.nnls, then normalise w to proportions summing to 1. This is
the required, default deconvolution: transparent, fast, dependency-light, and —
on the synthetic demo — recoverable against the planted ground truth. (The
optional cell2location / RCTD comparison is a separate v1.1 stage in a heavier
env; this NNLS path is the on-brand analogue of the hand-rolled DE/ORA in the
plant and AML pipelines.)

Per-gene weighting: a plain L2 fit in linear space is dominated by the handful
of very high-expression genes (immunoglobulins, etc.), so spiky immune
signatures absorb every spot — on real cross-platform data the fit collapses
onto a couple of cell types. We therefore divide each gene (both the signature
and the spot) by its mean signature level before solving, i.e. minimise
sum_g (1/mean_g^2) (S_g w - x_g)^2. This inverse-mean weighting puts genes on a
comparable scale so the fit reflects all cell types, not just the loudest genes.
It is still fully transparent (one deterministic factor per gene) and leaves the
matched synthetic demo — where every gene is already on the same scale —
unchanged within tolerance. The proper probabilistic correction for the
cross-platform shift is the v1.1 cell2location / RCTD bake-off.

Input:  spatial_norm.h5ad (X = linear normalised over HVGs) + signature.tsv
Output: proportions.tsv (spot x cell-type proportions, + dominant_type, coords)
Self-check: with --truth (planted proportions) it reports MAE / correlation and,
with --check, exits non-zero if the mean absolute error exceeds --tol.
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
