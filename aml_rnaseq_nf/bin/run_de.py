#!/usr/bin/env python3
"""Per-gene Welch's t-test on log2-CPM with Benjamini-Hochberg FDR.

Output columns: gene, mean_AML, mean_healthy, log2FC, t, pvalue, padj.
log2FC is AML - healthy (already in log space).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--lcpm", required=True, type=Path)
    p.add_argument("--meta", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    return p.parse_args()


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (q-values)."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    # Monotone cumulative minimum from the end.
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty_like(p)
    q[order] = np.clip(ranked, 0.0, 1.0)
    return q


def main() -> None:
    args = parse_args()
    lcpm = pd.read_csv(args.lcpm, sep="\t", index_col=0)
    meta = pd.read_csv(args.meta, sep="\t")

    aml_samples = meta.loc[meta["group"] == "AML", "sample_id"].tolist()
    hlt_samples = meta.loc[meta["group"] == "healthy", "sample_id"].tolist()
    if not aml_samples or not hlt_samples:
        raise SystemExit("Both 'AML' and 'healthy' samples are required.")

    x = lcpm[aml_samples].to_numpy()
    y = lcpm[hlt_samples].to_numpy()

    t_stat, p_val = stats.ttest_ind(x, y, axis=1, equal_var=False, nan_policy="omit")
    p_val = np.where(np.isnan(p_val), 1.0, p_val)
    padj = bh_fdr(p_val)

    out = pd.DataFrame({
        "gene": lcpm.index,
        "mean_AML": x.mean(axis=1),
        "mean_healthy": y.mean(axis=1),
        "log2FC": x.mean(axis=1) - y.mean(axis=1),
        "t": t_stat,
        "pvalue": p_val,
        "padj": padj,
    })
    out = out.sort_values("padj", kind="mergesort").reset_index(drop=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(
        f"wrote {args.out}: {len(out)} genes, "
        f"{(out['padj'] < 0.05).sum()} at FDR<0.05"
    )


if __name__ == "__main__":
    main()
