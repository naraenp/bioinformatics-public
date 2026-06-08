#!/usr/bin/env python3
"""Over-representation (GO / pathway) analysis on the significant DE genes.

A dependency-light enrichment step in the same spirit as the hand-rolled BH FDR
in aml_rnaseq_nf: read gene sets from a GMT file, take the DE genes that pass
the padj / |log2FC| thresholds as the foreground, use every tested gene as the
background universe, and score each gene set with a one-sided hypergeometric
test (over-representation). Per-term p-values are BH-adjusted.

This is what links the most differentially expressed genes back to *phenotype*:
the enriched terms name the biological processes (stress response, ion
transport, hormone signalling, …) that separate the tolerant and susceptible
lines.

GMT format (one gene set per line, tab-separated):
    term_id   term_name   gene1   gene2   ...

Output (enrichment.tsv): term, name, set_size, n_selected, k_overlap,
enrichment, pvalue, padj, overlap_genes.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--de", required=True, type=Path)
    p.add_argument("--gmt", required=True, type=Path)
    p.add_argument("--fdr", type=float, default=0.05, help="DE padj cutoff")
    p.add_argument("--lfc", type=float, default=1.0, help="DE |log2FC| cutoff")
    p.add_argument("--min-set", type=int, default=3)
    p.add_argument("--max-set", type=int, default=2000)
    p.add_argument("--out", required=True, type=Path)
    return p.parse_args()


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (matches aml_rnaseq_nf)."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty_like(p)
    q[order] = np.clip(ranked, 0.0, 1.0)
    return q


def read_gmt(path: Path) -> list[tuple[str, str, set[str]]]:
    sets = []
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            term, name, *genes = fields
            sets.append((term, name, {g for g in genes if g}))
    return sets


def main() -> None:
    args = parse_args()
    de = pd.read_csv(args.de, sep="\t")
    de = de.dropna(subset=["padj"])

    universe = set(de["gene"])
    sig_mask = (de["padj"] < args.fdr) & (de["log2FC"].abs() >= args.lfc)
    selected = set(de.loc[sig_mask, "gene"])
    M = len(universe)            # population size
    n = len(selected)           # successes in population
    if n == 0:
        raise SystemExit("no DE genes pass the thresholds; nothing to enrich")

    rows = []
    for term, name, genes in read_gmt(args.gmt):
        set_in_uni = genes & universe
        N = len(set_in_uni)     # draws
        if N < args.min_set or N > args.max_set:
            continue
        overlap = set_in_uni & selected
        k = len(overlap)
        if k == 0:
            continue
        # P(X >= k) with the hypergeometric survival function.
        pval = stats.hypergeom.sf(k - 1, M, n, N)
        expected = N * n / M
        rows.append({
            "term": term,
            "name": name,
            "set_size": N,
            "n_selected": n,
            "k_overlap": k,
            "enrichment": round(k / expected, 3) if expected else np.nan,
            "pvalue": pval,
            "overlap_genes": ",".join(sorted(overlap)),
        })

    if not rows:
        raise SystemExit("no gene set met the overlap / size criteria")

    out = pd.DataFrame(rows)
    out["padj"] = bh_fdr(out["pvalue"].to_numpy())
    out = out.sort_values("padj", kind="mergesort").reset_index(drop=True)
    out = out[["term", "name", "set_size", "n_selected", "k_overlap",
               "enrichment", "pvalue", "padj", "overlap_genes"]]
    out.to_csv(args.out, sep="\t", index=False)
    print(f"wrote {args.out}: {len(out)} terms tested, "
          f"{(out['padj'] < 0.05).sum()} at padj<0.05 "
          f"(foreground={n} DE genes / universe={M})")


if __name__ == "__main__":
    main()
