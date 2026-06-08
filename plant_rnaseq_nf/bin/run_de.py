#!/usr/bin/env python3
"""Differential expression with pydeseq2: tolerant vs. susceptible genotype.

Reads a raw gene-by-sample count matrix and a sample metadata table, fits a
DESeq2 model (median-of-ratios size factors, negative-binomial GLM, Wald test)
on the design ``~genotype``, and contrasts the tolerant genotype against the
susceptible one. Positive log2FC therefore means *higher in the tolerant line*.

Outputs:
  * de_results.tsv   gene, baseMean, log2FC, lfcSE, stat, pvalue, padj
  * norm_counts.tsv  log2(normalized counts + 1), gene x sample (for the heatmap)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--counts", required=True, type=Path,
                   help="raw counts, genes x samples, tab-separated")
    p.add_argument("--meta", required=True, type=Path,
                   help="metadata with columns sample_id, genotype")
    p.add_argument("--design", default="~genotype",
                   help="DESeq2 design formula; the contrasted factor "
                        "(genotype) must come last, e.g. '~condition + genotype'")
    p.add_argument("--tolerant", default="tolerant",
                   help="genotype label for the tolerant line")
    p.add_argument("--susceptible", default="susceptible",
                   help="genotype label for the susceptible line")
    p.add_argument("--out", required=True, type=Path, help="de_results.tsv")
    p.add_argument("--norm-out", required=True, type=Path, help="norm_counts.tsv")
    return p.parse_args()


def build_dds(counts: pd.DataFrame, meta: pd.DataFrame, design: str):
    """Construct and fit a DeseqDataSet, tolerant of the 0.4/0.5 API shift."""
    from pydeseq2.dds import DeseqDataSet

    try:                                    # pydeseq2 >= 0.4 formula API
        dds = DeseqDataSet(counts=counts, metadata=meta, design=design,
                           quiet=True)
    except TypeError:                       # older design_factors API
        factors = [t.strip() for t in design.lstrip("~").split("+") if t.strip()]
        dds = DeseqDataSet(counts=counts, metadata=meta,
                           design_factors=factors, quiet=True)
    dds.deseq2()
    return dds


def main() -> None:
    args = parse_args()
    from pydeseq2.ds import DeseqStats

    # DESeq2 wants integer counts as samples x genes.
    counts = pd.read_csv(args.counts, sep="\t", index_col=0)
    counts = counts.round().astype(int).T
    counts.index.name = "sample_id"

    meta = pd.read_csv(args.meta, sep="\t").set_index("sample_id")
    meta = meta.loc[counts.index]           # align order
    groups = {args.tolerant, args.susceptible}
    if not groups.issubset(set(meta["genotype"])):
        raise SystemExit(
            f"metadata 'genotype' must contain {sorted(groups)}; "
            f"found {sorted(set(meta['genotype']))}"
        )

    dds = build_dds(counts, meta, args.design)

    stat = DeseqStats(dds, contrast=["genotype", args.tolerant, args.susceptible],
                      quiet=True)
    stat.summary()
    res = stat.results_df.rename(columns={"log2FoldChange": "log2FC"}).copy()
    res.index.name = "gene"
    res = res.reset_index()
    res = res.sort_values("padj", kind="mergesort", na_position="last")
    res = res[["gene", "baseMean", "log2FC", "lfcSE", "stat", "pvalue", "padj"]]
    res.to_csv(args.out, sep="\t", index=False)

    # Normalized, log2 expression for the downstream heatmap (genes x samples).
    normed = dds.layers["normed_counts"]
    norm_df = pd.DataFrame(np.log2(normed + 1.0), index=counts.index,
                           columns=dds.var_names).T
    norm_df.index.name = "gene"
    norm_df.to_csv(args.norm_out, sep="\t")

    n_sig = int((res["padj"] < 0.05).sum())
    print(f"wrote {args.out}: {len(res)} genes, {n_sig} at padj<0.05")
    print(f"wrote {args.norm_out}: {norm_df.shape[0]} genes x {norm_df.shape[1]} samples")


if __name__ == "__main__":
    main()
