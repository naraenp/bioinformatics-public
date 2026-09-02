#!/usr/bin/env python3
"""Build a raw count matrix from recount3 (TCGA-LAML vs. GTEx whole blood).

Parses the GENCODE v26 GTF for an Ensembl -> HGNC map, loads both count
matrices, subsamples each cohort to --n-per-group, inner-joins on Ensembl ID,
maps to symbols, and keeps genes with CPM >= 1 in >= 25% of samples (an
edgeR-style filter on pooled expression, so it stays independent of the group
contrast). Writes counts_raw.tsv and metadata.tsv.

Caveat: GTEx has no bone marrow, so whole blood is the closest large healthy
comparator. Cohort is therefore confounded with disease, tissue and collection
protocol at once -- see docs/REPORT.md before reading anything biological into
the direction of a fold change.
"""
from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd

GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')
GENE_NAME_RE = re.compile(r'gene_name "([^"]+)"')


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--aml-counts", required=True, type=Path)
    p.add_argument("--healthy-counts", required=True, type=Path)
    p.add_argument("--gtf", required=True, type=Path)
    p.add_argument("--n-per-group", type=int, default=50)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--counts", required=True, type=Path)
    p.add_argument("--meta", required=True, type=Path)
    return p.parse_args()


def build_gene_map(gtf_path: Path) -> dict[str, str]:
    """Map Ensembl gene_id (no version) -> HGNC symbol from a GENCODE GTF."""
    mapping: dict[str, str] = {}
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.split("\t", 8)
            if len(cols) < 9 or cols[2] != "gene":
                continue
            gid_m = GENE_ID_RE.search(cols[8])
            gname_m = GENE_NAME_RE.search(cols[8])
            if gid_m and gname_m:
                mapping[gid_m.group(1).split(".", 1)[0]] = gname_m.group(1)
    return mapping


def load_recount3_counts(counts_path: Path) -> pd.DataFrame:
    """Load a recount3 gene_sums.gz TSV; strip the Ensembl version suffix."""
    df = pd.read_csv(counts_path, sep="\t", comment="#", index_col=0)
    df.index = df.index.str.split(".").str[0]
    df.index.name = "gene_id"
    return df


def cpm_filter(counts: pd.DataFrame, min_cpm: float = 1.0, min_frac: float = 0.25) -> pd.DataFrame:
    """Keep genes with CPM >= min_cpm in at least min_frac of samples."""
    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.divide(lib, axis=1) * 1e6
    keep = (cpm >= min_cpm).sum(axis=1) >= int(min_frac * counts.shape[1])
    return counts.loc[keep]


def main() -> None:
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    n = args.n_per_group

    gene_map = build_gene_map(args.gtf)
    print(f"[1/4] mapped {len(gene_map):,} Ensembl gene IDs to HGNC symbols")

    aml = load_recount3_counts(args.aml_counts)
    hlt = load_recount3_counts(args.healthy_counts)
    print(f"[2/4] loaded AML {aml.shape} and healthy {hlt.shape} matrices")
    if aml.shape[1] < n or hlt.shape[1] < n:
        raise SystemExit(
            f"need >= {n} samples per group; have {aml.shape[1]} AML and "
            f"{hlt.shape[1]} healthy"
        )

    # Subsample each cohort to a balanced n, then inner-join on Ensembl ID.
    aml = aml[list(rng.choice(aml.columns, size=n, replace=False))]
    hlt = hlt[list(rng.choice(hlt.columns, size=n, replace=False))]
    common = aml.index.intersection(hlt.index)
    combined = pd.concat([aml.loc[common], hlt.loc[common]], axis=1)

    # Map Ensembl -> HGNC symbol; drop unmapped, sum any duplicate symbols.
    symbols = combined.index.map(gene_map.get)
    combined = combined.loc[symbols.notna()]
    combined.index = symbols[symbols.notna()]
    combined.index.name = "gene"
    if combined.index.has_duplicates:
        combined = combined.groupby(level=0).sum()

    combined = cpm_filter(combined)
    print(f"[3/4] combined matrix after symbol map + CPM filter: {combined.shape}")

    meta = pd.DataFrame({
        "sample_id": list(combined.columns),
        "group": ["AML"] * n + ["healthy"] * n,
        "source": ["TCGA-LAML"] * n + ["GTEx-BLOOD"] * n,
    })

    combined.to_csv(args.counts, sep="\t")
    meta.to_csv(args.meta, sep="\t", index=False)
    print(f"[4/4] wrote {args.counts} {combined.shape} and {args.meta} {meta.shape}")


if __name__ == "__main__":
    main()
