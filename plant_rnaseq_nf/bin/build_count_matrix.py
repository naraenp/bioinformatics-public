#!/usr/bin/env python3
"""Tidy a featureCounts table into a clean gene-by-sample matrix + metadata.

featureCounts emits a wide table with a leading comment line, six annotation
columns (Geneid, Chr, Start, End, Strand, Length) and one count column per BAM
named by its file path. This collapses that to:

  * counts_raw.tsv   genes x samples, columns are sample_ids
  * metadata.tsv     sample_id, genotype  (carried from the samplesheet)

BAM count columns are matched back to sample_ids by filename stem
(``<sample_id>.sorted.bam`` -> ``<sample_id>``).
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--featurecounts", required=True, type=Path)
    p.add_argument("--samplesheet", required=True, type=Path)
    p.add_argument("--counts-out", required=True, type=Path)
    p.add_argument("--meta-out", required=True, type=Path)
    return p.parse_args()


def stem(col: str) -> str:
    """``/work/ab/cd/S1.sorted.bam`` -> ``S1``."""
    name = Path(col).name
    return re.sub(r"\.(sorted\.)?bam$", "", name)


def main() -> None:
    args = parse_args()
    sheet = pd.read_csv(args.samplesheet)
    if not {"sample_id", "genotype"}.issubset(sheet.columns):
        raise SystemExit("samplesheet needs 'sample_id' and 'genotype' columns")
    # Metadata carries every samplesheet column except the FASTQ paths, so extra
    # design covariates (e.g. 'condition') reach pydeseq2.
    meta_cols = [c for c in sheet.columns if c not in ("fastq_1", "fastq_2")]

    fc = pd.read_csv(args.featurecounts, sep="\t", comment="#")
    annotation = ["Chr", "Start", "End", "Strand", "Length"]
    fc = fc.drop(columns=[c for c in annotation if c in fc.columns])
    fc = fc.rename(columns={"Geneid": "gene"}).set_index("gene")
    fc.columns = [stem(c) for c in fc.columns]

    order = sheet["sample_id"].astype(str).tolist()
    missing = [s for s in order if s not in fc.columns]
    if missing:
        raise SystemExit(f"samples missing from featureCounts table: {missing}")
    counts = fc[order]
    counts.to_csv(args.counts_out, sep="\t")

    meta = sheet[meta_cols].copy()
    meta.to_csv(args.meta_out, sep="\t", index=False)

    print(f"wrote {args.counts_out}: {counts.shape[0]} genes x {counts.shape[1]} samples")
    print(f"wrote {args.meta_out}: {len(meta)} samples")


if __name__ == "__main__":
    main()
