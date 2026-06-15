#!/usr/bin/env python3
"""Self-check: did the pipeline recover the planted variants?

Compares the pipeline's called cohort VCF against the synthetic truth VCF that
make_demo_data.py planted. Both should be normalized + split to biallelic with
`bcftools norm -f genome.fa -m-` first, so a variant is keyed simply by
(CHROM, POS, REF, ALT). Reports site-level recall + precision and per-sample
genotype concordance, and exits non-zero if any threshold is missed — this is
the analog of the planted-DE-gene / planted-proportion checks in the sibling
pipelines, and it is what `run_local.sh --demo` and CI assert on.

Only PASS (or unfiltered ".") calls count toward precision; truth sites with at
least one non-reference sample count toward recall.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--called", required=True, type=Path)
    p.add_argument("--truth", required=True, type=Path)
    p.add_argument("--min-recall", type=float, default=0.90)
    p.add_argument("--min-precision", type=float, default=0.90)
    p.add_argument("--min-gt-concordance", type=float, default=0.85)
    return p.parse_args()


def dosage(gt: str) -> int | None:
    sep = "|" if "|" in gt else "/"
    alleles = gt.split(":")[0].split(sep)
    if any(a in (".", "") for a in alleles):
        return None
    return sum(1 for a in alleles if a != "0")


def load_vcf(path: Path, pass_only: bool = False) -> tuple[list[str], dict]:
    """Return (sample_ids, {(chrom,pos,ref,alt): {sample: dosage}})."""
    samples: list[str] = []
    sites: dict[tuple[str, int, str, str], dict[str, int | None]] = {}
    for line in path.read_text().splitlines():
        if line.startswith("##"):
            continue
        f = line.split("\t")
        if line.startswith("#CHROM"):
            samples = f[9:]
            continue
        if len(f) < 10:
            continue
        chrom, pos, ref, alt, filt, fmt = f[0], int(f[1]), f[3], f[4], f[6], f[8]
        if pass_only and filt not in (".", "PASS"):
            continue
        gt_idx = fmt.split(":").index("GT")
        key = (chrom, pos, ref, alt)
        per = {}
        for sid, s in zip(samples, f[9:]):
            sub = s.split(":")
            per[sid] = dosage(sub[gt_idx]) if gt_idx < len(sub) else None
        sites[key] = per
    return samples, sites


def has_alt(per: dict[str, int | None]) -> bool:
    return any(d is not None and d > 0 for d in per.values())


def main() -> None:
    args = parse_args()
    _, truth = load_vcf(args.truth)
    _, called = load_vcf(args.called, pass_only=True)

    truth_sites = {k: v for k, v in truth.items() if has_alt(v)}
    called_sites = {k: v for k, v in called.items() if has_alt(v)}

    tp = [k for k in truth_sites if k in called_sites]
    recall = len(tp) / len(truth_sites) if truth_sites else 0.0
    precision = len(tp) / len(called_sites) if called_sites else 0.0

    # Genotype concordance over the recovered sites (samples present in both).
    match = total = 0
    for k in tp:
        for sid, td in truth_sites[k].items():
            cd = called_sites[k].get(sid)
            if td is None or cd is None:
                continue
            total += 1
            match += int(td == cd)
    gt_conc = match / total if total else 0.0

    missed = sorted(set(truth_sites) - set(called_sites))
    print(f"truth sites (>=1 alt): {len(truth_sites)}")
    print(f"called PASS sites    : {len(called_sites)}")
    print(f"recovered (TP)       : {len(tp)}")
    print(f"recall    = {recall:.3f}  (>= {args.min_recall})")
    print(f"precision = {precision:.3f}  (>= {args.min_precision})")
    print(f"GT concord= {gt_conc:.3f}  (>= {args.min_gt_concordance})  "
          f"[{match}/{total}]")
    if missed:
        print(f"missed sites: {missed[:10]}{' ...' if len(missed) > 10 else ''}")

    ok = (recall >= args.min_recall and precision >= args.min_precision
          and gt_conc >= args.min_gt_concordance)
    if not ok:
        print("SELF-CHECK FAILED", file=sys.stderr)
        sys.exit(1)
    print("SELF-CHECK PASSED")


if __name__ == "__main__":
    main()
