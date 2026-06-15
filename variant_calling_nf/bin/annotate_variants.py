#!/usr/bin/env python3
"""Annotate a (normalized, biallelic) cohort VCF with a genic consequence.

Deliberately a small, transparent, hand-rolled annotator rather than a heavy
database tool (snpEff/VEP): it intersects each variant against the gene/exon
intervals in a GTF and classifies it as exonic / intronic / intergenic, tags
SNVs as transitions or transversions, and writes a tidy per-variant table plus
a per-sample dosage matrix for the downstream plots. The pure logic
(classification, Ts/Tv, GT parsing) is unit-tested in tests/.

Inputs:  --vcf cohort.filtered.vcf  --gtf genes.gtf
Outputs: --out annotated.tsv   (one row per variant: site, kind, ts/tv, region,
                                 filter, then one dosage column per sample)
"""
from __future__ import annotations

import argparse
from pathlib import Path

PURINES = {"A", "G"}
PYRIMIDINES = {"C", "T"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", required=True, type=Path)
    p.add_argument("--gtf", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    return p.parse_args()


def variant_kind(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    if len(alt) > len(ref):
        return "ins"
    if len(ref) > len(alt):
        return "del"
    return "mnv"


def ts_or_tv(ref: str, alt: str) -> str:
    """Transition (purine<->purine / pyrimidine<->pyrimidine) vs transversion."""
    if len(ref) != 1 or len(alt) != 1:
        return ""
    if (ref in PURINES and alt in PURINES) or (ref in PYRIMIDINES and alt in PYRIMIDINES):
        return "ts"
    return "tv"


def gt_dosage(gt: str) -> int | None:
    """Alt-allele dosage from a GT field; None for missing."""
    sep = "|" if "|" in gt else "/"
    alleles = gt.split(":")[0].split(sep)
    if any(a in (".", "") for a in alleles):
        return None
    return sum(1 for a in alleles if a != "0")


def read_intervals(gtf: Path) -> tuple[dict[str, list[tuple[int, int]]],
                                       dict[str, list[tuple[int, int]]]]:
    """Return (exon_intervals, gene_spans) per chrom, 1-based inclusive.

    Exon intervals come straight from `exon` features; gene spans are the
    min-start/max-end envelope of each gene_id's exons (so a variant inside the
    transcript footprint but outside an exon reads as intronic).
    """
    exons: dict[str, list[tuple[int, int]]] = {}
    gene_bounds: dict[tuple[str, str], list[int]] = {}
    for line in gtf.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        chrom, start, end, attr = f[0], int(f[3]), int(f[4]), f[8]
        exons.setdefault(chrom, []).append((start, end))
        gid = ""
        for field in attr.split(";"):
            field = field.strip()
            if field.startswith("gene_id"):
                gid = field.split('"')[1] if '"' in field else field.split()[-1]
                break
        key = (chrom, gid)
        if key not in gene_bounds:
            gene_bounds[key] = [start, end]
        else:
            gene_bounds[key][0] = min(gene_bounds[key][0], start)
            gene_bounds[key][1] = max(gene_bounds[key][1], end)
    genes: dict[str, list[tuple[int, int]]] = {}
    for (chrom, _gid), (s, e) in gene_bounds.items():
        genes.setdefault(chrom, []).append((s, e))
    return exons, genes


def classify_region(chrom: str, pos: int,
                    exons: dict[str, list[tuple[int, int]]],
                    genes: dict[str, list[tuple[int, int]]]) -> str:
    if any(s <= pos <= e for s, e in exons.get(chrom, [])):
        return "exonic"
    if any(s <= pos <= e for s, e in genes.get(chrom, [])):
        return "intronic"
    return "intergenic"


def parse_vcf(vcf: Path):
    """Yield (chrom, pos, ref, alt, filt, sample_ids, [dosages]) per record.

    Assumes a biallelic, normalized VCF (one ALT per line); multiallelic sites
    should be split upstream with `bcftools norm -m-`.
    """
    samples: list[str] = []
    for line in vcf.read_text().splitlines():
        if line.startswith("##"):
            continue
        f = line.split("\t")
        if line.startswith("#CHROM"):
            samples = f[9:]
            continue
        if len(f) < 10:
            continue
        chrom, pos, ref, alt, filt, fmt = f[0], int(f[1]), f[3], f[4], f[6], f[8]
        gt_idx = fmt.split(":").index("GT")
        dosages = []
        for s in f[9:]:
            sub = s.split(":")
            dosages.append(gt_dosage(sub[gt_idx]) if gt_idx < len(sub) else None)
        yield chrom, pos, ref, alt, filt, samples, dosages


def main() -> None:
    args = parse_args()
    exons, genes = read_intervals(args.gtf)

    rows, samples = [], []
    for chrom, pos, ref, alt, filt, samples, dosages in parse_vcf(args.vcf):
        kind = variant_kind(ref, alt)
        rows.append(dict(
            chrom=chrom, pos=pos, ref=ref, alt=alt, kind=kind,
            tstv=ts_or_tv(ref, alt),
            region=classify_region(chrom, pos, exons, genes),
            filter=filt if filt else ".",
            dosages=["" if d is None else str(d) for d in dosages],
        ))

    header = ["chrom", "pos", "ref", "alt", "kind", "tstv", "region", "filter",
              *samples]
    with open(args.out, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join([r["chrom"], str(r["pos"]), r["ref"], r["alt"],
                                r["kind"], r["tstv"], r["region"], r["filter"],
                                *r["dosages"]]) + "\n")

    n_pass = sum(1 for r in rows if r["filter"] in (".", "PASS"))
    ts = sum(1 for r in rows if r["tstv"] == "ts")
    tv = sum(1 for r in rows if r["tstv"] == "tv")
    by_region: dict[str, int] = {}
    for r in rows:
        by_region[r["region"]] = by_region.get(r["region"], 0) + 1
    print(f"wrote {args.out}: {len(rows)} variants ({n_pass} PASS), "
          f"{len(samples)} samples")
    print(f"  regions: {by_region}")
    print(f"  Ts/Tv = {ts}/{tv} = {ts / tv:.2f}" if tv else f"  Ts/Tv = {ts}/0")


if __name__ == "__main__":
    main()
