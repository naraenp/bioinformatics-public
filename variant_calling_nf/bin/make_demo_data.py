#!/usr/bin/env python3
"""Generate a tiny self-contained dataset that exercises the whole pipeline.

The real pipeline runs on a human GIAB sample and a slice of GRCh38 (see
fetch_real_data.sh), which is multi-GB and slow. For CI, local iteration, and
the smoke test we instead synthesize a *toy diploid genome* plus paired-end
reads carrying a known, planted set of germline variants. Reads are sampled as
exact substrings of each sample's mutated haplotypes, so BWA-MEM maps them
cleanly and the downstream stages (MarkDuplicates -> HaplotypeCaller -> joint
genotyping -> hard-filter -> annotation -> plots) run exactly as on real data.

The planted variants come in three flavours, giving the cohort a recoverable
group structure (the analog of the planted DE genes / cell-type proportions in
the sibling pipelines):
    shared          present (mostly het) in every sample
    groupA-specific alt only in the group-A samples
    groupB-specific alt only in the group-B samples

Writes into --outdir:
    genome.fa            toy reference (unmutated; what reads are called against)
    genes.gtf            gene/exon layout so variants get a genic consequence
    truth.vcf            planted variants in REFERENCE coords, per-sample GT
    samplesheet.csv      sample_id, group, fastq_1, fastq_2
    reads/<id>_R{1,2}.fastq.gz

`run_local.sh --demo` (and CI's parity check) call check_truth.py to confirm the
pipeline recovers truth.vcf within tolerance.
"""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np

BASES = np.array(list("ACGT"))
COMPLEMENT = str.maketrans("ACGT", "TGCA")
TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}  # purine<->purine etc.


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--chroms", type=int, default=2)
    p.add_argument("--chrom-len", type=int, default=12000)
    p.add_argument("--gene-len", type=int, default=1000)
    p.add_argument("--intergenic", type=int, default=500)
    p.add_argument("--reps", type=int, default=3, help="samples per group")
    p.add_argument("--read-len", type=int, default=100)
    p.add_argument("--frag-len", type=int, default=300)
    p.add_argument("--depth", type=int, default=40, help="approx per-haplotype coverage")
    p.add_argument("--n-shared", type=int, default=12, help="variants in all samples")
    p.add_argument("--n-group", type=int, default=6, help="per-group-specific variants")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def revcomp(s: str) -> str:
    return s.translate(COMPLEMENT)[::-1]


def make_reference(rng, n_chroms: int, length: int) -> dict[str, str]:
    """Random ACGT chromosomes (no N runs, so reads map uniquely)."""
    return {f"chr{c+1}": "".join(rng.choice(BASES, size=length))
            for c in range(n_chroms)}


def layout_genes(length: int, gene_len: int, gap: int) -> list[tuple[int, int]]:
    """Evenly tile (start, end) gene spans (0-based, end-exclusive) along a chrom."""
    spans, pos = [], gap
    while pos + gene_len <= length - gap:
        spans.append((pos, pos + gene_len))
        pos += gene_len + gap
    return spans


def plant_variants(rng, ref: dict[str, str], genes: dict[str, list[tuple[int, int]]],
                   n_shared: int, n_group: int) -> list[dict]:
    """Choose well-separated variant sites and assign each a category + alleles.

    Returns a list of dicts with reference-coordinate (chrom, pos0, ref, alt),
    a variant `kind` (snv/ins/del) and a `category` (shared/groupA/groupB). We
    keep sites >=200 bp apart so BWA-MEM / HaplotypeCaller resolve each cleanly
    and indels never overlap an adjacent variant.
    """
    chroms = list(ref)
    categories = (["shared"] * n_shared
                  + ["groupA"] * n_group + ["groupB"] * n_group)
    rng.shuffle(categories)

    # Bias most variants into exons so the consequence annotation has signal,
    # but leave a few intergenic.
    variants: list[dict] = []
    used: dict[str, list[int]] = {c: [] for c in chroms}

    def far_enough(chrom: str, pos: int) -> bool:
        return all(abs(pos - q) >= 200 for q in used[chrom])

    for i, cat in enumerate(categories):
        chrom = chroms[i % len(chroms)]
        seq = ref[chrom]
        exonic = (i % 4 != 0)            # ~75% exonic
        for _ in range(200):             # rejection-sample a clean site
            if exonic and genes[chrom]:
                gs, ge = genes[chrom][rng.integers(0, len(genes[chrom]))]
                pos = int(rng.integers(gs + 10, ge - 10))
            else:
                pos = int(rng.integers(50, len(seq) - 50))
            if far_enough(chrom, pos):
                break
        else:
            continue
        used[chrom].append(pos)

        roll = rng.random()
        ref_base = seq[pos]
        if roll < 0.7:                   # SNV
            # alternate transitions/transversions so Ts/Tv is interesting
            if rng.random() < 0.65:
                alt = TRANSITION[ref_base]
            else:
                alt = str(rng.choice([b for b in "ACGT"
                                      if b not in (ref_base, TRANSITION[ref_base])]))
            kind, vref, valt = "snv", ref_base, alt
        elif roll < 0.85:                # insertion (anchor + inserted bases)
            ins = "".join(rng.choice(BASES, size=int(rng.integers(1, 4))))
            kind, vref, valt = "ins", ref_base, ref_base + ins
        else:                            # deletion (anchor + deleted bases)
            ndel = int(rng.integers(1, 4))
            kind, vref, valt = "del", seq[pos:pos + 1 + ndel], ref_base
        variants.append(dict(chrom=chrom, pos0=pos, ref=vref, alt=valt,
                             kind=kind, category=cat))
    variants.sort(key=lambda v: (v["chrom"], v["pos0"]))
    return variants


def genotypes_for(variants: list[dict], samples: list[tuple[str, str]], rng
                  ) -> dict[str, list[int]]:
    """Per-variant genotype (0=hom-ref, 1=het, 2=hom-alt) for each sample."""
    gts: dict[str, list[int]] = {}
    for v in variants:
        row = []
        for _sid, grp in samples:
            if v["category"] == "shared":
                # mostly het, ~25% hom-alt — present in everyone
                g = 2 if rng.random() < 0.25 else 1
            elif v["category"] == grp:
                g = 1                    # group-specific: het in the owning group
            else:
                g = 0                    # absent in the other group
            row.append(g)
        gts[id(v)] = row
    return gts


def build_haplotype(ref_seq: str, applied: list[dict]) -> str:
    """Apply a sorted, non-overlapping list of variants to one reference seq."""
    out, i = [], 0
    for v in sorted(applied, key=lambda x: x["pos0"]):
        out.append(ref_seq[i:v["pos0"]])
        out.append(v["alt"])
        i = v["pos0"] + len(v["ref"])
    out.append(ref_seq[i:])
    return "".join(out)


def simulate_reads(rng, hap_seq: str, depth: int, read_len: int, frag_len: int,
                   sid: str, hap_id: int, f1, f2) -> None:
    """Tile paired-end reads across one haplotype to ~`depth` coverage."""
    if len(hap_seq) < frag_len:
        return
    n_pairs = max(1, depth * len(hap_seq) // (2 * read_len))
    q = "I" * read_len
    for k in range(n_pairs):
        s = int(rng.integers(0, len(hap_seq) - frag_len + 1))
        frag = hap_seq[s:s + frag_len]
        r1 = frag[:read_len]
        r2 = revcomp(frag[-read_len:])
        name = f"{sid}:h{hap_id}:{k}"
        f1.write(f"@{name}/1\n{r1}\n+\n{q}\n")
        f2.write(f"@{name}/2\n{r2}\n+\n{q}\n")


def write_truth_vcf(path: Path, ref: dict[str, str], variants: list[dict],
                    samples: list[tuple[str, str]], gts: dict[str, list[int]]) -> None:
    sample_ids = [s for s, _ in samples]
    lines = ["##fileformat=VCFv4.2"]
    for chrom, seq in ref.items():
        lines.append(f"##contig=<ID={chrom},length={len(seq)}>")
    lines.append('##INFO=<ID=CAT,Number=1,Type=String,Description="planted category">')
    lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
              "FORMAT", *sample_ids]
    lines.append("\t".join(header))
    code = {0: "0/0", 1: "0/1", 2: "1/1"}
    for v in variants:
        row = [v["chrom"], str(v["pos0"] + 1), ".", v["ref"], v["alt"], "100",
               "PASS", f"CAT={v['category']}", "GT",
               *[code[g] for g in gts[id(v)]]]
        lines.append("\t".join(row))
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    out = args.outdir
    (out / "reads").mkdir(parents=True, exist_ok=True)

    ref = make_reference(rng, args.chroms, args.chrom_len)
    genes = {c: layout_genes(len(seq), args.gene_len, args.intergenic)
             for c, seq in ref.items()}

    # genome.fa (the UNMUTATED reference)
    with open(out / "genome.fa", "w") as fh:
        for chrom, seq in ref.items():
            fh.write(f">{chrom}\n")
            for j in range(0, len(seq), 70):
                fh.write(seq[j:j + 70] + "\n")

    # genes.gtf — one gene == one exon spanning the gene (toy); enough for the
    # exonic/intronic/intergenic consequence call.
    with open(out / "genes.gtf", "w") as fh:
        gi = 0
        for chrom, spans in genes.items():
            for (gs, ge) in spans:
                gi += 1
                gid = f"GENE{gi:04d}"
                attr = f'gene_id "{gid}"; transcript_id "{gid}.1";'
                # GTF is 1-based, end-inclusive
                fh.write(f"{chrom}\tdemo\texon\t{gs+1}\t{ge}\t.\t+\t.\t{attr}\n")

    variants = plant_variants(rng, ref, genes, args.n_shared, args.n_group)

    samples = ([(f"A_s{r+1}", "groupA") for r in range(args.reps)]
               + [(f"B_s{r+1}", "groupB") for r in range(args.reps)])
    gts = genotypes_for(variants, samples, rng)

    sheet = ["sample_id,group,fastq_1,fastq_2"]
    for si, (sid, grp) in enumerate(samples):
        # Build the two haplotypes for this sample. het -> alt on hap1 only;
        # hom-alt -> alt on both; hom-ref -> neither.
        hap_apply: list[list[dict]] = [[], []]
        for v in variants:
            g = gts[id(v)][si]
            if g >= 1:
                hap_apply[0].append(v)
            if g == 2:
                hap_apply[1].append(v)
        r1_path = out / "reads" / f"{sid}_R1.fastq.gz"
        r2_path = out / "reads" / f"{sid}_R2.fastq.gz"
        with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
            for chrom, seq in ref.items():
                for hap_id in (0, 1):
                    applied = [v for v in hap_apply[hap_id] if v["chrom"] == chrom]
                    hap_seq = build_haplotype(seq, applied)
                    simulate_reads(rng, hap_seq, args.depth // 2, args.read_len,
                                   args.frag_len, sid, hap_id, f1, f2)
        sheet.append(f"{sid},{grp},reads/{r1_path.name},reads/{r2_path.name}")

    (out / "samplesheet.csv").write_text("\n".join(sheet) + "\n")
    write_truth_vcf(out / "truth.vcf", ref, variants, samples, gts)

    n_snv = sum(1 for v in variants if v["kind"] == "snv")
    n_indel = len(variants) - n_snv
    print(f"demo data -> {out}")
    print(f"  {args.chroms} chroms x {args.chrom_len} bp, {len(samples)} samples "
          f"({args.reps} groupA + {args.reps} groupB)")
    print(f"  {len(variants)} planted variants ({n_snv} SNV, {n_indel} indel); "
          f"{args.n_shared} shared, {args.n_group}+{args.n_group} group-specific")


if __name__ == "__main__":
    main()
