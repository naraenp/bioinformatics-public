#!/usr/bin/env python3
"""Generate a tiny self-contained dataset that exercises the whole pipeline.

The real pipeline runs on a rice genome and SRA FASTQ (see fetch_real_data.sh),
which is multi-GB and slow. For CI, local iteration, and the smoke test we
instead synthesize a *toy genome* plus paired-end reads with a known, planted
tolerant-vs-susceptible signal. The reads are exact substrings of the genome,
so HISAT2 aligns them cleanly and the downstream stages (featureCounts ->
pydeseq2 -> ORA -> plots) run exactly as they would on real data.

Writes into --outdir:
    genome.fa            toy reference (a few chromosomes of toy genes)
    genes.gtf            one exon per gene  (featureCounts -t exon -g gene_id)
    go_sets.gmt          toy GO sets; the planted-DE genes cluster in a few
    samplesheet.csv      sample_id, genotype, fastq_1, fastq_2
    reads/<id>_R{1,2}.fastq.gz
"""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np

BASES = np.array(list("ACGT"))
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--genes", type=int, default=60)
    p.add_argument("--gene-len", type=int, default=800)
    p.add_argument("--chroms", type=int, default=2)
    p.add_argument("--reps", type=int, default=3, help="replicates per genotype")
    p.add_argument("--read-len", type=int, default=75)
    p.add_argument("--base-depth", type=int, default=60,
                   help="mean read pairs per gene per sample")
    p.add_argument("--n-de", type=int, default=18, help="planted DE genes")
    p.add_argument("--lfc", type=float, default=2.0, help="planted log2 fold change")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def revcomp(s: str) -> str:
    return s.translate(COMPLEMENT)[::-1]


def main() -> None:
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    out = args.outdir
    (out / "reads").mkdir(parents=True, exist_ok=True)

    gene_ids = [f"PLANT{i:04d}" for i in range(args.genes)]
    gene_seq = {g: "".join(rng.choice(BASES, size=args.gene_len)) for g in gene_ids}

    # Lay genes out across chromosomes with fixed intergenic spacers.
    spacer = "N" * 100
    layout: dict[str, list[tuple[str, int, int]]] = {}   # chrom -> (gene, start, end) 1-based
    chrom_seq: dict[str, list[str]] = {f"chr{c+1}": [] for c in range(args.chroms)}
    pos = {f"chr{c+1}": 1 for c in range(args.chroms)}
    for i, g in enumerate(gene_ids):
        chrom = f"chr{i % args.chroms + 1}"
        chrom_seq[chrom].append(gene_seq[g])
        start = pos[chrom]
        end = start + args.gene_len - 1
        layout.setdefault(chrom, []).append((g, start, end))
        chrom_seq[chrom].append(spacer)
        pos[chrom] = end + len(spacer) + 1

    with open(out / "genome.fa", "w") as fh:
        for chrom, parts in chrom_seq.items():
            seq = "".join(parts)
            fh.write(f">{chrom}\n")
            for j in range(0, len(seq), 70):
                fh.write(seq[j:j + 70] + "\n")

    with open(out / "genes.gtf", "w") as fh:
        for chrom, entries in layout.items():
            for g, start, end in entries:
                attr = f'gene_id "{g}"; transcript_id "{g}.1";'
                fh.write(f"{chrom}\tdemo\texon\t{start}\t{end}\t.\t+\t.\t{attr}\n")

    # Planted DE: half up in tolerant, half down. log2FC sign per gene.
    de_genes = list(rng.choice(gene_ids, size=args.n_de, replace=False))
    direction = {g: (1.0 if k < args.n_de // 2 else -1.0)
                 for k, g in enumerate(de_genes)}
    base = {g: rng.uniform(0.6, 1.4) for g in gene_ids}   # per-gene baseline factor

    samples = ([(f"tol_r{r+1}", "tolerant") for r in range(args.reps)]
               + [(f"sus_r{r+1}", "susceptible") for r in range(args.reps)])

    sheet = ["sample_id,genotype,fastq_1,fastq_2"]
    rl = args.read_len
    for sid, geno in samples:
        r1_path = out / "reads" / f"{sid}_R1.fastq.gz"
        r2_path = out / "reads" / f"{sid}_R2.fastq.gz"
        with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
            for g in gene_ids:
                mu = args.base_depth * base[g]
                if g in de_genes and geno == "tolerant":
                    mu *= 2.0 ** (direction[g] * args.lfc)
                n_pairs = int(rng.poisson(mu))
                seq = gene_seq[g]
                frag = min(200, len(seq))
                for k in range(n_pairs):
                    s = int(rng.integers(0, len(seq) - frag + 1))
                    fragment = seq[s:s + frag]
                    r1 = fragment[:rl]
                    r2 = revcomp(fragment[-rl:])
                    name = f"{sid}:{g}:{k}"
                    q = "I" * rl
                    f1.write(f"@{name}/1\n{r1}\n+\n{q}\n")
                    f2.write(f"@{name}/2\n{r2}\n+\n{q}\n")
        sheet.append(f"{sid},{geno},reads/{r1_path.name},reads/{r2_path.name}")

    (out / "samplesheet.csv").write_text("\n".join(sheet) + "\n")

    # Toy GO sets. Concentrate the up-in-tolerant DE genes into one process so
    # the enrichment step has a real, recoverable signal.
    up = [g for g in de_genes if direction[g] > 0]
    down = [g for g in de_genes if direction[g] < 0]
    others = [g for g in gene_ids if g not in de_genes]
    half = len(others) // 2
    gmt = [
        ("GO:0009414", "response to water deprivation", up + others[:4]),
        ("GO:0006970", "response to osmotic stress", down + others[4:8]),
        ("GO:0015979", "photosynthesis", others[8:8 + half]),
        ("GO:0006811", "ion transport", others[8 + half:]),
        ("GO:0008150", "biological process (all genes)", gene_ids),
    ]
    with open(out / "go_sets.gmt", "w") as fh:
        for term, name, genes in gmt:
            fh.write("\t".join([term, name, *genes]) + "\n")

    print(f"demo data -> {out}")
    print(f"  {args.genes} genes on {args.chroms} chroms, "
          f"{len(samples)} samples ({args.reps}+{args.reps})")
    print(f"  {args.n_de} planted DE genes (|log2FC|={args.lfc}); "
          f"up-in-tolerant set enriches GO:0009414")


if __name__ == "__main__":
    main()
