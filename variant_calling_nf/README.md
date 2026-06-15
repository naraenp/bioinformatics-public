# variant_calling_nf

A **Nextflow DSL2** pipeline that takes **raw Illumina short reads** through the
**GATK germline short-variant best-practice path** to a filtered, annotated
**cohort VCF** and two interactive figures. It is the **DNA-variant** member of
the portfolio — where the sibling pipelines quantify *expression*
([`aml_rnaseq_nf`](../aml_rnaseq_nf), [`plant_rnaseq_nf`](../plant_rnaseq_nf)) or
map cell types in space ([`spatial_visium_nf`](../spatial_visium_nf)), this one
genotypes **SNVs and short indels** from the genome itself.

The default real target is the **Genome in a Bottle** sample **HG002** over a
slice of GRCh38 chr20, benchmarked against the **NIST v4.2.1** high-confidence
truth set — but nothing is HG002-specific: point it at any genome, GTF, and
samplesheet.

## Stages

| #  | Process            | Tool / script                  | Output                          |
|----|--------------------|--------------------------------|---------------------------------|
| 1  | `SUBSAMPLE`        | `bin/subsample_reads.sh` (seqtk) | downsampled FASTQ             |
| 2  | `QC_TRIM`          | `fastp`                        | trimmed FASTQ, QC reports       |
| 3  | `PREP_REFERENCE`   | `samtools faidx` + GATK dict + `bwa index` | indexed reference   |
| 4  | `ALIGN`            | `bwa mem` (+ read groups) + `samtools` | sorted, indexed BAM     |
| 5  | `MARK_DUPLICATES`  | `gatk MarkDuplicates`          | deduplicated BAM                |
| 6  | `CALL_VARIANTS`    | `gatk HaplotypeCaller -ERC GVCF` | per-sample GVCF               |
| 7  | `JOINT_GENOTYPE`   | `gatk CombineGVCFs` + `GenotypeGVCFs` | `cohort.vcf.gz`          |
| 8  | `FILTER_VARIANTS`  | `gatk VariantFiltration` (SNP + indel hard filters) | `cohort.filtered.vcf.gz` |
| 9  | `NORMALIZE`        | `bcftools norm -m-`            | `cohort.norm.vcf`               |
| 10 | `ANNOTATE`         | `bin/annotate_variants.py`     | `annotated.tsv` (region + Ts/Tv)|
| 11 | `MAKE_OVERVIEW`    | `bin/make_overview.py`         | `variant_overview.html`         |
| 12 | `MAKE_GENO_HEATMAP`| `bin/make_geno_heatmap.py`     | `genotype_heatmap.html`         |

All outputs land in `results/` via Nextflow `publishDir 'copy'`.

## Run

```bash
cd variant_calling_nf
mamba env create -f envs/variant_calling_env.yml   # or: conda env create -f ...
conda activate variant_calling_env

# Offline demo on a synthetic toy genome + planted variants (under a minute):
bash run_local.sh --demo
#   or with Nextflow:
nextflow run main.nf -profile conda \
    --samplesheet data/demo/samplesheet.csv --genome data/demo/genome.fa \
    --gtf data/demo/genes.gtf
```

`run_local.sh --demo` ends with a **self-check**: it normalizes the called
cohort VCF and the planted `truth.vcf` and asserts the pipeline recovered the
planted variants (recall / precision / genotype concordance) within tolerance.

For a **real run**, fetch the HG002 / GRCh38-chr20 inputs once and point the
pipeline at them:

```bash
bash fetch_real_data.sh          # GRCh38 chr20 + GIAB HG002 region reads + NIST truth
nextflow run main.nf -profile conda \
    --samplesheet data/real/samplesheet.csv --genome data/real/genome.fa \
    --gtf data/real/genes.gtf --intervals data/real/intervals.bed
```

Reads are **range-sliced** straight out of GIAB's chr20 300x BAM over HTTPS
(only the bytes for the target region are fetched), so a real germline-calling +
benchmarking run stays laptop-scale. `fetch_real_data.sh` prints the
`bcftools isec` commands that score the callset against the NIST truth.

### Containers

Besides `-profile conda`, the pipeline ships `docker` and `singularity`
profiles. The `Dockerfile` bakes `envs/variant_calling_env.yml` (the alignment +
GATK toolchain) into a
[micromamba](https://hub.docker.com/r/mambaorg/micromamba) image, so the
container matches the conda profile exactly:

```bash
docker build -t variant_calling_nf:0.1.0 .       # one-time image build
nextflow run main.nf -profile docker \           # or: -profile singularity
    --samplesheet data/demo/samplesheet.csv --genome data/demo/genome.fa \
    --gtf data/demo/genes.gtf
```

## Parameters

Override at the Nextflow CLI (`--param value`) or in `nextflow.config`.

| Parameter      | Default | Meaning                                            |
|----------------|---------|----------------------------------------------------|
| `samplesheet`  | `data/samplesheet.csv` | CSV: `sample_id,group,fastq_1,fastq_2` |
| `genome`       | `data/genome.fa`       | reference FASTA                     |
| `gtf`          | `data/genes.gtf`       | gene annotation (for the consequence call) |
| `intervals`    | *(empty)* | optional BED restricting where variants are called |
| `subsample`    | `0`     | read pairs/sample (0 = all)                        |
| `seed`         | `42`    | subsampling RNG seed                               |
| `ploidy`       | `2`     | sample ploidy for HaplotypeCaller                  |
| `topn`         | `40`    | variants shown in the genotype heatmap             |

## Methods

- **QC / trim.** `fastp` removes adapters and low-quality tails.
- **Alignment.** `bwa mem` to an indexed reference with per-sample read groups
  (`@RG … SM:`, required by GATK), sorted and indexed with `samtools`.
- **Duplicate marking.** `gatk MarkDuplicates` flags optical/PCR duplicates.
- **Calling.** `gatk HaplotypeCaller` in **GVCF mode** per sample (local
  de-novo assembly of haplotypes), then **joint genotyping** across the cohort
  with `CombineGVCFs` + `GenotypeGVCFs` — the GATK germline best practice.
- **Filtering.** Hard filters applied to **separate SNP and indel tracks**
  (`QD`, `FS`, `MQ`, `SOR` for SNPs; `QD`, `FS`, `SOR` for indels), the
  documented substitute for VQSR on small cohorts.
- **Normalization.** `bcftools norm -m-` splits multiallelics and left-aligns
  indels so variants are keyed canonically for annotation and benchmarking.
- **Annotation.** A small, transparent, hand-rolled annotator
  (`annotate_variants.py`) intersects each variant against the GTF
  (exonic / intronic / intergenic) and tags SNVs as transitions or
  transversions — the cohort **Ts/Tv** is a standard callset-QC read-out.
- **Visualization.** A per-sample variant-overview bar chart and a genotype
  heatmap of the most discriminating variants — both interactive Plotly in the
  naraen.net "Aequorea" abyssal-marine palette.

See [`docs/REPORT.md`](docs/REPORT.md) for an end-to-end run.

## Layout

```text
variant_calling_nf/
├── main.nf                  # Nextflow DSL2 workflow (12 processes)
├── nextflow.config          # manifest, params, profiles
├── fetch_real_data.sh       # GRCh38 chr20 + GIAB HG002 + NIST truth -> data/real/
├── run_local.sh             # mirrors main.nf without Nextflow (--demo offline)
├── bin/
│   ├── subsample_reads.sh
│   ├── make_demo_data.py    # toy genome + planted-variant reads + truth.vcf
│   ├── annotate_variants.py # genic consequence + Ts/Tv (hand-rolled, tested)
│   ├── check_truth.py       # demo self-check: recall / precision / GT concordance
│   ├── make_overview.py
│   └── make_geno_heatmap.py
├── envs/variant_calling_env.yml
├── tests/test_pipeline.py   # fast unit tests (no aligner / GATK needed)
├── docs/REPORT.md
├── data/                    # inputs (gitignored)
└── results/                 # outputs (gitignored)
```

## Reproducibility

- Deterministic subsampling (single `seed`); all paths project-relative.
- `nextflow.config` sets `PYTHONNOUSERSITE=1` so processes don't leak user
  site-packages.
- Pinned conda env in `envs/variant_calling_env.yml`.
- `run_local.sh --demo` runs the entire DAG offline on a synthetic genome and
  self-checks variant recovery — used as the validation path and to keep
  `main.nf` and `run_local.sh` in sync. The fast unit suite in `tests/` is the
  CI gate.

## Notes / next steps

- The demo plants both group-shared and group-specific variants, so the
  genotype heatmap shows recoverable group structure (the variant-calling analog
  of the planted DE genes in the RNA-seq pipelines).
- Hard-filtering stands in for VQSR/CNN filtering, which need large cohorts or
  pretrained models; both are drop-in upgrades.
- Somatic calling (`Mutect2`, tumor/normal) is a natural sibling of this
  germline path.
