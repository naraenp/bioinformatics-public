# bioinformatics-public

[![CI](https://github.com/naraenp/bioinformatics-public/actions/workflows/ci.yml/badge.svg)](https://github.com/naraenp/bioinformatics-public/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Portfolio of end-to-end bioinformatics projects covering single-cell RNA-seq,
TCR repertoire analysis, and translational interpretation workflows.

## What's inside

This repository demonstrates production-style scientific computation across six
independent subprojects: two analysis workflows (TCR repertoire, AML scRNA-seq)
and four Nextflow pipelines (two bulk RNA-seq differential expression, counts-in
and reads-in; one spatial transcriptomics deconvolution; and one germline DNA
short-variant caller):

- Real biological data handled with reproducible, version-pinned environments.
- Robust analysis methods (Seurat integration, anchor-based label transfer,
  tcrdist3 reference matching, clonality metrics).
- Outputs translated into interpretable figures, dashboards, and reports.

Each subproject is self-contained: its own conda environment, data layout, and
runtime instructions.

## Repository layout

```text
bioinformatics-public/
├── tcr_analysis/
│   ├── TCR_analysis.ipynb
│   ├── datasets/
│   ├── envs/
│   └── README.md
├── preleukemia_analysis/
│   ├── preleuk_analysis.R
│   ├── preleuk_analysis1.ipynb
│   ├── preleuk_analysis2.ipynb
│   ├── preleuk_dashboard/
│   ├── envs/
│   └── README.md
├── aml_rnaseq_nf/          # Nextflow bulk RNA-seq DE pipeline
│   ├── main.nf
│   ├── fetch_real_data.sh
│   ├── run_local.sh
│   ├── bin/
│   ├── docs/
│   └── envs/
├── plant_rnaseq_nf/        # Nextflow plant RNA-seq DE pipeline (reads-in)
│   ├── main.nf
│   ├── fetch_real_data.sh
│   ├── run_local.sh
│   ├── bin/
│   ├── docs/
│   └── envs/
├── spatial_visium_nf/      # Nextflow spatial transcriptomics deconvolution pipeline
│   ├── main.nf
│   ├── fetch_real_data.sh
│   ├── run_local.sh
│   ├── bin/
│   ├── docs/
│   └── envs/
├── variant_calling_nf/     # Nextflow germline DNA short-variant calling pipeline
│   ├── main.nf
│   ├── fetch_real_data.sh
│   ├── run_local.sh
│   ├── bin/
│   ├── docs/
│   └── envs/
└── .gitignore
```

## Project summaries

### TCR clonality analysis

Bulk TCR-β repertoire analysis comparing clonality and antigen-specific
clonotype frequencies in TB progressors vs. controllers. Standardizes V/J/CDR3
nomenclature to IMGT format and cross-references experimental repertoires
against curated IEDB and VDJdb reference sets. Environment definitions live in
`tcr_analysis/envs/`.

### AML scRNA-seq analysis

Reproducible scRNA-seq pipeline on 38 public AML patient samples, performing
QC, anchor-based integration, and reference-guided cell-type annotation to
characterize pre-leukemic populations. Extended with pseudotime, fate mapping,
and survival analysis of PLPS/Stem11 signatures against NCI clinical data, and
deployed via an interactive R Shiny dashboard. Environment definitions live in
`preleukemia_analysis/envs/`.

### AML bulk RNA-seq DE pipeline (Nextflow)

Small Nextflow DSL2 pipeline (`aml_rnaseq_nf/`) for a bulk RNA-seq differential
expression comparison (AML vs. healthy). It runs on real public RNA-seq data —
TCGA-LAML and GTEx whole blood, both pulled from
[recount3](https://rna.recount.bio/) for uniform alignment/quantification.
Four stages: `LOAD_COUNTS → NORMALIZE_COUNTS → RUN_DE → MAKE_VOLCANO`, with
`fetch_real_data.sh` for the one-time download and `run_local.sh` mirroring the
steps without Nextflow. See `aml_rnaseq_nf/docs/REPORT.md` for a run report.

### Plant RNA-seq DE pipeline (Nextflow)

Nextflow DSL2 pipeline (`plant_rnaseq_nf/`) taking **raw Illumina reads** from a
drought-tolerant (Apo) and drought-susceptible (IR64) rice cultivar through QC,
spliced **HISAT2 alignment**, `featureCounts` quantification, **pydeseq2**
differential expression (`~condition + genotype`), and hypergeometric **GO
enrichment** to interactive heatmap + enriched-process charts. The reads-in
companion to `aml_rnaseq_nf`: reads are stream-subsampled from ENA and the
reference comes from Ensembl Plants. `fetch_real_data.sh` pulls the data;
`run_local.sh --demo` runs the whole DAG offline on a toy genome as the CI smoke
test. See `plant_rnaseq_nf/docs/REPORT.md` for a run report.

### Spatial transcriptomics deconvolution pipeline (Nextflow)

Nextflow DSL2 pipeline (`spatial_visium_nf/`) taking a **10x Visium** tissue
section plus a matched **scRNA-seq reference** and mapping cell types back into
space: per-spot **NNLS deconvolution** against a reference signature and
**spatially variable genes** by Moran's I, rendered as interactive Plotly figures.
The real target is **10x Visium Human Breast Cancer** with the **Wu et al. 2021**
atlas (GSE176078) as the reference. The tissue-in companion to the bulk pipelines:
`fetch_real_data.sh` pulls the data, and `run_local.sh --demo` runs the whole DAG
offline on a synthetic section with **planted cell-type proportions** that the
deconvolution self-checks against, as the CI smoke test. See
`spatial_visium_nf/docs/REPORT.md` for a run report.

### Germline variant calling pipeline (Nextflow)

Nextflow DSL2 pipeline (`variant_calling_nf/`) taking **raw Illumina short
reads** through the **GATK germline short-variant best-practice path**:
**BWA-MEM** alignment with read groups, `MarkDuplicates`, **HaplotypeCaller** in
GVCF mode, **joint genotyping** (`CombineGVCFs` + `GenotypeGVCFs`), hard-filtering
on separate SNP/indel tracks, `bcftools` normalization, and a hand-rolled genic
consequence + Ts/Tv annotation, to interactive genotype-heatmap + variant-overview
charts. The DNA-variant companion to the expression pipelines: the real target is
**GIAB HG002** over a GRCh38 chr20 slice (reads range-sliced from a remote BAM,
benchmarked against the **NIST v4.2.1** truth via `bcftools isec`). `run_local.sh
--demo` runs the whole DAG offline on a synthetic diploid cohort with **planted
variants** that `check_truth.py` self-checks for recovery (recall / precision /
genotype concordance). See `variant_calling_nf/docs/REPORT.md` for a run report.

## Reproducibility notes

- No user-specific absolute paths. Scripts derive a project-relative working
  directory so collaborators can run them without edits.
- Dependency versions are pinned in each `envs/` definition.
- Raw data provenance is documented in each subproject README; large or
  experimental data files are excluded from version control via `.gitignore`.

## Getting started

See each subproject README for data expectations, environment creation steps,
run instructions, and output interpretation.
