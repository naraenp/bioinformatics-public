# bioinformatics-public

Portfolio of end-to-end bioinformatics projects covering single-cell RNA-seq,
TCR repertoire analysis, and translational interpretation workflows.

## What's inside

This repository demonstrates production-style scientific computation across two
independent subprojects:

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

## Reproducibility notes

- No user-specific absolute paths. Scripts derive a project-relative working
  directory so collaborators can run them without edits.
- Dependency versions are pinned in each `envs/` definition.
- Raw data provenance is documented in each subproject README; large or
  experimental data files are excluded from version control via `.gitignore`.

## Getting started

See each subproject README for data expectations, environment creation steps,
run instructions, and output interpretation.
