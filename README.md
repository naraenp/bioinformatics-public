# bioinformatics-public

Portfolio repository for bioinformatics projects spanning TCR repertoire work,
single-cell RNA-seq analysis, and translational interpretation workflows.

## Motivation

This repository is meant to show end-to-end scientific computation:

- handling real biological data,
- applying robust analysis methods,
- translating outputs into interpretable figures and reports.

It now also aims to improve engineering legibility by documenting setup,
dependencies, and expected project structure.

## Repository Layout

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

## Project Summaries

### 1) TCR Analysis

- compares T-cell receptor repertoire characteristics across cohorts,
- includes clonality/diversity analysis and reference dataset comparisons,
- primary environment definitions are in `tcr_analysis/envs/`.

### 2) Preleukemia scRNA-seq Analysis

- performs Seurat-based integration/annotation in R,
- includes downstream analytical steps and a Shiny dashboard for exploration,
- environment definitions are in `preleukemia_analysis/envs/`.

## Reproducibility Notes

- do not rely on user-specific absolute paths,
- prefer environment-driven or project-relative paths,
- pin dependency versions where possible for reruns,
- document raw data provenance when sharing publicly.

## Getting Started

See each subproject README for:

- data expectations,
- environment creation steps,
- run instructions,
- output interpretation.
