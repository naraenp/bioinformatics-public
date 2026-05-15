# AML scRNA-seq analysis

Reproducible scRNA-seq pipeline on **38 public AML patient samples**, performing
QC, anchor-based integration, and reference-guided cell-type annotation to
characterize pre-leukemic populations. Extended with pseudotime, fate mapping,
and survival analysis of PLPS/Stem11 signatures against NCI clinical data, and
deployed as an interactive R Shiny dashboard. Source data is drawn from [Zeng
et al., *Cell Genomics* (2023)](https://doi.org/10.1016/j.xgen.2023.100426) and
AML clinical data from NCI.

## Why it matters

Early characterization of pre-leukemic signatures supports risk stratification
and translational hypothesis generation in AML.

## Main contents

- `preleuk_analysis.R`: primary Seurat integration and annotation workflow.
- `preleuk_analysis1.ipynb`, `preleuk_analysis2.ipynb`: downstream Python
  analysis (pseudotime, fate mapping, metabolic activity, survival).
- `preleuk_dashboard/`: R Shiny app for interactive result exploration.
- `envs/`: conda environment definitions for each notebook.
- `packages.R`: R dependency installer for the Seurat workflow.

## Setup

Create the relevant conda environment from `envs/`:

```bash
conda env create -f envs/preleuk_env1.yml
conda env create -f envs/preleuk_env2.yml
```

Install R dependencies, then run the Seurat workflow:

```bash
Rscript packages.R
Rscript preleuk_analysis.R
```

## Reproducibility notes

The R workflow derives a project-relative working directory to avoid
user-specific path assumptions. Pin dependency versions when adding to env
files.

## Suggested next engineering step

Extract stable analysis blocks from the notebooks into versioned scripts with
parameterized CLI inputs for repeatable pipeline runs.
