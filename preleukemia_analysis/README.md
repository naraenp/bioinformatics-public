# Preleukemia scRNA-seq analysis

Reproducible scRNA-seq analysis of **38 mouse bone-marrow HSPC samples**
spanning eight preleukemic mutation models, performing QC, anchor-based
integration, and reference-guided cell-type annotation to characterize
preleukemic populations. Extended with pseudotime, fate mapping, and a survival
analysis of the paper's PLPS and Stem11 signatures in the TCGA-LAML cohort, and
deployed as an interactive R Shiny dashboard. Source data is drawn from [Isobe
et al., *Cell Genomics* (2023)](https://doi.org/10.1016/j.xgen.2023.100426);
the human survival arm uses 163 TCGA-LAML patients with both expression and
overall-survival data (NCI clinical data via cBioPortal).

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
