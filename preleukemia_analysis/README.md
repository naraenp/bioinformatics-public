## Preleukemia scRNA-seq Analysis

### Purpose

This project analyzes public single-cell RNA-seq data and related metadata to
investigate preleukemic cell populations and their potential disease trajectories.

### Why it matters

Early characterization of preleukemic signatures can support downstream risk
stratification and translational hypothesis generation.

### Main contents

- `preleuk_analysis.R`: main Seurat integration/annotation workflow.
- `preleuk_analysis1.ipynb` and `preleuk_analysis2.ipynb`: additional Python-based analysis.
- `preleuk_dashboard/`: Shiny UI for interactive result exploration.
- `envs/`: conda environment definitions.

### Setup

Create the relevant conda environment from `envs/`:

```bash
conda env create -f preleuk_env1.yml
conda env create -f preleuk_env2.yml
```

### Reproducibility note

The R workflow now derives a project-relative working directory to avoid
user-specific path assumptions.

### Suggested next engineering step

Extract stable analysis blocks from notebooks into versioned scripts with
parameterized CLI inputs for repeatable pipeline runs.
