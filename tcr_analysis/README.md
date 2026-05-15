# TCR clonality analysis

Bulk TCR-β repertoire analysis comparing **TB progressors vs. controllers**,
testing whether repertoire clonality and *M. tuberculosis*-specific clonotype
frequency differ between cohorts. Source data is drawn from [Musvosvi et al.,
*Nature Medicine* (2022)](https://doi.org/10.1038/s41591-022-02110-9);
methodology mirrors prior work published in *Frontiers in Immunology*.

## Why it matters

TCR diversity and clonality provide interpretable immune signatures linked to
disease progression and vaccine response. Antigen-specific clonotype tracking
adds a second axis beyond global repertoire shape.

## Main contents

- `TCR_analysis.ipynb`: end-to-end exploratory and analytical workflow.
- `datasets/`: experimental sample data and metadata.
- `envs/tcr_env.yml`: pinned conda environment.
- `envs/setup_env.sh`: helper script for environment creation and package setup.

## Setup

From `tcr_analysis/envs/`:

```bash
mamba env create -f tcr_env.yml
conda activate tcr_env
```

If needed, run the optional package setup:

```bash
bash setup_env.sh
```

## Engineering workflow

1. Keep exploratory work in notebooks.
2. Promote stable analysis logic into versioned scripts or functions.
3. Document parameter choices and output files for reproducibility.
