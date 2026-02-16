## TCR Analysis Project

### Purpose

This project analyzes T-cell receptor (TCR) repertoire characteristics to study
immune response patterns across participant groups.

### Why it matters

TCR diversity and clonality can provide interpretable immune signatures linked
to disease progression and vaccine response.

### Main contents

- `TCR_analysis.ipynb`: exploratory and analytical workflow.
- `datasets/`: experimental sample data and metadata.
- `envs/tcr_env.yml`: environment specification.
- `envs/setup_env.sh`: helper script for environment creation and package setup.

### Setup

From `tcr_analysis/envs/`:

```bash
mamba env create -f tcr_env.yml
conda activate tcr_env
```

If needed, run:

```bash
bash setup_env.sh
```

### Suggested engineering workflow

1. keep exploratory work in notebooks,
2. promote stable analysis logic into scripts/functions,
3. document parameter choices and output files for reproducibility.
