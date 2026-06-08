# aml_rnaseq_nf

A small **Nextflow DSL2** pipeline for a bulk RNA-seq differential-expression
comparison: **AML vs. healthy**. It runs on real public RNA-seq cohorts —
**TCGA-LAML** (AML) and **GTEx whole blood** (healthy) — pulled from the
[recount3 project](https://rna.recount.bio/), which re-aligns and re-quantifies
both sources through one uniform [Monorail](https://github.com/langmead-lab/monorail-external)
pipeline (STAR → GENCODE v26 gene sums) so the counts are directly comparable.

It's a compact workflow-engineering exercise: channels, processes,
`publishDir`, and profile-driven config on top of a transparent, dependency-light
biology layer (library-size CPM, Welch t-test, hand-rolled BH FDR, interactive
volcano).

## Stages

| # | Process            | Script                     | Output                          |
|---|--------------------|----------------------------|---------------------------------|
| 1 | `LOAD_COUNTS`      | `bin/load_counts.py`       | `counts_raw.tsv`, `metadata.tsv` |
| 2 | `NORMALIZE_COUNTS` | `bin/normalize_counts.py`  | `counts_lcpm.tsv`               |
| 3 | `RUN_DE`           | `bin/run_de.py`            | `de_results.tsv`                |
| 4 | `MAKE_VOLCANO`     | `bin/make_volcano.py`      | `volcano.html`                  |

All outputs land in `results/` via Nextflow `publishDir 'copy'`.

## Run

```bash
cd aml_rnaseq_nf
mamba env create -f envs/aml_rnaseq_env.yml   # or: conda env create -f ...
conda activate aml_rnaseq_env

# One-time data fetch (~130 MB into data/real/, gitignored):
bash fetch_real_data.sh

nextflow run main.nf -profile local           # or: -profile conda
```

`fetch_real_data.sh` downloads the two recount3 gene_sums files and the
GENCODE v26 basic annotation GTF (used to map Ensembl gene IDs to HGNC
symbols). To run the same four steps without Nextflow — handy for quick
iteration or a CI smoke test — use `bash run_local.sh`.

## Parameters

Override at the Nextflow CLI (`--param value`) or in `nextflow.config`.

| Parameter      | Default | Meaning                                       |
|----------------|---------|-----------------------------------------------|
| `outdir`       | `${projectDir}/results` | Per-run output directory      |
| `n_per_group`  | `50`    | Samples drawn per cohort (AML / healthy)      |
| `seed`         | `42`    | RNG seed for the balanced subsample           |
| `fdr`          | `0.05`  | BH FDR threshold for the volcano              |
| `lfc`          | `1.0`   | \|log2FC\| threshold for the volcano          |

## Methods

- **Loading.** Inner-join TCGA-LAML and GTEx-BLOOD gene sums on Ensembl ID,
  map to HGNC symbols via GENCODE v26, subsample `n_per_group` per cohort,
  and keep genes with `CPM ≥ 1` in ≥ 25% of samples (edgeR-style filter).
- **Normalization.** Per-sample CPM on library sizes, then `log2(CPM + 1)`.
- **DE test.** Per-gene Welch's t-test on log2-CPM
  (`scipy.stats.ttest_ind`, `equal_var=False`); p-values adjusted by
  Benjamini–Hochberg (hand-rolled to keep the dependency surface small).
- **Volcano.** Plotly `Scattergl` colored by significance, with canonical
  AML marker genes (FLT3, KIT, MEIS1, HOXA9, MPO, CD34, …) labeled.

See [`docs/REPORT.md`](docs/REPORT.md) for an end-to-end run on
50 AML vs. 50 healthy samples — provenance, the embedded volcano, and a
runtime profile.

> **Comparator caveat.** GTEx has no bone-marrow tissue, so whole peripheral
> blood is the closest large healthy bulk-RNA-seq comparator. Canonical AML
> markers are still cleanly recovered, but progenitor-associated genes can
> read as "up in AML" simply because mature blood lacks progenitor
> populations — see the report. Swapping in a healthy bone-marrow cohort
> (e.g. BLUEPRINT, GSE74246) is the obvious next step.

## Layout

```text
aml_rnaseq_nf/
├── main.nf                  # Nextflow DSL2 workflow
├── nextflow.config          # manifest, params, profiles
├── fetch_real_data.sh       # one-time recount3 + GENCODE download
├── run_local.sh             # mirrors main.nf without Nextflow
├── bin/
│   ├── load_counts.py       # TCGA-LAML + GTEx loader (recount3)
│   ├── normalize_counts.py
│   ├── run_de.py
│   └── make_volcano.py
├── data/real/               # populated by fetch_real_data.sh (gitignored)
├── docs/
│   ├── REPORT.md            # narrative run report
│   └── volcano.png          # static volcano (latest run)
├── envs/
│   └── aml_rnaseq_env.yml   # pinned conda env
├── tests/
│   └── test_de.py           # fast, data-free unit tests
└── results/                 # populated at runtime (gitignored)
```

## Reproducibility

- A single RNG seed controls the balanced subsample, so runs are
  deterministic on the same env + data files.
- No user-specific absolute paths; everything is project-relative.
- `nextflow.config` sets `PYTHONNOUSERSITE=1` so process environments don't
  leak the user site-packages.
- Pinned conda env in `envs/aml_rnaseq_env.yml`.
