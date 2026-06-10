# spatial_visium_nf

A **Nextflow DSL2** pipeline that takes a **10x Visium** spatial transcriptomics
section plus a **matched scRNA-seq reference** and maps cell types back into
tissue space: per-spot **cell-type deconvolution** (NNLS against a reference
signature matrix) and **spatially variable genes** (Moran's I), rendered as
interactive Plotly figures. It is the **spatial companion** to the bulk pipelines
[`aml_rnaseq_nf`](../aml_rnaseq_nf) (counts-in) and
[`plant_rnaseq_nf`](../plant_rnaseq_nf) (reads-in): this one is
*tissue-in → cell-types-in-space-out*.

The default real target is **10x Visium Human Breast Cancer** + the **Wu et al.
2021** breast-cancer scRNA-seq atlas (GSE176078), but nothing is breast-specific:
point it at any Visium sample + annotated reference with a cell-type column.

## Stages

| #  | Process              | Tool / script                  | Output                         |
|----|----------------------|--------------------------------|--------------------------------|
| 1  | `LOAD_SPATIAL`       | `bin/load_spatial.py` (scanpy) | `spatial.h5ad`                 |
| 2  | `LOAD_REFERENCE`     | `bin/load_reference.py`        | `reference.h5ad`               |
| 3  | `QC_SPATIAL`         | `bin/qc_spatial.py`            | `spatial_qc.h5ad`, `qc/*.json` |
| 4  | `NORMALIZE`          | `bin/normalize.py`             | normalized h5ads, `hvgs.txt`   |
| 5  | `BUILD_SIGNATURE`    | `bin/build_signature.py`       | `signature.tsv`                |
| 6  | `DECONVOLVE`         | `bin/deconvolve_nnls.py` (scipy NNLS) | `proportions.tsv`       |
| 7  | `SVG`                | `bin/svg_moran.py` (Moran's I) | `svg.tsv`                      |
| 8  | `MAKE_CELLTYPE_PLOT` | `bin/make_celltype_plot.py`    | `spatial_celltypes.html`       |
| 9  | `MAKE_SVG_PLOT`      | `bin/make_svg_plot.py`         | `spatial_svg.html`             |

All outputs land in `results/` via Nextflow `publishDir 'copy'`.

> **Optional v1.1 follow-up — `DECONVOLVE_REF`:** a cell2location / RCTD
> "method bake-off" against the NNLS default, kept in a *separate heavy env* so
> the core pipeline and CI stay light. Not in v1.

## Run

```bash
cd spatial_visium_nf
mamba env create -f envs/spatial_visium_env.yml   # or: conda env create -f ...
conda activate spatial_visium_env

# Offline demo on synthetic Visium + reference with PLANTED proportions (seconds):
bash run_local.sh --demo
#   or with Nextflow (defaults point at the generated demo data):
nextflow run main.nf -profile conda
```

`run_local.sh --demo` synthesizes a toy Visium grid + scRNA reference whose
per-spot cell-type proportions are **planted** and spatially structured, runs the
whole DAG, and **self-checks** that the NNLS deconvolution recovers the planted
proportions within tolerance — the spatial analogue of the planted-DE check in
`plant_rnaseq_nf` / `aml_rnaseq_nf`, and the CI gate.

For a **real run**, fetch the breast-cancer Visium section + Wu et al. atlas once:

```bash
bash fetch_real_data.sh          # downloads + arranges data/real/ (ref is ~0.56 GB)
# the script prints the exact invocation; in brief:
nextflow run main.nf -profile conda \
    --spatial_h5  data/real/spatial/filtered_feature_bc_matrix.h5 \
    --spatial_pos data/real/spatial/spatial/tissue_positions_list.csv \
    --ref_dir     data/real/reference --ref_meta data/real/reference/metadata.csv \
    --celltype_col celltype_major --truth ''
```

`--truth ''` disables the planted self-check (real data has no ground truth).

## Parameters

Override at the Nextflow CLI (`--param value`) or in `nextflow.config`.

| Parameter      | Default (demo)                       | Meaning                                  |
|----------------|--------------------------------------|------------------------------------------|
| `spatial_h5`   | `data/demo/.../filtered_feature_bc_matrix.h5` | Visium filtered matrix          |
| `spatial_pos`  | `data/demo/.../tissue_positions_list.csv` | spot positions                      |
| `ref_dir`      | `data/demo/reference`                | 10x-mtx reference directory              |
| `ref_meta`     | `data/demo/reference/metadata.csv`   | per-cell metadata (carries cell type)    |
| `celltype_col` | `cell_type`                          | metadata column with cell-type labels    |
| `truth`        | `data/demo/truth/proportions.csv`    | planted proportions; `''` disables the self-check |
| `n_hvg`        | `2000`                               | highly variable genes (shared space)     |
| `knn`          | `6`                                  | neighbours in the spatial graph (SVG)    |
| `nperm`        | `100`                                | permutations for Moran's I p-values      |
| `topn`         | `6`                                  | top SVGs shown in the SVG plot           |
| `seed`         | `42`                                 | RNG seed                                 |
| `tol`          | `0.12`                               | max MAE allowed in the deconvolution self-check |

## Methods

- **Load.** `scanpy` reads the 10x Visium filtered matrix + tissue positions
  (in-tissue spots, array coordinates in `obsm['spatial']`) and a CellRanger-style
  reference bundle + metadata carrying the cell-type labels — an identical code
  path on toy and real data.
- **QC.** Per-spot total counts, genes/spot and mitochondrial fraction
  (name-prefix), with threshold filtering and a JSON summary.
- **Normalize.** Library-size `normalize_total` + `log1p`; highly variable genes
  selected on the reference and intersected with the spatial genes, putting both
  matrices in one shared, comparable feature space.
- **Deconvolution (default).** Per spot, **non-negative least squares**
  (`scipy.optimize.nnls`) of its linear-normalized expression onto the
  per-cell-type reference signature matrix, normalized to proportions. Each gene
  is first divided by its mean signature level (inverse-mean weighting) so the fit
  is not dominated by a few very high-expression genes — without it, plain NNLS on
  real cross-platform data collapses onto one or two loud immune signatures. Still
  one deterministic factor per gene: transparent and dependency-light, the on-brand
  analogue of the hand-rolled DE/ORA in the bulk pipelines.
- **Spatially variable genes.** Moran's I on a symmetric kNN spot graph, computed
  from first principles and vectorized across genes, with a seeded permutation
  test. (Squidpy's `spatial_autocorr` is the production equivalent; the hand-rolled
  version keeps the maths legible and CI light.)
- **Visualization.** Two interactive Plotly figures in the naraen.net "Phalaena
  Automata" mauve palette: a spot map (dominant cell type, or any cell type's
  proportion) and the spatial expression of the top spatially variable genes —
  both iframe-ready for the portfolio.

See [`docs/REPORT.md`](docs/REPORT.md) for an end-to-end run.

## Layout

```text
spatial_visium_nf/
├── main.nf                  # Nextflow DSL2 workflow (9 processes)
├── nextflow.config          # manifest, params, profiles
├── fetch_real_data.sh       # Visium breast-cancer + Wu et al. atlas -> data/real/
├── run_local.sh             # mirrors main.nf without Nextflow (--demo offline)
├── bin/
│   ├── make_demo_data.py    # synthetic Visium + reference with PLANTED proportions
│   ├── load_spatial.py
│   ├── load_reference.py
│   ├── qc_spatial.py
│   ├── normalize.py
│   ├── build_signature.py
│   ├── deconvolve_nnls.py   # scipy NNLS + planted-proportions self-check
│   ├── svg_moran.py         # Moran's I on a kNN graph (hand-rolled)
│   ├── make_celltype_plot.py
│   └── make_svg_plot.py
├── assets/NO_TRUTH          # sentinel that disables the self-check on real data
├── envs/spatial_visium_env.yml
├── tests/test_pipeline.py   # fast unit tests (no scanpy stack needed)
├── docs/REPORT.md
├── data/                    # inputs / synthetic demo (gitignored)
└── results/                 # outputs (gitignored)
```

## Reproducibility

- Deterministic (single `seed`); all paths project-relative.
- `nextflow.config` sets `PYTHONNOUSERSITE=1` so processes don't leak user
  site-packages.
- Pinned conda env in `envs/spatial_visium_env.yml`.
- `run_local.sh --demo` runs the entire DAG offline on synthetic data — the CI
  smoke test — and `main.nf` ↔ `run_local.sh` produce **identical** proportions
  (a parity check), keeping the two paths in sync.

## Notes / next steps

- **`DECONVOLVE_REF` (v1.1):** cell2location / RCTD comparison in a separate heavy
  env — the "method bake-off" a standards-focused unit values.
- **Viewer (Part B):** an R Shiny app (SpatialExperiment) over the exported
  outputs — pick a cell type → spatial proportion map; pick a gene → spatial
  expression — deployed to shinyapps.io.
- The reference need not be patient-matched to the Visium section (it is a
  cell-type *reference*); the cross-platform shift this introduces is exactly what
  RCTD / cell2location correct for, motivating the bake-off.
