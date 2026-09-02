# spatial_visium_nf

A **Nextflow DSL2** pipeline that maps cell types back into tissue space. It takes
a **10x Visium** section plus a **scRNA-seq reference** and produces per-spot
**cell-type deconvolution** (NNLS) and **spatially variable genes** (Moran's I) as
interactive Plotly figures. It is the spatial companion to
[`aml_rnaseq_nf`](../aml_rnaseq_nf).

The default target is **10x Visium Human Breast Cancer** with the **Wu et al. 2021**
atlas (GSE176078) as the reference, but nothing is breast-specific: point it at any
Visium sample and annotated reference.

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

> **Planned:** a `DECONVOLVE_REF` stage comparing NNLS against cell2location and
> RCTD. It is kept out for now because those need a much heavier environment than
> the rest of the pipeline.

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
proportions within tolerance. This is the CI gate: if a change to the method stops
it recovering the planted values, the build fails.

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

### Containers

Besides `-profile conda`, the pipeline ships `docker` and `singularity`
profiles. The `Dockerfile` bakes `envs/spatial_visium_env.yml` (the
scanpy/squidpy/anndata stack) into a
[micromamba](https://hub.docker.com/r/mambaorg/micromamba) image, so the
container matches the conda profile exactly:

```bash
docker build -t spatial_visium_nf:0.1.0 .      # one-time image build
nextflow run main.nf -profile docker           # or: -profile singularity
```

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

- **Load.** `scanpy` reads the Visium matrix + tissue positions and a
  CellRanger-style reference bundle with cell-type labels. Toy and real data take
  the same code path.
- **QC.** Per-spot counts, genes/spot and mitochondrial fraction, with threshold
  filtering and a JSON summary.
- **Normalize.** `normalize_total` + `log1p`; highly variable genes are chosen on
  the reference and intersected with the spatial genes, putting both matrices in a
  shared feature space.
- **Deconvolution.** Non-negative least squares of each spot against the cell-type
  signature matrix, normalized to proportions. Genes are weighted by the inverse of
  their mean signature level, without which a few very high-expression genes
  dominate the fit. See `bin/deconvolve_nnls.py` for why.
- **Spatially variable genes.** Moran's I on a kNN spot graph with a seeded
  permutation test. Squidpy's `spatial_autocorr` is the production equivalent.
- **Visualization.** Two interactive Plotly figures: a spot map of cell-type
  proportions, and the spatial expression of the top spatially variable genes.

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

- **Compare against a probabilistic method.** NNLS gives a point estimate with no
  uncertainty. cell2location and RCTD model the counts directly and return a
  posterior, which is the principled way to handle the next point.
- **The reference is not patient-matched** to this section; it is a cell-type
  reference, so a cross-platform shift remains in the absolute proportions.
- **A Shiny viewer** over the exported outputs: pick a cell type for its spatial
  proportion map, or a gene for its spatial expression.
