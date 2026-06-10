# Run report — spatial_visium_nf

End-to-end run of the spatial deconvolution pipeline: the offline synthetic
validation (the CI gate) and the real breast-cancer target it is built for.

## Datasets

| Input | Source | Notes |
|-------|--------|-------|
| Spatial | 10x Genomics **Visium Human Breast Cancer (Block A, Section 1)**, Space Ranger 1.1.0 (CC BY 4.0) | filtered feature-barcode matrix + `spatial/` tissue positions |
| Reference | **Wu et al. 2021**, *A single-cell and spatially resolved atlas of human breast cancers*, Nat Genet — GEO **GSE176078** | scRNA-seq atlas; `metadata.csv` cell-type labels (`celltype_major`) |
| Demo | synthetic Visium grid + scRNA reference with **planted proportions** | no download; `bin/make_demo_data.py` |

All real URLs are verified and hardcoded in `fetch_real_data.sh`. The reference is
a breast-cancer cell-type *reference*, not patient-matched to the Visium section —
the realistic deconvolution setting (and the cross-platform case RCTD /
cell2location are designed for; that comparison is the v1.1 `DECONVOLVE_REF`
follow-up).

## Pipeline

`LOAD_SPATIAL → LOAD_REFERENCE → QC_SPATIAL → NORMALIZE → BUILD_SIGNATURE →
DECONVOLVE (NNLS) → SVG (Moran's I) → MAKE_CELLTYPE_PLOT / MAKE_SVG_PLOT.`

Deconvolution models each spot's library-size-normalized (linear) expression as a
non-negative mixture of the reference cell-type signatures and solves
`min ||S w − x||₂, w ≥ 0` with `scipy.optimize.nnls`, normalizing `w` to
proportions. Each gene (in both the signature and the spot) is first divided by
its mean signature level, i.e. the fit minimizes `Σ_g (1/mean_g²)(S_g w − x_g)²`.
Without this inverse-mean weighting a plain L2 fit in linear space is dominated by
a handful of very high-expression genes (immunoglobulins and the like), and on
real cross-platform data the deconvolution collapses onto one or two loud immune
signatures; the weighting is one deterministic factor per gene, so it stays fully
transparent, and it leaves the matched synthetic demo unchanged within tolerance.
SVGs are scored by Moran's I on a symmetric kNN spot graph with a seeded
permutation test.

## Validation (offline synthetic, the CI gate)

`bash run_local.sh --demo` synthesizes a reference of **5 cell types** (Tcell,
Bcell, Myeloid, Epithelial, Stroma; 300 cells each) and a **14×14 = 196-spot**
Visium grid whose per-spot proportions are planted with spatial structure (each
type dominates a region), then runs the whole DAG. With `seed=42`:

- **Deconvolution self-check** — recovered vs. planted proportions:
  **MAE = 0.0074, correlation = 0.999** (per-type MAE 0.007–0.009); well within
  the `tol = 0.12` gate. Dominant-type spots: Stroma 58, Bcell 46, Myeloid 45,
  Epithelial 26, Tcell 21.
- **Spatially variable genes** — because each cell type is spatially clustered,
  its marker genes are spatially autocorrelated: the **top SVGs are all planted
  markers** (e.g. `MK_Bcell_*` at Moran's I ≈ 0.87–0.89) while housekeeping genes
  sit near I ≈ 0.
- **Parity** — the Nextflow and `run_local.sh` paths produce **identical**
  `proportions.tsv`.
- **Unit tests** — `pytest spatial_visium_nf/tests/` → **7/7 pass** (NNLS
  recovery, kNN-graph + Moran's I, signature means, demo generators).

The figures are interactive Plotly HTML: `results/spatial_celltypes.html`
(dominant type, or any cell type's proportion map) and `results/spatial_svg.html`
(top-SVG spatial expression) — iframe-ready for the portfolio entry. Static PNG
thumbnails (`make_*_plot.py --png`) are best-effort via kaleido and are produced
in the portfolio-integration step; the interactive HTML is the canonical output.

## Results (real breast-cancer run)

`fetch_real_data.sh` pulls the **10x Visium Human Breast Cancer (Block A,
Section 1)** section (3,798 in-tissue spots, 36,601 genes) and the **Wu et al.
2021** atlas (100,064 cells, 9 `celltype_major` classes), and the pipeline runs
end-to-end on `data/real` with the self-check disabled (`--truth ''`). Spatial and
reference share **20,309 genes**; HVG selection on the reference reduces this to
the 2,000-gene signature space.

**Deconvolution.** The recovered architecture is the expected breast-tumor
composition. Mean per-spot proportions:

| Cell type | Mean proportion | Spots where dominant |
|-----------|----------------:|---------------------:|
| Cancer Epithelial | 0.38 | 2,201 |
| Normal Epithelial | 0.21 | 806 |
| CAFs | 0.11 | 284 |
| Myeloid | 0.09 | 168 |
| Endothelial | 0.08 | 132 |
| PVL (perivascular-like) | 0.04 | 46 |
| B-cells | 0.04 | 79 |
| T-cells | 0.04 | 46 |
| Plasmablasts | 0.02 | 36 |

Cancer epithelium dominates the section (58% of spots) with stromal (CAF, PVL,
endothelial) and immune (myeloid, lymphoid, plasmablast) compartments mapping to
coherent regions rather than scattering, which is the signal a spot map should
show. The interactive `spatial_celltypes.html` lets you switch between the
dominant-type view and any single cell type's proportion map.

**Spatially variable genes.** Of 1,632 expressed, sufficiently variable genes,
**949 are significant** (Moran's I permutation `p < 0.05`). The top SVGs read as a
tour of the tissue's compartments: **MUC1** (epithelial / breast-cancer mucin),
**MGP** and **IGFBP5** (stromal / matrix), **CD74** (antigen presentation),
**IGHG3** / **IGKC** (immunoglobulins, marking plasma-cell / lymphoid aggregates),
and **FTH1** (ferritin, myeloid). `spatial_svg.html` shows the top six in space.

| Rank | Gene | Moran's I | Compartment |
|-----:|------|----------:|-------------|
| 1 | MGP | 0.78 | stroma |
| 2 | IGFBP5 | 0.69 | stroma |
| 3 | MUC1 | 0.64 | epithelium / tumor |
| 4 | CD74 | 0.62 | immune (MHC-II) |
| 5 | MALAT1 | 0.60 | ubiquitous lncRNA |
| 6 | FTH1 | 0.60 | myeloid |
| 7 | IGHG3 | 0.60 | plasma cells |
| 8 | IGKC | 0.59 | plasma cells |

**Caveat (and the motivation for v1.1).** The Wu et al. reference is not patient-
matched to this section; it is a cross-platform cell-type reference. NNLS is the
transparent baseline and the inverse-mean gene weighting is what keeps it from
collapsing onto the loudest immune genes, but absolute proportions still carry the
platform shift. Quantifying that against a probabilistic method built for it
(`DECONVOLVE_REF`: cell2location / RCTD) is the v1.1 follow-up.

## What this project is / is not

- **Is:** a reproducible, readable Visium deconvolution + SVG workflow —
  tissue-in to cell-types-in-space — that runs offline on a laptop and self-checks
  against a known ground truth.
- **Is not:** a deconvolution benchmark. NNLS is the transparent default; the
  cell2location / RCTD comparison (`DECONVOLVE_REF`) is a flagged v1.1 stage in a
  separate heavy env. Moran's I is hand-rolled for legibility; squidpy's
  `spatial_autocorr` is the drop-in production equivalent.
