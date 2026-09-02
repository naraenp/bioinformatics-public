# Preleukemia scRNA-seq analysis (revised)

Reanalysis of 38 mouse bone-marrow HSPC samples across eight preleukemic mutation
models (Calr, Dnmt3a, Ezh2, Flt3-ITD, Idh1, Jak2, Npm1c, Utx) from
[Isobe et al., *Cell Genomics* 2023](https://doi.org/10.1016/j.xgen.2023.100426)
(GEO GSE227026), with a trajectory analysis, a TCGA-LAML survival arm, and a Shiny
dashboard.

This is a statistical revision of the original pipeline. Each stage is a Quarto
document under `analysis/`; the rendered HTML files are the analysis record and
explain each method choice in place.

## Statistical revisions vs. the original pipeline

| Decision | Original | Revised | Why |
| --- | --- | --- | --- |
| Cell calling | implicit (`min.features = 200`) | `emptyDrops` (FDR ≤ 0.001, knee retained) | tests barcodes against the ambient profile; recovers small cells, controls FDR |
| QC thresholds | one fixed cutoff for all 38 samples | per-sample 3-MAD adaptive (log scale) + 200-gene floor | libraries differ; fixed cutoffs conflate quality with biology |
| Doublets | none | `scDblFinder` per sample | undetected doublets fake "intermediate" states downstream |
| Integration | CCA anchors, arbitrary reference sample | Harmony on 50 PCs, sample as batch | scalable, no arbitrary reference, condition not given as nuisance |
| Annotation | anchor label transfer, no confidence | SingleR vs. Dahlin 2018 atlas, pruned scores + marker verification | integration-independent, calibrated confidence |
| Composition | chi-square on pooled cells | propeller: logit proportions, limma with `~ model + condition` | the mouse, not the cell, is the replicate (pseudoreplication) |
| DE | cell-level tests | pseudobulk + edgeR QL, `~ model + condition` | same replication logic; calibrated FDR (Squair et al. 2021) |
| Trajectory | ad-hoc notebook | DPT + CellRank (pseudotime kernel); sample-level fate summaries | no spliced counts ⇒ no velocity (stated limitation) |
| Survival | median split + log-rank | age-adjusted Cox PH on continuous score; KM for display; PH tested | dichotomization discards information; age is the dominant confounder |

## Headline results

| Stage | Result |
| --- | --- |
| 02 QC | 276,294 barcodes called as cells; 230,684 remain after adaptive QC and doublet removal (25,302 doublets). Per-sample gene-count thresholds ranged from 573 to 4,146 and mito thresholds from 3.3% to 4.7%. |
| 03 Integration | Harmony over 38 samples; per-cluster sample-mixing entropy 0.95 to 0.97 (1 = evenly mixed). |
| 04 Annotation | SingleR against the Dahlin 2018 atlas: 0.24% of cells low-confidence; label distribution consistent with an LK sort (4.6% HSCs). |
| 05 Composition | Chi-square on pooled cells: p < 10^-15. propeller (sample-level, model-adjusted): no cell type at FDR 0.05. All logit fold changes are negative, which means mutant samples are more heterogeneous than WT (each model expands a different lineage). |
| 05 DE | Pseudobulk edgeR: no shared mutant-versus-WT genes in HSCs and at most 13 in any cell type at FDR 0.05. A program shared across eight distinct mutations is not detectable at n = 38; per-model effects cannot be tested with 2 to 3 mice per arm. |
| 06 Trajectory | DPT from an *Hlf*-high HSC root orders lineages correctly (HSCs 0.02, late erythroid 0.74). CellRank fate probabilities toward five mature lineages, reported conditional on commitment. |
| 07 Survival | TCGA-LAML (n = 151): neither PLPS nor Stem11 is associated with overall survival after age adjustment (HR per SD 0.97 and 0.97; p = 0.75 and 0.81). The unadjusted result is also null. The proportional-hazards assumption holds. |

Once the animal or patient is the unit of inference and covariates are included,
most of the original pipeline's significant findings do not hold. That is the
correct result, not a disappointing one.

## Layout

```
analysis/    Quarto stages 01 to 07; rendered HTML files are the analysis record
bin/         helper functions and the per-sample QC worker
envs/        pinned conda envs: preleuk_r.yml, preleuk_py.yml
data/        raw/ and processed/ (gitignored; regenerated from stage 01 onward)
dashboard/   Shiny app; reads the small artifacts in dashboard/data/
```

## Reproduce

```bash
mamba env create -f envs/preleuk_r.yml
mamba env create -f envs/preleuk_py.yml
conda activate preleuk_r

# stages 01-05 (R); stage 02's per-sample QC is run in parallel first
quarto render analysis/01_data_acquisition.qmd
awk -F, 'NR>1{print $NF, $1}' data/raw/sample_metadata.csv |
  xargs -P 4 -n 2 sh -c 'Rscript bin/run_qc_one.R "$1" "$2" data/processed/02_qc' w
for s in 02_qc 03_integration 04_annotation 05_differential; do
  quarto render analysis/$s.qmd; done

# stages 06-07 (Python) + dashboard
QUARTO_PYTHON=$CONDA_PREFIX/../preleuk_py/bin/python \
  quarto render analysis/06_trajectory.qmd analysis/07_survival.qmd
Rscript bin/prepare_dashboard_data.R
R -e 'shiny::runApp("dashboard")'
```

Inputs: GEO h5 matrices are downloaded in stage 01 with an md5 manifest. The
annotation reference is the Dahlin et al. 2018 mouse HSPC atlas. The PLPS and
Stem11 signatures and the TCGA-LAML tables (cBioPortal) are under `data/external/`.

## Reproducibility notes

- Paths are anchored to the project root with `here` and the `.here` file; there
  are no user-specific absolute paths.
- All random steps are seeded. Parallel QC runs in separate processes, not forks.
- The env files are the only source of dependencies; versions are pinned.
- See `envs/INSTRUCTIONS.md` for environment setup, including a fix for the conda
  quarto package.
