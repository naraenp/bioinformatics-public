# Collects small dashboard inputs from data/processed/ into dashboard/data/.
# The app never loads the large objects; the biggest artifact is a 60k-cell UMAP sample.
# Run after stages 02-07: Rscript bin/prepare_dashboard_data.R

suppressPackageStartupMessages({
  library(here); library(tidyverse); library(Seurat)
})
set.seed(42)

dash_data <- here("dashboard", "data")
dir.create(dash_data, recursive = TRUE, showWarnings = FALSE)

# --- QC ---------------------------------------------------------------------
qc_dir <- here("data", "processed", "02_qc")
qc_summary <- read_csv(file.path(qc_dir, "qc_summary.csv"), show_col_types = FALSE)
saveRDS(qc_summary, file.path(dash_data, "qc_summary.rds"))

rank_curves <- list.files(qc_dir, pattern = "^meta_", full.names = TRUE) |>
  map(\(f) readRDS(f)$rank_curve) |> bind_rows() |>
  left_join(select(qc_summary, gsm, model), by = "gsm")
saveRDS(rank_curves, file.path(dash_data, "rank_curves.rds"))

# --- UMAP sample (annotated cells + trajectory overlay) ---------------------
seu <- readRDS(here("data", "processed", "04_annotation", "seurat_annotated.rds"))
umap_df <- as_tibble(Embeddings(seu, "umap"), rownames = "cell") |>
  rename(UMAP1 = 2, UMAP2 = 3) |>   # positional: col 1 is 'cell', 2-3 are the embedding
  bind_cols(seu@meta.data |> select(gsm, model, condition, celltype, seurat_clusters)) |>
  mutate(celltype = trimws(celltype))
rm(seu); invisible(gc())

traj <- read_csv(here("data", "processed", "06_trajectory", "trajectory_per_cell.csv"),
                 show_col_types = FALSE)
umap_df$dpt_pseudotime <- traj$dpt_pseudotime[match(umap_df$cell, traj[[1]])]

n_keep <- min(60000, nrow(umap_df))
umap_df <- umap_df |> slice_sample(n = n_keep)
saveRDS(umap_df, file.path(dash_data, "umap_sample.rds"))

# --- Composition ------------------------------------------------------------
ct_counts <- read_csv(here("data", "processed", "04_annotation", "celltype_counts.csv"),
                      show_col_types = FALSE) |>
  mutate(celltype = trimws(celltype)) |>          # atlas label "Middle " has a trailing space
  group_by(gsm) |> mutate(prop = n / sum(n)) |> ungroup()
saveRDS(ct_counts, file.path(dash_data, "celltype_props.rds"))
saveRDS(read_csv(here("data", "processed", "05_differential", "composition_propeller.csv"),
                 show_col_types = FALSE) |> mutate(celltype = trimws(celltype)),
        file.path(dash_data, "propeller.rds"))

# --- Differential expression ------------------------------------------------
saveRDS(read_csv(here("data", "processed", "05_differential", "pseudobulk_de_all.csv"),
                 show_col_types = FALSE),
        file.path(dash_data, "de_results.rds"))

# --- Survival ---------------------------------------------------------------
saveRDS(read_csv(here("data", "processed", "07_survival", "km_data.csv"),
                 show_col_types = FALSE),
        file.path(dash_data, "km_data.rds"))
saveRDS(read_csv(here("data", "processed", "07_survival", "cox_results.csv"),
                 show_col_types = FALSE),
        file.path(dash_data, "cox_results.rds"))

cat("dashboard/data ready:\n")
print(file.info(list.files(dash_data, full.names = TRUE))["size"] / 1e6)
