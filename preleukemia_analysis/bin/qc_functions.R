# QC helpers for stage 02, sourced by analysis/02_qc.qmd and bin/run_qc_one.R.

suppressPackageStartupMessages({
  library(DropletUtils)
  library(scater)      # brings scuttle: perCellQCMetrics, isOutlier
  library(scDblFinder)
  library(Matrix)
  library(Seurat)      # Read10X_h5: direct in-memory sparse read
  library(SingleCellExperiment)
})

#' Cell calling, adaptive QC, and doublet removal for one sample.
#'
#' 1. emptyDrops cell calling (FDR <= 0.001); barcodes above the knee are always kept.
#' 2. Per-sample QC: cells beyond 3 MADs (log scale) in library size, detected genes,
#'    or mito fraction are removed, plus a fixed floor of 200 detected genes.
#' 3. scDblFinder doublet removal on the QC-passing cells.
#'
#' Returns the cleaned SingleCellExperiment, a one-row summary, and a thinned
#' barcode-rank curve for plotting.
qc_one_sample <- function(h5_path, gsm, lower = 100, fdr_cutoff = 0.001,
                          nmads = 3, min_features = 200) {
  # Read10X_h5 gives an in-memory sparse matrix; read10xCounts returns a disk-backed
  # DelayedMatrix on .h5 input and is far slower for the steps below
  m <- Seurat::Read10X_h5(h5_path)
  if (is.list(m)) m <- m[["Gene Expression"]]

  # --- 1. cell calling ---------------------------------------------------
  br <- barcodeRanks(m, lower = lower)
  ed <- emptyDrops(m, lower = lower,
                   retain = metadata(br)$knee)      # barcodes above the knee are always cells
  called <- which(ed$FDR <= fdr_cutoff)             # ambient barcodes have NA FDR and drop out
  sce <- SingleCellExperiment(assays = list(counts = m[, called]))
  rowData(sce)$Symbol <- rownames(sce)              # downstream stages expect this column

  # --- 2. adaptive per-sample QC ----------------------------------------
  is_mito <- grepl("^mt-", rownames(sce))
  qc <- perCellQCMetrics(sce, subsets = list(mito = is_mito))
  low_lib   <- isOutlier(qc$sum,      log = TRUE, nmads = nmads, type = "lower")
  low_feat  <- isOutlier(qc$detected, log = TRUE, nmads = nmads, type = "lower")
  high_mito <- isOutlier(qc$subsets_mito_percent,  nmads = nmads, type = "higher")
  under_floor <- qc$detected < min_features
  discard <- low_lib | low_feat | high_mito | under_floor

  colData(sce) <- cbind(colData(sce), qc)
  thresholds <- c(
    lib   = attr(low_lib,   "thresholds")[["lower"]],
    feat  = attr(low_feat,  "thresholds")[["lower"]],
    mito  = attr(high_mito, "thresholds")[["higher"]]
  )
  sce <- sce[, !discard]

  # --- 3. doublets -------------------------------------------------------
  set.seed(1000)                                    # scDblFinder is stochastic
  sce <- scDblFinder(sce, verbose = FALSE)
  n_doublet <- sum(sce$scDblFinder.class == "doublet")
  sce <- sce[, sce$scDblFinder.class == "singlet"]

  summary_row <- data.frame(
    gsm            = gsm,
    n_barcodes     = ncol(m),
    knee           = metadata(br)$knee,
    n_called       = length(called),
    n_low_lib      = sum(low_lib),
    n_low_feat     = sum(low_feat),
    n_high_mito    = sum(high_mito),
    n_under_floor  = sum(under_floor & !(low_lib | low_feat | high_mito)),
    n_discard_qc   = sum(discard),
    thr_lib        = round(thresholds[["lib"]]),
    thr_feat       = round(thresholds[["feat"]]),
    thr_mito       = round(thresholds[["mito"]], 2),
    n_doublet      = n_doublet,
    n_final        = ncol(sce)
  )

  # thin the rank curve for plotting
  keep_pts <- unique(round(10^seq(0, log10(nrow(br)), length.out = 300)))
  rank_curve <- data.frame(gsm = gsm, rank = br$rank[keep_pts], total = br$total[keep_pts])

  list(sce = sce, summary = summary_row, rank_curve = rank_curve)
}
