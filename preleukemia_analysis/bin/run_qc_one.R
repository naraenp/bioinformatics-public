#!/usr/bin/env Rscript
# Runs stage-02 QC for one sample in its own R process.
# Usage: Rscript bin/run_qc_one.R <h5_path> <gsm> <out_dir>
# Separate processes are used instead of forking inside the qmd because HDF5 and
# BLAS are not fork-safe; parallelize with xargs -P at the shell level.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
h5_path <- args[1]; gsm <- args[2]; out_dir <- args[3]

out  <- file.path(out_dir, paste0("qc_", gsm, ".rds"))
meta <- file.path(out_dir, paste0("meta_", gsm, ".rds"))
if (file.exists(out) && file.exists(meta)) {
  message(gsm, ": already done, skipping"); quit(status = 0)
}

# locate qc_functions.R next to this script
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
source(file.path(script_dir, "qc_functions.R"))

set.seed(42)
message(gsm, ": starting")
t0 <- Sys.time()
res <- qc_one_sample(h5_path, gsm)
saveRDS(res$sce, out)
saveRDS(list(summary = res$summary, rank_curve = res$rank_curve), meta)
message(sprintf("%s: done in %.1f min (%d cells)", gsm,
                as.numeric(difftime(Sys.time(), t0, units = "mins")),
                res$summary$n_final))
