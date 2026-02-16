# Package bootstrap script for the R-based preleukemia workflow.
# Run this once before running preleuk_analysis.R.

required_packages <- c(
  "Seurat",
  "Matrix",
  "hdf5r",
  "patchwork",
  "ggplot2",
  "SeuratDisk",
  "shiny",
  "dplyr"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

cat("All required R packages are installed.\n")
