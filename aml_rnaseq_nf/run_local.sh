#!/usr/bin/env bash
# Run the pipeline end-to-end without Nextflow (handy for quick local
# iteration / CI smoke tests). Mirrors the four processes in main.nf.
#
# Requires the inputs in data/real/ — fetch them once with:
#   bash fetch_real_data.sh
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
BIN="$HERE/bin"
OUT="$HERE/results"
REAL="$HERE/data/real"
mkdir -p "$OUT"

N_PER_GROUP=50
SEED=42
FDR=0.05
LFC=1.0

python "$BIN/load_counts.py" \
    --aml-counts     "$REAL/tcga.gene_sums.LAML.G026.gz" \
    --healthy-counts "$REAL/gtex.gene_sums.BLOOD.G026.gz" \
    --gtf            "$REAL/gencode.v26.basic.annotation.gtf.gz" \
    --n-per-group "$N_PER_GROUP" --seed "$SEED" \
    --counts "$OUT/counts_raw.tsv" \
    --meta   "$OUT/metadata.tsv"

python "$BIN/normalize_counts.py" \
    --counts "$OUT/counts_raw.tsv" \
    --out    "$OUT/counts_lcpm.tsv"

python "$BIN/run_de.py" \
    --lcpm "$OUT/counts_lcpm.tsv" \
    --meta "$OUT/metadata.tsv" \
    --out  "$OUT/de_results.tsv"

python "$BIN/make_volcano.py" \
    --de "$OUT/de_results.tsv" \
    --fdr "$FDR" --lfc "$LFC" \
    --out "$OUT/volcano.html"

echo "done -> $OUT/volcano.html"
