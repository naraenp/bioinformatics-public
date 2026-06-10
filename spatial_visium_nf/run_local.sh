#!/usr/bin/env bash
# Run the spatial deconvolution pipeline end-to-end WITHOUT Nextflow — for quick
# iteration and as the CI smoke test. Mirrors the processes in main.nf one-for-one.
#
#   bash run_local.sh --demo     # synthesize a toy Visium + reference with PLANTED
#                                #   cell-type proportions, run, and self-check that
#                                #   NNLS recovers them (the parity/CI gate)
#   bash run_local.sh            # run on whatever is in data/real (fetched with
#                                #   fetch_real_data.sh)
#
# Override inputs with env vars: SPATIAL_H5, SPATIAL_POS, REF_DIR, REF_META,
# CELLTYPE_COL, TRUTH (planted proportions; enables the self-check).
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
BIN="$HERE/bin"
OUT="$HERE/results"
WORK="$OUT/_work"

# Use the env python explicitly when not already on PATH (keeps -profile local-ish
# behaviour); PYTHONNOUSERSITE mirrors the Nextflow env block.
export PYTHONNOUSERSITE=1
PY="${PYTHON:-python}"

NHVG="${NHVG:-2000}"
KNN="${KNN:-6}"
NPERM="${NPERM:-100}"
TOPN="${TOPN:-6}"
SEED="${SEED:-42}"
TOL="${TOL:-0.12}"

DATA="${DATA:-$HERE/data/real}"
if [[ "${1:-}" == "--demo" ]]; then
    DATA="$HERE/data/demo"
    echo ">> generating toy demo data in $DATA"
    "$PY" "$BIN/make_demo_data.py" --outdir "$DATA" --seed "$SEED"
    SPATIAL_H5="$DATA/spatial/filtered_feature_bc_matrix.h5"
    SPATIAL_POS="$DATA/spatial/spatial/tissue_positions_list.csv"
    REF_DIR="$DATA/reference"
    REF_META="$DATA/reference/metadata.csv"
    TRUTH="$DATA/truth/proportions.csv"
fi

SPATIAL_H5="${SPATIAL_H5:-$DATA/spatial/filtered_feature_bc_matrix.h5}"
SPATIAL_POS="${SPATIAL_POS:-$DATA/spatial/spatial/tissue_positions_list.csv}"
REF_DIR="${REF_DIR:-$DATA/reference}"
REF_META="${REF_META:-$DATA/reference/metadata.csv}"
CELLTYPE_COL="${CELLTYPE_COL:-cell_type}"
TRUTH="${TRUTH:-}"

mkdir -p "$OUT" "$WORK" "$OUT/qc"

echo ">> [LOAD_SPATIAL]"
"$PY" "$BIN/load_spatial.py" --h5 "$SPATIAL_H5" --positions "$SPATIAL_POS" \
    --out "$WORK/spatial.h5ad"

echo ">> [LOAD_REFERENCE]"
"$PY" "$BIN/load_reference.py" --dir "$REF_DIR" --metadata "$REF_META" \
    --celltype-col "$CELLTYPE_COL" --out "$WORK/reference.h5ad"

echo ">> [QC_SPATIAL]"
"$PY" "$BIN/qc_spatial.py" --in "$WORK/spatial.h5ad" \
    --out "$WORK/spatial_qc.h5ad" --qc-json "$OUT/qc/spatial_qc.json"

echo ">> [NORMALIZE]"
"$PY" "$BIN/normalize.py" --spatial "$WORK/spatial_qc.h5ad" \
    --reference "$WORK/reference.h5ad" --n-hvg "$NHVG" \
    --out-spatial "$WORK/spatial_norm.h5ad" \
    --out-reference "$WORK/reference_norm.h5ad" --hvg-out "$OUT/hvgs.txt"

echo ">> [BUILD_SIGNATURE]"
# load_reference.py canonicalises the chosen cell-type column onto obs['cell_type'],
# so build_signature reads that fixed name (not the source CELLTYPE_COL).
"$PY" "$BIN/build_signature.py" --reference "$WORK/reference_norm.h5ad" \
    --out "$OUT/signature.tsv"

echo ">> [DECONVOLVE] NNLS"
DECONV_ARGS=(--spatial "$WORK/spatial_norm.h5ad" --signature "$OUT/signature.tsv" \
    --out "$OUT/proportions.tsv")
if [[ -n "$TRUTH" && -f "$TRUTH" ]]; then
    DECONV_ARGS+=(--truth "$TRUTH" --check --tol "$TOL")
fi
"$PY" "$BIN/deconvolve_nnls.py" "${DECONV_ARGS[@]}"

echo ">> [SVG] Moran's I"
"$PY" "$BIN/svg_moran.py" --spatial "$WORK/spatial_norm.h5ad" \
    --k "$KNN" --n-perm "$NPERM" --seed "$SEED" --out "$OUT/svg.tsv"

echo ">> [MAKE_CELLTYPE_PLOT]"
"$PY" "$BIN/make_celltype_plot.py" --proportions "$OUT/proportions.tsv" \
    --out "$OUT/spatial_celltypes.html"

echo ">> [MAKE_SVG_PLOT]"
"$PY" "$BIN/make_svg_plot.py" --svg "$OUT/svg.tsv" \
    --spatial "$WORK/spatial_norm.h5ad" --topn "$TOPN" \
    --out "$OUT/spatial_svg.html"

echo "done -> $OUT/spatial_celltypes.html , $OUT/spatial_svg.html"
