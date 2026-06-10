#!/usr/bin/env bash
# One-time download of the REAL inputs into data/real/ (gitignored). Idempotent.
#
# Spatial:   10x Visium Human Breast Cancer (Block A, Section 1), Space Ranger
#            1.1.0 demo dataset (CC BY 4.0).
# Reference: Wu et al. 2021, "A single-cell and spatially resolved atlas of human
#            breast cancers" (Nat Genet), scRNA-seq atlas — GEO GSE176078. Its
#            metadata.csv carries the cell-type labels (celltype_major).
#
# The reference is not patient-matched to the Visium section — it is a
# breast-cancer cell-type reference, which is the realistic deconvolution setting
# (and exactly the cross-platform case RCTD / cell2location are built for; that
# comparison is the v1.1 DECONVOLVE_REF follow-up).
#
# Only curl/tar/gzip/awk are needed (no conda env). The reference tarball is
# ~0.56 GB — the one heavy download.
#
# All URLs below were verified (HTTP 206 + content-length) before being hardcoded.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REAL="$HERE/data/real"
mkdir -p "$REAL"

VISIUM_BASE="https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_1"
H5_URL="$VISIUM_BASE/V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5"
SPATIAL_URL="$VISIUM_BASE/V1_Breast_Cancer_Block_A_Section_1_spatial.tar.gz"
REF_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"

CELLTYPE_COL="${CELLTYPE_COL:-celltype_major}"

fetch() {  # url dest
    if [[ -s "$2" ]]; then echo "skip (exists): $(basename "$2")"; return; fi
    echo "downloading $(basename "$2")"
    curl -fsSL "$1" -o "$2"
}

# ---- spatial: filtered matrix + tissue positions ---------------------------
mkdir -p "$REAL/spatial"
fetch "$H5_URL" "$REAL/spatial/filtered_feature_bc_matrix.h5"
if [[ ! -s "$REAL/spatial/spatial/tissue_positions_list.csv" ]]; then
    echo "downloading + extracting spatial.tar.gz"
    curl -fsSL "$SPATIAL_URL" -o "$REAL/spatial/spatial.tar.gz"
    tar -xzf "$REAL/spatial/spatial.tar.gz" -C "$REAL/spatial"   # -> spatial/spatial/
    rm -f "$REAL/spatial/spatial.tar.gz"
fi

# ---- reference: Wu et al. atlas -> standard 10x mtx names -------------------
REF="$REAL/reference"
if [[ ! -s "$REF/matrix.mtx.gz" ]]; then
    echo "downloading + extracting reference (~0.56 GB, one-time)"
    mkdir -p "$REF"
    curl -fsSL "$REF_URL" -o "$REAL/ref.tar.gz"
    tar -xzf "$REAL/ref.tar.gz" -C "$REAL"
    SRC="$REAL/Wu_etal_2021_BRCA_scRNASeq"
    # CellRanger-style names so load_reference.py (sc.read_10x_mtx) reads it.
    gzip -c "$SRC/count_matrix_sparse.mtx"    > "$REF/matrix.mtx.gz"
    gzip -c "$SRC/count_matrix_barcodes.tsv"  > "$REF/barcodes.tsv.gz"
    # genes.tsv is single-column symbols; build a v3 features.tsv (id, name, type).
    awk 'BEGIN{FS=OFS="\t"} {print $1, $1, "Gene Expression"}' \
        "$SRC/count_matrix_genes.tsv" | gzip -c > "$REF/features.tsv.gz"
    cp "$SRC/metadata.csv" "$REF/metadata.csv"
    rm -rf "$SRC" "$REAL/ref.tar.gz"
fi

# Confirm the chosen cell-type column exists in the reference metadata.
if ! head -1 "$REF/metadata.csv" | tr ',' '\n' | grep -qx "$CELLTYPE_COL"; then
    echo "WARNING: '$CELLTYPE_COL' not found in $REF/metadata.csv header." >&2
    echo "         Available columns:" >&2
    head -1 "$REF/metadata.csv" | tr ',' '\n' | sed 's/^/           /' >&2
fi

cat <<EOF

done. Inputs in $REAL
Run the pipeline on the real breast-cancer data (no planted truth, so the
self-check is disabled with --truth ''):

  nextflow run main.nf -profile conda \\
      --spatial_h5  $REAL/spatial/filtered_feature_bc_matrix.h5 \\
      --spatial_pos $REAL/spatial/spatial/tissue_positions_list.csv \\
      --ref_dir     $REF \\
      --ref_meta    $REF/metadata.csv \\
      --celltype_col $CELLTYPE_COL \\
      --truth ''

  # or without Nextflow:
  DATA=$REAL CELLTYPE_COL=$CELLTYPE_COL TRUTH='' bash run_local.sh
EOF
