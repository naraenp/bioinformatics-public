#!/usr/bin/env bash
# One-time download of the pipeline inputs into data/real/ (gitignored, ~130 MB):
#   - TCGA-LAML and GTEx whole-blood gene_sums from recount3
#     (Monorail-aligned, GENCODE v26 gene-level counts).
#   - GENCODE v26 basic annotation GTF, for Ensembl -> HGNC symbol mapping.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
DEST="$HERE/data/real"
mkdir -p "$DEST"

base_recount=http://duffel.rail.bio/recount3/human/data_sources
base_gencode=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26

declare -a urls=(
  "$base_recount/tcga/gene_sums/ML/LAML/tcga.gene_sums.LAML.G026.gz"
  "$base_recount/gtex/gene_sums/OD/BLOOD/gtex.gene_sums.BLOOD.G026.gz"
  "$base_gencode/gencode.v26.basic.annotation.gtf.gz"
)

for url in "${urls[@]}"; do
    f="$DEST/$(basename "$url")"
    if [[ -s "$f" ]]; then
        echo "skip (exists): $f"
    else
        echo "fetching $url"
        curl -sSL -o "$f" "$url"
    fi
done

echo "done -> $DEST"
ls -lh "$DEST"
