#!/usr/bin/env bash
# One-time download of the REAL inputs for the rice drought tolerant-vs-
# susceptible comparison, into data/real/ (gitignored). Idempotent.
#
# Study: BioProject PRJNA338445 (Wilkins et al.) — paired-end Illumina RNA-seq
# of the drought-TOLERANT cultivar Apo and the drought-SUSCEPTIBLE cultivar
# IR64, each under control and drought stress. We contrast genotype while
# controlling for condition (design ~condition + genotype), so control + stress
# act as the two replicates per genotype.
#
# Reads are STREAM-SUBSAMPLED directly from ENA: only the first LITE_PAIRS read
# pairs of each run are pulled (curl stops once head closes the pipe), so the
# whole download is a few hundred MB rather than ~18 GB. Set LITE_PAIRS=0 for the
# full runs.
#
# Requires the conda env (hisat2): `conda activate plant_rnaseq_env`.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
BIN="$HERE/bin"
REAL="$HERE/data/real"
mkdir -p "$REAL/reads"

LITE_PAIRS="${LITE_PAIRS:-3000000}"
ENSEMBL_RELEASE="${ENSEMBL_RELEASE:-58}"
BASE="http://ftp.ensemblgenomes.org/pub/plants/release-${ENSEMBL_RELEASE}"
GENOME_URL="$BASE/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz"
GTF_URL="$BASE/gtf/oryza_sativa/Oryza_sativa.IRGSP-1.0.${ENSEMBL_RELEASE}.gtf.gz"

# run_accession  genotype  condition  sample_id
RUNS=$(cat <<'EOF'
SRR4017524	tolerant	control	apo_ctrl
SRR4017527	tolerant	stress	apo_stress
SRR4017522	susceptible	control	ir64_ctrl
SRR4017525	susceptible	stress	ir64_stress
EOF
)

fetch() {  # url dest
    if [[ -s "$2" ]]; then echo "skip (exists): $(basename "$2")"; return; fi
    echo "downloading $(basename "$2")"
    curl -fsSL "$1" -o "$2"
}

# Stream the first N read pairs from an ENA gzipped FASTQ without downloading it
# whole. pipefail is relaxed locally because curl is SIGPIPE'd once head is done.
stream_subsample() {  # url dest pairs
    local url="$1" dest="$2" pairs="$3"
    # >100 bytes guards against an empty-gzip stub from a previous failed stream.
    [[ -f "$dest" && "$(stat -c%s "$dest")" -gt 100 ]] && \
        { echo "skip (exists): $(basename "$dest")"; return; }
    echo "streaming $((pairs)) pairs -> $(basename "$dest")"
    ( set +o pipefail
      curl -fsS "$url" | zcat | head -n $((pairs * 4)) | gzip > "$dest" )
}

ena_fastq_urls() {  # run_accession -> prints "url1<TAB>url2"
    local run="$1"
    # The API always prepends a run_accession column, so fastq_ftp is field 2.
    curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=fastq_ftp&format=tsv" \
        | tail -n +2 | cut -f2 | tr ';' '\t' | sed 's#ftp.sra.ebi.ac.uk#http://ftp.sra.ebi.ac.uk#g'
}

# ---- reference genome + annotation -----------------------------------------
fetch "$GENOME_URL" "$REAL/genome.fa.gz"
fetch "$GTF_URL"    "$REAL/genes.gtf.gz"
[[ -s "$REAL/genome.fa" ]] || gunzip -k "$REAL/genome.fa.gz"
[[ -s "$REAL/genes.gtf" ]] || gunzip -k "$REAL/genes.gtf.gz"

# ---- reads + samplesheet ----------------------------------------------------
SHEET="$REAL/samplesheet.csv"
echo "sample_id,genotype,condition,fastq_1,fastq_2" > "$SHEET"
while IFS=$'\t' read -r run geno cond sid; do
    [[ -z "${run:-}" ]] && continue
    read -r u1 u2 < <(ena_fastq_urls "$run")
    r1="reads/${sid}_R1.fastq.gz"; r2="reads/${sid}_R2.fastq.gz"
    if [[ "$LITE_PAIRS" -gt 0 ]]; then
        stream_subsample "$u1" "$REAL/$r1" "$LITE_PAIRS"
        stream_subsample "$u2" "$REAL/$r2" "$LITE_PAIRS"
    else
        fetch "$u1" "$REAL/$r1"; fetch "$u2" "$REAL/$r2"
    fi
    echo "${sid},${geno},${cond},${r1},${r2}" >> "$SHEET"
done <<< "$RUNS"

# ---- GO gene sets (BioMart, IDs match the GTF) ------------------------------
if [[ ! -s "$REAL/rice_go.gmt" ]]; then
    python "$BIN/build_go_gmt.py" --out "$REAL/rice_go.gmt"
fi

# ---- HISAT2 index -----------------------------------------------------------
if [[ ! -s "$REAL/idx/genome.1.ht2" ]]; then
    echo "building HISAT2 index (slow step for a full genome)"
    mkdir -p "$REAL/idx"
    hisat2-build -p 4 "$REAL/genome.fa" "$REAL/idx/genome"
fi

cat <<EOF

done. Inputs in $REAL
Run the pipeline (design controls for drought condition):

  nextflow run main.nf -profile conda \\
      --samplesheet $SHEET --genome $REAL/genome.fa --gtf $REAL/genes.gtf \\
      --gmt $REAL/rice_go.gmt --design '~condition + genotype'

  # or without Nextflow:
  DATA=$REAL GMT=$REAL/rice_go.gmt DESIGN='~condition + genotype' \\
      TOLERANT=tolerant SUSCEPTIBLE=susceptible bash run_local.sh
EOF
