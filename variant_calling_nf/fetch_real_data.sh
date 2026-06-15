#!/usr/bin/env bash
# One-time download of the REAL inputs for a germline benchmarking run on the
# Genome in a Bottle (GIAB) sample HG002, into data/real/ (gitignored). Idempotent.
#
# Sample: GIAB / NIST HG002 (NA24385, Ashkenazi son) — the field-standard
# germline benchmark, shipping a high-confidence truth VCF + confident-region
# BED. We call variants in a small high-confidence slice of chr20 and benchmark
# the callset against the NIST v4.2.1 truth.
#
# To stay laptop-scale we do NOT download whole genomes:
#   * the reference is just GRCh38 chr20 (UCSC, ~20 MB gzipped);
#   * reads are RANGE-SLICED out of GIAB's chr20-only 300x BAM over HTTPS — only
#     the bytes for REGION are fetched (samtools uses the remote .bai), then
#     converted back to paired FASTQ (analogous to plant_rnaseq_nf's ENA stream).
#
# Requires the conda env (samtools/bcftools/bwa): `conda activate variant_calling_env`.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
REAL="$HERE/data/real"
mkdir -p "$REAL/reads"

# A ~2 Mb high-confidence window on chr20. Override REGION / SAMPLE_ID via env.
REGION="${REGION:-chr20:30000000-32000000}"
SAMPLE_ID="${SAMPLE_ID:-HG002}"
CHROM="${REGION%%:*}"; RANGE="${REGION#*:}"
START="${RANGE%-*}"; END="${RANGE#*-}"
START0=$((START - 1))   # BED is 0-based, half-open

GIAB="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab"
TRUTH_DIR="$GIAB/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38"
TRUTH_VCF_URL="$TRUTH_DIR/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED_URL="$TRUTH_DIR/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
BAM_URL="$GIAB/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x_chr20.bam"
REF_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gtf.gz"

fetch() {  # url dest
    if [[ -s "$2" ]]; then echo "skip (exists): $(basename "$2")"; return; fi
    echo "downloading $(basename "$2")"
    curl -fsSL "$1" -o "$2"
}

# ---- reference: GRCh38 chr20 (contig name "chr20" matches the BAM + truth) ---
fetch "$REF_URL" "$REAL/chr20.fa.gz"
[[ -s "$REAL/genome.fa" ]] || gunzip -kc "$REAL/chr20.fa.gz" > "$REAL/genome.fa"
samtools faidx "$REAL/genome.fa"

# ---- reads: range-slice REGION out of the remote chr20 BAM, back to FASTQ ----
R1="$REAL/reads/${SAMPLE_ID}_R1.fastq.gz"; R2="$REAL/reads/${SAMPLE_ID}_R2.fastq.gz"
if [[ ! -s "$R1" ]]; then
    echo "slicing $REGION reads from the remote GIAB 300x BAM (range requests)"
    samtools view -b "$BAM_URL" "$REGION" \
        | samtools collate -Ou - \
        | samtools fastq -1 "$R1" -2 "$R2" -0 /dev/null -s /dev/null -n 2>/dev/null
fi

# ---- annotation: GENCODE GTF, filtered to the calling chromosome -----------
fetch "$GTF_URL" "$REAL/gencode.gtf.gz"
if [[ ! -s "$REAL/genes.gtf" ]]; then
    zcat "$REAL/gencode.gtf.gz" | awk -v c="$CHROM" '$1==c' > "$REAL/genes.gtf"
fi

SHEET="$REAL/samplesheet.csv"
echo "sample_id,group,fastq_1,fastq_2" > "$SHEET"
echo "${SAMPLE_ID},${SAMPLE_ID},reads/${SAMPLE_ID}_R1.fastq.gz,reads/${SAMPLE_ID}_R2.fastq.gz" >> "$SHEET"

# ---- truth VCF + confident BED, restricted to REGION ------------------------
fetch "$TRUTH_VCF_URL"     "$REAL/truth.full.vcf.gz"
fetch "$TRUTH_VCF_URL.tbi" "$REAL/truth.full.vcf.gz.tbi"
fetch "$TRUTH_BED_URL"     "$REAL/confident.full.bed"

# confident regions intersected with REGION -> the pipeline's calling intervals
awk -v c="$CHROM" -v s="$START0" -v e="$END" 'BEGIN{OFS="\t"}
    $1==c && $3>s && $2<e {b=($2>s?$2:s); f=($3<e?$3:e); if(f>b) print c,b,f}' \
    "$REAL/confident.full.bed" > "$REAL/intervals.bed"
bcftools view -R "$REAL/intervals.bed" -O z -o "$REAL/truth.region.vcf.gz" \
    "$REAL/truth.full.vcf.gz"
bcftools index -t "$REAL/truth.region.vcf.gz"

cat <<EOF

done. Inputs in $REAL
Run the pipeline restricted to the high-confidence calling intervals:

  nextflow run main.nf -profile conda \\
      --samplesheet $SHEET --genome $REAL/genome.fa --gtf $REAL/genes.gtf \\
      --intervals $REAL/intervals.bed

  # or without Nextflow:
  DATA=$REAL INTERVALS=$REAL/intervals.bed bash run_local.sh

Benchmark the callset against the NIST truth (site-level precision / recall):

  bcftools norm -f $REAL/genome.fa -m- -O z -o results/called.norm.vcf.gz \\
      results/cohort.filtered.vcf.gz && bcftools index -t results/called.norm.vcf.gz
  bcftools isec -p results/bench -R $REAL/intervals.bed \\
      $REAL/truth.region.vcf.gz results/called.norm.vcf.gz
  #   results/bench/0000.vcf = false negatives (truth only)
  #   results/bench/0001.vcf = false positives (calls only)
  #   results/bench/0002.vcf = true positives  -> recall, precision
EOF
