#!/usr/bin/env bash
# Run the pipeline end-to-end WITHOUT Nextflow — for quick iteration and as the
# CI smoke test. Mirrors the processes in main.nf one-for-one.
#
#   bash run_local.sh --demo     # generate a toy genome + reads, then run
#   bash run_local.sh            # run on whatever is in data/ (e.g. real rice
#                                #   data fetched with fetch_real_data.sh)
#
# Override the input locations with env vars: DATA, SAMPLESHEET, GENOME, GTF, GMT.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
BIN="$HERE/bin"
OUT="$HERE/results"

SEED="${SEED:-42}"
SUBSAMPLE="${SUBSAMPLE:-0}"
FDR="${FDR:-0.05}"
LFC="${LFC:-1.0}"
TOPN="${TOPN:-40}"
DESIGN="${DESIGN:-~genotype}"
TOLERANT="${TOLERANT:-tolerant}"
SUSCEPTIBLE="${SUSCEPTIBLE:-susceptible}"
THREADS="${THREADS:-2}"

DATA="${DATA:-$HERE/data}"
if [[ "${1:-}" == "--demo" ]]; then
    DATA="$HERE/data/demo"
    echo ">> generating toy demo data in $DATA"
    python "$BIN/make_demo_data.py" --outdir "$DATA" --seed "$SEED"
fi

SAMPLESHEET="${SAMPLESHEET:-$DATA/samplesheet.csv}"
GENOME="${GENOME:-$DATA/genome.fa}"
GTF="${GTF:-$DATA/genes.gtf}"
GMT="${GMT:-$DATA/go_sets.gmt}"

WORK="$OUT/_work"
mkdir -p "$OUT" "$WORK"

mkdir -p "$WORK/idx"
if [[ -s "$WORK/idx/genome.1.ht2" ]]; then
    echo ">> [HISAT2_BUILD] reusing existing index"
else
    echo ">> [HISAT2_BUILD] indexing $(basename "$GENOME")"
    hisat2-build -p "$THREADS" "$GENOME" "$WORK/idx/genome" >/dev/null
fi

# Sample ids in samplesheet order (used to assemble the featureCounts BAM list).
mapfile -t SAMPLES < <(tail -n +2 "$SAMPLESHEET" | cut -d, -f1)

mkdir -p "$OUT/qc"
# Pull sample_id + fastq columns by HEADER NAME, so extra covariate columns
# (e.g. condition) don't shift positions. Emits "sid<TAB>fq1<TAB>fq2" per row.
awk -F, 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
         {print $h["sample_id"]"\t"$h["fastq_1"]"\t"$h["fastq_2"]}' \
    "$SAMPLESHEET" | while IFS=$'\t' read -r sid fq1 fq2; do
    [[ -z "$sid" ]] && continue
    echo ">> [$sid] subsample -> trim -> align"
    bash "$BIN/subsample_reads.sh" "$DATA/$fq1" "$DATA/$fq2" \
        "$WORK/${sid}_sub_R1.fastq.gz" "$WORK/${sid}_sub_R2.fastq.gz" \
        "$SUBSAMPLE" "$SEED"
    fastp -i "$WORK/${sid}_sub_R1.fastq.gz" -I "$WORK/${sid}_sub_R2.fastq.gz" \
        -o "$WORK/${sid}_trim_R1.fastq.gz" -O "$WORK/${sid}_trim_R2.fastq.gz" \
        --json "$OUT/qc/${sid}.fastp.json" --html "$OUT/qc/${sid}.fastp.html" \
        --thread "$THREADS" 2>/dev/null
    hisat2 -p "$THREADS" --no-unal -x "$WORK/idx/genome" \
        -1 "$WORK/${sid}_trim_R1.fastq.gz" -2 "$WORK/${sid}_trim_R2.fastq.gz" \
        2>"$WORK/${sid}.hisat2.log" \
        | samtools sort -@ "$THREADS" -o "$WORK/${sid}.sorted.bam" -
    samtools index "$WORK/${sid}.sorted.bam"
done

BAM_LIST=$(printf "%s\n" "${SAMPLES[@]}" | sed "s#^#$WORK/#; s#\$#.sorted.bam#" | tr '\n' ' ')

echo ">> [QUANTIFY] featureCounts"
featureCounts -p --countReadPairs -T "$THREADS" -t exon -g gene_id \
    -a "$GTF" -o "$WORK/featurecounts.txt" $BAM_LIST 2>"$WORK/featurecounts.log"
cp "$WORK/featurecounts.txt.summary" "$OUT/featurecounts.txt.summary"

echo ">> [BUILD_MATRIX]"
python "$BIN/build_count_matrix.py" \
    --featurecounts "$WORK/featurecounts.txt" --samplesheet "$SAMPLESHEET" \
    --counts-out "$OUT/counts_raw.tsv" --meta-out "$OUT/metadata.tsv"

echo ">> [RUN_DE] pydeseq2"
python "$BIN/run_de.py" \
    --counts "$OUT/counts_raw.tsv" --meta "$OUT/metadata.tsv" \
    --design "$DESIGN" --tolerant "$TOLERANT" --susceptible "$SUSCEPTIBLE" \
    --out "$OUT/de_results.tsv" --norm-out "$OUT/norm_counts.tsv"

echo ">> [ENRICH]"
python "$BIN/run_enrichment.py" \
    --de "$OUT/de_results.tsv" --gmt "$GMT" --fdr "$FDR" --lfc "$LFC" \
    --out "$OUT/enrichment.tsv"

echo ">> [MAKE_HEATMAP]"
python "$BIN/make_heatmap.py" \
    --norm "$OUT/norm_counts.tsv" --de "$OUT/de_results.tsv" \
    --meta "$OUT/metadata.tsv" --topn "$TOPN" --fdr "$FDR" --lfc "$LFC" \
    --out "$OUT/heatmap.html"

echo ">> [MAKE_ENRICH_PLT]"
python "$BIN/make_enrichment_plot.py" \
    --enrichment "$OUT/enrichment.tsv" --out "$OUT/enrichment.html"

echo "done -> $OUT/heatmap.html , $OUT/enrichment.html"
