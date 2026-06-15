#!/usr/bin/env bash
# Run the pipeline end-to-end WITHOUT Nextflow — for quick iteration and as the
# offline validation. Mirrors the processes in main.nf one-for-one.
#
#   bash run_local.sh --demo     # generate a toy genome + reads, run, self-check
#   bash run_local.sh            # run on whatever is in data/ (e.g. real GIAB
#                                #   data fetched with fetch_real_data.sh)
#
# Override input locations with env vars: DATA, SAMPLESHEET, GENOME, GTF, INTERVALS.
# Requires the conda env on PATH: `conda activate variant_calling_env`.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
BIN="$HERE/bin"
OUT="$HERE/results"

SEED="${SEED:-42}"
SUBSAMPLE="${SUBSAMPLE:-0}"
PLOIDY="${PLOIDY:-2}"
TOPN="${TOPN:-40}"
THREADS="${THREADS:-2}"
INTERVALS="${INTERVALS:-}"          # optional BED to restrict calling

DEMO=0
DATA="${DATA:-$HERE/data}"
if [[ "${1:-}" == "--demo" ]]; then
    DEMO=1
    DATA="$HERE/data/demo"
    echo ">> generating toy demo data in $DATA"
    python "$BIN/make_demo_data.py" --outdir "$DATA" --seed "$SEED"
fi

SAMPLESHEET="${SAMPLESHEET:-$DATA/samplesheet.csv}"
GENOME="${GENOME:-$DATA/genome.fa}"
GTF="${GTF:-$DATA/genes.gtf}"

WORK="$OUT/_work"
mkdir -p "$OUT" "$WORK"

# Calling intervals: '-L file.bed' when provided, else empty.
L_ARG=()
[[ -n "$INTERVALS" ]] && L_ARG=(-L "$INTERVALS")

# ---- reference prep (faidx + sequence dictionary + BWA index) ---------------
REF="$WORK/$(basename "$GENOME")"
cp -f "$GENOME" "$REF"
echo ">> [PREP_REFERENCE] faidx + dict + bwa index"
samtools faidx "$REF"
rm -f "${REF%.fa}.dict" "${REF%.*}.dict"
gatk CreateSequenceDictionary -R "$REF" >/dev/null 2>"$WORK/dict.log"
if [[ ! -s "$REF.bwt" ]]; then
    bwa index "$REF" 2>"$WORK/bwaindex.log"
fi

# Sample ids in samplesheet order (used to assemble the GVCF list).
mapfile -t SAMPLES < <(tail -n +2 "$SAMPLESHEET" | cut -d, -f1)

mkdir -p "$OUT/qc"
GVCFS=()
# Pull sample_id + fastq columns by HEADER NAME so extra covariate columns
# (e.g. group) don't shift positions.
awk -F, 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
         {print $h["sample_id"]"\t"$h["fastq_1"]"\t"$h["fastq_2"]}' \
    "$SAMPLESHEET" | while IFS=$'\t' read -r sid fq1 fq2; do
    [[ -z "$sid" ]] && continue
    echo ">> [$sid] subsample -> trim -> align -> markdup -> call"
    bash "$BIN/subsample_reads.sh" "$DATA/$fq1" "$DATA/$fq2" \
        "$WORK/${sid}_sub_R1.fastq.gz" "$WORK/${sid}_sub_R2.fastq.gz" \
        "$SUBSAMPLE" "$SEED"
    fastp -i "$WORK/${sid}_sub_R1.fastq.gz" -I "$WORK/${sid}_sub_R2.fastq.gz" \
        -o "$WORK/${sid}_trim_R1.fastq.gz" -O "$WORK/${sid}_trim_R2.fastq.gz" \
        --json "$OUT/qc/${sid}.fastp.json" --html "$OUT/qc/${sid}.fastp.html" \
        --thread "$THREADS" 2>/dev/null
    bwa mem -t "$THREADS" \
        -R "@RG\tID:${sid}\tSM:${sid}\tPL:ILLUMINA\tLB:${sid}" \
        "$REF" "$WORK/${sid}_trim_R1.fastq.gz" "$WORK/${sid}_trim_R2.fastq.gz" \
        2>"$WORK/${sid}.bwa.log" \
        | samtools sort -@ "$THREADS" -o "$WORK/${sid}.sorted.bam" -
    samtools index "$WORK/${sid}.sorted.bam"
    gatk MarkDuplicates -I "$WORK/${sid}.sorted.bam" \
        -O "$WORK/${sid}.dedup.bam" -M "$WORK/${sid}.dupmetrics.txt" \
        --CREATE_INDEX true >/dev/null 2>"$WORK/${sid}.markdup.log"
    gatk --java-options "-Xmx4g" HaplotypeCaller -R "$REF" \
        -I "$WORK/${sid}.dedup.bam" -O "$WORK/${sid}.g.vcf.gz" \
        -ERC GVCF --sample-ploidy "$PLOIDY" "${L_ARG[@]}" \
        >/dev/null 2>"$WORK/${sid}.hc.log"
done

# (subshell above can't export the array; rebuild the GVCF list here)
for sid in "${SAMPLES[@]}"; do GVCFS+=(-V "$WORK/${sid}.g.vcf.gz"); done

echo ">> [JOINT_GENOTYPE] CombineGVCFs + GenotypeGVCFs"
gatk CombineGVCFs -R "$REF" "${GVCFS[@]}" -O "$WORK/cohort.g.vcf.gz" \
    >/dev/null 2>"$WORK/combine.log"
gatk GenotypeGVCFs -R "$REF" -V "$WORK/cohort.g.vcf.gz" \
    -O "$WORK/cohort.vcf.gz" "${L_ARG[@]}" >/dev/null 2>"$WORK/genotype.log"

echo ">> [FILTER_VARIANTS] GATK hard filters (SNP + indel tracks)"
gatk SelectVariants -R "$REF" -V "$WORK/cohort.vcf.gz" --select-type-to-include SNP \
    -O "$WORK/snps.vcf.gz" >/dev/null 2>"$WORK/selsnp.log"
gatk VariantFiltration -R "$REF" -V "$WORK/snps.vcf.gz" \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "SOR > 3.0" --filter-name "SOR3" \
    -O "$WORK/snps.filt.vcf.gz" >/dev/null 2>"$WORK/filtsnp.log"
gatk SelectVariants -R "$REF" -V "$WORK/cohort.vcf.gz" --select-type-to-include INDEL \
    -O "$WORK/indels.vcf.gz" >/dev/null 2>"$WORK/selindel.log"
gatk VariantFiltration -R "$REF" -V "$WORK/indels.vcf.gz" \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 200.0" --filter-name "FS200" \
    --filter-expression "SOR > 10.0" --filter-name "SOR10" \
    -O "$WORK/indels.filt.vcf.gz" >/dev/null 2>"$WORK/filtindel.log"
gatk MergeVcfs -I "$WORK/snps.filt.vcf.gz" -I "$WORK/indels.filt.vcf.gz" \
    -O "$OUT/cohort.filtered.vcf.gz" >/dev/null 2>"$WORK/merge.log"

echo ">> [NORMALIZE] bcftools norm (split multiallelic, left-align)"
bcftools norm -f "$REF" -m- -O v -o "$OUT/cohort.norm.vcf" \
    "$OUT/cohort.filtered.vcf.gz" 2>"$WORK/norm.log"

echo ">> [ANNOTATE] genic consequence + Ts/Tv"
python "$BIN/annotate_variants.py" \
    --vcf "$OUT/cohort.norm.vcf" --gtf "$GTF" --out "$OUT/annotated.tsv"

echo ">> [MAKE_OVERVIEW]"
python "$BIN/make_overview.py" --annotated "$OUT/annotated.tsv" \
    --out "$OUT/variant_overview.html"

echo ">> [MAKE_GENO_HEATMAP]"
python "$BIN/make_geno_heatmap.py" --annotated "$OUT/annotated.tsv" \
    --meta "$SAMPLESHEET" --topn "$TOPN" --out "$OUT/genotype_heatmap.html"

if [[ "$DEMO" -eq 1 ]]; then
    echo ">> [SELF-CHECK] recover planted variants from $DATA/truth.vcf"
    bcftools norm -f "$REF" -m- -O v -o "$WORK/truth.norm.vcf" "$DATA/truth.vcf" \
        2>"$WORK/truthnorm.log"
    python "$BIN/check_truth.py" \
        --called "$OUT/cohort.norm.vcf" --truth "$WORK/truth.norm.vcf"
fi

echo "done -> $OUT/variant_overview.html , $OUT/genotype_heatmap.html"
