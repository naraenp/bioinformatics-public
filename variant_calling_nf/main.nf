#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * variant_calling_nf
 * A Nextflow DSL2 pipeline that takes raw Illumina short reads through the GATK
 * germline short-variant best-practice path to a filtered, annotated cohort VCF
 * and two interactive figures.
 *
 *   SUBSAMPLE        -> seqtk downsample (reproducible, laptop-friendly)
 *   QC_TRIM          -> fastp adapter/quality trim
 *   PREP_REFERENCE   -> samtools faidx + GATK sequence dictionary + BWA index
 *   ALIGN            -> BWA-MEM (+ read groups) -> sorted, indexed BAM
 *   MARK_DUPLICATES  -> GATK MarkDuplicates
 *   CALL_VARIANTS    -> GATK HaplotypeCaller in GVCF mode (per sample)
 *   JOINT_GENOTYPE   -> CombineGVCFs + GenotypeGVCFs -> cohort VCF
 *   FILTER_VARIANTS  -> GATK hard filters (separate SNP / indel tracks)
 *   NORMALIZE        -> bcftools norm (split multiallelic, left-align)
 *   ANNOTATE         -> genic consequence (exonic/intronic/intergenic) + Ts/Tv
 *   MAKE_OVERVIEW    -> interactive per-sample variant overview
 *   MAKE_GENO_HEATMAP-> interactive genotype heatmap of discriminating variants
 *
 * Real run: fetch a GRCh38 slice + GIAB HG002 reads with `bash fetch_real_data.sh`.
 * Offline demo / CI: `bash run_local.sh --demo` builds a toy genome + reads with
 * planted variants and self-checks recovery against the truth set.
 */

params.outdir      = "${projectDir}/results"
params.samplesheet = "${projectDir}/data/samplesheet.csv"
params.genome      = "${projectDir}/data/genome.fa"
params.gtf         = "${projectDir}/data/genes.gtf"
params.intervals   = ""          // optional BED to restrict calling (real run)
params.subsample   = 0           // read pairs per sample; 0 = use all reads
params.seed        = 42
params.ploidy      = 2
params.topn        = 40

// Optional interval restriction, threaded into the GATK steps when set.
ivl = params.intervals ? "-L ${params.intervals}" : ""


process SUBSAMPLE {
    tag "${sid}"

    input:
    tuple val(sid), val(grp), path(r1), path(r2)

    output:
    tuple val(sid), val(grp), path("${sid}_sub_R1.fastq.gz"), path("${sid}_sub_R2.fastq.gz")

    script:
    """
    subsample_reads.sh ${r1} ${r2} ${sid}_sub_R1.fastq.gz ${sid}_sub_R2.fastq.gz \\
        ${params.subsample} ${params.seed}
    """
}

process QC_TRIM {
    tag "${sid}"
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: "*.{html,json}"

    input:
    tuple val(sid), val(grp), path(r1), path(r2)

    output:
    tuple val(sid), val(grp), path("${sid}_trim_R1.fastq.gz"), path("${sid}_trim_R2.fastq.gz"), emit: trimmed
    path "${sid}.fastp.{html,json}"

    script:
    """
    fastp -i ${r1} -I ${r2} \\
        -o ${sid}_trim_R1.fastq.gz -O ${sid}_trim_R2.fastq.gz \\
        --json ${sid}.fastp.json --html ${sid}.fastp.html \\
        --thread ${task.cpus}
    """
}

process PREP_REFERENCE {
    tag "${genome.name}"

    input:
    path genome

    output:
    path "ref"

    script:
    """
    mkdir ref
    cp ${genome} ref/genome.fa
    samtools faidx ref/genome.fa
    gatk CreateSequenceDictionary -R ref/genome.fa
    bwa index ref/genome.fa
    """
}

process ALIGN {
    tag "${sid}"

    input:
    tuple val(sid), val(grp), path(r1), path(r2)
    path ref

    output:
    tuple val(sid), path("${sid}.sorted.bam"), path("${sid}.sorted.bam.bai")

    script:
    """
    bwa mem -t ${task.cpus} \\
        -R "@RG\\tID:${sid}\\tSM:${sid}\\tPL:ILLUMINA\\tLB:${sid}" \\
        ${ref}/genome.fa ${r1} ${r2} \\
        | samtools sort -@ ${task.cpus} -o ${sid}.sorted.bam -
    samtools index ${sid}.sorted.bam
    """
}

process MARK_DUPLICATES {
    tag "${sid}"

    input:
    tuple val(sid), path(bam), path(bai)

    output:
    tuple val(sid), path("${sid}.dedup.bam"), path("${sid}.dedup.bai")

    script:
    """
    gatk MarkDuplicates -I ${bam} -O ${sid}.dedup.bam \\
        -M ${sid}.dupmetrics.txt --CREATE_INDEX true
    """
}

process CALL_VARIANTS {
    tag "${sid}"

    input:
    tuple val(sid), path(bam), path(bai)
    path ref

    output:
    tuple path("${sid}.g.vcf.gz"), path("${sid}.g.vcf.gz.tbi")

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" HaplotypeCaller \\
        -R ${ref}/genome.fa -I ${bam} -O ${sid}.g.vcf.gz \\
        -ERC GVCF --sample-ploidy ${params.ploidy} ${ivl}
    """
}

process JOINT_GENOTYPE {
    publishDir "${params.outdir}", mode: 'copy', pattern: "cohort.vcf.gz*"

    input:
    path gvcfs
    path tbis
    path ref

    output:
    tuple path("cohort.vcf.gz"), path("cohort.vcf.gz.tbi")

    script:
    """
    args=""
    for v in *.g.vcf.gz; do args="\$args -V \$v"; done
    gatk CombineGVCFs -R ${ref}/genome.fa \$args -O combined.g.vcf.gz
    gatk GenotypeGVCFs -R ${ref}/genome.fa -V combined.g.vcf.gz \\
        -O cohort.vcf.gz ${ivl}
    """
}

process FILTER_VARIANTS {
    publishDir "${params.outdir}", mode: 'copy', pattern: "cohort.filtered.vcf.gz*"

    input:
    tuple path(vcf), path(tbi)
    path ref

    output:
    path "cohort.filtered.vcf.gz", emit: vcf
    path "cohort.filtered.vcf.gz.tbi"

    script:
    """
    gatk SelectVariants -R ${ref}/genome.fa -V ${vcf} --select-type-to-include SNP -O snps.vcf.gz
    gatk VariantFiltration -R ${ref}/genome.fa -V snps.vcf.gz \\
        --filter-expression "QD < 2.0"   --filter-name QD2 \\
        --filter-expression "FS > 60.0"  --filter-name FS60 \\
        --filter-expression "MQ < 40.0"  --filter-name MQ40 \\
        --filter-expression "SOR > 3.0"  --filter-name SOR3 \\
        -O snps.filt.vcf.gz
    gatk SelectVariants -R ${ref}/genome.fa -V ${vcf} --select-type-to-include INDEL -O indels.vcf.gz
    gatk VariantFiltration -R ${ref}/genome.fa -V indels.vcf.gz \\
        --filter-expression "QD < 2.0"   --filter-name QD2 \\
        --filter-expression "FS > 200.0" --filter-name FS200 \\
        --filter-expression "SOR > 10.0" --filter-name SOR10 \\
        -O indels.filt.vcf.gz
    gatk MergeVcfs -I snps.filt.vcf.gz -I indels.filt.vcf.gz -O cohort.filtered.vcf.gz
    """
}

process NORMALIZE {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path vcf
    path ref

    output:
    path "cohort.norm.vcf"

    script:
    """
    bcftools norm -f ${ref}/genome.fa -m- -O v -o cohort.norm.vcf ${vcf}
    """
}

process ANNOTATE {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path norm
    path gtf

    output:
    path "annotated.tsv"

    script:
    """
    annotate_variants.py --vcf ${norm} --gtf ${gtf} --out annotated.tsv
    """
}

process MAKE_OVERVIEW {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path annotated

    output:
    path "variant_overview.html"

    script:
    """
    make_overview.py --annotated ${annotated} --out variant_overview.html
    """
}

process MAKE_GENO_HEATMAP {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path annotated
    path samplesheet

    output:
    path "genotype_heatmap.html"

    script:
    """
    make_geno_heatmap.py --annotated ${annotated} --meta ${samplesheet} \\
        --topn ${params.topn} --out genotype_heatmap.html
    """
}

workflow {
    sheet = file(params.samplesheet, checkIfExists: true)
    base  = sheet.parent

    reads_ch = Channel.fromPath(sheet)
        .splitCsv(header: true)
        .map { row -> tuple(
            row.sample_id,
            row.group,
            file("${base}/${row.fastq_1}", checkIfExists: true),
            file("${base}/${row.fastq_2}", checkIfExists: true),
        ) }

    genome_ch = Channel.fromPath(params.genome, checkIfExists: true)
    gtf_ch    = Channel.fromPath(params.gtf,    checkIfExists: true)

    // .first() turns the single reference bundle into a value channel so it is
    // reused for every sample (otherwise it pairs with just one read set).
    ref = PREP_REFERENCE(genome_ch).first()

    SUBSAMPLE(reads_ch)
    QC_TRIM(SUBSAMPLE.out)
    ALIGN(QC_TRIM.out.trimmed, ref)
    MARK_DUPLICATES(ALIGN.out)
    CALL_VARIANTS(MARK_DUPLICATES.out, ref)

    // Collect per-sample GVCFs (+ their indexes) for joint genotyping.
    gvcfs = CALL_VARIANTS.out.map { it[0] }.collect()
    tbis  = CALL_VARIANTS.out.map { it[1] }.collect()

    JOINT_GENOTYPE(gvcfs, tbis, ref)
    FILTER_VARIANTS(JOINT_GENOTYPE.out, ref)
    NORMALIZE(FILTER_VARIANTS.out.vcf, ref)
    ANNOTATE(NORMALIZE.out, gtf_ch)
    MAKE_OVERVIEW(ANNOTATE.out)
    MAKE_GENO_HEATMAP(ANNOTATE.out, sheet)
}
