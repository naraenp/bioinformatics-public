#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * plant_rnaseq_nf
 * A Nextflow DSL2 pipeline that takes raw Illumina short reads from two plant
 * genotypes — a stress-TOLERANT and a stress-SUSCEPTIBLE line — all the way to
 * differential expression and a phenotype-facing functional summary.
 *
 *   SUBSAMPLE       -> seqtk downsample (reproducible, laptop-friendly)
 *   QC_TRIM         -> fastp adapter/quality trim + FastQC
 *   HISAT2_BUILD    -> spliced-aligner index from the reference genome
 *   ALIGN           -> HISAT2 -> sorted, indexed BAM
 *   QUANTIFY        -> featureCounts gene-level counts
 *   BUILD_MATRIX    -> tidy counts_raw.tsv + metadata.tsv
 *   RUN_DE          -> pydeseq2 tolerant vs. susceptible -> de_results.tsv
 *   ENRICH          -> hypergeometric GO/pathway ORA -> enrichment.tsv
 *   MAKE_HEATMAP    -> interactive heatmap of the top DE genes
 *   MAKE_ENRICH_PLT -> interactive enriched-process bar chart
 *
 * Real run: fetch a rice genome + SRA FASTQ with `bash fetch_real_data.sh`.
 * Offline demo / CI: `bash run_local.sh --demo` builds a toy genome + reads.
 */

params.outdir      = "${projectDir}/results"
params.samplesheet = "${projectDir}/data/samplesheet.csv"
params.genome      = "${projectDir}/data/genome.fa"
params.gtf         = "${projectDir}/data/genes.gtf"
params.gmt         = "${projectDir}/data/go_sets.gmt"

params.subsample   = 0          // read pairs per sample; 0 = use all reads
params.seed        = 42
params.fdr         = 0.05
params.lfc         = 1.0
params.topn        = 40
params.design      = "~genotype"
params.tolerant    = "tolerant"
params.susceptible = "susceptible"


process SUBSAMPLE {
    tag "${sid}"

    input:
    tuple val(sid), val(geno), path(r1), path(r2)

    output:
    tuple val(sid), val(geno), path("${sid}_sub_R1.fastq.gz"), path("${sid}_sub_R2.fastq.gz")

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
    tuple val(sid), val(geno), path(r1), path(r2)

    output:
    tuple val(sid), val(geno), path("${sid}_trim_R1.fastq.gz"), path("${sid}_trim_R2.fastq.gz"), emit: trimmed
    path "${sid}.fastp.{html,json}"

    script:
    """
    fastp -i ${r1} -I ${r2} \\
        -o ${sid}_trim_R1.fastq.gz -O ${sid}_trim_R2.fastq.gz \\
        --json ${sid}.fastp.json --html ${sid}.fastp.html \\
        --thread ${task.cpus}
    """
}

process HISAT2_BUILD {
    tag "${genome.name}"

    input:
    path genome

    output:
    path "idx"

    script:
    """
    mkdir idx
    hisat2-build -p ${task.cpus} ${genome} idx/genome
    """
}

process ALIGN {
    tag "${sid}"

    input:
    tuple val(sid), val(geno), path(r1), path(r2)
    path index

    output:
    path "${sid}.sorted.bam"

    script:
    """
    hisat2 -p ${task.cpus} --no-unal -x ${index}/genome \\
        -1 ${r1} -2 ${r2} \\
        | samtools sort -@ ${task.cpus} -o ${sid}.sorted.bam -
    samtools index ${sid}.sorted.bam
    """
}

process QUANTIFY {
    publishDir "${params.outdir}", mode: 'copy', pattern: "featurecounts.txt.summary"

    input:
    path bams
    path gtf

    output:
    path "featurecounts.txt", emit: counts
    path "featurecounts.txt.summary"

    script:
    """
    featureCounts -p --countReadPairs -T ${task.cpus} \\
        -t exon -g gene_id -a ${gtf} \\
        -o featurecounts.txt ${bams}
    """
}

process BUILD_MATRIX {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path featurecounts
    path samplesheet

    output:
    path "counts_raw.tsv", emit: counts
    path "metadata.tsv",   emit: meta

    script:
    """
    build_count_matrix.py \\
        --featurecounts ${featurecounts} \\
        --samplesheet ${samplesheet} \\
        --counts-out counts_raw.tsv \\
        --meta-out metadata.tsv
    """
}

process RUN_DE {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path counts
    path meta

    output:
    path "de_results.tsv", emit: de
    path "norm_counts.tsv", emit: norm

    script:
    """
    run_de.py \\
        --counts ${counts} --meta ${meta} \\
        --design '${params.design}' \\
        --tolerant ${params.tolerant} --susceptible ${params.susceptible} \\
        --out de_results.tsv --norm-out norm_counts.tsv
    """
}

process ENRICH {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path de
    path gmt

    output:
    path "enrichment.tsv"

    script:
    """
    run_enrichment.py \\
        --de ${de} --gmt ${gmt} \\
        --fdr ${params.fdr} --lfc ${params.lfc} \\
        --out enrichment.tsv
    """
}

process MAKE_HEATMAP {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path norm
    path de
    path meta

    output:
    path "heatmap.html"

    script:
    """
    make_heatmap.py \\
        --norm ${norm} --de ${de} --meta ${meta} \\
        --topn ${params.topn} --fdr ${params.fdr} --lfc ${params.lfc} \\
        --out heatmap.html
    """
}

process MAKE_ENRICH_PLT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path enrichment

    output:
    path "enrichment.html"

    script:
    """
    make_enrichment_plot.py --enrichment ${enrichment} --out enrichment.html
    """
}

workflow {
    sheet = file(params.samplesheet, checkIfExists: true)
    base  = sheet.parent

    reads_ch = Channel.fromPath(sheet)
        .splitCsv(header: true)
        .map { row -> tuple(
            row.sample_id,
            row.genotype,
            file("${base}/${row.fastq_1}", checkIfExists: true),
            file("${base}/${row.fastq_2}", checkIfExists: true),
        ) }

    genome_ch = Channel.fromPath(params.genome, checkIfExists: true)
    gtf_ch    = Channel.fromPath(params.gtf,    checkIfExists: true)
    gmt_ch    = Channel.fromPath(params.gmt,    checkIfExists: true)

    // .first() turns the single index emission into a value channel so it is
    // reused for every sample (otherwise it pairs with just one read set).
    index = HISAT2_BUILD(genome_ch).first()

    SUBSAMPLE(reads_ch)
    QC_TRIM(SUBSAMPLE.out)
    ALIGN(QC_TRIM.out.trimmed, index)

    QUANTIFY(ALIGN.out.collect(), gtf_ch)
    BUILD_MATRIX(QUANTIFY.out.counts, sheet)
    RUN_DE(BUILD_MATRIX.out.counts, BUILD_MATRIX.out.meta)
    ENRICH(RUN_DE.out.de, gmt_ch)
    MAKE_HEATMAP(RUN_DE.out.norm, RUN_DE.out.de, BUILD_MATRIX.out.meta)
    MAKE_ENRICH_PLT(ENRICH.out)
}
