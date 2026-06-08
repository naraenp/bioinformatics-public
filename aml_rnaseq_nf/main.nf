#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * aml_rnaseq_nf
 * A small Nextflow DSL2 pipeline for a bulk RNA-seq differential-expression
 * comparison: AML vs. healthy. It runs on real public RNA-seq cohorts —
 * TCGA-LAML (AML) and GTEx whole blood (healthy) — pulled from recount3,
 * which aligns and quantifies both sources through one uniform pipeline so
 * their gene-level counts are directly comparable.
 *
 *   LOAD_COUNTS       -> join + symbol-map + filter -> counts_raw.tsv, metadata.tsv
 *   NORMALIZE_COUNTS  -> log2(CPM + 1) matrix
 *   RUN_DE            -> per-gene Welch t-test on log2-CPM + BH FDR
 *   MAKE_VOLCANO      -> interactive Plotly volcano (AML markers labeled)
 *
 * Fetch the inputs once with `bash fetch_real_data.sh` (see README).
 */

params.outdir   = "${projectDir}/results"
params.n_per_group = 50
params.seed     = 42
params.fdr      = 0.05
params.lfc      = 1.0

// recount3 + GENCODE inputs, populated by fetch_real_data.sh (gitignored).
params.aml_counts     = "${projectDir}/data/real/tcga.gene_sums.LAML.G026.gz"
params.healthy_counts = "${projectDir}/data/real/gtex.gene_sums.BLOOD.G026.gz"
params.gtf            = "${projectDir}/data/real/gencode.v26.basic.annotation.gtf.gz"

process LOAD_COUNTS {
    tag "TCGA-LAML + GTEx BLOOD, n_per_group=${params.n_per_group}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path aml_counts
    path healthy_counts
    path gtf

    output:
    path "counts_raw.tsv", emit: counts
    path "metadata.tsv",   emit: meta

    script:
    """
    load_counts.py \\
        --aml-counts ${aml_counts} \\
        --healthy-counts ${healthy_counts} \\
        --gtf ${gtf} \\
        --n-per-group ${params.n_per_group} \\
        --seed ${params.seed} \\
        --counts counts_raw.tsv \\
        --meta metadata.tsv
    """
}

process NORMALIZE_COUNTS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path counts

    output:
    path "counts_lcpm.tsv", emit: lcpm

    script:
    """
    normalize_counts.py --counts ${counts} --out counts_lcpm.tsv
    """
}

process RUN_DE {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path lcpm
    path meta

    output:
    path "de_results.tsv", emit: de

    script:
    """
    run_de.py --lcpm ${lcpm} --meta ${meta} --out de_results.tsv
    """
}

process MAKE_VOLCANO {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path de

    output:
    path "volcano.html"

    script:
    """
    make_volcano.py \\
        --de ${de} \\
        --fdr ${params.fdr} \\
        --lfc ${params.lfc} \\
        --out volcano.html
    """
}

workflow {
    aml_ch     = Channel.fromPath(params.aml_counts,     checkIfExists: true)
    healthy_ch = Channel.fromPath(params.healthy_counts, checkIfExists: true)
    gtf_ch     = Channel.fromPath(params.gtf,            checkIfExists: true)

    LOAD_COUNTS(aml_ch, healthy_ch, gtf_ch)
    NORMALIZE_COUNTS(LOAD_COUNTS.out.counts)
    RUN_DE(NORMALIZE_COUNTS.out.lcpm, LOAD_COUNTS.out.meta)
    MAKE_VOLCANO(RUN_DE.out.de)
}
