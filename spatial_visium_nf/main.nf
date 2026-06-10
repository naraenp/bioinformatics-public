#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * spatial_visium_nf
 * A Nextflow DSL2 pipeline that takes a 10x Visium spatial sample plus a matched
 * scRNA-seq reference and maps cell types back into tissue space: per-spot
 * deconvolution (NNLS against a reference signature) and spatially variable genes.
 *
 *   LOAD_SPATIAL        -> read Visium filtered matrix + spot positions -> AnnData
 *   LOAD_REFERENCE      -> read annotated scRNA-seq reference (cell-type labels)
 *   QC_SPATIAL          -> per-spot QC (counts / genes / mito %) + filter
 *   NORMALIZE           -> library-size normalise + log1p + shared HVG space
 *   BUILD_SIGNATURE     -> per-cell-type reference signature matrix
 *   DECONVOLVE          -> NNLS per spot -> cell-type proportions (self-checks on demo)
 *   SVG                 -> spatially variable genes (Moran's I, kNN graph)
 *   MAKE_CELLTYPE_PLOT  -> interactive spatial map of the deconvolution
 *   MAKE_SVG_PLOT       -> interactive spatial expression of the top SVGs
 *
 * Real run: fetch a 10x Visium breast-cancer section + the Wu et al. 2021 atlas
 *           with `bash fetch_real_data.sh`, then override the params it prints.
 * Offline demo / CI: `bash run_local.sh --demo` synthesizes a toy Visium +
 *           reference with PLANTED proportions and checks NNLS recovers them.
 *           `main.nf` defaults point at that demo data for an out-of-box run.
 */

params.outdir       = "${projectDir}/results"
params.spatial_h5   = "${projectDir}/data/demo/spatial/filtered_feature_bc_matrix.h5"
params.spatial_pos  = "${projectDir}/data/demo/spatial/spatial/tissue_positions_list.csv"
params.ref_dir      = "${projectDir}/data/demo/reference"
params.ref_meta     = "${projectDir}/data/demo/reference/metadata.csv"
params.celltype_col = "cell_type"
params.truth        = "${projectDir}/data/demo/truth/proportions.csv"  // "" to disable

params.n_hvg = 2000
params.knn   = 6
params.nperm = 100
params.topn  = 6
params.seed  = 42
params.tol   = 0.12


process LOAD_SPATIAL {
    input:
    path h5
    path positions

    output:
    path "spatial.h5ad"

    script:
    """
    load_spatial.py --h5 ${h5} --positions ${positions} --out spatial.h5ad
    """
}

process LOAD_REFERENCE {
    input:
    path ref_dir
    path ref_meta

    output:
    path "reference.h5ad"

    script:
    """
    load_reference.py --dir ${ref_dir} --metadata ${ref_meta} \\
        --celltype-col ${params.celltype_col} --out reference.h5ad
    """
}

process QC_SPATIAL {
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: "*.json"

    input:
    path spatial

    output:
    path "spatial_qc.h5ad", emit: qc
    path "spatial_qc.json"

    script:
    """
    qc_spatial.py --in ${spatial} --out spatial_qc.h5ad --qc-json spatial_qc.json
    """
}

process NORMALIZE {
    publishDir "${params.outdir}", mode: 'copy', pattern: "hvgs.txt"

    input:
    path spatial_qc
    path reference

    output:
    path "spatial_norm.h5ad",   emit: spatial
    path "reference_norm.h5ad", emit: reference
    path "hvgs.txt"

    script:
    """
    normalize.py --spatial ${spatial_qc} --reference ${reference} \\
        --n-hvg ${params.n_hvg} \\
        --out-spatial spatial_norm.h5ad --out-reference reference_norm.h5ad \\
        --hvg-out hvgs.txt
    """
}

process BUILD_SIGNATURE {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path reference_norm

    output:
    path "signature.tsv"

    script:
    """
    build_signature.py --reference ${reference_norm} --out signature.tsv
    """
}

process DECONVOLVE {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path spatial_norm
    path signature
    path truth

    output:
    path "proportions.tsv"

    script:
    def check = truth.name != 'NO_TRUTH' ? "--truth ${truth} --check --tol ${params.tol}" : ""
    """
    deconvolve_nnls.py --spatial ${spatial_norm} --signature ${signature} \\
        --out proportions.tsv ${check}
    """
}

process SVG {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path spatial_norm

    output:
    path "svg.tsv"

    script:
    """
    svg_moran.py --spatial ${spatial_norm} --k ${params.knn} \\
        --n-perm ${params.nperm} --seed ${params.seed} --out svg.tsv
    """
}

process MAKE_CELLTYPE_PLOT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path proportions

    output:
    path "spatial_celltypes.html"

    script:
    """
    make_celltype_plot.py --proportions ${proportions} --out spatial_celltypes.html
    """
}

process MAKE_SVG_PLOT {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path svg
    path spatial_norm

    output:
    path "spatial_svg.html"

    script:
    """
    make_svg_plot.py --svg ${svg} --spatial ${spatial_norm} \\
        --topn ${params.topn} --out spatial_svg.html
    """
}

workflow {
    spatial_h5  = Channel.fromPath(params.spatial_h5,  checkIfExists: true)
    spatial_pos = Channel.fromPath(params.spatial_pos, checkIfExists: true)
    ref_dir     = Channel.fromPath(params.ref_dir,     checkIfExists: true, type: 'dir')
    ref_meta    = Channel.fromPath(params.ref_meta,    checkIfExists: true)

    // Optional planted-truth self-check: a sentinel disables it on real data.
    truth = params.truth ? file(params.truth, checkIfExists: true)
                         : file("${projectDir}/assets/NO_TRUTH")

    spatial = LOAD_SPATIAL(spatial_h5, spatial_pos)
    reference = LOAD_REFERENCE(ref_dir, ref_meta)

    QC_SPATIAL(spatial)
    NORMALIZE(QC_SPATIAL.out.qc, reference)

    BUILD_SIGNATURE(NORMALIZE.out.reference)
    DECONVOLVE(NORMALIZE.out.spatial, BUILD_SIGNATURE.out, truth)
    SVG(NORMALIZE.out.spatial)

    MAKE_CELLTYPE_PLOT(DECONVOLVE.out)
    MAKE_SVG_PLOT(SVG.out, NORMALIZE.out.spatial)
}
