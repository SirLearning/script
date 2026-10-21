nextflow.enable.dsl=2

/*
 * Gene expression resource: turn the TPM atlases under 00data/12cs_tpm into
 * per-gene summaries keyed by RefSeq v1.1 gene ID, which is what the VEP Gene
 * column carries. Downstream this lets burden and association runs restrict to
 * expressed genes and annotate hits with a tissue profile.
 *
 * Two input modes:
 *   --gene_expr_from_dir <dir>   every *tpm*.tsv[.gz] in the directory; a
 *                                default_metadata.txt beside them is picked up
 *   --gene_expr_tpm <files>      comma-separated matrices
 *
 * Launch: partial_router.nf --partial_task gene_expr --mod test_wap_refchr \
 *           --gene_expr_from_dir /data/home/tusr1/01projects/vmap4/00data/12cs_tpm/ramirez2018_refseq1
 */

include { gene_expression_summary } from '../../../modules/local/genotype/annotate/annotate_expression.nf'

workflow RUN_GENE_EXPR {
    main:
    if (!params.output_dir || !params.job) {
        error 'RUN_GENE_EXPR: params.output_dir and params.job are required.'
    }
    if (!params.user_dir) {
        error 'RUN_GENE_EXPR: params.user_dir is required (Conda stats env).'
    }
    if (!params.mod) {
        error 'RUN_GENE_EXPR: params.mod is required (publish folder).'
    }

    def no_metadata = file("${moduleDir}/../../../assets/NO_METADATA", checkIfExists: true)
    def gff = file(params.gene_expr_gff ?: params.ref_v1_gff, checkIfExists: true)

    def matrices = []
    def default_metadata = null
    if (params.gene_expr_from_dir) {
        def dir = file(params.gene_expr_from_dir, checkIfExists: true)
        matrices = file("${dir}/*tpm*.tsv*").sort()
        if (!matrices) {
            error "RUN_GENE_EXPR: no *tpm*.tsv[.gz] under ${dir}"
        }
        def candidate = file("${dir}/default_metadata.txt")
        if (candidate.exists()) {
            default_metadata = candidate
        }
    }
    else if (params.gene_expr_tpm) {
        matrices = params.gene_expr_tpm.toString().split(',').collect {
            file(it.trim(), checkIfExists: true)
        }
    }
    else {
        error 'RUN_GENE_EXPR: give --gene_expr_from_dir or --gene_expr_tpm.'
    }

    def metadata = params.gene_expr_metadata
        ? file(params.gene_expr_metadata, checkIfExists: true)
        : (default_metadata ?: no_metadata)

    def jobs = matrices.collect { matrix ->
        def prefix = matrix.name.replaceAll(/\.tsv(\.gz)?$/, '').replaceAll(/\./, '_')
        tuple(prefix, matrix, metadata)
    }
    log.info "${params.c_green}gene_expr: ${jobs.size()} matrices, metadata=${metadata.name}${params.c_reset}"

    gene_expression_summary(channel.from(jobs), gff)

    emit:
    table = gene_expression_summary.out.table
}
