nextflow.enable.dsl=2

/*
 * Gene expression atlases (00data/12cs_tpm) normalized onto RefSeq v1.1 gene
 * IDs, so the VEP Gene column joins to them without string surgery. The v1.1
 * GFF supplies the v1.0 -> v1.1 map through its previous_id attribute.
 */

process gene_expression_summary {
    tag "gene expression ${prefix}"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.gene_expr.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_1'

    input:
    tuple val(prefix), path(tpm_table), path(metadata)
    path gff

    output:
    tuple val(prefix), path("*.gene_expr.tsv"), emit: table
    path "*.th.tsv", emit: thresholds
    path "*.png", emit: plots, optional: true
    path "*.log", emit: logs

    script:
    def metadata_arg = metadata.name == 'NO_METADATA' ? 'None' : "\"${metadata}\""
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("${prefix}.gene_expr.log", "w")
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.gene.expression import run_gene_expression

    run_gene_expression(
        "${tpm_table}",
        "${gff}",
        "${prefix}",
        metadata_path=${metadata_arg},
        sample_col="${params.gene_expr_sample_col}",
        group_col="${params.gene_expr_group_col}",
        detect_min=float("${params.gene_expr_detect_min}"),
    )
    """
}
