nextflow.enable.dsl=2

/*
 * WAP phenotype matrix aligned to assoc.keep, ready for plink2 --pheno.
 *
 * No mixed-model adjustment here: BLUE/BLUP across replicates and environments
 * is a separate decision and runs before this alignment step.
 */

process wap_phenotype_matrix {
    tag "wap phenotype ${params.mod}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/phenotype", mode: 'copy', pattern: "*.pheno.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.pheno.qc.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_1'

    input:
    path phenotype
    path keep
    val prefix

    output:
    path "*.pheno.tsv", emit: pheno
    path "*.pheno.qc.tsv", emit: qc
    path "*.th.tsv", emit: thresholds
    path "*.png", emit: plots, optional: true
    path "*.log", emit: logs

    script:
    def traits = params.phenotype_wap_traits ? "\"${params.phenotype_wap_traits}\"" : 'None'
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("phenotype_wap.log", "w")
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.phenotype.wap_matrix import run_wap_phenotype_matrix

    run_wap_phenotype_matrix(
        "${phenotype}",
        "${keep}",
        "${prefix}",
        id_col="${params.phenotype_wap_id_col}",
        traits=${traits},
        aggregate_duplicates=${params.phenotype_wap_aggregate_duplicates ? 'True' : 'False'},
    )
    """
}
