nextflow.enable.dsl=2

process plink2_gwas_glm {
    tag "plink2 glm: ${id}"
    label 'cpus_8'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/compute", mode: 'copy', pattern: "*.{glm.linear,log}"

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)
    path phenotype
    val trait

    output:
    path "${id}.gwas.glm.linear", emit: glm
    val 'plink2', emit: source

    script:
    """
    set -euo pipefail
    exec > ${id}.plink2_glm.log 2>&1
    python3 -c "from genetics.genomics.plink.results_io import prepare_plink_phenotype_table; prepare_plink_phenotype_table('${phenotype}', '${trait}', '${id}.pheno.tsv')"
    plink2 --pfile ${prefix} \\
        --allow-extra-chr \\
        --pheno ${id}.pheno.tsv \\
        --pheno-name ${trait} \\
        --glm omit-ref no-x-sex hide-covar \\
        --threads ${task.cpus} \\
        --out ${id}.gwas
    """
}

process gcta_gwas {
    tag "gcta gwas: ${id}"
    label 'cpus_8'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/integrated/${params.mod}/compute", mode: 'copy', pattern: "*.{fastGWA,gwa,log}"

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)
    path phenotype
    val trait

    output:
    path "${id}.gcta.fastGWA", emit: glm
    val 'gcta', emit: source

    script:
    """
    set -euo pipefail
    exec > ${id}.gcta_gwas.log 2>&1
    python3 -c "from genetics.genomics.plink.results_io import prepare_plink_phenotype_table; prepare_plink_phenotype_table('${phenotype}', '${trait}', '${id}.pheno.txt')"
    plink2 --pfile ${prefix} --allow-extra-chr --make-bed --threads ${task.cpus} --out ${id}.bed
    ${params.wheat_gcta_bin} \\
        --bfile ${id}.bed \\
        --pheno ${id}.pheno.txt \\
        --mpheno 1 \\
        --out ${id}.gcta \\
        --thread-num ${task.cpus}
    """
}
