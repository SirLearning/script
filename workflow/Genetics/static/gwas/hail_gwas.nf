nextflow.enable.dsl=2

process HAIL_GWAS {
    tag "hail gwas ${trait}"
    publishDir "${params.outdir}/gwas/hail", mode: 'copy'

    input:
    tuple val(meta), path(vcf)
    path pheno
    path covar
    val trait

    output:
    path "*.gwas.tsv", emit: results
    path "*.log"

    script:
    def prefix = "${meta.id}.${trait}"
    def covar_arg = covar.name != 'NO_FILE' ? "--covar ${covar}" : ""
    def covar_cols = params.covar_names ? "--covar_cols ${params.covar_names.join(',')}" : ""
    def ref = params.reference_genome ?: ""
    def ref_arg = ref ? "--reference ${ref}" : ""
    
    """
    python ${params.src_dir}/python/genetics/hail/gwas.py \
        --vcf ${vcf} \
        --pheno ${pheno} \
        --out ${prefix} \
        --response ${trait} \
        --id_col ${params.pheno_id_col ?: 'IID'} \
        ${covar_arg} \
        ${covar_cols} \
        ${ref_arg} > ${prefix}.log 2>&1
    """
}
