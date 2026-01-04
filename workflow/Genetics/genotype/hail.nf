nextflow.enable.dsl=2

process vcf_to_mt_hail {
    tag { meta.id ? "hail vcf2mt ${meta.id}" : 'hail vcf2mt' }
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy'

    input:
    tuple val(meta), val(prefix), path(vcf)

    output:
    tuple val(meta), path("${prefix}.mt"), emit: mt

    script:
    def ref = params.reference_genome ?: ""
    def ref_arg = ref ? "--reference ${ref}" : ""
    """
    python ${params.src_dir}/python/genetics/hail/vcf_to_mt.py \\
        --vcf ${vcf} \\
        --out ${prefix}.mt \\
        ${ref_arg}
    """
}

process filter_hail {
    tag { meta.id ? "hail filter ${meta.id}" : 'hail filter' }
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy'

    input:
    tuple val(meta), path(vcf), val(job_config)

    output:
    tuple val(meta), val("${meta.id}.hail.filtered"), path("${meta.id}.hail.filtered.vcf.gz"), emit: vcf

    script:
    def id = meta.id
    def ref = params.reference_genome ?: ""
    def ref_arg = ref ? "--reference ${ref}" : ""
    """
    python ${params.src_dir}/python/genetics/hail/filter.py \
        --vcf ${vcf} \
        --out ${id}.hail.filtered \
        --maf ${maf} \
        --mac ${mac} \
        --min_alleles ${min_alleles} \
        --max_alleles ${max_alleles} \
        --max_missing ${max_missing} \
        ${ref_arg}
    
    mv ${id}.hail.filtered.vcf.bgz ${id}.hail.filtered.vcf.gz
    """
}


process HAIL_QC {
    tag "hail qc ${meta.id}"
    container "${params.hail_container}"
    publishDir "${params.outdir}/genotype/stats/hail", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.variant_qc.tsv"), path("*.sample_qc.tsv"), emit: qc_stats
    path "*.log"

    script:
    def prefix = meta.id
    def ref = params.reference_genome ?: ""
    def ref_arg = ref ? "--reference ${ref}" : ""
    """
    python ${params.src_dir}/python/genetics/hail/qc.py \
        --vcf ${vcf} \
        --out ${prefix} \
        ${ref_arg} > ${prefix}.log 2>&1
    """
}

process HAIL_KINSHIP {
    tag "hail kinship ${meta.id}"
    publishDir "${params.outdir}/genotype/kinship/hail", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.kinship.tsv"), emit: kinship
    path "*.log"

    script:
    def prefix = meta.id
    def ref = params.reference_genome ?: ""
    def ref_arg = ref ? "--reference ${ref}" : ""
    def k = params.pca_k ?: 10
    """
    python ${params.src_dir}/python/genetics/hail/kinship.py \
        --vcf ${vcf} \
        --out ${prefix} \
        --method pc_relate \
        --k ${k} \
        ${ref_arg} > ${prefix}.log 2>&1
    """
}

process HAIL_PCA {
    tag "hail pca ${meta.id}"
    publishDir "${params.outdir}/genotype/population_structure/hail", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.scores.tsv"), emit: scores
    tuple val(meta), path("*.loadings.tsv"), emit: loadings
    tuple val(meta), path("*.eigenvalues.txt"), emit: eigenvalues
    path "*.log"

    script:
    def prefix = meta.id
    def k = params.pca_k ?: 10
    def ref = params.reference_genome ?: ""
    def ref_arg = ref ? "--reference ${ref}" : ""
    """
    python ${params.src_dir}/python/genetics/hail/pca.py \
        --vcf ${vcf} \
        --out ${prefix} \
        --k ${k} \
        ${ref_arg} > ${prefix}.log 2>&1
    """
}

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

