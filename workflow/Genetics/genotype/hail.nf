nextflow.enable.dsl=2

workflow HAIL {
    take:
    // Expect a channel: [ val(id), path(vcf), val(job_config) ]
    vcf_in

    main:
    // 1. Convert VCF to Hail MatrixTable
    mt = vcf_to_mt_hail(vcf_in.map{ id, vcf, job_config -> tuple(id, id, vcf) })

    // 2. Filter VCF using Hail
    filtered_vcf = filter_hail(mt.mt.map{ id, mtx -> tuple(id, id + '.hail', mtx.path) }, vcf_in.map{ id, vcf, job_config -> job_config })

    // 3. Hail QC
    qc_stats = hail_qc(filtered_vcf.vcf)

    // 4. Hail Kinship
    kinship = hail_kinship(filtered_vcf.vcf)

    // 5. Hail PCA
    pca_results = hail_pca(filtered_vcf.vcf)

    emit:
    vcf = filtered_vcf.vcf
    qc_stats = qc_stats.qc_stats
    kinship = kinship.kinship
    pca_scores = pca_results.scores
    pca_loadings = pca_results.loadings
    pca_eigenvalues = pca_results.eigenvalues
}

process vcf_to_mt_hail {
    tag { id ? "hail vcf2mt ${id}" : 'hail vcf2mt' }
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy'

    input:
    tuple val(id), val(prefix), path(vcf)

    output:
    tuple val(id), path("${prefix}.mt"), emit: mt

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
    tag { id ? "hail filter ${id}" : 'hail filter' }
    publishDir "${params.output_dir}/${params.job}/process", mode: 'copy'

    input:
    tuple val(id), path(vcf), val(job_config)

    output:
    tuple val(id), val("${id}.hail.filtered"), path("${id}.hail.filtered.vcf.gz"), emit: vcf

    script:
    def ref = params.reference_genome ?: ""
    def ref_arg = ref ? "--reference ${ref}" : ""
    // --- Filter Parameters ---
    // standard filter parameters (std)
    def maf_val = params.maf ?: 0.05
    def mac_val = params.mac ?: 2
    def min_alleles = params.min_alleles ?: 2
    def max_alleles = params.max_alleles ?: 2
    def max_missing = params.max_missing ?: 0.05
    // additional parameters can be added here
    def qual = params.qual ?: 30
    def hwe_pval =  params.hwe_pval ?: 1e-6
    """
    python ${params.src_dir}/python/genetics/hail/filter.py \
        --vcf ${vcf} \
        --out ${id}.hail.filtered \
        --maf ${maf_val} \
        --mac ${mac_val} \
        --min_alleles ${min_alleles} \
        --max_alleles ${max_alleles} \
        --max_missing ${max_missing} \
        ${ref_arg}
    
    mv ${id}.hail.filtered.vcf.bgz ${id}.hail.filtered.vcf.gz
    """
}


process hail_qc {
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

process hail_kinship {
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

process hail_pca {
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

process hail_normal_gwas {
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

