nextflow.enable.dsl=2

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
