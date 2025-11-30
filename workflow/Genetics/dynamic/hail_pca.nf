nextflow.enable.dsl=2

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
