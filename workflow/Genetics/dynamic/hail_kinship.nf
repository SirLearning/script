nextflow.enable.dsl=2

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
