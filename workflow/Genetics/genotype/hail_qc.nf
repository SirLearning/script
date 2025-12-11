nextflow.enable.dsl=2

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
