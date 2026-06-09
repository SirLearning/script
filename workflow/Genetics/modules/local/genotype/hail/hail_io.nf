nextflow.enable.dsl=2

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
