nextflow.enable.dsl=2

process build_hail_from_vcf {
    tag "build_hail_from_vcf"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.mt") into ch_hail_mt

    script:
    """
    hail_script_build_from_vcf.py \\
        --input_vcf ${vcf} \\
        --output_mt ${meta}.mt
    """
}