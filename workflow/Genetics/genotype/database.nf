nextflow.enable.dsl=2

workflow database {
    take:
    // Expect a channel: [ val(meta), path(vcf) ]
    vcf_in

    main:
    // Build Hail MatrixTable from VCF
    hail_mt = build_hail_from_vcf(vcf_in)

    emit:
    mt = hail_mt.mt
}

process build_hail_from_vcf {
    tag "build_hail_from_vcf"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.mt"), emit: mt

    script:
    """
    hail_script_build_from_vcf.py \\
        --input_vcf ${vcf} \\
        --output_mt ${meta}.mt
    """
}