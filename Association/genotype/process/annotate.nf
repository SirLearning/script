nextflow.enable.dsl=2
// Using bcftools annotate to add functional annotations to VCF files

workflow annotate {
    take:
    // Expect a combined channel: [ val(meta), path(vcf), val(job_config) ]
    vcf_in

    main:
    // Annotate VCF using bcftools
    annotated_vcf = annotate_vcf_bcftools(vcf_in)

    emit:
    vcf = annotated_vcf.vcf
}
