nextflow.enable.dsl=2

include { build_hail_from_vcf } from './database_build.nf'
include {
    annotate_vcf_bcftools
    prepare_tab_file
    generate_gene_bed
} from './database_annotate.nf'

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
