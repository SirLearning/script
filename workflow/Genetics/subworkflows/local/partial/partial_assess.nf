nextflow.enable.dsl=2

/*
 * Partial assess tracks (compute in processor_assess; plots in stats_assess).
 * Wired from partial_router.nf (--partial_task assess_plink | assess_vcf).
 */

include { plink2_assess_debug_slice } from '../../../modules/local/genotype/processor/processor_assess.nf'
include { assess_plink_debug_plots } from '../../../modules/local/genotype/stats/stats_assess.nf'
include { quick_count; bcftools_qc_assess } from '../../../modules/local/genotype/processor/processor_assess.nf'
include { dumpnice_vcf_qc_assess } from '../../../modules/local/genotype/stats/stats_assess.nf'

workflow RUN_ASSESS_PLINK_DEBUG {
    main:
    if (!(params.mod in ['test_thin', 'test_common_thin'])) {
        error "RUN_ASSESS_PLINK_DEBUG: params.mod must be test_thin or test_common_thin (got: ${params.mod})."
    }
    if (!params.output_dir || !params.job) {
        error 'params.output_dir and params.job are required.'
    }

    def sub = ['A', 'B', 'D', 'Others']
    def chrPick = [
        'A'     : '1',
        'B'     : '3',
        'D'     : '5',
        'Others': '0',
    ]
    def ch_in = channel.from(sub.collect { sg -> tuple(sg, "${sg}_test.plink2", chrPick[sg]) })
    plink2_assess_debug_slice(ch_in)
    def b = plink2_assess_debug_slice.out.bundle
    def plot_in = b.map { id, ac, vm, co, mh, gq -> tuple(id, ac, vm, co, mh, gq) }
    assess_plink_debug_plots(plot_in)
}

workflow RUN_ASSESS_VCF_DEBUG {
    main:
    if (!params.mod) {
        error 'params.mod is required.'
    }
    if (!params.output_dir || !params.job) {
        error 'params.output_dir and params.job are required.'
    }

    def sub = ['A', 'B', 'D', 'Others']
    def vcfDir = params.assess_vcf_dir ?: "${params.output_dir}/${params.job}/process/${params.mod}/export"
    def ch_vcf = channel.from(
        sub.collect { sg -> tuple(sg, file("${vcfDir}/${sg}.debug.vcf.gz", checkIfExists: true)) }
    )

    quick_count(ch_vcf)
    bcftools_qc_assess(ch_vcf)
    dumpnice_vcf_qc_assess(ch_vcf)
}
