nextflow.enable.dsl=2

/*
 * Join existing WGS-core and WAP PLINK2 afreq tables (same variant IDs).
 *
 *   nextflow run .../subworkflows/tmp/dual_af_from_afreq.nf \
 *     -c .../nextflow.config \
 *     --user_dir ... --output_dir ... --job test_plink --mod dual_af \
 *     --wgs_afreq /path/to/wgs.afreq --wap_afreq /path/to/wap.afreq
 */

include { dual_af_from_afreq } from '../../modules/local/genotype/stats/stats_cohort.nf'

workflow {
    if (!params.output_dir || !params.job || !params.user_dir) {
        error 'dual_af_from_afreq: --output_dir, --job, and --user_dir are required.'
    }
    if (!params.mod) {
        error 'dual_af_from_afreq: --mod is required (publish folder).'
    }
    if (!params.wgs_afreq || !params.wap_afreq) {
        error 'dual_af_from_afreq: --wgs_afreq and --wap_afreq are required.'
    }
    dual_af_from_afreq(channel.value(tuple(
        params.dual_af_prefix ?: "${params.mod}.cohort",
        file(params.wgs_afreq, checkIfExists: true),
        file(params.wap_afreq, checkIfExists: true),
    )))
}
