#!/usr/bin/env nextflow
/*
 * test_thin four-step cascade on the n500 RelMean×RelCV reliable libraries.
 *
 * For keep 80% and keep 60% separately:
 *   PopDep reliable BED → excess het (P>=1e-6) → MAC≥5 → WAP MAC≥2.
 * Genome counts are test_thin counts / thin_rate (default 0.001).
 * Does not overwrite the older mean-vs-SD density cascade tables.
 *
 * Launch from a vmap4 run folder:
 *   nextflow run .../subworkflows/tmp/filter_cascade_cv_reliable.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_thin \
 *     --popdep_cv_reliable_dir /data1/.../process/_popdep_n500_cv_reliable
 */

nextflow.enable.dsl=2

include { variant_filter_cascade_cv_reliable } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = 'test_thin'
params.filter_cascade_het_mod = 'test_common_only'
params.filter_cascade_wap_mod = 'test_wap'
params.filter_cascade_wap_keep = null
params.filter_cascade_het_alpha = '1e-6'
params.popdep_cv_reliable_dir = null

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'filter_cascade_cv_reliable: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir) {
        error 'filter_cascade_cv_reliable: --user_dir is required (conda stats)'
    }
    if (!params.src_dir) {
        error 'filter_cascade_cv_reliable: --src_dir is required'
    }
    def bed_root = params.popdep_cv_reliable_dir
    if (!bed_root) {
        error 'filter_cascade_cv_reliable: --popdep_cv_reliable_dir is required'
    }
    file("${bed_root}/keep80/chr001.popdep.cv_reliable.keep80.bed.gz", checkIfExists: true)
    file("${bed_root}/keep60/chr001.popdep.cv_reliable.keep60.bed.gz", checkIfExists: true)

    def wap_mod = params.filter_cascade_wap_mod ?: 'test_wap'
    def keep_name = params.filter_cascade_wap_keep ?: "${wap_mod}.keep"
    file("${params.output_dir}/${params.job}/process/${wap_mod}/${keep_name}", checkIfExists: true)

    log.info "${params.c_green}test_thin cv-reliable cascade keep 80 and 60 from ${bed_root}${params.c_reset}"
    variant_filter_cascade_cv_reliable(Channel.value(true))
}
