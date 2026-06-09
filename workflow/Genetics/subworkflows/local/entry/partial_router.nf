nextflow.enable.dsl=2

/*
 * Single entry for partial reruns (assess, stats redraw, wheat-from-plink).
 * Launch: nextflow run .../subworkflows/local/entry/partial_router.nf -c .../nextflow.config --partial_task <name> ...
 *
 * Tasks: assess_plink | assess_vcf | ld_redraw | mac_stats | mac_dist_redraw |
 *        chr_counts | chr_compare | rebuild_lib_stats | wheat_from_plink
 */

include { RUN_ASSESS_PLINK_DEBUG; RUN_ASSESS_VCF_DEBUG } from '../partial/partial_assess.nf'
include {
    RUN_LD_PLOTS_REDRAW
    RUN_MAC_STATS_FROM_GCOUNT
    RUN_MAC_DIST_LOG_REDRAW
    RUN_CHR_VARIANT_COUNTS
    RUN_CHR_VARIANT_COMPARE
    RUN_REBUILLD_LIB_STATS
} from '../partial/partial_stats.nf'
include { RUN_WHEAT_FROM_PLINK } from '../wheat/wheat_integrated_study.nf'

workflow {
    if (!params.partial_task) {
        error """partial_router.nf: --partial_task is required. Choices:
          assess_plink, assess_vcf, ld_redraw, mac_stats, mac_dist_redraw,
          chr_counts, chr_compare, rebuild_lib_stats, wheat_from_plink"""
    }

    if (params.partial_task == 'assess_plink') {
        RUN_ASSESS_PLINK_DEBUG()
    } else if (params.partial_task == 'assess_vcf') {
        RUN_ASSESS_VCF_DEBUG()
    } else if (params.partial_task == 'ld_redraw') {
        RUN_LD_PLOTS_REDRAW()
    } else if (params.partial_task == 'mac_stats') {
        RUN_MAC_STATS_FROM_GCOUNT()
    } else if (params.partial_task == 'mac_dist_redraw') {
        RUN_MAC_DIST_LOG_REDRAW()
    } else if (params.partial_task == 'chr_counts') {
        RUN_CHR_VARIANT_COUNTS()
    } else if (params.partial_task == 'chr_compare') {
        RUN_CHR_VARIANT_COMPARE()
    } else if (params.partial_task == 'rebuild_lib_stats') {
        RUN_REBUILLD_LIB_STATS()
    } else if (params.partial_task == 'wheat_from_plink') {
        RUN_WHEAT_FROM_PLINK()
    } else {
        error "Unknown partial_task: ${params.partial_task}"
    }
}
