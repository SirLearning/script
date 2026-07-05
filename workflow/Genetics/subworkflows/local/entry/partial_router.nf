nextflow.enable.dsl=2

/*
 * Single entry for partial reruns (assess, stats redraw, wheat-from-plink).
 * Launch: nextflow run .../subworkflows/local/entry/partial_router.nf -c .../nextflow.config --partial_task <name> ...
 *
 * Tasks: assess_plink | assess_vcf | ld_redraw | mac_stats | mac_miss_bin50_sample |
 *        mac_dist_redraw | mq_missing_reg | popdep_missing_reg | chr_counts | chr_compare | chr_pi_compare |
 *        thin_common_miss_compare | rebuild_lib_stats | wheat_from_plink | abstract_mq_50_bams | main_raw_popdepth
 */

include { RUN_ASSESS_PLINK_DEBUG; RUN_ASSESS_VCF_DEBUG } from '../partial/partial_assess.nf'
include {
    RUN_LD_PLOTS_REDRAW
    RUN_MAC_STATS_FROM_GCOUNT
    RUN_MAC_MISS_BIN50_SAMPLE
    RUN_MAC_DIST_LOG_REDRAW
    RUN_MQ_MISSING_REG
    RUN_POPDEP_MISSING_REG
    RUN_CHR_VARIANT_COUNTS
    RUN_CHR_VARIANT_COMPARE
    RUN_CHR_PI_COMPARE
    RUN_THIN_COMMON_MISSING_COMPARE
    RUN_REBUILLD_LIB_STATS
} from '../partial/partial_stats.nf'
include { RUN_ABSTRACT_MQ_50_BAMS } from '../partial/partial_abstract_mq.nf'
include { RUN_MAIN_RAW_POPDEPTH } from '../partial/partial_main_raw_popdepth.nf'
include { RUN_WHEAT_FROM_PLINK } from '../wheat/wheat_integrated_study.nf'

workflow {
    if (!params.partial_task) {
        error """partial_router.nf: --partial_task is required. Choices:
          assess_plink, assess_vcf, ld_redraw, mac_stats, mac_miss_bin50_sample,
          mac_dist_redraw, mq_missing_reg, popdep_missing_reg, chr_counts, chr_compare, chr_pi_compare,
          thin_common_miss_compare, rebuild_lib_stats, wheat_from_plink, abstract_mq_50_bams, main_raw_popdepth"""
    }

    if (params.partial_task == 'assess_plink') {
        RUN_ASSESS_PLINK_DEBUG()
    } else if (params.partial_task == 'assess_vcf') {
        RUN_ASSESS_VCF_DEBUG()
    } else if (params.partial_task == 'ld_redraw') {
        RUN_LD_PLOTS_REDRAW()
    } else if (params.partial_task == 'mac_stats') {
        RUN_MAC_STATS_FROM_GCOUNT()
    } else if (params.partial_task == 'mac_miss_bin50_sample') {
        RUN_MAC_MISS_BIN50_SAMPLE()
    } else if (params.partial_task == 'mac_dist_redraw') {
        RUN_MAC_DIST_LOG_REDRAW()
    } else if (params.partial_task == 'mq_missing_reg') {
        RUN_MQ_MISSING_REG()
    } else if (params.partial_task == 'popdep_missing_reg') {
        RUN_POPDEP_MISSING_REG()
    } else if (params.partial_task == 'chr_counts') {
        RUN_CHR_VARIANT_COUNTS()
    } else if (params.partial_task == 'chr_compare') {
        RUN_CHR_VARIANT_COMPARE()
    } else if (params.partial_task == 'chr_pi_compare') {
        RUN_CHR_PI_COMPARE()
    } else if (params.partial_task == 'thin_common_miss_compare') {
        RUN_THIN_COMMON_MISSING_COMPARE()
    } else if (params.partial_task == 'rebuild_lib_stats') {
        RUN_REBUILLD_LIB_STATS()
    } else if (params.partial_task == 'wheat_from_plink') {
        RUN_WHEAT_FROM_PLINK()
    } else if (params.partial_task == 'abstract_mq_50_bams') {
        RUN_ABSTRACT_MQ_50_BAMS()
    } else if (params.partial_task == 'main_raw_popdepth') {
        RUN_MAIN_RAW_POPDEPTH()
    } else {
        error "Unknown partial_task: ${params.partial_task}"
    }
}
