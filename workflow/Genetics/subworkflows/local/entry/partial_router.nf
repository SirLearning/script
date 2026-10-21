nextflow.enable.dsl=2

/*
 * Single entry for partial reruns (assess, stats redraw, wheat-from-plink).
 * Launch: nextflow run .../subworkflows/local/entry/partial_router.nf -c .../nextflow.config --partial_task <name> ...
 *
 * Tasks: assess_plink | assess_vcf | ld_redraw | ld_rare_common | mac_stats | mac_het_excess | mac_miss_bin50_sample |
 *        mac_dist_redraw | mq_missing_reg | popdep_missing_reg | chr_counts | chr_compare | chr_pi_compare |
 *        thin_common_miss_compare | rebuild_lib_stats | wheat_from_plink | abstract_mq_50_bams |
 *        aneuploidy_mapq | aneuploidy_breachr | aneuploidy_breachr_cohort | main_raw_popdepth |
 *        sample_manifest | cohort_freq | annotate_vep | annotate_snpeff | site_anno | gene_expr | variant_sets |
 *        assoc_wap | phenotype_wap
 *
 * Cohort-to-association chain (each step consumes the previous step's publish):
 *   sample_manifest -> cohort_freq -> annotate_vep -> variant_sets -> assoc_wap
 *   phenotype_wap feeds assoc_wap in parallel.
 *   site_anno and gene_expr are optional inputs to variant_sets: they add the
 *   conservation, regulatory and expressed-gene sets on top of the coding ones.
 */

include { RUN_ASSESS_PLINK_DEBUG; RUN_ASSESS_VCF_DEBUG } from '../partial/partial_assess.nf'
include {
    RUN_LD_PLOTS_REDRAW
    RUN_LD_RARE_COMMON
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
include { RUN_MAC_HET_EXCESS } from '../partial/partial_mac_het_excess.nf'
include { RUN_ANEUPLOIDY_MAPQ } from '../partial/partial_aneuploidy_mapq.nf'
include { RUN_BREACHR_FLAGGED } from '../partial/partial_breachr_flagged.nf'
include { RUN_BREACHR_COHORT } from '../partial/partial_breachr_cohort.nf'
include { RUN_MAIN_RAW_POPDEPTH } from '../partial/partial_main_raw_popdepth.nf'
include { RUN_WHEAT_FROM_PLINK } from '../wheat/wheat_integrated_study.nf'
include { RUN_SAMPLE_MANIFEST } from '../partial/partial_sample_manifest.nf'
include { RUN_COHORT_FREQ } from '../partial/partial_cohort_freq.nf'
include { RUN_ANNOTATE_VEP } from '../partial/partial_annotate_vep.nf'
include { RUN_ANNOTATE_SNPEFF } from '../partial/partial_annotate_snpeff.nf'
include { RUN_SITE_ANNO } from '../partial/partial_site_anno.nf'
include { RUN_GENE_EXPR } from '../partial/partial_gene_expr.nf'
include { RUN_VARIANT_SETS } from '../partial/partial_variant_sets.nf'
include { RUN_ASSOC_WAP } from '../partial/partial_assoc_wap.nf'
include { RUN_PHENOTYPE_WAP } from '../partial/partial_phenotype_wap.nf'

workflow {
    if (!params.partial_task) {
        error """partial_router.nf: --partial_task is required. Choices:
          assess_plink, assess_vcf, ld_redraw, ld_rare_common, mac_stats, mac_miss_bin50_sample,
          mac_dist_redraw, mq_missing_reg, popdep_missing_reg, chr_counts, chr_compare, chr_pi_compare,
          thin_common_miss_compare, rebuild_lib_stats, wheat_from_plink, abstract_mq_50_bams,
          aneuploidy_mapq, aneuploidy_breachr, aneuploidy_breachr_cohort, main_raw_popdepth,
          mac_het_excess, sample_manifest, cohort_freq, annotate_vep, annotate_snpeff,
          site_anno, gene_expr, variant_sets, assoc_wap, phenotype_wap"""
    }

    if (params.partial_task == 'assess_plink') {
        RUN_ASSESS_PLINK_DEBUG()
    } else if (params.partial_task == 'assess_vcf') {
        RUN_ASSESS_VCF_DEBUG()
    } else if (params.partial_task == 'ld_redraw') {
        RUN_LD_PLOTS_REDRAW()
    } else if (params.partial_task == 'ld_rare_common') {
        RUN_LD_RARE_COMMON()
    } else if (params.partial_task == 'mac_stats') {
        RUN_MAC_STATS_FROM_GCOUNT()
    } else if (params.partial_task == 'mac_het_excess') {
        RUN_MAC_HET_EXCESS()
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
    } else if (params.partial_task == 'aneuploidy_mapq') {
        RUN_ANEUPLOIDY_MAPQ()
    } else if (params.partial_task == 'aneuploidy_breachr') {
        RUN_BREACHR_FLAGGED()
    } else if (params.partial_task == 'aneuploidy_breachr_cohort') {
        RUN_BREACHR_COHORT()
    } else if (params.partial_task == 'main_raw_popdepth') {
        RUN_MAIN_RAW_POPDEPTH()
    } else if (params.partial_task == 'sample_manifest') {
        RUN_SAMPLE_MANIFEST()
    } else if (params.partial_task == 'cohort_freq') {
        RUN_COHORT_FREQ()
    } else if (params.partial_task == 'annotate_vep') {
        RUN_ANNOTATE_VEP()
    } else if (params.partial_task == 'annotate_snpeff') {
        RUN_ANNOTATE_SNPEFF()
    } else if (params.partial_task == 'site_anno') {
        RUN_SITE_ANNO()
    } else if (params.partial_task == 'gene_expr') {
        RUN_GENE_EXPR()
    } else if (params.partial_task == 'variant_sets') {
        RUN_VARIANT_SETS()
    } else if (params.partial_task == 'assoc_wap') {
        RUN_ASSOC_WAP()
    } else if (params.partial_task == 'phenotype_wap') {
        RUN_PHENOTYPE_WAP()
    } else {
        error "Unknown partial_task: ${params.partial_task}"
    }
}
