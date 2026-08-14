#!/usr/bin/env nextflow
/*
 * Outlier-review plots for PopDep LOO influence vs sample BAM/depth QC.
 *
 * Prerequisite:
 *   {popdep_loo_publish_root}/stats/n500_chr{chr}_{max_bp}_loo.info.tsv
 *   N500 panel QC table (default: 00data/05taxaBamMap/popdep/tb.N500_gam_credible_deep.gam_credible_deep_top500.info.tsv)
 *
 * Publishes under {popdep_loo_publish_root}/stats/{plots,info,logs}
 *
 * Launch from a vmap4 run folder. Example:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/10stats.genome/62run_popdep_n500_loo_outlier_review/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/10stats.genome/62run_popdep_n500_loo_outlier_review
 *   nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/tmp/popdep_loo_outlier_review.nf \
 *     -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink
 */

nextflow.enable.dsl=2

include { popdep_loo_outlier_review } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'

workflow {
    if (!params.home_dir || !params.user_dir || !params.src_dir) {
        error 'popdep_loo_outlier_review: --home_dir, --user_dir, --src_dir are required'
    }
    if (!params.output_dir || !params.job) {
        error 'popdep_loo_outlier_review: --output_dir and --job are required'
    }

    def publish_root = params.popdep_loo_publish_root
        ?: "${params.output_dir}/${params.job}/process/${params.popdep_loo_publish_slug}"
    def out_prefix = params.popdep_loo_summary_prefix
        ?: "n500_chr${params.popdep_loo_chr}_${params.popdep_loo_bench_max_bp}_loo"
    def info_tsv = params.popdep_loo_info_tsv
        ?: "${publish_root}/stats/${out_prefix}.info.tsv"
    def sample_qc_tsv = params.popdep_loo_sample_qc_tsv
        ?: "${params.home_dir}/00data/05taxaBamMap/popdep/tb.N500_gam_credible_deep.gam_credible_deep_top500.info.tsv"

    if (!file(info_tsv).exists()) {
        error "popdep_loo_outlier_review: missing influence summary ${info_tsv}"
    }
    if (!file(sample_qc_tsv).exists()) {
        error "popdep_loo_outlier_review: missing sample QC table ${sample_qc_tsv}"
    }

    log.info "popdep_loo_outlier_review: info_tsv=${info_tsv}"
    log.info "popdep_loo_outlier_review: sample_qc_tsv=${sample_qc_tsv}"
    log.info "popdep_loo_outlier_review: publish_root=${publish_root}"

    popdep_loo_outlier_review(
        tuple(
            file(info_tsv, checkIfExists: true),
            file(sample_qc_tsv, checkIfExists: true),
        )
    )
}
