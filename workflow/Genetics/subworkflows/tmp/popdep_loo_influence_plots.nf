#!/usr/bin/env nextflow
/*
 * Plot leave-one-out PopDep influence summary from an existing *.info.tsv.
 *
 * Prerequisite: popdep_n500_loo_chr1.nf (or popdep_loo_influence_summary) has written
 *   {popdep_loo_publish_root}/stats/n500_chr{chr}_{max_bp}_loo.info.tsv
 *
 * Publishes under {popdep_loo_publish_root}/stats/{plots,info,logs}
 *
 * Launch from a vmap4 run folder. Example:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/10stats.genome/61run_popdep_n500_loo_influence_plots/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/10stats.genome/61run_popdep_n500_loo_influence_plots
 *   nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/tmp/popdep_loo_influence_plots.nf \
 *     -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink
 */

nextflow.enable.dsl=2

include { popdep_loo_influence_plots } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'

workflow {
    if (!params.home_dir || !params.user_dir || !params.src_dir) {
        error 'popdep_loo_influence_plots: --home_dir, --user_dir, --src_dir are required'
    }
    if (!params.output_dir || !params.job) {
        error 'popdep_loo_influence_plots: --output_dir and --job are required'
    }

    def publish_root = params.popdep_loo_publish_root
        ?: "${params.output_dir}/${params.job}/process/${params.popdep_loo_publish_slug}"
    def out_prefix = params.popdep_loo_summary_prefix
        ?: "n500_chr${params.popdep_loo_chr}_${params.popdep_loo_bench_max_bp}_loo"
    def info_tsv = params.popdep_loo_info_tsv
        ?: "${publish_root}/stats/${out_prefix}.info.tsv"

    if (!file(info_tsv).exists()) {
        error "popdep_loo_influence_plots: missing influence summary ${info_tsv}"
    }

    log.info "popdep_loo_influence_plots: info_tsv=${info_tsv}"
    log.info "popdep_loo_influence_plots: publish_root=${publish_root}"

    popdep_loo_influence_plots(file(info_tsv, checkIfExists: true))
}
