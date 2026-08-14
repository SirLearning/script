#!/usr/bin/env nextflow
/*
 * Plot aneuploidy same-chromosome LOO influence (per-chr breakdown + chr1 vs samechr comparison).
 *
 * Prerequisite:
 *   - popdep_n500_loo_aneuploidy_samechr.nf completed
 *   - popdep_loo_aneuploidy_samechr_compare.nf (or existing compare TSV)
 *   - chr1 full LOO info TSV (n500_chr1_10000000_loo.info.tsv)
 *
 * Publishes under _popdep_n500_loo_aneuploidy_samechr/stats/{plots,info,logs}
 * and per-chr plots under chr{N}/plots/
 *
 * Launch:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/10stats.genome/69run_popdep_n500_loo_aneuploidy_plots/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/10stats.genome/69run_popdep_n500_loo_aneuploidy_plots
 *   nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/tmp/popdep_loo_aneuploidy_samechr_plots.nf \
 *     -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink
 */

nextflow.enable.dsl=2

include { popdep_loo_aneuploidy_samechr_plots } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'

workflow {
    if (!params.home_dir || !params.user_dir || !params.src_dir) {
        error 'popdep_loo_aneuploidy_samechr_plots: --home_dir, --user_dir, --src_dir are required'
    }

    def publish_root = params.popdep_loo_aneuploidy_compare_publish_root
        ?: "${params.output_dir}/${params.job}/process/_popdep_n500_loo_aneuploidy_samechr"
    def compare_tsv = params.popdep_loo_aneuploidy_compare_tsv
        ?: "${publish_root}/stats/n500_aneuploidy_samechr_vs_chr1.info.tsv"
    def chr1_info = params.popdep_loo_chr1_info_tsv
        ?: "${params.output_dir}/${params.job}/process/_popdep_n500_loo_chr1_10mb/stats/n500_chr1_10000000_loo.info.tsv"
    def manifest = params.popdep_loo_aneuploidy_manifest
        ?: "${params.home_dir}/00data/popdep_bench/n500_aneuploidy_samechr_loo.manifest.tsv"

    if (!file(manifest).exists()) {
        error "popdep_loo_aneuploidy_samechr_plots: missing manifest ${manifest}"
    }
    if (!file(chr1_info).exists()) {
        error "popdep_loo_aneuploidy_samechr_plots: missing chr1 LOO summary ${chr1_info}"
    }

    log.info "popdep_loo_aneuploidy_samechr_plots: publish_root=${publish_root}"
    log.info "popdep_loo_aneuploidy_samechr_plots: manifest=${manifest}"

    popdep_loo_aneuploidy_samechr_plots(
        compare_tsv,
        publish_root,
        file(chr1_info, checkIfExists: true),
        file(manifest, checkIfExists: true),
    )
}
