#!/usr/bin/env nextflow
/*
 * Pattern influence × sample quality quadrant with aneuploidy highlights.
 *
 * 1) N500 chr1 full LOO — pattern quadrant, manifest samples as gold ★
 * 2) Aneuploidy same-chr LOO — combined + per-chr pattern quadrants
 *
 * Launch:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/10stats.genome/73run_popdep_loo_pattern_quadrant_aneuploidy/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/10stats.genome/73run_popdep_loo_pattern_quadrant_aneuploidy
 *   nextflow run .../popdep_loo_pattern_quadrant_aneuploidy.nf -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink
 */

nextflow.enable.dsl=2

include { popdep_loo_quadrant_review; popdep_loo_aneuploidy_samechr_pattern_quadrant } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.popdep_loo_quadrant_influence_col = 'influence_pattern'

workflow {
    if (!params.home_dir || !params.user_dir || !params.src_dir) {
        error 'popdep_loo_pattern_quadrant_aneuploidy: --home_dir, --user_dir, --src_dir required'
    }

    def chr1_root = params.popdep_loo_publish_root
        ?: "${params.output_dir}/${params.job}/process/${params.popdep_loo_publish_slug ?: '_popdep_n500_loo_chr1_10mb'}"
    def out_prefix = params.popdep_loo_summary_prefix
        ?: "n500_chr${params.popdep_loo_chr ?: 1}_${params.popdep_loo_bench_max_bp ?: 10000000}_loo"
    def info_tsv = params.popdep_loo_info_tsv
        ?: "${chr1_root}/stats/${out_prefix}.info.tsv"
    def sample_qc = params.popdep_loo_sample_qc_tsv
        ?: "${params.home_dir}/00data/05taxaBamMap/popdep/tb.N500_gam_credible_deep.gam_credible_deep_top500.info.tsv"
    def manifest = params.popdep_loo_aneuploidy_manifest
        ?: "${params.home_dir}/00data/popdep_bench/n500_aneuploidy_samechr_loo.manifest.tsv"
    def aneu_root = params.popdep_loo_aneuploidy_compare_publish_root
        ?: "${params.output_dir}/${params.job}/process/_popdep_n500_loo_aneuploidy_samechr"

    // Reuse chr1 N500 thresholds from prior 72run when available (optional)
    def pat_thr = params.popdep_loo_quadrant_influence_threshold
    def qual_thr = params.popdep_loo_quadrant_quality_threshold

    popdep_loo_quadrant_review(
        tuple(
            file(info_tsv, checkIfExists: true),
            file(sample_qc, checkIfExists: true),
        )
    )

    popdep_loo_aneuploidy_samechr_pattern_quadrant(
        file(manifest, checkIfExists: true),
        aneu_root,
        file(sample_qc, checkIfExists: true),
        pat_thr != null ? pat_thr : 'null',
        qual_thr != null ? qual_thr : 'null',
    )
}
