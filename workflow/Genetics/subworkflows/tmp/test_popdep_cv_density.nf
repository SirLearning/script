#!/usr/bin/env nextflow
/*
 * Linear RelativeDepth mean-vs-CV density on the test_popdep site panel
 * (all ~100M sites). Mean range is [0, 6]; CV uses the observed min/max
 * among in-range positive-mean sites (no CV clip). 90_000 equal-linear
 * cells (300×300). Keep the densest non-empty cells until cumulative
 * in-range sites reach 60, 70, 80, 85, 90, 95, 97, 98 and 99%. Writes the
 * implied cell-top fraction. Heatmap counts in-range sites; keep/drop
 * JointGrids use the published scatter (point alpha 0.2, mean 0–6, CV
 * free). Also writes RelDepth CV vs mean by genomic class (mean 0–6).
 * --popdep_cv_counts_npz reuses a previous histogram.
 * --popdep_cv_plot_keep_drop false writes the curve and class plot
 * without keep/drop JointGrids.
 *
 * Launch from a vmap4 run folder:
 *   nextflow run .../subworkflows/tmp/test_popdep_cv_density.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_popdep
 */

nextflow.enable.dsl=2

include { test_popdep_cv_density } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = 'test_popdep'
params.test_popdep_sites_dir = null
params.test_popdep_scatter_tsv = null
params.popdep_cv_n_cells = 90000
params.popdep_cv_counts_npz = null
params.popdep_cv_plot_keep_drop = true

workflow {
    if (!params.output_dir || !params.job) {
        error 'test_popdep_cv_density: --output_dir and --job are required'
    }
    if (!params.user_dir) {
        error 'test_popdep_cv_density: --user_dir is required (conda stats)'
    }
    if (!params.src_dir) {
        error 'test_popdep_cv_density: --src_dir is required'
    }

    log.info "${params.c_green}test_popdep RelMean vs CV density n_cells=${params.popdep_cv_n_cells} mean=[0, 6] CV=observed${params.c_reset}"
    test_popdep_cv_density(Channel.value(true))
}
