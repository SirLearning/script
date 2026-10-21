#!/usr/bin/env nextflow
/*
 * log10 RelativeDepth mean-vs-SD density on the test_popdep site panel
 * (all ~100M sites). 90_000 equal-log cells (300×300). Keep the densest
 * non-empty cells until cumulative sites reach 60%, 80% and 90% (not those
 * fractions of cells). Writes the implied cell-top fraction. Heatmap counts
 * every site; keep/drop JointGrids use the published scatter (point alpha 0.2).
 * Linear RelDepth keep 60/80/90 scatters use axes [0, 10]; log JointGrids keep
 * the observed log10 range.
 * --popdep_log_counts_npz reuses a previous histogram and skips the 100M scan.
 * --popdep_log_plot_keep_drop false writes the site-coverage curve without
 * replotting keep/drop JointGrids.
 *
 * Launch from a vmap4 run folder:
 *   nextflow run .../subworkflows/tmp/test_popdep_log_density.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_popdep
 */

nextflow.enable.dsl=2

include { test_popdep_log_density } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = 'test_popdep'
params.test_popdep_sites_dir = null
params.test_popdep_scatter_tsv = null
params.popdep_log_n_cells = 90000
params.popdep_log_counts_npz = null
params.popdep_log_plot_keep_drop = true

workflow {
    if (!params.output_dir || !params.job) {
        error 'test_popdep_log_density: --output_dir and --job are required'
    }
    if (!params.user_dir) {
        error 'test_popdep_log_density: --user_dir is required (conda stats)'
    }
    if (!params.src_dir) {
        error 'test_popdep_log_density: --src_dir is required'
    }

    log.info "${params.c_green}test_popdep log10 RelDepth density n_cells=${params.popdep_log_n_cells}${params.c_reset}"
    test_popdep_log_density(Channel.value(true))
}
