#!/usr/bin/env nextflow
/*
 * Colour published test_popdep scatter points by genome-wide PopDep
 * RelativeDepth 60%/80% trusted cells. Kept = blue, dropped = grey.
 * Mean-vs-SD JointGrids (Depth and RelativeDepth; linear and log) with
 * stacked x/y marginals. Does not rebuild sites or overwrite class plots.
 *
 * Launch from a vmap4 run folder:
 *   nextflow run .../subworkflows/tmp/test_popdep_density_keep_replot.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_popdep
 */

nextflow.enable.dsl=2

include { test_popdep_density_keep_plots } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = 'test_popdep'
params.test_popdep_scatter_tsv = null
params.popdep_density_masks = null

workflow {
    if (!params.output_dir || !params.job) {
        error 'test_popdep_density_keep_replot: --output_dir and --job are required'
    }
    if (!params.user_dir) {
        error 'test_popdep_density_keep_replot: --user_dir is required (conda stats)'
    }
    if (!params.src_dir) {
        error 'test_popdep_density_keep_replot: --src_dir is required'
    }

    def scatter = params.test_popdep_scatter_tsv
        ?: "${params.output_dir}/${params.job}/stats/${params.mod}/info/test_popdep.popdep.scatter.sites.tsv"
    def masks = params.popdep_density_masks
        ?: "${params.output_dir}/${params.job}/stats/test_thin/info/all.variant.popdep_density.masks.npz"

    log.info "${params.c_green}test_popdep density keep/drop from ${scatter}${params.c_reset}"
    test_popdep_density_keep_plots(
        tuple(file(scatter, checkIfExists: true), file(masks, checkIfExists: true))
    )
}
