#!/usr/bin/env nextflow
/*
 * GAM te residual outlier diagnostics for existing sg_cov ibs_depth_miss info tables.
 *
 * Fits te(IBS, log10 Depth), flags samples with residual ≥ upper-tail cutoff (observed missing
 * above GAM prediction only; default top gam_residual_outlier_frac by signed residual),
 * publishes group bar charts and highlighted scatter plots under stats/{mod}/.
 */

nextflow.enable.dsl=2

include { sample_gam_residual_outlier_plots } from '../../modules/local/genotype/stats/stats_sample.nf'

params.job = 'test_plink'
params.stats_mod = 'test_thin'
params.sample_topic = 'sg_cov'
params.gam_residual_outlier_frac = 0.01

workflow {
    stats_base = "${params.output_dir}/${params.job}/stats/${params.stats_mod}"
    topic = params.sample_topic

    residual_in = Channel.from(['A', 'B', 'D', 'Others'])
        .map { mod ->
            tuple(
                mod,
                topic,
                params.stats_mod,
                "${stats_base}/info/${mod}.sample.${topic}.ibs_depth_miss.info.tsv",
            )
        }

    sample_gam_residual_outlier_plots(residual_in)
}
