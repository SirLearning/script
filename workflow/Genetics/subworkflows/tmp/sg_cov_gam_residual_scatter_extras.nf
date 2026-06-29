#!/usr/bin/env nextflow
/*
 * Replot GAM residual scatter extras from existing *.gam_residual.info.tsv (no GAM refit).
 *
 * Emits four PNGs per subgenome:
 *   *.gam_residual.dist.png — signed residual histogram
 *   *.gam_residual.dist.logy.png — same, log-scaled Y axis
 *   *.gam_residual_outliers_vs_{logdepth,ibs}.residual_cmap.png  — color = signed residual (gray/transparent at 0)
 *   *.gam_residual_outliers_vs_{logdepth,ibs}.outlier_size.png   — blue/red + size = |residual|
 */

nextflow.enable.dsl=2

include { sample_gam_residual_scatter_extras } from '../../modules/local/genotype/stats/stats_sample.nf'

params.job = 'test_plink'
params.stats_mod = 'test_thin'
params.sample_topic = 'sg_cov'

workflow {
    stats_base = "${params.output_dir}/${params.job}/stats/${params.stats_mod}"
    topic = params.sample_topic

    extras_in = Channel.from(['A', 'B', 'D', 'Others'])
        .map { mod ->
            tuple(
                mod,
                topic,
                params.stats_mod,
                "${stats_base}/info/${mod}.sample.${topic}.gam_residual.info.tsv",
            )
        }

    sample_gam_residual_scatter_extras(extras_in)
}
