#!/usr/bin/env nextflow
/*
 * Hybrid heatmap: F_MISS from test_thin sg_cov, IBS from test_common_thin sg_cov,
 * subgenome mosdepth depth on Y. Publishes under test_thin stats (plots + info).
 *
 * Requires existing {mod}.sample.sg_cov.ibs_depth_miss.info.tsv in both job stats trees.
 * Pass --thin_stats_mod / --common_stats_mod (or rely on defaults below) so publishDir resolves.
 */

nextflow.enable.dsl=2

include { sample_hybrid_thinmiss_commonibs_heatmap } from '../../modules/local/genotype/stats/stats_sample.nf'

params.job = 'test_plink'
params.thin_stats_mod = 'test_thin'
params.common_stats_mod = 'test_common_thin'
params.sample_topic = 'sg_cov'

workflow {
    stats_base = "${params.output_dir}/${params.job}/stats"
    topic = params.sample_topic

    hybrid_in = Channel.from(['A', 'B', 'D', 'Others'])
        .map { mod ->
            tuple(
                mod,
                topic,
                "${stats_base}/${params.thin_stats_mod}/info/${mod}.sample.${topic}.ibs_depth_miss.info.tsv",
                "${stats_base}/${params.common_stats_mod}/info/${mod}.sample.${topic}.ibs_depth_miss.info.tsv",
            )
        }

    sample_hybrid_thinmiss_commonibs_heatmap(hybrid_in)
}
