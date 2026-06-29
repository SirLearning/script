#!/usr/bin/env nextflow
/*
 * Re-publish sample.sg_cov stats (mosdepth subgenome mean depth) for an existing PLINK job.
 *
 * Reads smiss/scount from params.process_dir; writes via publishDir to
 * params.output_dir/params.job/stats/params.mod/{info,plots,thresholds,logs}.
 *
 * Launch from a vmap4 run folder (see doc/NF_CMD.md). Example:
 *   nextflow run .../subworkflows/tmp/sg_cov_replot.nf -c .../nextflow.config \
 *     --mod test_thin \
 *     --process_dir /data1/.../test_plink/process/test_thin
 */

nextflow.enable.dsl=2

include {
    sample_sg_coverage_stats
    plot_subgenome_gam_ibs_depth_compare_sg
} from '../../modules/local/genotype/stats/stats_sample.nf'

params.job = 'test_plink'
params.mod = null
params.process_dir = null

workflow {
    if (!params.mod) {
        error 'sg_cov_replot: --mod is required (e.g. test_thin, test_common_thin)'
    }
    if (!params.process_dir) {
        error 'sg_cov_replot: --process_dir is required (e.g. .../test_plink/process/test_thin)'
    }

    group_file = file("${params.output_dir}/meta_data/sample_group.txt")
    sample_dir = file("${params.process_dir}/sample")

    smiss = Channel.from(['A', 'B', 'D', 'Others'])
        .map { id -> tuple(id, "sub${id}", file("${sample_dir}/${id}.info.smiss")) }

    scount = Channel.from(['A', 'B', 'D', 'Others'])
        .map { id -> tuple(id, "sub${id}", file("${sample_dir}/${id}.info.scount")) }

    sg_out = sample_sg_coverage_stats(group_file, smiss, scount)

    plot_subgenome_gam_ibs_depth_compare_sg(
        sg_out.info
            .map { id, chr, infos ->
                def files = infos instanceof List ? infos : [infos]
                tuple(id, files.find { f -> f.name.endsWith('.sample.sg_cov.ibs_depth_miss.info.tsv') })
            }
            .filter { id, info -> id in ['A', 'B', 'D'] && info != null }
            .map { id, info -> info }
            .collect()
    )
}
