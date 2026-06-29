#!/usr/bin/env nextflow
/*
 * Plot per-chromosome rel_depth (mosdepth) vs F_MISS (test_thin sample_chr_miss).
 *
 * Publish:
 *   {output_dir}/{job}/stats/{mod}/plots/{chr}.rel_depth.vs_fmiss.png
 *   {output_dir}/{job}/stats/{mod}/info/sample_chr_miss.depth_miss_merged.tsv
 *
 * Example:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/13sample_chr_miss/02run_depth_miss_plot
 *   cd /data/home/tusr1/01projects/vmap4/13sample_chr_miss/02run_depth_miss_plot
 *   nextflow run .../subworkflows/tmp/sample_chr_miss_depth_plot.nf \\
 *     -c .../nextflow.config \\
 *     --user_dir /data/home/tusr1 \\
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \\
 *     --job benchmark \\
 *     --mod sample_chr_miss \\
 *     --aneuploidy_mod aneuploidy_check
 */

nextflow.enable.dsl=2

include { plot_rel_depth_vs_fmiss } from '../../modules/local/genotype/stats/stats_sample_chr_miss.nf'

params.job = 'benchmark'
params.mod = 'sample_chr_miss'
params.aneuploidy_mod = 'aneuploidy_check'
params.group_file = "${params.output_dir}/meta_data/sample_group.txt"

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'sample_chr_miss_depth_plot: --output_dir, --job, and --mod are required'
    }

    def chr_depth_dir = "${params.output_dir}/${params.job}/stats/${params.aneuploidy_mod}/info/chr_depth"
    def smiss_long = "${params.output_dir}/${params.job}/process/${params.mod}/info/sample_chr_miss.long.tsv"
    def group_file = file(params.group_file)

    plot_rel_depth_vs_fmiss(chr_depth_dir, smiss_long, group_file)
}
