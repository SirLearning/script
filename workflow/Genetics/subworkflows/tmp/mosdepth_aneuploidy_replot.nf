#!/usr/bin/env nextflow
/*
 * Replot mosdepth aneuploidy rel_depth distributions from published chr_depth TSVs.
 *
 * Launch from a vmap4 run folder (see doc/project_knowledge/workspace/vmap4_11aneuploidy.yaml):
 *   mkdir -p /data/home/tusr1/01projects/vmap4/11aneuploidy/02run_replot/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/11aneuploidy/02run_replot
 *   nextflow run .../subworkflows/tmp/mosdepth_aneuploidy_replot.nf \\
 *     -c .../nextflow.config \\
 *     --home_dir /data/home/tusr1/01projects/vmap4 \\
 *     --user_dir /data/home/tusr1 \\
 *     --src_dir /data/home/tusr1/git/script/src \\
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \\
 *     --job benchmark \\
 *     --mod aneuploidy_check
 */

nextflow.enable.dsl=2

include { mosdepth_aneuploidy_replot } from '../../modules/local/germplasm/mosdepth_aneuploidy.nf'

params.job = 'benchmark'
params.mod = 'aneuploidy_check'
params.group_file = "${params.output_dir}/meta_data/sample_group.txt"

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'mosdepth_aneuploidy_replot: --output_dir, --job, and --mod are required'
    }

    def stats_base = "${params.output_dir}/${params.job}/stats/${params.mod}"
    def chr_depth_glob = "${stats_base}/info/chr_depth/*.depth.tsv"
    def group_file = file(params.group_file)
    def flagged_file = file("${stats_base}/info/flagged_samples.tsv")

    chr_tsvs = Channel.fromPath(chr_depth_glob).collect()
    mosdepth_aneuploidy_replot(chr_tsvs, group_file, flagged_file)
}
