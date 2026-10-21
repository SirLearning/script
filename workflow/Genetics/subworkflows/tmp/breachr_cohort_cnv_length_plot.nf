#!/usr/bin/env nextflow
/*
 * CNV length distribution plots (deletion / duplication) from published cohort segments.
 *
 * Launch from a vmap4 run folder:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/11aneuploidy/20run_breachr_cohort_cnv_length_plot/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/11aneuploidy/20run_breachr_cohort_cnv_length_plot
 *   nextflow run .../subworkflows/tmp/breachr_cohort_cnv_length_plot.nf \\
 *     -c .../nextflow.config \\
 *     -c .../conf/breachr_cohort.config \\
 *     --home_dir /data/home/tusr1/01projects/vmap4 \\
 *     --user_dir /data/home/tusr1 \\
 *     --src_dir /data/home/tusr1/git/script/src \\
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \\
 *     --job benchmark \\
 *     --mod aneuploidy_check
 */

nextflow.enable.dsl=2

include { breachr_cnv_length_plots } from '../../modules/local/germplasm/breachr_segment.nf'

params.job = 'benchmark'
params.mod = 'aneuploidy_check'
params.breachr_info_subdir = 'breachr_cohort'

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'breachr_cohort_cnv_length_plot: --output_dir, --job, and --mod are required'
    }

    def stats_base = "${params.output_dir}/${params.job}/stats/${params.mod}"
    def segments_all = file(
        "${stats_base}/info/${params.breachr_info_subdir ?: 'breachr_cohort'}/breachr_segments_all.tsv",
        checkIfExists: true,
    )

    breachr_cnv_length_plots(segments_all)
}
