#!/usr/bin/env nextflow
/*
 * Replot breachr mosdepth reconciliation summary charts from published TSVs.
 *
 * Launch from a vmap4 run folder (see doc/project_knowledge/workspace/vmap4_11aneuploidy.yaml):
 *   mkdir -p /data/home/tusr1/01projects/vmap4/11aneuploidy/18run_breachr_reconcile_replot/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/11aneuploidy/18run_breachr_reconcile_replot
 *   nextflow run .../subworkflows/tmp/breachr_reconcile_replot.nf \\
 *     -c .../nextflow.config \\
 *     --home_dir /data/home/tusr1/01projects/vmap4 \\
 *     --user_dir /data/home/tusr1 \\
 *     --src_dir /data/home/tusr1/git/script/src \\
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \\
 *     --job benchmark \\
 *     --mod aneuploidy_check
 */

nextflow.enable.dsl=2

include { breachr_reconcile_plots } from '../../modules/local/germplasm/breachr_segment.nf'

params.job = 'benchmark'
params.mod = 'aneuploidy_check'

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'breachr_reconcile_replot: --output_dir, --job, and --mod are required'
    }

    def stats_base = "${params.output_dir}/${params.job}/stats/${params.mod}"
    def reconcile_all = file("${stats_base}/info/breachr_check/breachr_reconcile_all.tsv", checkIfExists: true)
    def segments_all = file("${stats_base}/info/breachr_check/breachr_segments_all.tsv", checkIfExists: true)
    def flagged_file = file("${stats_base}/info/flagged_samples.tsv", checkIfExists: true)

    breachr_reconcile_plots(reconcile_all, flagged_file, segments_all)
}
