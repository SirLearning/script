#!/usr/bin/env nextflow
/*
 * Plot confirmed whole-chromosome aneuploidy counts by chromosome from BreachR reconciliation TSVs.
 *
 * Launch from a vmap4 run folder:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/08stats.genome/45run_aneuploidy_chromosome_plot/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/08stats.genome/45run_aneuploidy_chromosome_plot
 *   nextflow run .../subworkflows/tmp/aneuploidy_chromosome_plot.nf \\
 *     -c .../nextflow.config \\
 *     --user_dir /data/home/tusr1 \\
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \\
 *     --job benchmark \\
 *     --mod aneuploidy_check \\
 *     --aneuploidy_breachr_input_dir /data1/dazheng_tusr1/vmap4.VCF.v1/benchmark/stats/aneuploidy_check/info/breachr_check
 */

nextflow.enable.dsl=2

include { plot_confirmed_aneuploidy_chromosome } from '../../modules/local/genotype/stats/stats_aneuploidy.nf'

params.job = 'benchmark'
params.mod = 'aneuploidy_check'

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'aneuploidy_chromosome_plot: --output_dir, --job, and --mod are required'
    }

    def breachr_dir = params.aneuploidy_breachr_input_dir
        ?: "${params.output_dir}/${params.job}/stats/${params.mod}/info/breachr_check"

    plot_confirmed_aneuploidy_chromosome(breachr_dir)
}
