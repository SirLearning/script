#!/usr/bin/env nextflow
/*
 * Cohort 1 Mb population depth grid + QC scatters from published mapq_check TSVs.
 *
 * Publishes (default unfiltered depth, ``-c conf/mapq_bin_popdep_all.config``):
 *   {output_dir}/{job}/stats/{mod}/info/mapq_popdep/allmapq_1mb*.popdep_bins.tsv
 *   {output_dir}/{job}/stats/{mod}/plots/mapq_popdep/allmapq_1mb*.png
 *
 * Legacy MQ30+ approximation: use ``-c conf/mapq_bin_popdep.config`` (prefix mapq30_1mb).
 *
 * Launch from a vmap4 run folder:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/11aneuploidy/22run_mapq30_1mb_popdep/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/11aneuploidy/22run_mapq30_1mb_popdep
 *   nextflow run .../subworkflows/tmp/mapq_bin_popdep.nf \\
 *     -c .../nextflow.config \\
 *     -c .../conf/mapq_bin_popdep.config \\
 *     --home_dir /data/home/tusr1/01projects/vmap4 \\
 *     --user_dir /data/home/tusr1 \\
 *     --src_dir /data/home/tusr1/git/script/src \\
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \\
 *     --job benchmark \\
 *     --mod aneuploidy_check
 */

nextflow.enable.dsl=2

include { mapq_bin_popdep_build; mapq_bin_popdep_replot } from '../../modules/local/germplasm/mapq_bin_popdep.nf'

params.job = 'benchmark'
params.mod = 'aneuploidy_check'
params.mapq_popdep_prefix = 'allmapq_1mb'
params.mapq_popdep_depth_mode = 'all'
params.mapq_popdep_engine = 'polars'
params.mapq_popdep_cpus = 32
params.mapq_popdep_min_reads = 1
params.mapq_popdep_min_samples = 3
params.mapq_check_dir = null
params.mapq_popdep_plot_only = false

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'mapq_bin_popdep: --output_dir, --job, and --mod are required'
    }

    def prefix = params.mapq_popdep_prefix ?: 'mapq30_1mb'
    def popdep_info = "${params.output_dir}/${params.job}/stats/${params.mod}/info/mapq_popdep"

    if (params.mapq_popdep_plot_only) {
        def grid_tsv = file("${popdep_info}/${prefix}.popdep_bins.tsv", checkIfExists: true)
        mapq_bin_popdep_replot(grid_tsv)
    } else {
        def mapq_dir = params.mapq_check_dir
            ?: "${params.output_dir}/${params.job}/stats/${params.mod}/info/mapq_check"
        mapq_bin_popdep_build(mapq_dir)
    }
}
