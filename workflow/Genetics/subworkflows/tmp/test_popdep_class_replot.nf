#!/usr/bin/env nextflow
/*
 * Redraw test_popdep class plots from existing chromosome site tables.
 * Does not rescale or rebuild sites. Reads params.test_popdep_sites_dir
 * (default {output_dir}/{job}/process/{mod}/sites) and republishes
 * stats/{mod}/{plots,info,logs}.
 *
 * Launch from a vmap4 run folder:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/14stats.genome/24run_test_popdep_reldepth_plots/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/14stats.genome/24run_test_popdep_reldepth_plots
 *   nextflow run .../subworkflows/tmp/test_popdep_class_replot.nf -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink \
 *     --mod test_popdep \
 *     --test_popdep_sites_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/test_popdep/sites
 */

nextflow.enable.dsl=2

include { test_popdep_class_plots } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = 'test_popdep'
params.test_popdep_sites_dir = null
params.test_popdep_scatter_tsv = null

workflow {
    if (!params.mod) {
        error 'test_popdep_class_replot: --mod test_popdep is required.'
    }
    if (!params.output_dir || !params.job) {
        error 'test_popdep_class_replot: --output_dir and --job are required.'
    }

    test_popdep_class_plots(Channel.value(true))
}
