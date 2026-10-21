#!/usr/bin/env nextflow
/*
 * Correct test_popdep site tables: invert June CrossChr zero-padded nTaxa
 * (A=7675, B=7594, D=7738, Others=8285) back to n=500. chr032 (PopDepFull)
 * is copied unchanged. Then redraw class-coloured popdep plots.
 *
 * Reads existing {output_dir}/{job}/process/test_popdep/sites/*.popdep.sites.tsv
 * and republishes over the same process/stats trees.
 *
 * Launch from a vmap4 run folder:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/14stats.genome/23run_test_popdep_ntaxa500_rescale/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/14stats.genome/23run_test_popdep_ntaxa500_rescale
 *   nextflow run .../subworkflows/tmp/test_popdep_ntaxa_rescale.nf -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink \
 *     --mod test_popdep
 */

nextflow.enable.dsl=2

include { rescale_test_popdep_chr_sites } from '../../modules/local/genotype/processor/processor_test.nf'
include { test_popdep_class_plots } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = 'test_popdep'
params.test_popdep_sites_dir = null

workflow {
    if (!params.mod) {
        error 'test_popdep_ntaxa_rescale: --mod test_popdep is required.'
    }
    if (!params.output_dir || !params.job) {
        error 'test_popdep_ntaxa_rescale: --output_dir and --job are required.'
    }

    ch_chr = Channel.fromList((0..44).collect { String.format('chr%03d', it) })
    ch_sites = rescale_test_popdep_chr_sites(ch_chr)
    test_popdep_class_plots(ch_sites.sites.collect().map { true })
}
