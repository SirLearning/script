#!/usr/bin/env nextflow
/*
 * Replot variant MAF/miss distributions and MAF-vs-miss regression with sample_group coloring.
 *
 * Prerequisite: merged {A,B,D,Others}_test.plink2 and variant/*.info.{afreq,vmiss} under process/{mod}/.
 *
 * Launch from a vmap4 run folder. Example:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/08stats.genome/101run_variant_maf_miss_group_replot/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/08stats.genome/101run_variant_maf_miss_group_replot
 *   nextflow run .../subworkflows/tmp/variant_maf_miss_group_replot.nf -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 --job test_plink --mod test_thin
 */

nextflow.enable.dsl=2

include { variant_maf_miss_group_stats } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = null

workflow {
    if (!params.mod) {
        error 'variant_maf_miss_group_replot: --mod is required (e.g. test_thin)'
    }
    if (!params.output_dir || !params.job) {
        error 'variant_maf_miss_group_replot: --output_dir and --job are required'
    }

    def proc = "${params.output_dir}/${params.job}/process/${params.mod}"
    def subgenomes = ['A', 'B', 'D', 'Others']

    ch = channel.from(
        subgenomes.collect { sg ->
            tuple(
                sg,
                "sub_${sg}",
                file("${proc}/variant/${sg}.info.afreq"),
                file("${proc}/variant/${sg}.info.vmiss"),
            )
        }
    )
    variant_maf_miss_group_stats(ch)
}
