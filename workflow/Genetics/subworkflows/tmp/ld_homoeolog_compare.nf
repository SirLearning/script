#!/usr/bin/env nextflow
/*
 * Homoeolog LD compare: three vertical axes (A/B/D), one colour per chromosome
 * number. Reads published *.variant.ld_decay.by_chr.info.tsv (no vcor / PLINK2).
 *
 *   nextflow run .../subworkflows/tmp/ld_homoeolog_compare.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_wap_common_refchr
 */

nextflow.enable.dsl=2

include { variant_ld_homoeolog_compare_plot } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = null

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'ld_homoeolog_compare: --output_dir, --job, and --mod are required'
    }
    def info_dir = "${params.output_dir}/${params.job}/stats/${params.mod}/info"
    def a_by_chr = file("${info_dir}/A.variant.ld_decay.by_chr.info.tsv")
    def b_by_chr = file("${info_dir}/B.variant.ld_decay.by_chr.info.tsv")
    def d_by_chr = file("${info_dir}/D.variant.ld_decay.by_chr.info.tsv")
    if (!a_by_chr.exists() || !b_by_chr.exists() || !d_by_chr.exists()) {
        error "ld_homoeolog_compare: missing A/B/D by_chr TSV under ${info_dir}"
    }
    variant_ld_homoeolog_compare_plot(a_by_chr, b_by_chr, d_by_chr)
}
