#!/usr/bin/env nextflow
/*
 * Preliminary tables/plots from published SINGER ARG texts (no new SINGER run).
 *
 *   nextflow run .../subworkflows/tmp/singer_arg_preview.nf \
 *     -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_common_inter_2mb_wap_china_refchr
 */

nextflow.enable.dsl=2

include { singer_arg_preview } from '../../modules/local/genotype/stats/stats_singer.nf'

params.job = 'test_plink'
params.mod = 'test_common_inter_2mb_wap_china_refchr'

workflow {
    if (!params.output_dir || !params.job || !params.mod || !params.user_dir) {
        error 'singer_arg_preview: --output_dir, --job, --mod, and --user_dir are required'
    }
    def pub = "${params.output_dir}/${params.job}/process/${params.mod}"
    def windows_tsv = file("${pub}/singer_windows.tsv")
    def singer_dir = file("${pub}/singer")
    if (!windows_tsv.exists()) {
        error "singer_arg_preview: missing ${windows_tsv}"
    }
    if (!singer_dir.exists()) {
        error "singer_arg_preview: missing ${singer_dir}"
    }
    log.info "${params.c_green}SINGER ARG preview${params.c_reset}: ${windows_tsv} + ${singer_dir}"
    singer_arg_preview(windows_tsv, singer_dir.toString())
}
