#!/usr/bin/env nextflow
/*
 * Local relatedness ranking from published SINGER ARG texts (no new sampling).
 *
 *   nextflow run .../subworkflows/tmp/singer_arg_relatedness.nf \
 *     -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_common_inter_2mb_wap_china_refchr
 */

nextflow.enable.dsl=2

include { singer_arg_relatedness } from '../../modules/local/genotype/stats/stats_singer.nf'

params.job = 'test_plink'
params.mod = 'test_common_inter_2mb_wap_china_refchr'
params.singer_relate_n_positions = 11
params.singer_relate_top_k = 5
params.singer_relate_top_pairs = 200

workflow {
    if (!params.output_dir || !params.job || !params.mod || !params.user_dir) {
        error 'singer_arg_relatedness: --output_dir, --job, --mod, and --user_dir are required'
    }
    def pub = "${params.output_dir}/${params.job}/process/${params.mod}"
    def windows_tsv = file("${pub}/singer_windows.tsv")
    def keep_file = file("${pub}/WAP_China.keep")
    def singer_dir = file("${pub}/singer")
    if (!windows_tsv.exists()) {
        error "singer_arg_relatedness: missing ${windows_tsv}"
    }
    if (!keep_file.exists()) {
        error "singer_arg_relatedness: missing ${keep_file}"
    }
    if (!singer_dir.exists()) {
        error "singer_arg_relatedness: missing ${singer_dir}"
    }
    log.info "${params.c_green}SINGER ARG relatedness${params.c_reset}: ${singer_dir} x ${params.singer_relate_n_positions} positions"
    singer_arg_relatedness(windows_tsv, keep_file, singer_dir.toString())
}
