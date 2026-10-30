#!/usr/bin/env nextflow
/*
 * Publish D/WAP BAM QC distribution figures and summary tables.
 *
 * The figures and tables are produced in a vmap4 run folder by
 * plot_d_wap_qc.py and plot_wap_multiregion.py. This publish-only workflow
 * copies them into the standard stats tree:
 *   {output_dir}/{job}/stats/{mod}/plots/*.png
 *   {output_dir}/{job}/stats/{mod}/info/*.tsv
 *
 * Launch from the producing run folder:
 *   nextflow run .../subworkflows/tmp/d_wap_bam_qc_dist.nf \
 *     -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod bam_qc_dist \
 *     --qc_dir <run folder>
 */

nextflow.enable.dsl=2

params.job = 'test_plink'
params.mod = 'bam_qc_dist'
params.qc_dir = null
params.run_tag = '1'

process publish_d_wap_bam_qc_dist {
    tag "publish D/WAP BAM QC"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.tsv"

    input:
    val qc_dir
    val run_tag

    output:
    path "*.png", emit: plots
    path "*.tsv", emit: info

    script:
    """
    echo "publish run ${run_tag}"
    cp ${qc_dir}/*.png .
    cp ${qc_dir}/*.tsv .
    """
}

workflow {
    if (!params.output_dir || !params.job || !params.mod || !params.qc_dir) {
        error 'd_wap_bam_qc_dist: --output_dir, --job, --mod, and --qc_dir are required'
    }
    publish_d_wap_bam_qc_dist(params.qc_dir, params.run_tag)
}
