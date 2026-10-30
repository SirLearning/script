#!/usr/bin/env nextflow
/*
 * Publish already-rendered IBS histograms into the stats plots tree.
 *
 * Used for analysis figures produced in a run folder (the selfing-model
 * filtered combine-vs-random comparison); Nextflow writes the publish tree so
 * the run, the naming and the provenance stay auditable.
 *
 *   nextflow run .../subworkflows/tmp/publish_ibs_hist.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_common_inter \
 *     --src_glob '<run>/plots/*.hist.png'
 */

nextflow.enable.dsl=2

params.job = 'test_plink'
params.mod = 'test_common_inter'
params.src_glob = null

process publish_ibs_hist {
    tag "publish IBS hist: ${name}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 1
    memory '2.GB'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"

    input:
    tuple val(name), path(file, stageAs: 'src.png')

    output:
    path("*.png")

    script:
    """
    cp src.png ${name}
    """
}

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'publish_ibs_hist: --output_dir, --job, and --mod are required'
    }
    if (!params.src_glob) {
        error 'publish_ibs_hist: --src_glob is required'
    }
    publish_ibs_hist(
        Channel.fromPath(params.src_glob).map { file -> tuple(file.name, file) }
    )
}
