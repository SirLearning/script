#!/usr/bin/env nextflow
/*
 * Run ct_depth_with_mosdepth for one or more absolute BAM paths.
 *
 * Launch from a vmap4 run folder (level-1 module, not under 00data). Example:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/12depth/01run_D_0132_mosdepth
 *   cd /data/home/tusr1/01projects/vmap4/12depth/01run_D_0132_mosdepth
 *   nextflow run .../subworkflows/tmp/mosdepth_single_bam.nf \\
 *     -c .../nextflow.config -c .../conf/mosdepth_single_bam.config \\
 *     --user_dir /data/home/tusr1 \\
 *     --home_dir /data/home/tusr1/01projects/vmap4 \\
 *     --depth_publish_dir /data/home/tusr1/01projects/vmap4/00data/04depth/04D \\
 *     --bam_list /path/to/bam_list.tsv
 *
 * bam_list TSV (no header): prefix<TAB>absolute_bam_path
 * BAI is resolved as {bam_path}.bai and staged alongside the BAM.
 * Use prefix D_0132.bam to match existing 04depth/04D mosdepth file names.
 */

nextflow.enable.dsl=2

params.depth_publish_dir = "${params.home_dir}/00data/04depth/04D"
params.bam_list = null

process ct_depth_with_mosdepth {
    tag "ct_depth_${prefix}"

    conda "${params.user_dir}/miniconda3/envs/stats"

    publishDir "${params.depth_publish_dir}", mode: 'copy', pattern: '*.mosdepth.*'
    publishDir "${params.depth_publish_dir}", mode: 'copy', pattern: 'ct_depth_*.log'

    input:
    tuple val(prefix), path(input_bam), path(input_bai)

    output:
    tuple val(prefix), path("${prefix}.mosdepth.summary.txt"), emit: summary

    script:
    """
    set -euo pipefail
    exec > ct_depth_${prefix}.log 2>&1

    echo "Calculating coverage depth for ${input_bam} with mosdepth..."
    mosdepth \\
        -t ${task.cpus} \\
        -n ${prefix} \\
        ${input_bam}
    echo "Coverage depth calculated and saved as ${prefix}.mosdepth.*"
    """
}

workflow {
    if (!params.bam_list) {
        error 'mosdepth_single_bam: --bam_list is required (TSV: prefix<TAB>absolute_bam_path)'
    }

    Channel
        .fromPath(params.bam_list)
        .splitCsv(sep: '\t', header: false)
        .map { row ->
            def prefix = row[0] as String
            def bam_path = row[1] as String
            def bam = file(bam_path)
            def bai = file("${bam_path}.bai")
            tuple(prefix, bam, bai)
        }
        | ct_depth_with_mosdepth
}
