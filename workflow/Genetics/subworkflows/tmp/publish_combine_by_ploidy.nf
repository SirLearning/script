#!/usr/bin/env nextflow
/*
 * Republish existing combine_ibs by-ploidy TSVs.
 *
 * The first all-ploidy runs produced these tables in their work directories
 * before the process output declaration included them. This tiny entry stages
 * the existing files and publishes them through publishDir, without rewriting
 * the analysis.
 */

nextflow.enable.dsl=2

params.job = 'test_plink'
params.mod = 'test_common_inter'
params.combine_by_ploidy_files = null

process publish_combine_by_ploidy {
    tag "publish ${tsv.name}"
    cpus 1
    memory '1 GB'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "published/*.tsv", saveAs: { filename -> file(filename).name }

    input:
    path tsv

    output:
    path("published/*.tsv"), emit: tsvs

    script:
    """
    mkdir -p published
    cp ${tsv} published/
    """
}

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'publish_combine_by_ploidy: --output_dir, --job, and --mod are required'
    }
    if (!params.combine_by_ploidy_files) {
        error 'publish_combine_by_ploidy: --combine_by_ploidy_files is required'
    }
    def files = params.combine_by_ploidy_files.toString()
        .split(',')
        .collect { it.trim() }
        .findAll { it }
    def paths = files.collect { file(it, checkIfExists: true) }
    publish_combine_by_ploidy(Channel.fromList(paths))
}
