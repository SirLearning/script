#!/usr/bin/env nextflow
/*
 * Grouped-bar histogram of PCoA edge lengths, pass vs fail, one figure per
 * cohort x ploidy group. Same visual language as the by-ploidy IBS histograms
 * (plot_grouped_histogram: side-by-side bars, density per series, dpi 300).
 *
 * Reads the published IBS-PCoA coordinates in process/<mod>/pcoa_ibs/.
 *
 *   nextflow run .../subworkflows/tmp/combine_ibs_edge_distance_hist.nf \
 *     -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_common_inter
 */

nextflow.enable.dsl=2

params.job = 'test_plink'
params.mod = 'test_common_inter'
params.ploidies = 'AA,AABB,AABBDD,DD'

process combine_ibs_edge_distance_hist {
    tag "combine IBS PCoA edge hist: ${prefix} ${group}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 2
    memory '16.GB'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(prefix), val(threshold), val(group), path(pairs), path(eigenvec)

    output:
    path("*.png"), emit: plots, optional: true
    path("${prefix}.pcoa_edge_distance.${group}.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${prefix}.pcoa_edge_distance.${group}.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.sample.combine_ibs_validate import (
        plot_combine_ibs_edge_distance,
    )

    plot_combine_ibs_edge_distance(
        "${pairs}",
        "${eigenvec}",
        "${prefix}",
        subgenome="${group}",
        threshold=${threshold},
    )
    """
}

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'combine_ibs_edge_distance_hist: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir || !params.src_dir) {
        error 'combine_ibs_edge_distance_hist: --user_dir and --src_dir are required'
    }
    def info = "${params.output_dir}/${params.job}/stats/${params.mod}/info"
    def process_dir = "${params.output_dir}/${params.job}/process/${params.mod}"

    def requested = params.ploidies.toString().split(',').collect { it.trim() }
    def eigenvecs = Channel
        .fromPath("${process_dir}/pcoa_ibs/*.pcoa.eigenvec")
        .map { file -> tuple(file.name.split('\\.')[0], file) }
        .filter { group, file -> requested.contains(group) }

    // Cutoff per cohort, plus the populations that cohort actually spans.
    def cohorts = Channel.of(
        tuple('all.sample.combine_ibs', 0.995 as double, 'AA,AABB,AABBDD,DD'),
        tuple('all.sample.combine_ibs_ipk_hzau_watkins', 0.99 as double, 'AABBDD'),
    )
    def specs = eigenvecs
        .combine(cohorts)
        .filter { item ->
            // item = (group, eigenvec, prefix, threshold, allowed)
            item[4].toString().split(',').collect { it.trim() }.contains(item[0])
        }
        .map { group, eigenvec, prefix, threshold, allowed ->
            tuple(
                prefix,
                threshold,
                group,
                file("${info}/${prefix}.pairs.tsv", checkIfExists: true),
                eigenvec,
            )
        }
    combine_ibs_edge_distance_hist(specs)
}
