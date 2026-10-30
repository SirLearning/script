#!/usr/bin/env nextflow
/*
 * Coverage distribution before vs after the final combine_map merge.
 *
 * Merges come from the authoritative final.ALL.taxaMergeMap.txt: "after" is
 * the sum of each line's source coverages in final.ALL.taxaBamMap.txt, the
 * value TIGER uses for coverage = .. "Before" is every source sample's
 * Coverage-Of-All-Bams. By default only multi-member lines are used; set
 * --include_singletons to keep every line. --sources_pattern restricts the
 * lines to one cohort, with combine_map.tsv used only to label that cohort.
 * Grouped histogram bars match the existing *.vs_random.hist.png style.
 *
 * Publishes under {output_dir}/{job}/stats/{mod}/{plots,logs}:
 *   combine_map.merge_coverage.before_after.png       (KDE curves)
 *   combine_map.merge_coverage.before_after.hist.png  (bar histogram)
 *   combine_map.merge_coverage.before_after.log
 *
 * Launch from a vmap4 run folder:
 *   nextflow run .../subworkflows/tmp/combine_map_coverage_hist.nf \
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
params.merge_map = '/data/home/tusr1/git/script/src/resource/sample_combine/final.ALL.taxaMergeMap.txt'
params.combine_map = '/data/home/tusr1/git/script/src/resource/sample_combine/combine_map.tsv'
params.taxa_map = '/data/home/tusr1/01projects/vmap4/00data/05taxaBamMap/final.ALL.taxaBamMap.txt'
params.output_prefix = 'combine_map.merge_coverage.before_after'
params.bins = 100
params.x_max = 30.0
params.sources_pattern = null
params.cohort_label = null
params.include_singletons = false

process combine_map_coverage_hist {
    tag "combine_map coverage hist"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 2
    memory '8.GB'
    errorStrategy 'terminate'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple path(merge_map), path(taxa_map), path(combine_map)

    output:
    path("*.png"), emit: plots
    path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import logging
    import sys
    sys.stdout = open("${params.output_prefix}.log", "w", buffering=1)
    sys.stderr = sys.stdout
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.sample.combine_ibs_validate import (
        plot_merge_coverage_hist,
    )

    plot_merge_coverage_hist(
        "${merge_map}",
        "${taxa_map}",
        "${params.output_prefix}",
        combine_map_path="${combine_map}",
        bins=int("${params.bins}"),
        x_max=float("${params.x_max}"),
        sources_pattern=${params.sources_pattern ? '"' + params.sources_pattern + '"' : 'None'},
        cohort_label=${params.cohort_label ? '"' + params.cohort_label + '"' : 'None'},
        include_singletons=${params.include_singletons.toString().toLowerCase() in ['true', '1', 'yes'] ? 'True' : 'False'},
    )
    """
}

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'combine_map_coverage_hist: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir || !params.src_dir) {
        error 'combine_map_coverage_hist: --user_dir and --src_dir are required'
    }
    combine_map_coverage_hist(
        Channel.value([
            file(params.merge_map, checkIfExists: true),
            file(params.taxa_map, checkIfExists: true),
            file(params.combine_map, checkIfExists: true),
        ])
    )
}
