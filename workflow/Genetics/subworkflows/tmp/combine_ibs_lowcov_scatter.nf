#!/usr/bin/env nextflow
/*
 * Scatter the lower-coverage member's depth against hom-SNP IBS for trusted
 * combine pairs.
 *
 * The Watkins pairs in the combine-group table are treated as the ground-truth
 * merge set; for each pair ``low_coverage`` is the smaller of ``coverage1`` and
 * ``coverage2``. The figures mark the 0.99 merge cutoff, a 0.999 review line,
 * and the 3x coverage floor. Pairs from combine groups with
 * ``--combine_highlight_union_n`` members (default 48, i.e. WATDE0278 and
 * WATDE0305) are drawn in a separate colour.
 *
 * The pair table's depth columns come from the dedup snapshot, which lacks
 * some genotyped Watkins samples. ``--combine_coverage_map`` (default
 * final.ALL.taxaBamMap.txt) fills those depths by sample ID so every pair
 * gets a low_coverage.
 *
 * ``--combine_highlight_pairs`` names another pair table (for example the
 * test_common_inter table) and ``--combine_highlight_min_ibs`` selects the
 * pairs whose metric exceeds that value there; they are drawn in a third
 * colour under ``--combine_highlight_label``. This lets the same pairs be
 * followed across test sets.
 *
 * Publishes under {output_dir}/{job}/stats/{mod}/{plots,logs}:
 *   all.sample.combine_ibs_ipk_hzau_watkins.watkins_lowcov_vs_homsnp.png / .linear.png
 *   all.sample.combine_ibs_ipk_hzau_watkins.watkins_lowcov_vs_homsnp.high_missing_rate.png
 *   all.sample.combine_ibs_ipk_hzau_watkins.watkins_lowcov_vs_homsnp.high_missing_rate.linear.png
 *   all.sample.combine_ibs_ipk_hzau_watkins.watkins_lowcov_vs_homsnp.missing_sum.png
 *   all.sample.combine_ibs_ipk_hzau_watkins.watkins_lowcov_vs_homsnp.missing_sum.linear.png
 *   all.sample.combine_ibs_ipk_hzau_watkins.watkins_lowcov_vs_homsnp.log
 *
 * Launch from a vmap4 run folder:
 *   nextflow run .../subworkflows/tmp/combine_ibs_lowcov_scatter.nf \
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
params.combine_pairs = null
params.combine_coverage_map = '/data/home/tusr1/01projects/vmap4/00data/05taxaBamMap/final.ALL.taxaBamMap.txt'
params.combine_pair_type = 'within_watkins'
params.combine_metric = 'hom_ibs'
params.combine_cutoff = 0.99
params.combine_coverage_floor = 3.0
params.combine_highlight_union_n = 48
params.combine_highlight_pairs = null
params.combine_highlight_min_ibs = null
params.combine_highlight_label = 'selected pairs'
params.combine_label = null
params.combine_y_min = null
params.combine_output_prefix = 'all.sample.combine_ibs_ipk_hzau_watkins.watkins_lowcov_vs_homsnp'

process combine_ibs_lowcov_scatter {
    tag "combine IBS low-coverage scatter: ${params.combine_pair_type}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 2
    memory '8.GB'
    errorStrategy 'terminate'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple path(pairs), path(coverage_map), path(miss_a), path(miss_b), path(miss_d), path(highlight_pairs, stageAs: 'highlight_pairs.tsv')

    output:
    path("*.png"), emit: plots
    path("${params.combine_output_prefix}.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import logging
    import sys
    sys.stdout = open("${params.combine_output_prefix}.log", "w", buffering=1)
    sys.stderr = sys.stdout
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.sample.combine_ibs_validate import (
        plot_combine_ibs_low_coverage,
    )

    plot_combine_ibs_low_coverage(
        "${pairs}",
        "${params.combine_output_prefix}",
        pair_type="${params.combine_pair_type}",
        metric="${params.combine_metric}",
        cutoff=float("${params.combine_cutoff}"),
        coverage_floor=float("${params.combine_coverage_floor}"),
        coverage_map_path="${coverage_map}",
        highlight_union_n=${params.combine_highlight_union_n ? params.combine_highlight_union_n : 'None'},
        highlight_pairs_path=${highlight_pairs ? "'" + highlight_pairs + "'" : "None"},
        highlight_min_ibs=${params.combine_highlight_min_ibs ? params.combine_highlight_min_ibs : "None"},
        highlight_label="${params.combine_highlight_label}",
        missing_rate_paths=["${miss_a}", "${miss_b}", "${miss_d}"],
        label=${params.combine_label ? '"' + params.combine_label + '"' : 'None'},
        y_min=${params.combine_y_min ? params.combine_y_min : 'None'},
    )
    """
}

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'combine_ibs_lowcov_scatter: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir || !params.src_dir) {
        error 'combine_ibs_lowcov_scatter: --user_dir and --src_dir are required'
    }
    def info = "${params.output_dir}/${params.job}/stats/${params.mod}/info"
    def pairs = params.combine_pairs
        ? file(params.combine_pairs, checkIfExists: true)
        : file("${info}/all.sample.combine_ibs_ipk_hzau_watkins.homsnp.pairs.tsv", checkIfExists: true)
    def coverage_map = file(params.combine_coverage_map, checkIfExists: true)
    def miss_a = file("${info}/A.sample.miss.rate.info.tsv", checkIfExists: true)
    def miss_b = file("${info}/B.sample.miss.rate.info.tsv", checkIfExists: true)
    def miss_d = file("${info}/D.sample.miss.rate.info.tsv", checkIfExists: true)
    def highlight = params.combine_highlight_pairs
        ? file(params.combine_highlight_pairs, checkIfExists: true)
        : []
    combine_ibs_lowcov_scatter(
        Channel.value([pairs, coverage_map, miss_a, miss_b, miss_d, highlight])
    )
}
