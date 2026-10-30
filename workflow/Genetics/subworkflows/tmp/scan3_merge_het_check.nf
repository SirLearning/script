#!/usr/bin/env nextflow
/*
 * Merged-germplasm heterozygosity check for one FastCall3 scan3 chromosome.
 *
 * scan3 sums the raw AD of a germplasm's source taxa, so a wrongly merged group
 * would show up as excess heterozygosity. This compares the PLINK2 --het table
 * of the chromosome VCF against the sub-genome taxaMergeMap used to build it.
 * Merged samples also genotype more sites (several sources summed), so the
 * comparison is repeated inside call-rate bins: raw distributions alone are
 * dominated by low-coverage singletons carrying error-driven heterozygotes.
 *
 * Publishes under {output_dir}/{job}/stats/{mod}/{info,plots,logs}:
 *   <prefix>.sample.het.info.tsv
 *   <prefix>.het.summary.info.tsv
 *   <prefix>.het.call_rate_bins.info.tsv
 *   <prefix>.sample.het.distribution.png   merged vs unmerged histogram
 *   <prefix>.sample.het.vs_call_rate.png   median heterozygosity per call-rate bin
 *   <prefix>.merge_het_check.log
 *
 * Launch from a vmap4 run folder:
 *   nextflow run .../subworkflows/tmp/scan3_merge_het_check.nf \
 *     -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_merge_chr6 \
 *     --het <chr006_het.het> --merge_map <final.D.taxaMergeMap.txt>
 */

nextflow.enable.dsl=2

params.job = 'test_plink'
params.mod = 'test_merge_chr6'
params.het = null
params.merge_map = null
params.output_prefix = 'chr006'
params.title_prefix = ''
params.highlight_prefix = null
params.exclude_prefix = null

process scan3_merge_het_check {
    tag "merged heterozygosity check"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 2
    memory '8.GB'
    errorStrategy 'terminate'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.info.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple path(het), path(merge_map)

    output:
    path("*.info.tsv"), emit: info
    path("*.png"), emit: plots
    path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import logging
    import sys
    sys.stdout = open("${params.output_prefix}.merge_het_check.log", "w", buffering=1)
    sys.stderr = sys.stdout
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.sample.merge_het_check import run_merge_het_check

    run_merge_het_check(
        "${het}",
        "${merge_map}",
        "${params.output_prefix}",
        info_dir=".",
        plot_dir=".",
        title_prefix="${params.title_prefix}",
        highlight_prefix=${params.highlight_prefix ? '"' + params.highlight_prefix + '"' : 'None'},
        exclude_prefix=${params.exclude_prefix ? '"' + params.exclude_prefix + '"' : 'None'},
    )
    """
}

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'scan3_merge_het_check: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir || !params.src_dir) {
        error 'scan3_merge_het_check: --user_dir and --src_dir are required'
    }
    if (!params.het || !params.merge_map) {
        error 'scan3_merge_het_check: --het and --merge_map are required'
    }
    scan3_merge_het_check(
        Channel.value([
            file(params.het, checkIfExists: true),
            file(params.merge_map, checkIfExists: true),
        ])
    )
}
