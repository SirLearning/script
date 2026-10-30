#!/usr/bin/env nextflow
/*
 * Add histogram-bar distribution plots from the published combine_ibs pair
 * tables. Existing KDE plots are left untouched; this entry only writes the
 * *.hist.png counterparts and the per-ploidy histogram panels.
 *
 * Launch from a vmap4 run folder:
 *   nextflow run .../subworkflows/tmp/combine_ibs_hist_replot.nf \
 *     -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_common_inter
 *
 * Optional: --hist_bins N overrides the bar count of every histogram
 * (defaults kept in the Python replot helper: 80 AABBDD / 60 by-ploidy).
 * --ibs_thresholds A,B draws each IBS value as a dashed reference line on the
 * hom-SNP panels only (default 0.99,0.995; pass null to omit them). The merge
 * criterion is defined on hom_ibs, so PLINK IBS panels stay unmarked.
 */

nextflow.enable.dsl=2

params.job = 'test_plink'
params.mod = 'test_common_inter'
params.random_pairs = null
params.hist_bins = null
params.ibs_thresholds = '0.99,0.995'

process combine_ibs_hist_replot {
    tag "combine IBS histogram: ${prefix}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    cpus 2
    memory '16.GB'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(prefix), val(metric), path(pairs), path(random_pairs)

    output:
    path("*.png"), emit: plots
    path("${prefix}.hist.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${prefix}.hist.log", "w", buffering=1)
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.sample.combine_ibs_validate import replot_combine_ibs_hist

    replot_combine_ibs_hist(
        "${pairs}",
        "${random_pairs}",
        "${prefix}",
        metric="${metric}",
        bins=${params.hist_bins ? params.hist_bins : 'None'},
        ibs_threshold=${metric == 'hom_ibs' && params.ibs_thresholds ? '[' + params.ibs_thresholds.toString().split(',').collect { it.trim() }.join(',') + ']' : 'None'},
    )
    """
}

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'combine_ibs_hist_replot: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir || !params.src_dir) {
        error 'combine_ibs_hist_replot: --user_dir and --src_dir are required'
    }
    def info = "${params.output_dir}/${params.job}/stats/${params.mod}/info"
    def random_pairs = params.random_pairs
        ? file(params.random_pairs, checkIfExists: true)
        : file("${info}/all.sample.identity.ploidy_group_ibs.pairs.tsv", checkIfExists: true)
    def specs = [
        ['all.sample.combine_ibs', 'plink_ibs'],
        ['all.sample.combine_ibs_ipk_hzau_watkins', 'plink_ibs'],
        ['all.sample.combine_ibs.homsnp', 'hom_ibs'],
        ['all.sample.combine_ibs_ipk_hzau_watkins.homsnp', 'hom_ibs'],
    ].collect { prefix, metric ->
        tuple(prefix, metric, file("${info}/${prefix}.pairs.tsv", checkIfExists: true), random_pairs)
    }
    combine_ibs_hist_replot(Channel.fromList(specs))
}
