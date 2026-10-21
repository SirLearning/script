#!/usr/bin/env nextflow
/*
 * Overall P+ vs q from published *.maf_scatter.info.tsv, restricted to
 * pairs with distance ≤ params.pplus_vs_q_max_kb (default 1000 kb = 1 Mb).
 * Does not recompute dosages or redraw the default 5 Mb P+ vs q figure.
 *
 * Launch from a vmap4 run folder; log in doc/progress/YYYY-MM-DD.md (#### NF replay).
 *   nextflow run .../subworkflows/tmp/ld_pplus_vs_q_overall_replot.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_wap_thin_refchr
 */

nextflow.enable.dsl=2

include { variant_ld_pplus_vs_q_overall_replot } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = null
params.pplus_vs_q_max_kb = 1000
params.pplus_vs_q_suffix = '1mb'

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'ld_pplus_vs_q_overall_replot: --output_dir, --job, and --mod are required'
    }

    def info_dir = "${params.output_dir}/${params.job}/stats/${params.mod}/info"
    Channel
        .from(['A', 'B', 'D', 'Others'])
        .map { id ->
            tuple(
                id,
                file("${info_dir}/${id}.variant.ld_decay.rare_common.maf_scatter.info.tsv"),
            )
        }
        .filter { id, info -> info.exists() }
        .set { ch_info }
    variant_ld_pplus_vs_q_overall_replot(ch_info)
}
