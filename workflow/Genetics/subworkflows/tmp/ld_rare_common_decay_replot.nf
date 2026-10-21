#!/usr/bin/env nextflow
/*
 * Replot *.variant.ld_decay.rare_common.png from published
 * *.variant.ld_decay.rare_common.info.tsv (no vcor). Rare-common uses
 * alpha; rare-rare stays opaque so overlap stays visible.
 *
 * Launch from a vmap4 run folder; log in doc/progress/YYYY-MM-DD.md (#### NF replay).
 *   nextflow run .../subworkflows/tmp/ld_rare_common_decay_replot.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_wap_thin_refchr
 */

nextflow.enable.dsl=2

include { variant_ld_rare_common_decay_replot } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = null

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'ld_rare_common_decay_replot: --output_dir, --job, and --mod are required'
    }

    def info_dir = "${params.output_dir}/${params.job}/stats/${params.mod}/info"
    Channel
        .from(['A', 'B', 'D', 'Others'])
        .map { id ->
            tuple(
                id,
                file("${info_dir}/${id}.variant.ld_decay.rare_common.info.tsv"),
            )
        }
        .filter { id, info -> info.exists() }
        .set { ch_info }
    variant_ld_rare_common_decay_replot(ch_info)
}
