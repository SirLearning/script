#!/usr/bin/env nextflow
/*
 * Replot LD MAF-vs-MAF scatters from published *.maf_scatter.info.tsv with
 * r_max^2 = 0.8 / 0.6 / 0.2 contour pairs (six curves), rewrite the pair table
 * with sign-aware |D'| from pairwise pfile dosages (not r2 / r_max^2), and
 * plot |D'| decay, a two-colour |D'| vs distance scatter, and fine-bin
 * |D'| density lines (frequency polygons; log panel is log10(density + 1))
 * and signed-D' class proportions / |D'|=1 stack / interior density.
 *
 * Launch from a vmap4 run folder; log in doc/progress/YYYY-MM-DD.md (#### NF replay).
 *   nextflow run .../subworkflows/tmp/ld_maf_scatter_r2max_replot.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_thin
 */

nextflow.enable.dsl=2

include { variant_ld_maf_scatter_r2max_replot } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = null

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'ld_maf_scatter_r2max_replot: --output_dir, --job, and --mod are required'
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
    variant_ld_maf_scatter_r2max_replot(ch_info)
}
