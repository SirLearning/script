#!/usr/bin/env nextflow
/*
 * Within-bin MAF LD decay from published *.info.ld.vcor + *.info.afreq.
 * Each curve uses only pairs whose two sites fall in the same MAF interval:
 * (0, 0.01], (0.01, 0.05], (0.05, 0.1], (0.1, 0.2], (0.2, 0.5].
 * Does not overwrite *.variant.ld_decay.rare_common.png.
 * Also writes pooled A+B+D all.variant.ld_decay.maf_bin.png.
 *
 * Launch from a vmap4 run folder; log in doc/progress/YYYY-MM-DD.md (#### NF replay).
 *   nextflow run .../subworkflows/tmp/ld_maf_bin_decay.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_wap_thin_refchr
 *
 * Replot from published *.maf_bin.info.tsv (no vcor):
 *     --ld_maf_bin_plot_only true
 */

nextflow.enable.dsl=2

include {
    variant_ld_decay_maf_bin_plot
    variant_ld_maf_bin_decay_replot
    variant_ld_maf_bin_decay_all
} from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = null
params.process_dir = null
params.ld_maf_bin_plot_only = false

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'ld_maf_bin_decay: --output_dir, --job, and --mod are required'
    }

    def process_dir = params.process_dir
        ? params.process_dir
        : "${params.output_dir}/${params.job}/process/${params.mod}"
    def subgenomes = ['A', 'B', 'D']

    if (params.ld_maf_bin_plot_only) {
        def info_dir = "${params.output_dir}/${params.job}/stats/${params.mod}/info"
        Channel
            .from(subgenomes)
            .map { sg ->
                tuple(
                    sg,
                    file("${info_dir}/${sg}.variant.ld_decay.maf_bin.info.tsv"),
                )
            }
            .filter { id, info -> info.exists() }
            .set { ch_info }
        variant_ld_maf_bin_decay_replot(ch_info)
        variant_ld_maf_bin_decay_all(ch_info.map { _id, info -> info }.collect())
        return
    }

    Channel
        .from(subgenomes)
        .map { sg ->
            tuple(
                sg,
                "sub_${sg}",
                file("${process_dir}/variant/${sg}.info.ld.vcor"),
                file("${process_dir}/variant/${sg}.info.afreq"),
            )
        }
        .filter { id, chr, vcor, afreq -> vcor.exists() && vcor.size() > 1024 && afreq.exists() }
        .set { ch_ld }
    variant_ld_decay_maf_bin_plot(ch_ld)
    variant_ld_maf_bin_decay_all(
        variant_ld_decay_maf_bin_plot.out.info.map { _id, _chr, info -> info }.collect()
    )
}