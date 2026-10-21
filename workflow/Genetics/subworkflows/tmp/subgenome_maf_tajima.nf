#!/usr/bin/env nextflow
/*
 * Stacked A/B/D MAF histograms (linear + log y) from published *.info.afreq,
 * and stacked Tajima's D histograms from PLINK2 VCF 4.2 export + vcftools --TajimaD.
 * Intermediate VCFs stay in the work dir (not published).
 *
 * Launch from a vmap4 run folder; log in doc/progress/YYYY-MM-DD.md (#### NF replay).
 *   nextflow run .../subworkflows/tmp/subgenome_maf_tajima.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_wap_thin_refchr
 *
 * Replot Tajima's D from published *.Tajima.D (no VCF export):
 *     --tajima_d_plot_only true
 */

nextflow.enable.dsl=2

include { plink2_vcftools_tajima_d } from '../../modules/local/genotype/processor/processor_plink2.nf'
include {
    variant_subgenome_maf_plot
    variant_subgenome_tajima_plot
} from '../../modules/local/genotype/stats/stats_variant.nf'
include {
    listMergedSubgenomeTestPfileTuples
} from '../../modules/local/infra/infra_plink_reuse.nf'

params.job = 'test_plink'
params.mod = null
params.process_dir = null
params.tajima_d_window_bp = 100000
params.tajima_d_plot_only = false

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'subgenome_maf_tajima: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir) {
        error 'subgenome_maf_tajima: --user_dir is required (conda stats)'
    }

    def process_dir = params.process_dir
        ? params.process_dir
        : "${params.output_dir}/${params.job}/process/${params.mod}"
    def window_bp = (params.tajima_d_window_bp ?: 100000) as int
    def window_tag = "${window_bp.intdiv(1000)}kb"
    def subgenomes = ['A', 'B', 'D']

    Channel
        .from(subgenomes)
        .map { sg ->
            def afreq = file("${process_dir}/variant/${sg}.info.afreq")
            if (!afreq.exists()) {
                error "subgenome_maf_tajima: missing ${afreq}"
            }
            afreq
        }
        .collect()
        .set { ch_afreq }
    variant_subgenome_maf_plot(ch_afreq)

    if (params.tajima_d_plot_only) {
        def info_dir = "${params.output_dir}/${params.job}/stats/${params.mod}/info"
        Channel
            .from(subgenomes)
            .map { sg ->
                def tajima = file("${info_dir}/${sg}.variant.tajima_d.${window_tag}.Tajima.D")
                if (!tajima.exists()) {
                    error "subgenome_maf_tajima plot-only: missing ${tajima}"
                }
                tajima
            }
            .collect()
            .set { ch_tajima }
        variant_subgenome_tajima_plot(ch_tajima)
        return
    }

    Channel
        .from(listMergedSubgenomeTestPfileTuples(process_dir))
        .filter { id, chr, prefix, pgen, psam, pvar -> id in subgenomes }
        .set { ch_pfile }
    plink2_vcftools_tajima_d(ch_pfile)
    variant_subgenome_tajima_plot(plink2_vcftools_tajima_d.out.tajima.collect())
}