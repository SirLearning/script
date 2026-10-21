#!/usr/bin/env nextflow
/*
 * All-pairs unphased LD (5 Mb window) plus LD-decay / cross-chr plots from merged
 * *_test.plink2. Use when a test mod skipped TEST_PLINK_STATS (e.g. test_wap_common_refchr).
 *
 *   nextflow run .../subworkflows/tmp/ld_decay_from_pfile.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_wap_common_refchr \
 *     --process_dir /data1/.../test_plink/process/test_wap_common_refchr
 *
 * Replot from existing `variant/*.info.ld.vcor` (no PLINK2):
 *     --ld_decay_plot_only true
 */

nextflow.enable.dsl=2

include { calc_plink_ld_unphased; calc_plink_ld_crosschr_random } from '../../modules/local/genotype/processor/processor_plink2.nf'
include { variant_ld_decay_plot; variant_ld_crosschr_plot } from '../../modules/local/genotype/stats/stats_variant.nf'
include {
    hasMergedSubgenomeTestPfiles
    listMergedSubgenomeTestPfileTuples
} from '../../modules/local/infra/infra_plink_reuse.nf'

params.job = 'test_plink'
params.mod = null
params.process_dir = null
params.ld_decay_plot_only = false

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'ld_decay_from_pfile: --output_dir, --job, and --mod are required'
    }
    def process_dir = params.process_dir
        ? params.process_dir
        : "${params.output_dir}/${params.job}/process/${params.mod}"

    if (params.ld_decay_plot_only) {
        def subgenomes = ['A', 'B', 'D', 'Others']
        Channel
            .from(subgenomes)
            .map { sg ->
                def vcor = file("${process_dir}/variant/${sg}.info.ld.vcor")
                if (!vcor.exists()) {
                    error "ld_decay_from_pfile plot-only: missing ${vcor}"
                }
                tuple(sg, "sub_${sg}", vcor)
            }
            .set { ch_ld }
        variant_ld_decay_plot(ch_ld)
        return
    }

    if (!hasMergedSubgenomeTestPfiles(process_dir)) {
        error "ld_decay_from_pfile: missing merged *_test.plink2 under ${process_dir}"
    }

    Channel
        .from(listMergedSubgenomeTestPfileTuples(process_dir))
        .set { ch_pfile }

    calc_plink_ld_unphased(ch_pfile)
    calc_plink_ld_crosschr_random(ch_pfile)
    variant_ld_decay_plot(calc_plink_ld_unphased.out.ld)
    variant_ld_crosschr_plot(calc_plink_ld_crosschr_random.out.ld)
}
