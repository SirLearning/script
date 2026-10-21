#!/usr/bin/env nextflow
/*
 * Fine-bin all-pairs mean r2 and the per-chromosome distance where that mean
 * first falls to --ld_decay_r2_threshold (default 0.5). Reads existing
 * variant/{A,B,D}.info.ld.vcor. Tagged outputs (ld_decay.r2_0p5.1kb) do not
 * overwrite the default 50 kb ld_decay tables.
 *
 *   nextflow run .../subworkflows/tmp/ld_decay_r2_threshold.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_wap_common_refchr \
 *     --ld_decay_r2_threshold 0.5 --ld_decay_fine_bin_kb 1 --ld_decay_fine_max_kb 100
 */

nextflow.enable.dsl=2

include {
    variant_ld_decay_r2_threshold_plot
} from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = null
params.process_dir = null

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'ld_decay_r2_threshold: --output_dir, --job, and --mod are required'
    }
    def process_dir = params.process_dir
        ? params.process_dir
        : "${params.output_dir}/${params.job}/process/${params.mod}"

    Channel
        .from(['A', 'B', 'D'])
        .map { sg ->
            def vcor = file("${process_dir}/variant/${sg}.info.ld.vcor")
            if (!vcor.exists()) {
                error "ld_decay_r2_threshold: missing ${vcor}"
            }
            tuple(sg, "sub_${sg}", vcor)
        }
        .set { ch_ld }
    variant_ld_decay_r2_threshold_plot(ch_ld)
}
