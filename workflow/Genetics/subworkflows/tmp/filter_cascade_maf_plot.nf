#!/usr/bin/env nextflow
/*
 * Publish het → MAC≥5 → WAP MAC≥2 filter cascade counts + MAF rare/common plot.
 *
 * Prerequisite: merged {A,B,D,Others}_test.plink2 + variant/*.info.gcount under process/{mod}/;
 * F_panel het files under process/{filter_cascade_het_mod}/sample/*.info.het;
 * WAP keep list under process/{filter_cascade_wap_mod}/.
 *
 * Launch from a vmap4 run folder. Example:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/12stats.genome/37run_test_thin_filter_cascade_publish
 *   cd /data/home/tusr1/01projects/vmap4/12stats.genome/37run_test_thin_filter_cascade_publish
 *   nextflow run .../subworkflows/tmp/filter_cascade_maf_plot.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_thin
 */

nextflow.enable.dsl=2

include { variant_filter_cascade_maf_plot } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = null
params.filter_cascade_het_mod = 'test_common_only'
params.filter_cascade_wap_mod = 'test_wap'
params.filter_cascade_wap_keep = null
params.filter_cascade_het_alpha = '1e-6'

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'filter_cascade_maf_plot: --output_dir, --job, and --mod are required'
    }

    def wap_mod = params.filter_cascade_wap_mod ?: 'test_wap'
    def keep_name = params.filter_cascade_wap_keep ?: "${wap_mod}.keep"
    def wap_keep = file("${params.output_dir}/${params.job}/process/${wap_mod}/${keep_name}")

    variant_filter_cascade_maf_plot(Channel.of(wap_keep))
}
