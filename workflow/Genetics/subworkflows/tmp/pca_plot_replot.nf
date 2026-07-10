#!/usr/bin/env nextflow
/*
 * Replot PCA (+ PC1–6 variance bar) from existing process/{mod}/pca/*.pca.eigenvec|eigenval.
 *
 * Example:
 *   cd /data/home/tusr1/01projects/vmap4/08stats.genome/118run_pca_pc6_replot_test_thin
 *   nextflow run .../subworkflows/tmp/pca_plot_replot.nf -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 --job test_plink --mod test_thin
 */

nextflow.enable.dsl=2

include { listMergedSubgenomeTestPcaTuples } from '../../modules/local/infra/infra_plink_reuse.nf'
include { sample_plink1_pca_plot } from '../../modules/local/genotype/stats/stats_sample.nf'

params.job = 'test_plink'
params.mod = null
params.pca_dir = null

workflow {
    if (!params.mod) {
        error 'pca_plot_replot: --mod is required'
    }
    if (!(params.mod in ['test_thin', 'test_common_thin', 'test_common_only', 'test_common_inter'])) {
        error "pca_plot_replot: unsupported mod ${params.mod}"
    }
    if (!params.output_dir || !params.job) {
        error 'pca_plot_replot: --output_dir and --job are required'
    }

    def pca_dir = params.pca_dir ?: "${params.output_dir}/${params.job}/process/${params.mod}/pca"
    def group_file = file("${params.output_dir}/meta_data/sample_group.txt")
    def ch_pca = channel.from(listMergedSubgenomeTestPcaTuples(pca_dir))
    sample_plink1_pca_plot(group_file, ch_pca)
}
