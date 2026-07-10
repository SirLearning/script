#!/usr/bin/env nextflow
/*
 * PLINK2 approx PCA + plot on merged *_test.plink2 (PLINK1 --pca returns all-zero
 * eigenvectors at ~7700 samples on this cohort; IBS/MDS handled separately).
 *
 * Example:
 *   cd /data/home/tusr1/01projects/vmap4/08stats.genome/111run_plink1_pca_replot_test_thin
 *   nextflow run .../subworkflows/tmp/plink1_pca_replot.nf -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 --job test_plink --mod test_thin
 */

nextflow.enable.dsl=2

include { listMergedSubgenomeTestPfileTuples } from '../../modules/local/infra/infra_plink_reuse.nf'
include { plink1_pca; sample_plink1_pca_plot } from '../../modules/local/genotype/stats/stats_sample.nf'

params.job = 'test_plink'
params.mod = null

workflow {
    if (!params.mod) {
        error 'plink1_pca_replot: --mod is required (e.g. test_thin)'
    }
    if (!(params.mod in ['test_thin', 'test_common_thin', 'test_common_only', 'test_common_inter'])) {
        error "plink1_pca_replot: unsupported mod ${params.mod}"
    }
    if (!params.output_dir || !params.job) {
        error 'plink1_pca_replot: --output_dir and --job are required'
    }

    def proc = params.process_dir ?: "${params.output_dir}/${params.job}/process/${params.mod}"
    def group_file = file("${params.output_dir}/meta_data/sample_group.txt")

    def ch_pfile = channel.from(listMergedSubgenomeTestPfileTuples(proc))
    def pca_out = plink1_pca(ch_pfile)
    sample_plink1_pca_plot(group_file, pca_out.pca)
}
