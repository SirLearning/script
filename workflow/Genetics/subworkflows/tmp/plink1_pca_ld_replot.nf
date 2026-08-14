#!/usr/bin/env nextflow
/*
 * LD-pruned PLINK2 approx PCA + plot on merged *_test.plink2.
 *
 * Pruning: plink2 --indep-pairwise pca_ld_window pca_ld_step pca_ld_r2 (defaults 50 5 0.2),
 * then --extract *.prune.in before --pca approx.
 *
 * Example:
 *   cd /data/home/tusr1/01projects/vmap4/10stats.genome/51run_pca_ld_test_thin
 *   nextflow run .../subworkflows/tmp/plink1_pca_ld_replot.nf -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 --job test_plink --mod test_thin
 */

nextflow.enable.dsl=2

include { listMergedSubgenomeTestPfileTuples } from '../../modules/local/infra/infra_plink_reuse.nf'
include { plink1_pca_ld_prune; sample_plink1_pca_ld_plot } from '../../modules/local/genotype/stats/stats_sample.nf'

params.job = 'test_plink'
params.mod = null

workflow {
    if (!params.mod) {
        error 'plink1_pca_ld_replot: --mod is required (e.g. test_thin)'
    }
    if (!(params.mod in ['test_thin', 'test_common_thin', 'test_common_only', 'test_common_inter', 'test_rare_only'])) {
        error "plink1_pca_ld_replot: unsupported mod ${params.mod}"
    }
    if (!params.output_dir || !params.job) {
        error 'plink1_pca_ld_replot: --output_dir and --job are required'
    }

    def proc = params.process_dir ?: "${params.output_dir}/${params.job}/process/${params.mod}"
    def group_file = file("${params.output_dir}/meta_data/sample_group.txt")

    // Others is mostly chr0; LD prune leaves too few non-chr0 sites for --pca approx.
    def ch_pfile = channel.from(listMergedSubgenomeTestPfileTuples(proc))
        .filter { id, chr, prefix, pgen, psam, pvar -> id != 'Others' }
    def pca_ld_out = plink1_pca_ld_prune(ch_pfile)
    sample_plink1_pca_ld_plot(group_file, pca_ld_out.pca)
}
