#!/usr/bin/env nextflow
/*
 * PLINK1 genotype PCA + IBS-matrix MDS for merged test plink bfiles.
 *
 * Prerequisite: merged {A,B,D,Others}_test.plink.{bed,bim,fam} under process/{mod}/.
 *
 * Launch from a vmap4 run folder. Example:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/08stats.genome/107run_ibs_pca_replot/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/08stats.genome/107run_ibs_pca_replot
 *   nextflow run .../subworkflows/tmp/ibs_pca_replot.nf -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 --job test_plink --mod test_thin
 */

nextflow.enable.dsl=2

include { listMergedSubgenomeTestBfileTuples; listMergedSubgenomeTestPfileTuples } from '../../modules/local/infra/infra_plink_reuse.nf'
include {
    plink1_ibs_square
    plink1_pca
    sample_ibs_mds_plot
    sample_plink1_pca_plot
} from '../../modules/local/genotype/stats/stats_sample.nf'

params.job = 'test_plink'
params.mod = null

workflow {
    if (!params.mod) {
        error 'ibs_pca_replot: --mod is required (e.g. test_thin, test_common_thin, test_common_only)'
    }
    if (!(params.mod in ['test_thin', 'test_common_thin', 'test_common_only', 'test_common_inter'])) {
        error "ibs_pca_replot: unsupported mod ${params.mod}"
    }
    if (!params.output_dir || !params.job) {
        error 'ibs_pca_replot: --output_dir and --job are required'
    }

    def proc = params.process_dir ?: "${params.output_dir}/${params.job}/process/${params.mod}"
    def group_file = file("${params.output_dir}/meta_data/sample_group.txt")

    def ch_bfile = channel.from(listMergedSubgenomeTestBfileTuples(proc))
    def ch_pfile = channel.from(listMergedSubgenomeTestPfileTuples(proc))
    def ibs_out = plink1_ibs_square(ch_bfile)
    def pca_out = plink1_pca(ch_pfile)
    sample_ibs_mds_plot(group_file, ibs_out.ibs)
    sample_plink1_pca_plot(group_file, pca_out.pca)
}
