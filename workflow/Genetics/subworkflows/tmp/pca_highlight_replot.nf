#!/usr/bin/env nextflow
/*
 * Replot PCA scatter with highlighted sample IDs from existing eigenvec/eigenval.
 *
 * Example:
 *   cd /data/home/tusr1/01projects/vmap4/10stats.genome/56run_pca_highlight_ab0043_test_common_inter
 *   nextflow run .../subworkflows/tmp/pca_highlight_replot.nf -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 --job test_plink --mod test_common_inter \
 *     --highlight_samples AB_0043 --pca_kinds pca,pca_ld
 */

nextflow.enable.dsl=2

include {
    listMergedSubgenomeTestPcaTuples
    listMergedSubgenomeTestPcaLdTuples
} from '../../modules/local/infra/infra_plink_reuse.nf'
include { sample_plink1_pca_highlight_plot } from '../../modules/local/genotype/stats/stats_sample.nf'

params.job = 'test_plink'
params.mod = null
params.highlight_samples = null
params.pca_kinds = 'pca,pca_ld'
params.pca_dir = null
params.pca_ld_dir = null

workflow {
    if (!params.mod) {
        error 'pca_highlight_replot: --mod is required'
    }
    if (!params.highlight_samples) {
        error 'pca_highlight_replot: --highlight_samples is required (comma-separated IDs)'
    }
    if (!params.output_dir || !params.job) {
        error 'pca_highlight_replot: --output_dir and --job are required'
    }

    def kinds = params.pca_kinds.toString().split(',').collect { it.trim() }.findAll { it }
    def group_file = file("${params.output_dir}/meta_data/sample_group.txt")
    def pca_dir = params.pca_dir ?: "${params.output_dir}/${params.job}/process/${params.mod}/pca"
    def pca_ld_dir = params.pca_ld_dir ?: "${params.output_dir}/${params.job}/process/${params.mod}/pca_ld"

    def ch_plot = channel.empty()
    if ('pca' in kinds) {
        ch_plot = ch_plot.mix(
            channel.from(listMergedSubgenomeTestPcaTuples(pca_dir))
                .map { id, chr, eigenvec, eigenval -> tuple(id, chr, eigenvec, eigenval, 'pca') }
        )
    }
    if ('pca_ld' in kinds) {
        ch_plot = ch_plot.mix(
            channel.from(listMergedSubgenomeTestPcaLdTuples(pca_ld_dir))
                .map { id, chr, eigenvec, eigenval -> tuple(id, chr, eigenvec, eigenval, 'pca_ld') }
        )
    }
    if (kinds.any { !(it in ['pca', 'pca_ld']) }) {
        error "pca_highlight_replot: unsupported pca_kinds ${kinds} (use pca and/or pca_ld)"
    }

    sample_plink1_pca_highlight_plot(
        ch_plot.map { id, chr, eigenvec, eigenval, kind -> group_file },
        ch_plot.map { id, chr, eigenvec, eigenval, kind -> tuple(id, chr, eigenvec, eigenval) },
        ch_plot.map { id, chr, eigenvec, eigenval, kind -> kind },
        ch_plot.map { id, chr, eigenvec, eigenval, kind -> params.highlight_samples },
    )
}
