#!/usr/bin/env nextflow
/*
 * Replot the VEP atlas Circos from published window counts, adding:
 * outermost IWGSC cytobands (arm / pericentromere / centromere; no cytoband legend),
 * LD blocks from --vep_atlas_ld_blocks (default test_wap_thin_refchr 1 Mb coverage),
 * then gene counts (yellow), predicted LoF (purple), common (red), rare (blue).
 * Legend sits in the centre hole.
 *
 * Launch from a vmap4 run folder; log in doc/progress/YYYY-MM-DD.md (#### NF replay).
 *   nextflow run .../subworkflows/tmp/vep_atlas_circos_replot.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_wap_refchr \
 *     --vep_atlas_ld_blocks .../test_wap_thin_refchr/info/all.variant.circos_ld_map.blocks_1mb.info.tsv \
 *     --vep_atlas_gene_gff .../IWGSC_v1.1_HC_20170706.gff3
 */

nextflow.enable.dsl=2

include { variant_vep_atlas_circos_replot } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = 'test_wap_refchr'
params.vep_atlas_bin_bp = 10000000
params.vep_atlas_peri_flank_bp = 50000000
params.vep_atlas_ld_blocks = null
params.vep_atlas_gene_gff = null

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'vep_atlas_circos_replot: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir) {
        error 'vep_atlas_circos_replot: --user_dir is required (conda stats)'
    }
    if (!params.src_dir) {
        error 'vep_atlas_circos_replot: --src_dir is required'
    }
    if (!params.vep_atlas_ld_blocks) {
        error 'vep_atlas_circos_replot: --vep_atlas_ld_blocks is required'
    }
    if (!params.vep_atlas_gene_gff) {
        error 'vep_atlas_circos_replot: --vep_atlas_gene_gff is required'
    }

    def info_dir = "${params.output_dir}/${params.job}/stats/${params.mod}/info"
    def windows = file("${info_dir}/all.variant.vep_atlas.windows.tsv", checkIfExists: true)
    def ld_blocks = file(params.vep_atlas_ld_blocks, checkIfExists: true)
    def gene_gff = file(params.vep_atlas_gene_gff, checkIfExists: true)

    log.info "${params.c_green}vep_atlas Circos replot ${params.mod}${params.c_reset}"
    variant_vep_atlas_circos_replot(tuple(windows, ld_blocks, gene_gff))
}
