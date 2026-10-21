#!/usr/bin/env nextflow
/*
 * Gene-level Manhattan / QQ from published *.burden.tsv (no genotype export).
 *
 * Launch from a vmap4 run folder; log in doc/progress/YYYY-MM-DD.md (#### NF replay).
 *   nextflow run .../subworkflows/tmp/burden_replot.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_wap_refchr
 */

nextflow.enable.dsl=2

include { plot_wap_burden } from '../../modules/local/static/gwas/gwas_plot.nf'

params.job = 'test_plink'
params.mod = null
params.burden_replot_dir = null

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'burden_replot: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir) {
        error 'burden_replot: --user_dir is required (Conda stats env).'
    }

    def assoc_dir = params.burden_replot_dir
        ?: "${params.output_dir}/${params.job}/process/${params.mod}/assoc"
    def gff = file(params.assoc_wap_burden_gff ?: params.ref_v1_gff_hc, checkIfExists: true)
    def tables = file("${assoc_dir}/*.burden.tsv")
    if (!tables) {
        error "burden_replot: no *.burden.tsv under ${assoc_dir}"
    }

    log.info "${params.c_green}burden replot ${tables.size()} tables from ${assoc_dir}${params.c_reset}"
    plot_wap_burden(
        channel.from(tables).map { table ->
            tuple(table, table.name.replaceAll(/\.burden\.tsv$/, ''))
        },
        gff,
    )
}
