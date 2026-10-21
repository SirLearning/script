#!/usr/bin/env nextflow
/*
 * Publish LoF-vs-gene-length (and LoF ancestry class) plots into
 * stats/<mod>/plots from a LoF site table + Lulab GFF.
 *
 * Launch from a vmap4 run folder; log in doc/progress/YYYY-MM-DD.md (#### NF replay).
 *   nextflow run .../subworkflows/tmp/lof_vs_gene_length_replot.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_wap_exon \
 *     --lof_table .../test_wap_exon.lof_ancestry.tsv \
 *     --lof_ancestry_png .../test_wap_exon.lof_ancestry.class.png
 */

nextflow.enable.dsl=2

params.job = 'test_plink'
params.mod = 'test_wap_exon'
params.lof_table = null
params.lof_ancestry_png = null
params.lof_plot_prefix = 'all.gene'

process plot_lof_vs_gene_length {
    tag "lof vs length ${params.mod}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.lof_vs_gene_length.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"

    input:
    path lof_table
    path gff
    path ancestry_png

    output:
    path "${params.lof_plot_prefix}.lof_vs*.png", emit: plots
    path "${params.lof_plot_prefix}.lof_ancestry.class.png", emit: ancestry_plot
    path "${params.lof_plot_prefix}.lof_vs_gene_length.tsv", emit: info
    path "${params.lof_plot_prefix}.lof_vs_gene_length.th.tsv", emit: thresholds

    script:
    """
    #!/usr/bin/env python
    import shutil
    import sys

    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.lof_ancestry import run_lof_vs_gene_length

    run_lof_vs_gene_length("${lof_table}", "${gff}", "${params.lof_plot_prefix}")
    shutil.copy("${ancestry_png}", "${params.lof_plot_prefix}.lof_ancestry.class.png")
    """
}

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'lof_vs_gene_length_replot: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir) {
        error 'lof_vs_gene_length_replot: --user_dir is required (Conda stats env).'
    }
    if (!params.lof_table) {
        error 'lof_vs_gene_length_replot: --lof_table is required'
    }
    if (!params.lof_ancestry_png) {
        error 'lof_vs_gene_length_replot: --lof_ancestry_png is required'
    }

    def gff = file(params.test_wap_exon_gff ?: params.ref_v1_gff, checkIfExists: true)
    plot_lof_vs_gene_length(
        file(params.lof_table, checkIfExists: true),
        gff,
        file(params.lof_ancestry_png, checkIfExists: true),
    )
}
