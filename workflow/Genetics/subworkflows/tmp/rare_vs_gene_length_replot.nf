#!/usr/bin/env nextflow
/*
 * Count rare_wap (WAP MAF<1%) sites per gene from published variant_sets
 * tables, plot vs Lulab CDS / gene span, and publish into stats/<mod>.
 *
 * Launch from a vmap4 run folder; log in doc/progress/YYYY-MM-DD.md (#### NF replay).
 *   nextflow run .../subworkflows/tmp/rare_vs_gene_length_replot.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_wap_exon
 */

nextflow.enable.dsl=2

params.job = 'test_plink'
params.mod = 'test_wap_exon'
params.rare_plot_prefix = 'all.gene'
params.variant_sets_dir = null

process plot_rare_vs_gene_length {
    tag "rare vs length ${params.mod}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    stageInMode 'symlink'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.rare_vs_gene_length.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: "*.th.tsv"

    input:
    path sets_a
    path sets_b
    path sets_d
    path gff

    output:
    path "${params.rare_plot_prefix}.rare*.png", emit: plots
    path "${params.rare_plot_prefix}.rare_vs_gene_length.tsv", emit: info
    path "${params.rare_plot_prefix}.rare_vs_gene_length.th.tsv", emit: thresholds

    script:
    """
    #!/usr/bin/env python
    import sys

    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.lof_ancestry import (
        count_flagged_per_gene,
        run_lof_vs_gene_length,
    )

    counts = count_flagged_per_gene(
        ["${sets_a}", "${sets_b}", "${sets_d}"],
        flag_col="rare_wap",
        count_col="n_rare",
    )
    run_lof_vs_gene_length(
        counts,
        "${gff}",
        "${params.rare_plot_prefix}",
        count_col="n_rare",
        y_noun="Rare (WAP MAF<1%)",
        file_stem="rare",
    )
    """
}

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'rare_vs_gene_length_replot: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir) {
        error 'rare_vs_gene_length_replot: --user_dir is required (Conda stats env).'
    }

    def info_dir = params.variant_sets_dir ?: "${params.output_dir}/${params.job}/stats/${params.mod}/info"
    def gff = file(params.test_wap_exon_gff ?: params.ref_v1_gff, checkIfExists: true)
    plot_rare_vs_gene_length(
        file("${info_dir}/A.variant_sets.tsv", checkIfExists: true),
        file("${info_dir}/B.variant_sets.tsv", checkIfExists: true),
        file("${info_dir}/D.variant_sets.tsv", checkIfExists: true),
        gff,
    )
}
