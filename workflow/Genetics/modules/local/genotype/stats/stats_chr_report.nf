nextflow.enable.dsl=2

process report_plink_chr_variant_counts {
    tag "chr variant counts: ${params.mod}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    val process_root

    output:
    path("${params.mod}.chr_variant_counts.tsv"), emit: by_chr
    path("${params.mod}.chr_variant_counts.by_ref.tsv"), emit: by_ref
    path("${params.mod}.chr_variant_counts.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${params.mod}.chr_variant_counts.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.plink.results_io import summarize_chr_variant_counts

    print("Summarizing per-chromosome variant counts from ${process_root} ...")
    per_chr, by_ref = summarize_chr_variant_counts(
        "${process_root}",
        "${params.mod}",
    )
    n_rows = len(per_chr) - 1
    total = int(per_chr.loc[per_chr['plink_chr'] == 'total', 'n_variants'].iloc[0])
    print(f"Wrote {n_rows} PLINK chromosome rows; total variants={total}")
    print(by_ref.to_string(index=False))
    """
}

process plot_thin_common_chr_variant_compare {
    tag "thin vs common-thin chr variant plots"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/thin_common_compare/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/thin_common_compare/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/thin_common_compare/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(thin_process_dir), val(common_process_dir), path(thin_by_ref), path(common_by_ref), val(output_prefix)

    output:
    path("${output_prefix}.thin_common_chr_compare.info.tsv"), emit: info
    path("${output_prefix}.variant.genome_mb_density.info.tsv"), emit: density_info
    path("${output_prefix}.variant.genome_mb_density.line.png"), emit: line_plot
    path("${output_prefix}.variant.chr_distribution.line.png"), emit: chr_line_plot
    path("${output_prefix}.variant.common_fraction.bar.png"), emit: bar_plot
    path("${output_prefix}.thin_common_chr_compare.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${output_prefix}.thin_common_chr_compare.log", "w")
    sys.stderr = sys.stdout

    from infra.utils.graph import plot_bar_chart, plot_genome_binned_density_panels
    from infra.utils.io import save_df_to_tsv
    from genetics.genomics.plink.results_io import (
        build_genome_density_compare_table,
        compute_genome_bin_counts,
        merge_thin_common_chr_by_ref,
        plot_thin_common_by_ref_distribution_line,
    )

    BIN_SIZE_BP = 5_000_000

    print("Plotting test_thin vs test_common_thin variant comparison ...")
    merged = merge_thin_common_chr_by_ref("${thin_by_ref}", "${common_by_ref}")
    save_df_to_tsv(merged, "${output_prefix}.thin_common_chr_compare.info.tsv")

    plot_thin_common_by_ref_distribution_line(
        merged,
        filename="${output_prefix}.variant.chr_distribution.line.png",
        title='Variant count by chromosome (test_thin vs test_common_thin)',
        y_label='Variant count',
    )

    thin_bins = compute_genome_bin_counts("${thin_process_dir}", bin_size_bp=BIN_SIZE_BP)
    common_bins = compute_genome_bin_counts("${common_process_dir}", bin_size_bp=BIN_SIZE_BP)
    density = build_genome_density_compare_table(thin_bins, common_bins)
    save_df_to_tsv(density, "${output_prefix}.variant.genome_mb_density.info.tsv")

    plot_genome_binned_density_panels(
        density,
        x_col='bin_index',
        ref_name_col='ref_name',
        panel_specs=[
            {'y_col': 'test_thin', 'title': 'test_thin', 'color': '#4c72b0'},
            {'y_col': 'test_common_thin', 'title': 'test_common_thin', 'color': '#c44e52'},
        ],
        filename="${output_prefix}.variant.genome_mb_density.line.png",
        suptitle='Genome-wide variant density (5 Mb bins, vmap4 sequence order)',
        bin_size_mb=5,
    )

    plot_bar_chart(
        merged['ref_name'].tolist(),
        merged['common_fraction'].tolist(),
        title='Common-thin variant fraction of test_thin by chromosome',
        ylabel='Fraction (test_common_thin / test_thin)',
        filename="${output_prefix}.variant.common_fraction.bar.png",
        ylim=(0.0, 1.05),
        figure_size=(14, 5),
        rotate_xlabels=45,
    )

    print("Per-chromosome summary:")
    print(merged.to_string(index=False))
    print(f"Genome 5 Mb bins: {len(density)}")
    """
}

process plot_thin_common_chr_pi_compare {
    tag "thin vs common-thin chr pi plots"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/thin_common_pi_compare/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/thin_common_pi_compare/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/thin_common_pi_compare/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(thin_process_dir), val(common_process_dir), val(output_prefix)

    output:
    path("${output_prefix}.thin_common_pi_compare.info.tsv"), emit: info
    path("${output_prefix}.pi.genome_mb_density.info.tsv"), emit: density_info
    path("${output_prefix}.pi.genome_mb_density.line.png"), emit: line_plot
    path("${output_prefix}.pi.chr_distribution.line.png"), emit: chr_line_plot
    path("${output_prefix}.pi.common_fraction.bar.png"), emit: bar_plot
    path("${output_prefix}.thin_common_pi_compare.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${output_prefix}.thin_common_pi_compare.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.plink.nucleotide_diversity import plot_thin_common_pi_compare

    print("Computing SNP-only nucleotide diversity (pi) for test_thin vs test_common_thin ...")
    merged, _thin_bins, density = plot_thin_common_pi_compare(
        "${thin_process_dir}",
        "${common_process_dir}",
        "${output_prefix}",
        bin_size_bp=5_000_000,
    )
    print("Per-chromosome pi summary:")
    print(merged.to_string(index=False))
    print(f"Genome 5 Mb bins: {len(density)}")
    """
}
