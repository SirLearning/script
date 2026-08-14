nextflow.enable.dsl=2

def metaStatsPublishBase() {
    "${params.output_dir}/${params.job}/stats/${params.mod}"
}

process ld_block_size_plot {
    tag 'wu2025_ld_block'
    cpus 1
    memory '1 GB'
    publishDir "${metaStatsPublishBase()}/plots", mode: 'copy', pattern: '*.png'
    publishDir "${metaStatsPublishBase()}/info", mode: 'copy', pattern: '*.tsv'
    publishDir "${metaStatsPublishBase()}/logs", mode: 'copy', pattern: '*.log'
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path ld_block_csv

    output:
    path('wu2025_ld_block_size_distribution.png'), emit: distribution_plot
    path('wu2025_ld_block_size_distribution_logy.png'), emit: distribution_logy_plot
    path('wu2025_ld_block_size_vs_n_genes.png'), emit: size_vs_genes_plot
    path('wu2025_ld_block_sizes.tsv'), emit: blocks_tsv
    path('wu2025_ld_block_size_summary.tsv'), emit: summary_tsv
    path('ld_block_size_plot.log'), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("ld_block_size_plot.log", "w")
    sys.stderr = sys.stdout

    from genetics.germplasm.meta.ld_block_size_plot import run_ld_block_size_plot

    print("LD block size plot: csv=${ld_block_csv}")
    outputs = run_ld_block_size_plot(
        csv_path="${ld_block_csv}",
        output_dir=".",
    )
    for name, path in outputs.items():
        print(f"Wrote {name}: {path}")
    """
}
