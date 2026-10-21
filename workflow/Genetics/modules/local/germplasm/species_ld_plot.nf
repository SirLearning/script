nextflow.enable.dsl=2

def metaStatsPublishBase() {
    "${params.output_dir}/${params.job}/stats/${params.mod}"
}

process species_ld_plot {
    tag 'species_ld_vs_genome'
    cpus 1
    memory '1 GB'
    publishDir "${metaStatsPublishBase()}/plots", mode: 'copy', pattern: '*.png'
    publishDir "${metaStatsPublishBase()}/info", mode: 'copy', pattern: '*.tsv'
    publishDir "${metaStatsPublishBase()}/logs", mode: 'copy', pattern: '*.log'
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path species_ld_csv

    output:
    path('species_ld_vs_genome.png'), emit: scatter_plot
    path('species_ld_vs_genome.labeled.png'), emit: scatter_plot_labeled
    path('species_ld_vs_genome.linear_x.png'), emit: scatter_plot_linear_x
    path('species_ld_vs_genome.linear_x.labeled.png'), emit: scatter_plot_linear_x_labeled
    path('species_ld_vs_genome.tsv'), emit: table_tsv
    path('species_ld_plot.log'), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("species_ld_plot.log", "w")
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.germplasm.meta.species_ld_plot import run_species_ld_plot

    print("Species LD plot: csv=${species_ld_csv}")
    outputs = run_species_ld_plot(
        csv_path="${species_ld_csv}",
        output_dir=".",
    )
    for name, path in outputs.items():
        print(f"Wrote {name}: {path}")
    """
}
