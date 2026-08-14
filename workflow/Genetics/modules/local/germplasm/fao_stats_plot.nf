nextflow.enable.dsl=2

def faoStatsPublishBase() {
    "${params.output_dir}/${params.job}/stats/${params.mod}"
}

process fao_stats_plot {
    tag 'fao_wheat_stats'
    cpus 1
    memory '1 GB'
    publishDir "${faoStatsPublishBase()}/plots", mode: 'copy', pattern: '*.png'
    publishDir "${faoStatsPublishBase()}/info", mode: 'copy', pattern: '*.tsv'
    publishDir "${faoStatsPublishBase()}/logs", mode: 'copy', pattern: '*.log'
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path production_csv
    path supply_csv

    output:
    path('wheat_production_area.png'), emit: production_plot
    path('wheat_cereal_percentage.png'), emit: cereal_plot
    path('wheat_total_percentage.png'), emit: total_plot
    path('wheat_production_area.tsv'), emit: production_tsv
    path('wheat_supply_percentage.tsv'), emit: supply_tsv
    path('fao_stats_plot.log'), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("fao_stats_plot.log", "w")
    sys.stderr = sys.stdout

    from genetics.germplasm.fao.fao_stats_plot import run_fao_stats_plots

    print("FAO stats plot: production=${production_csv}, supply=${supply_csv}")
    outputs = run_fao_stats_plots(
        production_csv="${production_csv}",
        supply_csv="${supply_csv}",
        output_dir=".",
    )
    for name, path in outputs.items():
        print(f"Wrote {name}: {path}")
    """
}
