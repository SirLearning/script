nextflow.enable.dsl=2

def aneuploidyStatsPublishBase() {
    "${params.output_dir}/${params.job}/stats/${params.mod}"
}

process plot_confirmed_aneuploidy_chromosome {
    tag "confirmed aneuploidy chromosome plot"
    publishDir "${aneuploidyStatsPublishBase()}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${aneuploidyStatsPublishBase()}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${aneuploidyStatsPublishBase()}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    val(breachr_input_dir)

    output:
    path("*.png"), emit: plots
    path("*.events.tsv"), emit: events
    path("*.counts.tsv"), emit: counts
    path("*.aneuploidy_plot.log"), emit: log

    script:
    def plot_id = params.aneuploidy_plot_id ?: params.job
    def plot_title = params.aneuploidy_plot_title ? "\"${params.aneuploidy_plot_title}\"" : 'None'
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("${plot_id}.confirmed_aneuploidy.aneuploidy_plot.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample.aneuploidy_plot import plot_confirmed_aneuploidy_chromosome_distribution

    outputs = plot_confirmed_aneuploidy_chromosome_distribution(
        "${breachr_input_dir}",
        ".",
        plot_id="${plot_id}",
        plot_title=${plot_title},
    )
    print("Wrote plot:", outputs["plot"])
    """
}
