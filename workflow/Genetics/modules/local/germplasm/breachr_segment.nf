nextflow.enable.dsl=2

def breachrPublishBase() {
    "${params.output_dir}/${params.job}/stats/${params.mod}/info/breachr_check"
}

process breachr_segment_sample {
    tag "${sample}"
    cpus 1
    memory '1 GB'
    publishDir "${breachrPublishBase()}/segments", mode: 'copy', pattern: "*.breachr_segments.tsv"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(sample), path(mapq_tsv)

    output:
    tuple val(sample), path("${sample}.breachr_segments.tsv"), emit: segments

    script:
    """
    #!/usr/bin/env python
    from breachr import call_mapq_depth_bin_tsv

    calls = call_mapq_depth_bin_tsv("${mapq_tsv}", sample="${sample}")
    calls.to_csv("${sample}.breachr_segments.tsv", sep="\\t", index=False)
    """
}

process breachr_reconcile_flagged {
    tag "reconcile"
    cpus 1
    memory '2 GB'
    publishDir "${breachrPublishBase()}", mode: 'copy'
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path segment_files
    path flagged_samples

    output:
    path "breachr_segments_all.tsv", emit: segments_all
    path "breachr_reconcile_all.tsv", emit: reconcile_all
    path "cnv_misclassified_as_aneuploidy.tsv", emit: misclassified
    path "breachr_reconcile_summary.tsv", emit: summary

    script:
    """
    #!/usr/bin/env python
    from pathlib import Path

    from breachr.reconcile import ReconcileConfig, write_reconcile_outputs

    segment_paths = sorted(Path(".").glob("*.breachr_segments.tsv"))
    write_reconcile_outputs(
        segment_paths=segment_paths,
        flagged_path=Path("${flagged_samples}"),
        out_segments_all=Path("breachr_segments_all.tsv"),
        out_misclassified=Path("cnv_misclassified_as_aneuploidy.tsv"),
        out_summary=Path("breachr_reconcile_summary.tsv"),
        config=ReconcileConfig(
            min_cnv_length_bp=${params.breachr_min_cnv_length_bp ?: 2000000},
            min_cnv_bins=${params.breachr_min_cnv_bins ?: 2},
        ),
    )
    """
}

process breachr_mapq_depth_cap_plot {
    tag "${sample}"
    cpus 1
    memory '1 GB'
    publishDir "${breachrPublishBase()}/plots/mapq_depth_cap", mode: 'copy', pattern: "*.mapq_depth_cap.genome.png"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(sample), path(mapq_tsv), path(segments_tsv), path(reconcile_all), path(flagged_samples)

    output:
    path("${sample}.mapq_depth_cap.genome.png"), emit: png

    script:
    """
    #!/usr/bin/env python
    from pathlib import Path

    from breachr.mapq_depth_cap_plot import plot_breachr_flagged_mapq_depth_cap

    plot_breachr_flagged_mapq_depth_cap(
        mapq_tsv_path=Path("${mapq_tsv}"),
        sample="${sample}",
        flagged_path=Path("${flagged_samples}"),
        reconcile_all_path=Path("${reconcile_all}"),
        segments_path=Path("${segments_tsv}"),
        output_png=Path("${sample}.mapq_depth_cap.genome.png"),
    )
    """
}

process breachr_reconcile_plots {
    tag "plots"
    cpus 1
    memory '1 GB'
    publishDir "${breachrPublishBase()}/plots", mode: 'copy'
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path reconcile_all
    path flagged_samples
    path segments_all

    output:
    path "*.png", emit: plots
    path "breachr_reconcile_sample_summary.tsv", emit: sample_summary

    script:
    """
    #!/usr/bin/env python
    from pathlib import Path

    from breachr.reconcile_plot import plot_breachr_reconcile

    plot_breachr_reconcile(
        reconcile_all_path=Path("${reconcile_all}"),
        flagged_path=Path("${flagged_samples}"),
        segments_all_path=Path("${segments_all}"),
        plot_dir=Path("."),
        sample_summary_path=Path("breachr_reconcile_sample_summary.tsv"),
    )
    """
}
