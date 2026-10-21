nextflow.enable.dsl=2

process singer_arg_preview {
    tag 'SINGER ARG preview'
    cpus 1
    memory '4 GB'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: '*.info.tsv'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: '*.th.tsv'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: '*.png'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: '*.log'
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path windows_tsv
    val singer_dir

    output:
    path '*.info.tsv', emit: info
    path '*.th.tsv', emit: th
    path '*.png', emit: plots
    path '*.log', emit: log

    script:
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("all.variant.singer_arg.log", "w")
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.singer_arg import plot_singer_arg_preview

    print("SINGER ARG preview: windows=${windows_tsv} singer_dir=${singer_dir}")
    written = plot_singer_arg_preview(
        "${windows_tsv}",
        "${singer_dir}",
        "all.variant.singer_arg",
    )
    print(f"n_ok={written['n_ok']} n_fail={written['n_fail']}")
    """
}

process singer_arg_relatedness {
    tag 'SINGER ARG relatedness'
    cpus 1
    memory '16 GB'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: '*.info.tsv'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/thresholds", mode: 'copy', pattern: '*.th.tsv'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: '*.png'
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: '*.log'
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path windows_tsv
    path keep_file
    val singer_dir

    output:
    path '*.info.tsv', emit: info
    path '*.th.tsv', emit: th
    path '*.png', emit: plots
    path '*.log', emit: log

    script:
    def n_pos = params.singer_relate_n_positions ?: 11
    def top_k = params.singer_relate_top_k ?: 5
    def top_pairs = params.singer_relate_top_pairs ?: 200
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("all.variant.singer_arg.relatedness.log", "w")
    sys.stderr = sys.stdout
    if "${params.src_dir}" and "${params.src_dir}" != "null":
        sys.path.insert(0, "${params.src_dir}/python")

    from genetics.genomics.variant.singer_relatedness import plot_singer_arg_relatedness

    print("SINGER ARG relatedness: windows=${windows_tsv} keep=${keep_file} singer_dir=${singer_dir}")
    written = plot_singer_arg_relatedness(
        "${windows_tsv}",
        "${singer_dir}",
        "${keep_file}",
        "all.variant.singer_arg",
        n_positions=${n_pos},
        top_k=${top_k},
        top_pairs=${top_pairs},
    )
    print(f"n_samples={written['n_samples']} n_windows={written['n_windows']}")
    """
}
