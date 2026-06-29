nextflow.enable.dsl=2

def relDepthFmissPublishBase() {
    "${params.output_dir}/${params.job}/stats/${params.mod}"
}

process plot_rel_depth_vs_fmiss {
    tag "rel_depth vs F_MISS plots"
    publishDir "${relDepthFmissPublishBase()}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${relDepthFmissPublishBase()}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${relDepthFmissPublishBase()}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    val(chr_depth_dir)
    val(smiss_long_path)
    val(group_file)

    output:
    path("*.png"), emit: plots
    path("sample_chr_miss.depth_miss_merged.tsv"), emit: merged
    path("plot_rel_depth_vs_fmiss.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("plot_rel_depth_vs_fmiss.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.sample.sample_chr_miss import plot_rel_depth_vs_fmiss_from_tables

    paths = plot_rel_depth_vs_fmiss_from_tables(
        "${chr_depth_dir}",
        "${smiss_long_path}",
        ".",
        group_file="${group_file}",
        merged_out="sample_chr_miss.depth_miss_merged.tsv",
    )
    print("Wrote {} plots".format(len(paths)))
    """
}
