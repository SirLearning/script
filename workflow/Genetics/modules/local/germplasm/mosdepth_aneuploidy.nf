nextflow.enable.dsl=2

def aneuploidyPublishBase() {
    "${params.output_dir}/${params.job}/stats/${params.mod}"
}

process mosdepth_aneuploidy_sample {
    tag "${cohort}:${sample}"
    cpus 1
    memory '512 MB'
    publishDir "${aneuploidyPublishBase()}/logs", mode: 'copy', pattern: "*.aneuploidy.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(cohort), val(sample), path(summary)

    output:
    path("${sample}.${cohort}.chr_depth.tsv"), emit: chr
    path("${sample}.${cohort}.sample_summary.tsv"), emit: sample
    path("${sample}.${cohort}.aneuploidy.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("${sample}.${cohort}.aneuploidy.log", "w")
    sys.stderr = sys.stdout

    from genetics.germplasm.sample.mosdepth_aneuploidy import analyze_sample_summary

    print("mosdepth aneuploidy: cohort=${cohort}, sample=${sample}, summary=${summary}")
    chr_df, sample_df = analyze_sample_summary(
        "${summary}",
        sample="${sample}",
        cohort="${cohort}",
    )
    chr_df.to_csv("${sample}.${cohort}.chr_depth.tsv", sep="\\t", index=False)
    sample_df.to_csv("${sample}.${cohort}.sample_summary.tsv", sep="\\t", index=False)
    if not sample_df.empty:
        row = sample_df.iloc[0]
        print(
            f"baseline={row.get('baseline_depth', float('nan')):.4f}, "
            f"n_abnormal={row.get('n_abnormal_chr', 0)}, "
            f"flag={row.get('aneuploid_flag', False)}"
        )
    """
}

process mosdepth_aneuploidy_collect {
    tag "merge aneuploidy tables"
    publishDir "${aneuploidyPublishBase()}/info", mode: 'copy', pattern: "{sample_aneuploidy_summary,flagged_samples}.tsv"
    publishDir "${aneuploidyPublishBase()}/info/chr_depth", mode: 'copy', pattern: "chr_depth/**", saveAs: { filename -> new File(filename.toString()).name }
    publishDir "${aneuploidyPublishBase()}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path(chr_tsvs, stageAs: 'chr/*')
    path(sample_tsvs, stageAs: 'sample/*')

    output:
    path("sample_aneuploidy_summary.tsv"), emit: sample
    path("flagged_samples.tsv"), emit: flagged
    path("chr_depth/*.tsv"), emit: chr_depth
    path("mosdepth_aneuploidy_collect.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys
    from pathlib import Path

    sys.stdout = open("mosdepth_aneuploidy_collect.log", "w")
    sys.stderr = sys.stdout

    from genetics.germplasm.sample.mosdepth_aneuploidy import publish_aneuploidy_collect_outputs

    chr_paths = sorted(Path("chr").glob("*.tsv"))
    sample_paths = sorted(Path("sample").glob("*.tsv"))
    print("Collecting {} chr tables and {} sample summaries".format(len(chr_paths), len(sample_paths)))
    _, sample_all, flagged, chr_tables = publish_aneuploidy_collect_outputs(
        chr_paths,
        sample_paths,
        out_dir=".",
    )
    n_flag = len(flagged)
    n_all = len(sample_all)
    print("Merged samples={}, flagged={}, chr_tables={}".format(
        n_all, n_flag, len(chr_tables),
    ))
    """
}

process mosdepth_aneuploidy_replot {
    tag "replot aneuploidy figures"
    publishDir "${aneuploidyPublishBase()}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${aneuploidyPublishBase()}/info", mode: 'copy', pattern: "flagged_sample_chr_rel_depth.profile.tsv"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path(chr_tsvs, stageAs: 'chr_depth/*')
    val(group_file)
    val(flagged_samples_file)

    output:
    path("*.png"), emit: plots
    path("flagged_sample_chr_rel_depth.profile.tsv"), emit: profile_info, optional: true

    script:
    """
    #!/usr/bin/env python
    from genetics.germplasm.sample.mosdepth_aneuploidy import replot_rel_depth_distributions_from_chr_depth

    paths = replot_rel_depth_distributions_from_chr_depth(
        "chr_depth",
        ".",
        group_file="${group_file}",
        flagged_samples_path="${flagged_samples_file}",
    )
    print("Wrote {} plots".format(len(paths)))
    """
}
