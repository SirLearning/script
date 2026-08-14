nextflow.enable.dsl=2

def mapqPublishBase() {
    "${params.output_dir}/${params.job}/stats/${params.mod}/info/mapq_check"
}

process bam_mapq_depth_bin_sample {
    tag "${sample}"
    cpus 2
    memory '2 GB'
    publishDir "${mapqPublishBase()}", mode: 'copy', pattern: "*.mapq_depth_bin.tsv"
    publishDir "${mapqPublishBase()}", mode: 'copy', pattern: "*.mapq_depth.log"
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(sample), path(bam)

    output:
    path("${sample}.mapq_depth_bin.tsv"), emit: tsv
    path("${sample}.mapq_depth.log"), emit: log

    script:
    def make_plots = params.mapq_plot ? 'True' : 'False'
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("${sample}.mapq_depth.log", "w")
    sys.stderr = sys.stdout

    from genetics.germplasm.sample.bam_mapq_depth_bin import run_sample_mapq_depth_bins

    print("MAPQ-depth bin scan: sample=${sample}, bam=${bam}, plot=${make_plots}")
    tsv_path, plot_paths = run_sample_mapq_depth_bins(
        sample="${sample}",
        bam_path="${bam}",
        output_tsv="${sample}.mapq_depth_bin.tsv",
        plot_dir=".",
        bin_size_bp=${params.mapq_bin_size_bp ?: 1000000},
        samtools="${params.mapq_samtools ?: 'samtools'}",
        threads=${task.cpus},
        make_plots=${make_plots},
    )
    print(f"Wrote table {tsv_path} and {len(plot_paths)} plots")
    """
}
