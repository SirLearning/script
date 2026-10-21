nextflow.enable.dsl=2

def mapqPopdepInfoBase() {
    "${params.output_dir}/${params.job}/stats/${params.mod}/info/mapq_popdep"
}

def mapqPopdepPlotBase() {
    "${params.output_dir}/${params.job}/stats/${params.mod}/plots/mapq_popdep"
}

process mapq_bin_popdep_replot {
    tag 'mapq30_1mb_replot'
    cpus 1
    memory '8 GB'
    publishDir "${mapqPopdepPlotBase()}", mode: 'copy', pattern: '*.png'
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    path grid_tsv

    output:
    path '*.png', emit: plots

    script:
    def prefix = params.mapq_popdep_prefix ?: 'allmapq_1mb'
    def depth_mode = params.mapq_popdep_depth_mode ?: 'all'
    """
    #!/usr/bin/env python
    from genetics.germplasm.sample.mapq_bin_popdep import run_mapq_bin_popdep_replot

    plot_paths = run_mapq_bin_popdep_replot(
        grid_tsv="${grid_tsv}",
        output_dir=".",
        prefix="${prefix}",
        depth_mode="${depth_mode}",
    )
    print(f"plots={len(plot_paths)}")
    """
}

process mapq_bin_popdep_build {
    tag { params.mapq_popdep_prefix ?: 'allmapq_1mb' }
    cpus { params.mapq_popdep_cpus ?: 32 }
    memory { params.mapq_popdep_memory ?: '64 GB' }
    publishDir "${mapqPopdepInfoBase()}", mode: 'copy', pattern: '*.popdep_bins.tsv'
    publishDir "${mapqPopdepPlotBase()}", mode: 'copy', pattern: '*.png'
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    val mapq_dir

    output:
    path "${params.mapq_popdep_prefix ?: 'allmapq_1mb'}.popdep_bins.tsv", emit: grid
    path '*.popdep_bins.tsv', emit: all_tsv
    path '*.png', emit: plots

    script:
    def prefix = params.mapq_popdep_prefix ?: 'allmapq_1mb'
    def depth_mode = params.mapq_popdep_depth_mode ?: 'all'
    def workers = params.mapq_popdep_cpus ?: task.cpus
    def engine = params.mapq_popdep_engine ?: 'polars'
    """
    #!/usr/bin/env python
    from genetics.germplasm.sample.mapq_bin_popdep import run_mapq_bin_popdep_cohort

    grid, tsv_paths, plot_paths = run_mapq_bin_popdep_cohort(
        mapq_dir="${mapq_dir}",
        output_dir=".",
        prefix="${prefix}",
        max_workers=${workers},
        min_reads=${params.mapq_popdep_min_reads ?: 1},
        min_samples=${params.mapq_popdep_min_samples ?: 3},
        engine="${engine}",
        plot=True,
        depth_mode="${depth_mode}",
    )
    print(f"grid rows={len(grid)} tsvs={len(tsv_paths)} plots={len(plot_paths)}")
    """
}
