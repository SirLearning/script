nextflow.enable.dsl=2

/*
 * Selection-statistic plots (FST, pi, Tajima D, ideogram density).
 * Wire from upstream compute workflows when fst/pi/tajima/ideogram inputs exist.
 */

process plot_fst_boxplot {
    tag "plot fst"
    publishDir "${params.output_dir}/${params.job}/population_genetics/selection/plots/fst", mode: 'copy'

    input:
    path fst_files
    val output_name

    output:
    path "${output_name}"

    script:
    """
    Rscript ${params.src_dir}/r/genetics/dynamic/plot/plot_fst_boxplot.R \\
        --input_dir . \\
        --output ${output_name} \\
        --pattern ".fst"
    """
}

process plot_pi_distribution {
    tag "plot pi"
    publishDir "${params.output_dir}/${params.job}/population_genetics/selection/plots/pi", mode: 'copy'

    input:
    path pi_files
    val output_name

    output:
    path "${output_name}"

    script:
    """
    Rscript ${params.src_dir}/r/genetics/dynamic/plot/plot_pi_distribution.R \\
        --input_dir . \\
        --output ${output_name} \\
        --pattern ".pi"
    """
}

process plot_tajima_d_distribution {
    tag "plot tajima d"
    publishDir "${params.output_dir}/${params.job}/population_genetics/selection/plots/tajima_d", mode: 'copy'

    input:
    path tajima_files
    val output_name

    output:
    path "${output_name}"

    script:
    """
    Rscript ${params.src_dir}/r/genetics/dynamic/plot/plot_tajima_d.R \\
        --input_dir . \\
        --output ${output_name} \\
        --pattern ".Tajima.D"
    """
}

process plot_ideogram_density {
    tag "plot ideogram"
    publishDir "${params.output_dir}/${params.job}/population_genetics/selection/plots/ideogram", mode: 'copy'

    input:
    path karyotype
    path density
    val output_name

    output:
    path "${output_name}"

    script:
    """
    Rscript ${params.src_dir}/r/genetics/dynamic/plot/plot_ideogram_density.R \\
        --karyotype ${karyotype} \\
        --density ${density} \\
        --output ${output_name}
    """
}
