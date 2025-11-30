nextflow.enable.dsl=2

process PLOT_FST_BOXPLOT {
    tag "plot fst"
    publishDir "${params.outdir}/plots/fst", mode: 'copy'

    input:
    path fst_files // List of files or a directory containing them
    val output_name

    output:
    path "${output_name}"

    script:
    """
    # If input is a list of files, they are staged in current dir.
    # If input is a directory, we point to it.
    
    # We assume fst_files is a collection of files. 
    # We can pass "." as input dir since files are staged here.
    
    Rscript ${params.src_dir}/r/genetics/dynamic/plot/plot_fst_boxplot.R \\
        --input_dir . \\
        --output ${output_name} \\
        --pattern ".fst"
    """
}

process PLOT_PI_DISTRIBUTION {
    tag "plot pi"
    publishDir "${params.outdir}/plots/pi", mode: 'copy'

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

process PLOT_TAJIMA_D_DISTRIBUTION {
    tag "plot tajima d"
    publishDir "${params.outdir}/plots/tajima_d", mode: 'copy'

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

process PLOT_IDEOGRAM {
    tag "plot ideogram"
    publishDir "${params.outdir}/plots/ideogram", mode: 'copy'

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

workflow plots {
    take:
    ch_fst_files
    ch_pi_files
    ch_tajima_files
    
    main:
    if (ch_fst_files) {
        PLOT_FST_BOXPLOT(ch_fst_files.collect(), "fst_boxplot.pdf")
    }
    
    if (ch_pi_files) {
        PLOT_PI_DISTRIBUTION(ch_pi_files.collect(), "pi_distribution.pdf")
    }
    
    if (ch_tajima_files) {
        PLOT_TAJIMA_D_DISTRIBUTION(ch_tajima_files.collect(), "tajima_d_distribution.pdf")
    }
}
