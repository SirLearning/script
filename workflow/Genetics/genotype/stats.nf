nextflow.enable.dsl=2

// --- Genotype statistics workflow ---
workflow plink_stats {
    take:
    // Expect a channel: [ val(id), path(smiss) ]
    smiss
    vmiss
    gcount
    afreq
    hardy
    popdep

    main:
    def stats_out = sample_missing_stats(smiss)

    emit:
    bfiles = stats_out.bfiles
    missing = stats_out.missing
    freq = stats_out.freq
    het = stats_out.het
    pca = stats_out.pca
    plots = stats_out.plots
    pca_plot = stats_out.pca_plot
}

process sample_missing_stats {
    tag "compute missing rate threshold: ${id}"
    publishDir "${params.output_dir}/${params.job}/assess/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/assess/thresholds", mode: 'copy', pattern: "*.smiss_th.tsv"
    publishDir "${params.output_dir}/${params.job}/assess/logs", mode: 'copy', pattern: "*.log"
    conda 'stats'

    input:
    tuple val(id), val(chr), path(smiss)

    output:
    tuple val(id), path("*.smiss_th.tsv"), emit: smiss_th
    tuple val(id), path("*.png"), emit: plots
    tuple val(id), path("*.log"), emit: logs

    script:
    """
    #!/usr/bin/env python
    import sys
    
    sys.stdout = open("${id}.smiss_ana.log", "w")
    sys.stderr = sys.stdout

    from python_script.genomics.sample.smiss_ana import calculate_missing_threshold, plot_missing_dist

    print(f"Processing sample missing rate for ${id}...")
    plot_missing_dist("${smiss}", "${chr}.missing_dist")
    
    print(f"Calculating threshold...")
    calculate_missing_threshold("${smiss}", "${id}.smiss_th.tsv")
    """
}

process sample_het_stats {
    tag "compute heterozygosity rate: ${id}"
    publishDir "${params.output_dir}/${params.job}/assess/het", mode: 'copy', pattern: "*.het_stats.tsv"
    publishDir "${params.output_dir}/${params.job}/assess/logs", mode: 'copy', pattern: "*.log"
    conda 'stats'

    input:
    tuple val(id), val(chr), path(het)

    output:
    tuple val(id), path("${id}.het_stats.tsv"), emit: het_stats
    tuple val(id), path("${id}.het_ana.log"), emit: log

    script:
    """
    #!/usr/bin/env python
    import sys
    
    sys.stdout = open("${id}.het_ana.log", "w")
    sys.stderr = sys.stdout

    from python_script.genomics.sample.het_ana import compute_het_stats

    print(f"Processing heterozygosity rate for ${id}...")
    compute_het_stats("${het}", "${id}.het_stats.tsv")
    """
}

process plink_pca {
    tag "${id}"
    publishDir "${params.output_dir}/${params.job}/stats/pca", mode: 'copy'

    input:
    tuple val(id), path(bed), path(bim), path(fam)

    output:
    tuple val(id), path("${id}.eigenvec"), path("${id}.eigenval"), emit: pca

    script:
    """
    set -euo pipefail
    exec > ${id}.plink_pca.log 2>&1

    plink --bfile ${id} \\
        --pca ${params.pc_num} \\
        --out ${id} \\
        --allow-extra-chr \\
        --chr-set 42 \\
        --threads ${task.cpus}
    """
}

process PLOT_PLINK_PCA {
    tag "${id}"
    publishDir "${params.output_dir}/${params.job}/stats/plots", mode: 'copy'

    input:
    tuple val(id), path(eigenvec), path(eigenval)
    path metadata

    output:
    tuple val(id), path("${id}.pca_plot.pdf"), emit: pca_plot

    script:
    def md_arg = metadata.name != 'NO_METADATA' ? "--metadata ${metadata}" : ""
    """
    set -euo pipefail
    Rscript ${params.src_dir}/r/genetics/plot_pca.r \\
        --eigenvec ${eigenvec} \\
        --eigenval ${eigenval} \\
        ${md_arg} \\
        --output ${id}.pca_plot.pdf
    """
}

process vcftools_vcf_qc_r {
    tag "${id}" ? "plot qc ${id}" : 'plot qc'
    publishDir "${params.output_dir}/${params.job}/stats", mode: 'copy'

    input:
    tuple val(id), path(imiss), path(lmiss), path(het), path(frq), path(idepth), path(ldepth)

    output:
    path("${id}.qc_plots.pdf")

    script:
    """
    set -euo pipefail
    Rscript ${params.src_dir}/r/genetics/vcf_qc_plot.r \\
        --imiss ${imiss} \\
        --lmiss ${lmiss} \\
        --het ${het} \\
        --frq ${frq} \\
        --depth ${idepth} \\
        --site_depth ${ldepth} \\
        --output ${id}.qc_plots.png
    """
}

