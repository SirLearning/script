nextflow.enable.dsl=2

// --- Genotype statistics workflow ---

workflow stats {
    take:
        ch_vcf // tuple val(id), path(vcf)

    main:
        // 1. Convert VCF to PLINK binary format
        ch_bfiles = PLINK_FROM_VCF(ch_vcf)
        
        // 2. Calculate statistics
        ch_missing = PLINK_MISSING(ch_bfiles.bfiles)
        ch_freq = PLINK_FREQ(ch_bfiles.bfiles)
        ch_het = PLINK_HET(ch_bfiles.bfiles)
        ch_pca = PLINK_PCA(ch_bfiles.bfiles)

        ch_pca_plot = channel.empty()
        if (params.enable_pca_plot) {
            def md = params.sample_metadata ? file(params.sample_metadata) : file('NO_METADATA')
            ch_pca_plot = PLOT_PLINK_PCA(ch_pca.pca, md)
        }

        ch_plots = channel.empty()
        if (params.enable_simple_plots) {
            ch_plots = PLOT_PLINK_QC(ch_missing.missing, ch_freq.freq, ch_het.het)
        }

    emit:
        bfiles = ch_bfiles.bfiles
        missing = ch_missing.missing
        freq = ch_freq.freq
        het = ch_het.het
        pca = ch_pca.pca
        plots = ch_plots
        pca_plot = ch_pca_plot
}

workflow plink_stats {
    take:
        smiss
        vmiss
        gcount
        afreq
        hardy

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
    publishDir "${params.output_dir}/${params.job}/assess/threshold", mode: 'copy', pattern: "*.smiss_th.tsv"
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

process vcftools_vcf_qc {
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
        --output ${id}.qc_plots.pdf
    """
}

// --- old code ---

process PLINK_FROM_VCF {
    tag "${id}" 
    publishDir "${params.output_dir}/${params.job}/stats/plink", mode: 'copy'

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("${id}.bed"), path("${id}.bim"), path("${id}.fam"), emit: bfiles

    script:
    """
    set -euo pipefail
    plink --vcf ${vcf} \\
        --make-bed \\
        --out ${id} \\
        --allow-extra-chr \\
        --chr-set 42 \\
        --double-id \\
        --threads ${task.cpus}
    """
}

process PLINK_MISSING {
    tag "${id}" 
    publishDir "${params.output_dir}/${params.job}/stats/missing", mode: 'copy'

    input:
    tuple val(id), path(bed), path(bim), path(fam)

    output:
    tuple val(id), path("${id}.imiss"), path("${id}.lmiss"), emit: missing

    script:
    """
    set -euo pipefail
    plink --bfile ${id} --missing --out ${id} --allow-extra-chr --chr-set 42 --threads ${task.cpus}
    """
}

process PLINK_FREQ {
    tag "${id}" 
    publishDir "${params.output_dir}/${params.job}/stats/freq", mode: 'copy'

    input:
    tuple val(id), path(bed), path(bim), path(fam)

    output:
    tuple val(id), path("${id}.frq"), emit: freq

    script:
    """
    set -euo pipefail
    plink --bfile ${id} --freq --out ${id} --allow-extra-chr --chr-set 42 --threads ${task.cpus}
    """
}

process PLINK_HET {
    tag "${id}" 
    publishDir "${params.output_dir}/${params.job}/stats/het", mode: 'copy'

    input:
    tuple val(id), path(bed), path(bim), path(fam)

    output:
    tuple val(id), path("${id}.het"), emit: het

    script:
    """
    set -euo pipefail
    plink --bfile ${id} --het --out ${id} --allow-extra-chr --chr-set 42 --threads ${task.cpus}
    """
}

process PLINK_PCA {
    tag "${id}"
    publishDir "${params.output_dir}/${params.job}/stats/pca", mode: 'copy'

    input:
    tuple val(id), path(bed), path(bim), path(fam)

    output:
    tuple val(id), path("${id}.eigenvec"), path("${id}.eigenval"), emit: pca

    script:
    """
    set -euo pipefail
    plink --bfile ${id} --pca 5 --out ${id} --allow-extra-chr --chr-set 42 --threads ${task.cpus}
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

process PLOT_PLINK_QC {
    tag "${id}"
    publishDir "${params.output_dir}/${params.job}/stats/plots", mode: 'copy'

    input:
    tuple val(id), path(imiss), path(lmiss)
    tuple val(id2), path(frq)
    tuple val(id3), path(het)

    output:
    tuple val(id), path("${id}.qc_plots.pdf"), emit: qc_plots

    script:
    """
    set -euo pipefail
    Rscript ${params.src_dir}/r/genetics/plink_qc_plot.r \\
        --imiss ${imiss} \\
        --lmiss ${lmiss} \\
        --frq ${frq} \\
        --het ${het} \\
        --output ${id}.qc_plots.pdf
    """
}
