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
