nextflow.enable.dsl=2

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

// --- Assess debug plots (consume processor assess outputs; wired from tmp/assess_*.nf) ---

process assess_plink_debug_plots {
    tag "assess plots ${id}"
    label 'cpus_2'
    cpus 2
    memory '8.GB'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(id), path(acount), path(vmiss), path(counts), path(mac_hist), path(gq_summary)

    output:
    path("*.png"), emit: plots
    path("*.log"), emit: logs
    path("*.tsv"), emit: tables

    script:
    """
    #!/usr/bin/env python
    import sys
    sys.stdout = open("${id}.assess_plots.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.assess_slice import ana_assess_plink_debug_slice

    ana_assess_plink_debug_slice(
        acount_path="${acount}",
        vmiss_path="${vmiss}",
        mac_hist_path="${mac_hist}",
        counts_path="${counts}",
        output_prefix="${id}.assess",
    )
    """
}

process dumpnice_vcf_qc_assess {
    tag "dumpnice vcf qc assess: ${id}"
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/stats/${params.mod}/plots", mode: 'copy', pattern: "*.pdf"

    input:
    tuple val(id), path(vcf)

    output:
    tuple val(id), path("*.pdf"), emit: qc_plots

    script:
    def inputdir = params.dumpnice_inputdir
    """
    set -euo pipefail
    if [ -d "${inputdir}" ]; then
        rsync -a "${inputdir}/" ./ || true
    fi
    Rscript ${params.src_dir}/r/dumpnice/vcf/12_VcfQc.r || true
    ls -1 *.pdf >/dev/null 2>&1 || touch ${id}.vcf_qc_placeholder.pdf
    """
}
