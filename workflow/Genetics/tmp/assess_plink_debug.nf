nextflow.enable.dsl=2

/*
 * Debug assess track for test_thin / test_common_thin:
 * PLINK2 --freq counts / --missing on a narrow chromosome slice from existing *_test.plink2,
 * MAF-bin table from .acount (no bcftools), singleton / MAC summaries in Python, then
 * plots via infra.utils.graph.
 *
 * Per subgenome, one representative PLINK chromosome: A=1, B=3, D=5, Others=0.
 */
include { plink2_assess_debug_slice } from '../modules/local/genotype/assess.nf'

process assess_plink_debug_plots {
    tag "assess plots ${id}"
    label 'cpus_2'
    cpus 2
    memory '8.GB'
    conda "${params.user_dir}/miniconda3/envs/stats"
    publishDir "${params.output_dir}/${params.job}/assess/${params.mod}/plots", mode: 'copy', pattern: "*.png"
    publishDir "${params.output_dir}/${params.job}/assess/${params.mod}/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/assess/${params.mod}/logs", mode: 'copy', pattern: "*.log"

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

workflow {
    if (!(params.mod in ['test_thin', 'test_common_thin'])) {
        error "assess_plink_debug.nf: params.mod must be test_thin or test_common_thin (got: ${params.mod})."
    }
    if (!params.output_dir || !params.job) {
        error "params.output_dir and params.job are required."
    }

    def sub = ['A', 'B', 'D', 'Others']
    def chrPick = [
        'A'     : '1',
        'B'     : '3',
        'D'     : '5',
        'Others': '0',
    ]
    def ch_in = channel.from(sub.collect { sg -> tuple(sg, "${sg}_test.plink2", chrPick[sg]) })
    plink2_assess_debug_slice(ch_in)
    def b = plink2_assess_debug_slice.out.bundle
    def plot_in = b.map { id, ac, vm, co, mh, gq -> tuple(id, ac, vm, co, mh, gq) }
    assess_plink_debug_plots(plot_in)
}
