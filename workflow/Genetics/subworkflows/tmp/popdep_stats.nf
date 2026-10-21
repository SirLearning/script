#!/usr/bin/env nextflow
/*
 * Popdep QC stats from existing *.popdep.info.tsv + PLINK2 vmiss tables.
 *
 * Runs extended QC plots (depth + relative depth), and Mahalanobis (depth + relative depth).
 * Publishes under {output_dir}/{job}/stats/{mod}/{info,plots,logs}
 *
 * Prerequisite: popdep_annotate.nf has written
 *   {output_dir}/{job}/process/{mod}/variant/{A,B,D,Others}.popdep.info.tsv
 *
 * Launch from a vmap4 run folder. Example:
 *   mkdir -p /data/home/tusr1/01projects/vmap4/10stats.genome/40run_popdep_n500_stats_thin/run_logs
 *   cd /data/home/tusr1/01projects/vmap4/10stats.genome/40run_popdep_n500_stats_thin
 *   nextflow run .../subworkflows/tmp/popdep_stats.nf -c .../nextflow.config \
 *     --home_dir /data/home/tusr1/01projects/vmap4 \
 *     --user_dir /data/home/tusr1 \
 *     --src_dir /data/home/tusr1/git/script/src \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink \
 *     --mod test_thin
 */

nextflow.enable.dsl=2

include {
    variant_popdep_qc_plots
    variant_popdep_mahalanobis
} from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = null

workflow {
    if (!params.mod) {
        error 'popdep_stats: --mod is required (e.g. test_thin, test_common_thin)'
    }
    if (!params.output_dir || !params.job) {
        error 'popdep_stats: --output_dir and --job are required'
    }
    def vmiss_mod = params.popdep_vmiss_mod ?: params.mod
    if (!(vmiss_mod in ['test_thin', 'test_common_thin'])) {
        error "popdep_stats: vmiss mod must be test_thin or test_common_thin (got ${vmiss_mod})"
    }

    def proc_popdep = "${params.output_dir}/${params.job}/process/${params.mod}/variant"
    def proc_vmiss = "${params.output_dir}/${params.job}/process/${vmiss_mod}/variant"
    def subgenomes = ['A', 'B', 'D', 'Others']

    ch_popdep = channel.from(
        subgenomes.collect { sg ->
            tuple(sg, "sub_${sg}", file("${proc_popdep}/${sg}.popdep.info.tsv"))
        }
    )
    ch_vmis = channel.from(
        subgenomes.collect { sg ->
            tuple(sg, "sub_${sg}", file("${proc_vmiss}/${sg}.info.vmiss"))
        }
    )

    variant_popdep_qc_plots(ch_popdep.combine(ch_vmis, by: [0, 1]))
    variant_popdep_mahalanobis(ch_popdep)
}
