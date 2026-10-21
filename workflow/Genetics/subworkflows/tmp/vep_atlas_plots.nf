#!/usr/bin/env nextflow
/*
 * Composition Figure 3 plus window-count Circos / mirror / violin from published
 * VEP tables. Genome-wide sites: --mod (default test_wap_refchr). Exon-interval
 * sites: --vep_atlas_exon_mod (default test_wap_exon; GFF exon, including UTR).
 * y-axis is counts per 10 Mb window, not per callable Mb.
 * Common/rare tracks use --mod; rare predicted LoF uses --vep_atlas_exon_mod.
 *
 * Launch from a vmap4 run folder; log in doc/progress/YYYY-MM-DD.md (#### NF replay).
 *   nextflow run .../subworkflows/tmp/vep_atlas_plots.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_wap_refchr \
 *     --vep_atlas_exon_mod test_wap_exon
 */

nextflow.enable.dsl=2

include { variant_vep_atlas } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = 'test_wap_refchr'
params.vep_atlas_exon_mod = 'test_wap_exon'
params.vep_atlas_bin_bp = 10000000

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'vep_atlas_plots: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir) {
        error 'vep_atlas_plots: --user_dir is required (conda stats)'
    }
    if (!params.src_dir) {
        error 'vep_atlas_plots: --src_dir is required'
    }

    def genome_dir = "${params.output_dir}/${params.job}/process/${params.mod}"
    def exon_dir = "${params.output_dir}/${params.job}/process/${params.vep_atlas_exon_mod}"
    def ids = ['A', 'B', 'D']

    def genome_vep = ids.collect { sg -> file("${genome_dir}/variant/${sg}.vep.tsv", checkIfExists: true) }
    def exon_vep = ids.collect { sg -> file("${exon_dir}/variant/${sg}.vep.tsv", checkIfExists: true) }
    def genome_afreq = ids.collect { sg -> file("${genome_dir}/variant/${sg}.info.afreq", checkIfExists: true) }
    def genome_acount = ids.collect { sg -> file("${genome_dir}/variant/${sg}.assoc.acount", checkIfExists: true) }
    def genome_pvar = ids.collect { sg -> file("${genome_dir}/${sg}_test.plink2.pvar", checkIfExists: true) }
    def exon_afreq = ids.collect { sg -> file("${exon_dir}/variant/${sg}.info.afreq", checkIfExists: true) }
    def exon_acount = ids.collect { sg -> file("${exon_dir}/variant/${sg}.assoc.acount", checkIfExists: true) }
    def exon_pvar = ids.collect { sg -> file("${exon_dir}/${sg}_test.plink2.pvar", checkIfExists: true) }

    log.info "${params.c_green}vep_atlas ${params.mod} common/rare + exon LoF ${params.vep_atlas_exon_mod}${params.c_reset}"
    variant_vep_atlas(tuple(genome_vep, exon_vep, genome_afreq, genome_acount, genome_pvar, exon_afreq, exon_acount, exon_pvar))
}
