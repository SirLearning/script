#!/usr/bin/env nextflow
/*
 * First three cascade steps on test_wap_refchr and test_wap_exon, then the
 * VEP atlas from the sites that remain.
 *
 * Steps, for keep 80% and keep 60% separately:
 *   PopDep reliable BED → excess het (P>=1e-6) → MAC≥5.
 * Inbreeding F is plink2 --het on that mod's own pfile (work dir only).
 * Pfiles under process/ are not rewritten. Atlas PNGs use a cv_reliable
 * prefix so the existing all.variant.vep_atlas.* figures stay in place.
 *
 * Launch from a vmap4 run folder:
 *   nextflow run .../subworkflows/tmp/wap_cv_reliable_vep_atlas.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_wap_refchr \
 *     --popdep_cv_reliable_dir /data1/.../process/_popdep_n500_cv_reliable
 */

nextflow.enable.dsl=2

include { variant_wap_cv_reliable_filter; variant_vep_atlas_cv_reliable } from '../../modules/local/genotype/stats/stats_variant.nf'

params.job = 'test_plink'
params.mod = 'test_wap_refchr'
params.vep_atlas_exon_mod = 'test_wap_exon'
params.vep_atlas_bin_bp = 10000000
params.filter_cascade_het_alpha = '1e-6'
params.popdep_cv_reliable_dir = null
params.vep_atlas_ld_blocks = '/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/stats/test_wap_thin_refchr/info/all.variant.circos_ld_map.blocks_1mb.info.tsv'
params.vep_atlas_gene_gff = '/data/public_data/wheat/00_genome/reference/annotation/gene/download/iwgsc_refseqv1.1_genes_2017July06/IWGSC_v1.1_HC_20170706.gff3'

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'wap_cv_reliable_vep_atlas: --output_dir, --job, and --mod are required'
    }
    if (!params.user_dir) {
        error 'wap_cv_reliable_vep_atlas: --user_dir is required (conda stats)'
    }
    if (!params.src_dir) {
        error 'wap_cv_reliable_vep_atlas: --src_dir is required'
    }
    def bed_root = params.popdep_cv_reliable_dir
    if (!bed_root) {
        error 'wap_cv_reliable_vep_atlas: --popdep_cv_reliable_dir is required'
    }
    file("${bed_root}/keep80/chr001.popdep.cv_reliable.keep80.bed.gz", checkIfExists: true)
    file("${bed_root}/keep60/chr001.popdep.cv_reliable.keep60.bed.gz", checkIfExists: true)

    def genome_mod = params.mod
    def exon_mod = params.vep_atlas_exon_mod
    def genome_dir = "${params.output_dir}/${params.job}/process/${genome_mod}"
    def exon_dir = "${params.output_dir}/${params.job}/process/${exon_mod}"
    def ids = ['A', 'B', 'D']
    def genome_vep = ids.collect { sg -> file("${genome_dir}/variant/${sg}.vep.tsv", checkIfExists: true) }
    def exon_vep = ids.collect { sg -> file("${exon_dir}/variant/${sg}.vep.tsv", checkIfExists: true) }
    def genome_afreq = ids.collect { sg -> file("${genome_dir}/variant/${sg}.info.afreq", checkIfExists: true) }
    def genome_acount = ids.collect { sg -> file("${genome_dir}/variant/${sg}.assoc.acount", checkIfExists: true) }
    def genome_pvar = ids.collect { sg -> file("${genome_dir}/${sg}_test.plink2.pvar", checkIfExists: true) }
    def exon_afreq = ids.collect { sg -> file("${exon_dir}/variant/${sg}.info.afreq", checkIfExists: true) }
    def exon_acount = ids.collect { sg -> file("${exon_dir}/variant/${sg}.assoc.acount", checkIfExists: true) }
    def exon_pvar = ids.collect { sg -> file("${exon_dir}/${sg}_test.plink2.pvar", checkIfExists: true) }
    def ld_blocks = file(params.vep_atlas_ld_blocks, checkIfExists: true)
    def gene_gff = file(params.vep_atlas_gene_gff, checkIfExists: true)

    log.info "${params.c_green}cv-reliable three-step ${genome_mod} + ${exon_mod}, then VEP atlas${params.c_reset}"
    def filtered = variant_wap_cv_reliable_filter(Channel.of(genome_mod, exon_mod))
    def genome80 = filtered.keep80.filter { mod, _sites -> mod == genome_mod }.map { _mod, sites -> sites }
    def genome60 = filtered.keep60.filter { mod, _sites -> mod == genome_mod }.map { _mod, sites -> sites }
    def exon80 = filtered.keep80.filter { mod, _sites -> mod == exon_mod }.map { _mod, sites -> sites }
    def exon60 = filtered.keep60.filter { mod, _sites -> mod == exon_mod }.map { _mod, sites -> sites }
    variant_vep_atlas_cv_reliable(
        genome80
            .combine(genome60)
            .combine(exon80)
            .combine(exon60)
            .combine(Channel.value(ld_blocks))
            .combine(Channel.value(gene_gff))
    )
}
