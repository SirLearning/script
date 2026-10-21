#!/usr/bin/env nextflow
/*
 * Three-ring Circos map:
 * inner rare-variant density and middle common-variant density from --mod
 * (typically test_wap_thin_refchr); outer LD blocks from --circos_ld_mod
 * (typically test_wap_common_refchr; defaults to --mod).
 *
 * LD blocks: adjacent r2 >= 0.8, merge if junction Q1(|D'|) > 0.98.
 * PLINK 1.9 --r2 dprime on merged A/B/D bfiles (not the 40+ GB vcor).
 *
 *   nextflow run .../subworkflows/tmp/circos_ld_map.nf \
 *     -c .../nextflow.config \
 *     --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
 *     --job test_plink --mod test_wap_thin_refchr \
 *     --circos_ld_mod test_wap_common_refchr
 *
 * --circos_ld_map_plot_only true replots from published
 * stats/<mod>/info/all.variant.circos_ld_map.{density,blocks}.info.tsv.
 */

nextflow.enable.dsl=2

include { plink1_ld_dprime_window } from '../../modules/local/genotype/processor/processor_plink2.nf'
include {
    variant_circos_ld_map
    variant_circos_ld_map_replot
} from '../../modules/local/genotype/stats/stats_variant.nf'
include {
    listMergedSubgenomeTestPlink1BedTuples
} from '../../modules/local/infra/infra_plink_reuse.nf'

params.job = 'test_plink'
params.mod = 'test_wap_thin_refchr'
params.circos_ld_mod = null
params.process_dir = null
params.circos_ld_process_dir = null
params.circos_ld_map_plot_only = false

workflow {
    if (!params.output_dir || !params.job || !params.mod) {
        error 'circos_ld_map: --output_dir, --job, and --mod are required'
    }

    if (params.circos_ld_map_plot_only) {
        def info_dir = "${params.output_dir}/${params.job}/stats/${params.mod}/info"
        def density = file("${info_dir}/all.variant.circos_ld_map.density.info.tsv")
        def blocks = file("${info_dir}/all.variant.circos_ld_map.blocks.info.tsv")
        if (!density.exists() || !blocks.exists()) {
            error "circos_ld_map plot-only: missing ${density} or ${blocks}"
        }
        variant_circos_ld_map_replot(tuple(density, blocks))
        return
    }

    def density_dir = params.process_dir
        ? params.process_dir
        : "${params.output_dir}/${params.job}/process/${params.mod}"
    def ld_mod = params.circos_ld_mod ?: params.mod
    def ld_dir = params.circos_ld_process_dir
        ? params.circos_ld_process_dir
        : "${params.output_dir}/${params.job}/process/${ld_mod}"

    def beds = listMergedSubgenomeTestPlink1BedTuples(ld_dir)
    if (beds.size() < 3) {
        error "circos_ld_map: need A/B/D ${ld_dir}/*_test.plink.{bed,bim,fam}"
    }

    Channel
        .from(['A', 'B', 'D'])
        .map { sg ->
            def pvar = file("${density_dir}/${sg}_test.plink2.pvar")
            def afreq = file("${density_dir}/variant/${sg}.info.afreq")
            if (!pvar.exists() || !afreq.exists()) {
                error "circos_ld_map: missing density ${pvar} or ${afreq}"
            }
            tuple(sg, pvar, afreq)
        }
        .set { ch_density }

    Channel
        .from(['A', 'B', 'D'])
        .map { sg ->
            def pvar = file("${ld_dir}/${sg}_test.plink2.pvar")
            if (!pvar.exists()) {
                error "circos_ld_map: missing LD pvar ${pvar}"
            }
            tuple(sg, pvar)
        }
        .set { ch_ld_pvar }

    Channel
        .from(beds)
        .set { ch_bed }
    plink1_ld_dprime_window(ch_bed)
    ch_ld = plink1_ld_dprime_window.out.ld

    ch_density
        .join(ch_ld_pvar)
        .join(ch_ld)
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { rows ->
            tuple(
                rows.collect { it[1] },
                rows.collect { it[2] },
                rows.collect { it[3] },
                rows.collect { it[4] },
            )
        }
        .set { ch_plot }
    variant_circos_ld_map(ch_plot)
}
