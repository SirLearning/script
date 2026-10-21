nextflow.enable.dsl=2

/*
 * Wheat selfing-model excess-heterozygosity probability sweep.
 *
 * Estimates F_panel as the median of PLINK2 --het sample F on common SNPs
 * (params.het_source_mod, default test_common_only), then scores rare-variant
 * gcount tables under params.mod (e.g. test_rare_only). Does not write a
 * filtered pfile; publishes retain/drop counts across P_excess thresholds,
 * MAC heatmaps, popdep mean/SD scatters coloured by the het filter, and
 * optional R_H=H/MAC band diagnostics (`mac_het_rh_diag`).
 *
 * Launch: partial_router.nf --partial_task mac_het_excess --mod test_rare_only
 */

include { listMergedSubgenomeTestPfileTuples } from '../../../modules/local/infra/infra_plink_reuse.nf'
include {
    plink2_het_inbreeding
    plink2_het_extract
} from '../../../modules/local/genotype/processor/processor_plink2.nf'
include {
    variant_mac_het_excess_stats
    variant_mac_het_excess_combine
    variant_mac_hetfrac_heatmap_retained
    variant_popdep_het_filter_scatter
    variant_mac_het_rh_diagnostics
    variant_mac_het_rh_combine
    variant_mac_het_rh_sample_burden
} from '../../../modules/local/genotype/stats/stats_variant.nf'

workflow RUN_MAC_HET_EXCESS {
    main:
    if (!params.mod) {
        error 'RUN_MAC_HET_EXCESS: params.mod is required (e.g. test_rare_only).'
    }
    if (!params.output_dir || !params.job) {
        error 'RUN_MAC_HET_EXCESS: params.output_dir and params.job are required.'
    }

    def het_mod = params.het_source_mod ?: 'test_common_only'
    def gcount_mod = params.mac_het_variant_mod ?: params.mod
    def common_proc = "${params.output_dir}/${params.job}/process/${het_mod}"
    def rare_variant = "${params.output_dir}/${params.job}/process/${gcount_mod}/variant"
    def popdep_variant = "${params.output_dir}/${params.job}/process/${params.mod}/variant"
    def het_dir = "${common_proc}/sample"
    def subgenomes = ['A', 'B', 'D', 'Others']

    log.info "mac_het_excess: gcount ${rare_variant}; popdep ${popdep_variant}"

    def existing_het = []
    def missing_sg = []
    subgenomes.each { sg ->
        def hetf = file("${het_dir}/${sg}.info.het")
        if (hetf.exists()) {
            existing_het << tuple(sg, "sub_${sg}", hetf)
        } else {
            missing_sg << sg
        }
    }

    ch_existing = channel.fromList(existing_het)
    if (missing_sg) {
        log.info "mac_het_excess: PLINK2 --het for ${missing_sg} from ${common_proc}"
        def pfiles = listMergedSubgenomeTestPfileTuples(common_proc).findAll { row -> row[0] in missing_sg }
        pfiles.each { row ->
            if (!row[3].exists()) {
                error "RUN_MAC_HET_EXCESS: missing pfile ${row[3]} under ${common_proc}"
            }
        }
        ch_het = ch_existing.mix(plink2_het_inbreeding(channel.fromList(pfiles)).het)
    } else {
        log.info "mac_het_excess: reuse existing --het tables under ${het_dir}"
        ch_het = ch_existing
    }

    ch_gcount = channel.fromList(
        subgenomes.collect { sg ->
            tuple(sg, "sub_${sg}", file("${rare_variant}/${sg}.info.gcount", checkIfExists: true))
        }
    )
    ch_joined = ch_gcount.join(ch_het, by: [0, 1])

    if ("${params.mac_het_sweep}" != 'false') {
        stats_out = variant_mac_het_excess_stats(ch_joined)
        variant_mac_het_excess_combine(stats_out.sweep.map { _id, _chr, sweep -> sweep }.collect())
    }
    if ("${params.mac_het_heatmap}" != 'false') {
        variant_mac_hetfrac_heatmap_retained(ch_joined)
    }
    if ("${params.mac_het_popdep_scatter}" != 'false') {
        ch_popdep = channel.fromList(
            subgenomes.collect { sg ->
                tuple(sg, "sub_${sg}", file("${popdep_variant}/${sg}.popdep.info.tsv", checkIfExists: true))
            }
        )
        variant_popdep_het_filter_scatter(ch_joined.join(ch_popdep, by: [0, 1]))
    }
    if ("${params.mac_het_rh_diag}" != 'false') {
        def no_mapq30 = file("${moduleDir}/../../../assets/NO_POPDEP_MAPQ30", checkIfExists: true)
        def mq30_mod = params.mac_het_popdep_mapq30_mod
        def mq30_variant = mq30_mod ? "${params.output_dir}/${params.job}/process/${mq30_mod}/variant" : null
        ch_popdep_rh = channel.fromList(
            subgenomes.collect { sg ->
                tuple(sg, "sub_${sg}", file("${popdep_variant}/${sg}.popdep.info.tsv", checkIfExists: true))
            }
        )
        ch_mq30 = channel.fromList(
            subgenomes.collect { sg ->
                def mq30f = mq30_variant ? file("${mq30_variant}/${sg}.popdep.info.tsv", checkIfExists: true) : no_mapq30
                tuple(sg, "sub_${sg}", mq30f, mq30_variant ? true : false)
            }
        )
        ch_mq = channel.fromList(
            subgenomes.collect { sg ->
                tuple(sg, "sub_${sg}", file("${rare_variant}/${sg}.mq.info.tsv", checkIfExists: true))
            }
        )
        ch_vmiss = channel.fromList(
            subgenomes.collect { sg ->
                tuple(sg, "sub_${sg}", file("${rare_variant}/${sg}.info.vmiss", checkIfExists: true))
            }
        )
        ch_rh_in = ch_joined
            .join(ch_popdep_rh, by: [0, 1])
            .join(ch_mq30, by: [0, 1])
            .join(ch_mq, by: [0, 1])
            .join(ch_vmiss, by: [0, 1])
        rh_out = variant_mac_het_rh_diagnostics(ch_rh_in)
        variant_mac_het_rh_combine(rh_out.th.map { _id, _chr, th -> th }.collect())
        if ("${params.mac_het_rh_sample}" != 'false') {
            def rare_proc = "${params.output_dir}/${params.job}/process/${gcount_mod}"
            ch_pfile = channel.fromList(listMergedSubgenomeTestPfileTuples(rare_proc))
            ch_extract_pfile = ch_pfile.join(rh_out.extract, by: [0, 1])
            het_all = plink2_het_extract(ch_extract_pfile)
            variant_mac_het_rh_sample_burden(het_all.het)
        }
    }
}
