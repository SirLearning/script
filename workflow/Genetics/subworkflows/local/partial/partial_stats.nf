nextflow.enable.dsl=2

/*
 * Partial stats reruns from existing process/<mod>/ or stats/<mod>/ artefacts.
 * Wired from partial_router.nf (--partial_task ld_redraw | mac_stats | mac_dist_redraw | chr_counts | chr_compare | rebuild_lib_stats).
 */

include { variant_ld_decay_plot; variant_ld_crosschr_plot; variant_mac_stats; variant_mac_maf_reg; variant_mac_missing_reg; variant_mac_missing_reg_bin50_sample; variant_mac_dist_log_redraw } from '../../../modules/local/genotype/stats/stats_variant.nf'
include { report_plink_chr_variant_counts; plot_thin_common_chr_variant_compare } from '../../../modules/local/genotype/stats/stats_chr_report.nf'
include { test_plink_stats as TEST_PLINK_STATS } from '../../../modules/local/genotype/stats/stats.nf'

workflow RUN_LD_PLOTS_REDRAW {
    main:
    if (!params.mod) {
        error 'params.mod is required (e.g. test_thin, test_common_thin).'
    }
    if (!params.output_dir || !params.job) {
        error 'params.output_dir and params.job are required.'
    }

    def subgenomes = ['A', 'B', 'D', 'Others']
    def base = "${params.output_dir}/${params.job}/process/${params.mod}/variant"

    def ch_ld = channel.from(
        subgenomes.collect { sg ->
            tuple(sg, "sub_${sg}", file("${base}/${sg}.info.ld.vcor", checkIfExists: true))
        }
    )
    def ch_ld_cross = channel.from(
        subgenomes.collect { sg ->
            tuple(sg, "sub_${sg}", file("${base}/${sg}.info.ld.crosschr.vcor", checkIfExists: true))
        }
    )

    variant_ld_decay_plot(ch_ld)
    variant_ld_crosschr_plot(ch_ld_cross)
}

workflow RUN_MAC_STATS_FROM_GCOUNT {
    main:
    if (!params.mod) {
        error 'params.mod is required (e.g. test_thin, test_common_thin, test_rebulld_lib).'
    }
    if (!params.output_dir || !params.job) {
        error 'params.output_dir and params.job are required.'
    }

    def proc = "${params.output_dir}/${params.job}/process/${params.mod}/variant"

    if (params.mod in ['test_thin', 'test_common_thin']) {
        def subgenomes = ['A', 'B', 'D', 'Others']
        def ch_gcount = channel.from(
            subgenomes.collect { sg ->
                tuple(sg, "sub_${sg}", file("${proc}/${sg}.info.gcount"))
            }
        )
        def ch_vmis = channel.from(
            subgenomes.collect { sg ->
                tuple(sg, "sub_${sg}", file("${proc}/${sg}.info.vmiss"))
            }
        )
        variant_mac_stats(ch_gcount)
        variant_mac_maf_reg(ch_gcount)
        variant_mac_missing_reg(ch_gcount.combine(ch_vmis, by: [0, 1]))
    } else if (params.mod == 'test_rebulld_lib') {
        def ch_gcount = channel.of(tuple('chr002', '2', file("${proc}/chr002.info.gcount")))
        def ch_vmis = channel.of(tuple('chr002', '2', file("${proc}/chr002.info.vmiss")))
        variant_mac_stats(ch_gcount)
        variant_mac_maf_reg(ch_gcount)
        variant_mac_missing_reg(ch_gcount.combine(ch_vmis, by: [0, 1]))
    } else {
        error "RUN_MAC_STATS_FROM_GCOUNT: unsupported mod ${params.mod}."
    }
}

workflow RUN_MAC_MISS_BIN50_SAMPLE {
    main:
    if (!params.mod) {
        error 'params.mod is required (e.g. test_thin, test_common_thin, test_rebulld_lib).'
    }
    if (!params.output_dir || !params.job) {
        error 'params.output_dir and params.job are required.'
    }

    def proc = "${params.output_dir}/${params.job}/process/${params.mod}/variant"

    if (params.mod in ['test_thin', 'test_common_thin']) {
        def subgenomes = ['A', 'B', 'D', 'Others']
        def ch_gcount = channel.from(
            subgenomes.collect { sg ->
                tuple(sg, "sub_${sg}", file("${proc}/${sg}.info.gcount"))
            }
        )
        def ch_vmis = channel.from(
            subgenomes.collect { sg ->
                tuple(sg, "sub_${sg}", file("${proc}/${sg}.info.vmiss"))
            }
        )
        variant_mac_missing_reg_bin50_sample(ch_gcount.combine(ch_vmis, by: [0, 1]))
    } else if (params.mod == 'test_rebulld_lib') {
        def ch_gcount = channel.of(tuple('chr002', '2', file("${proc}/chr002.info.gcount")))
        def ch_vmis = channel.of(tuple('chr002', '2', file("${proc}/chr002.info.vmiss")))
        variant_mac_missing_reg_bin50_sample(ch_gcount.combine(ch_vmis, by: [0, 1]))
    } else {
        error "RUN_MAC_MISS_BIN50_SAMPLE: unsupported mod ${params.mod}."
    }
}

workflow RUN_MAC_DIST_LOG_REDRAW {
    main:
    if (!params.output_dir || !params.job) {
        error 'params.output_dir and params.job are required.'
    }

    def base = "${params.output_dir}/${params.job}/stats"
    def entries = [
        ['test_thin', 'A', 'sub_A'],
        ['test_thin', 'B', 'sub_B'],
        ['test_thin', 'D', 'sub_D'],
        ['test_thin', 'Others', 'sub_Others'],
        ['test_common_thin', 'A', 'sub_A'],
        ['test_common_thin', 'B', 'sub_B'],
        ['test_common_thin', 'D', 'sub_D'],
        ['test_common_thin', 'Others', 'sub_Others'],
        ['test_rebulld_lib', 'chr002', '2'],
    ]

    def ch = channel.from(
        entries.collect { plink_mod, id, chr ->
            tuple(
                plink_mod,
                id,
                chr,
                file("${base}/${plink_mod}/info/${id}.variant.mac.info.tsv", checkIfExists: true),
            )
        }
    )

    variant_mac_dist_log_redraw(ch)
}

workflow RUN_CHR_VARIANT_COUNTS {
    main:
    if (!params.output_dir || !params.job || !params.mod) {
        error 'params.output_dir, params.job, and params.mod are required.'
    }

    def process_root = params.process_dir ?: "${params.output_dir}/${params.job}/process/${params.mod}"
    log.info "Counting variants under process dir: ${process_root}"

    report_plink_chr_variant_counts(channel.of(process_root))
}

workflow RUN_CHR_VARIANT_COMPARE {
    main:
    if (!params.output_dir || !params.job) {
        error 'params.output_dir and params.job are required.'
    }

    def thin_root = params.thin_process_dir ?: "${params.output_dir}/${params.job}/process/test_thin"
    def common_root = params.common_process_dir ?: "${params.output_dir}/${params.job}/process/test_common_thin"
    def thin_by_ref = file("${thin_root}/info/test_thin.chr_variant_counts.by_ref.tsv")
    def common_by_ref = file("${common_root}/info/test_common_thin.chr_variant_counts.by_ref.tsv")
    def output_prefix = params.output_prefix ?: 'test_thin_vs_test_common_thin'

    if (!thin_by_ref.exists()) {
        error "Missing ${thin_by_ref}; run RUN_CHR_VARIANT_COUNTS for test_thin first."
    }
    if (!common_by_ref.exists()) {
        error "Missing ${common_by_ref}; run RUN_CHR_VARIANT_COUNTS for test_common_thin first."
    }

    log.info "Thin process:   ${thin_root}"
    log.info "Common process: ${common_root}"

    plot_thin_common_chr_variant_compare(
        channel.of(tuple(thin_root, common_root, thin_by_ref, common_by_ref, output_prefix))
    )
}

workflow RUN_REBUILLD_LIB_STATS {
    main:
    if (params.mod != 'test_rebulld_lib') {
        error "RUN_REBUILLD_LIB_STATS: params.mod must be test_rebulld_lib (got: ${params.mod})."
    }
    if (!params.output_dir || !params.job) {
        error 'params.output_dir and params.job are required.'
    }

    def proc = "${params.output_dir}/${params.job}/process/${params.mod}"
    def id = 'chr002'
    def chr = '2'

    TEST_PLINK_STATS(
        channel.empty(),
        channel.empty(),
        channel.of(tuple(id, chr, file("${proc}/sample/${id}.info.smiss"))),
        channel.of(tuple(id, chr, file("${proc}/sample/${id}.info.scount"))),
        channel.of(tuple(id, chr, file("${proc}/variant/${id}.info.vmiss"))),
        channel.of(tuple(id, chr, file("${proc}/variant/${id}.info.gcount"))),
        channel.of(tuple(id, chr, file("${proc}/variant/${id}.info.afreq"))),
        channel.of(tuple(id, chr, file("${proc}/variant/${id}.info.hardy"))),
    )
}
