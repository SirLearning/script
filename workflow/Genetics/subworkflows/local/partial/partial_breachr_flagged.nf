nextflow.enable.dsl=2

/*
 * Run breachr segment calling on mosdepth aneuploidy-flagged samples that already
 * have published MAPQ-depth bin TSVs, then reconcile vs mosdepth chr flags.
 *
 * Launch via partial_router.nf --partial_task aneuploidy_breachr
 */

include { breachr_segment_sample; breachr_reconcile_flagged; breachr_reconcile_plots; breachr_mapq_depth_cap_plot } from '../../../modules/local/germplasm/breachr_segment.nf'

workflow RUN_BREACHR_FLAGGED {
    main:
    if (params.mod != 'aneuploidy_check') {
        error "RUN_BREACHR_FLAGGED: params.mod must be aneuploidy_check (got: ${params.mod})."
    }
    if (!params.output_dir || !params.job) {
        error 'RUN_BREACHR_FLAGGED: params.output_dir and params.job are required.'
    }

    def stats_base = "${params.output_dir}/${params.job}/stats/${params.mod}"
    def mapq_dir = "${stats_base}/info/mapq_check"
    def flagged_path = file("${stats_base}/info/flagged_samples.tsv", checkIfExists: true)

    def mapq_inputs = flaggedMapqInputChannel(flagged_path, mapq_dir)

    breachr_segment_sample(mapq_inputs)

    def flagged_file = channel.value(flagged_path)
    breachr_reconcile_flagged(
        breachr_segment_sample.out.segments.map { _sample, seg -> seg }.collect(),
        flagged_file,
    )
    breachr_reconcile_plots(
        breachr_reconcile_flagged.out.reconcile_all,
        flagged_file,
        breachr_reconcile_flagged.out.segments_all,
    )

    def reconcile_tsv = breachr_reconcile_flagged.out.reconcile_all
    flaggedMapqInputChannel(flagged_path, mapq_dir)
        .join(breachr_segment_sample.out.segments)
        .combine(reconcile_tsv)
        .combine(flagged_file)
        .map { sample, mapq_tsv, segments_tsv, reconcile_all, flagged_samples ->
            tuple(sample, mapq_tsv, segments_tsv, reconcile_all, flagged_samples)
        }
        .set { cap_plot_inputs }

    breachr_mapq_depth_cap_plot(cap_plot_inputs)
}

def flaggedMapqInputChannel(flaggedPath, mapqDir) {
    channel
        .fromPath(flaggedPath)
        .splitCsv(sep: '\t', header: true)
        .filter { row ->
            def flag = row.aneuploid_flag?.toString()?.trim()?.toLowerCase()
            flag in ['true', '1', 'yes']
        }
        .map { row ->
            def sample = row.sample.toString().trim()
            def mapq_tsv = resolveMapqDepthBinTsv(sample, mapqDir)
            tuple(sample, mapq_tsv)
        }
        .filter { _sample, mapq_tsv -> mapq_tsv.exists() }
}

/*
 * Resolve published MAPQ TSV path for a flagged sample ID.
 * flagged_samples.tsv may use mosdepth names with ``_deduped`` suffix while
 * mapq_check/ publishes taxaBamMap sample IDs without that suffix.
 */
def resolveMapqDepthBinTsv(String sample, String mapqDir) {
    def direct = file("${mapqDir}/${sample}.mapq_depth_bin.tsv")
    if (direct.exists()) {
        return direct
    }
    if (sample.endsWith('_deduped')) {
        def base = sample.replaceAll(/_deduped$/, '')
        def alias = file("${mapqDir}/${base}.mapq_depth_bin.tsv")
        if (alias.exists()) {
            return alias
        }
    }
    return direct
}
