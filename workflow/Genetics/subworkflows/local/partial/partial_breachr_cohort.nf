nextflow.enable.dsl=2

/*
 * Full-cohort breachr segment calling on all published MAPQ-depth bin TSVs.
 *
 * Publishes per-sample segments under info/breachr_cohort/segments/ and a merged
 * breachr_segments_all.tsv. Reconcile / cap plots for flagged samples remain
 * on partial_task aneuploidy_breachr (breachr_check/).
 *
 * Launch via partial_router.nf --partial_task aneuploidy_breachr_cohort
 */

include {
    breachr_segment_sample
    breachr_merge_segments_cohort
} from '../../../modules/local/germplasm/breachr_segment.nf'

workflow RUN_BREACHR_COHORT {
    main:
    if (params.mod != 'aneuploidy_check') {
        error "RUN_BREACHR_COHORT: params.mod must be aneuploidy_check (got: ${params.mod})."
    }
    if (!params.output_dir || !params.job) {
        error 'RUN_BREACHR_COHORT: params.output_dir and params.job are required.'
    }

    def stats_base = "${params.output_dir}/${params.job}/stats/${params.mod}"
    def mapq_dir = "${stats_base}/info/mapq_check"
    def cohort_seg_dir = "${stats_base}/info/${params.breachr_info_subdir ?: 'breachr_cohort'}/segments"
    def skip_existing = params.breachr_skip_existing != false

    def existing_samples = [] as Set
    if (skip_existing) {
        def publish_path = file(cohort_seg_dir)
        if (publish_path.exists()) {
            publish_path.eachFileMatch(~/^.+\.breachr_segments\.tsv$/) { f ->
                existing_samples << f.name.replace('.breachr_segments.tsv', '')
            }
        }
        log.info "aneuploidy_breachr_cohort: skip_existing under ${cohort_seg_dir} (${existing_samples.size()} segment TSVs)"
    }

    def mapq_inputs = channel
        .fromPath("${mapq_dir}/*.mapq_depth_bin.tsv")
        .map { mapq_tsv ->
            def sample = mapq_tsv.name.replace('.mapq_depth_bin.tsv', '')
            tuple(sample, mapq_tsv)
        }
        .filter { sample, _mapq_tsv -> !skip_existing || !existing_samples.contains(sample) }

    breachr_segment_sample(mapq_inputs)

    breachr_merge_segments_cohort(
        channel.value(cohort_seg_dir),
        breachr_segment_sample.out.segments.map { _sample, seg -> seg }.collect(),
    )
}
