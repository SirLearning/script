nextflow.enable.dsl=2

/*
 * Full-cohort MAPQ + per-bin depth scan for aneuploidy_check.
 *
 * Reads final.ALL.taxaBamMap.txt, scans each BAM once, publishes:
 *   {output_dir}/{job}/stats/{mod}/info/mapq_check/{sample}.mapq_depth_bin.tsv
 *
 * Launch via partial_router.nf --partial_task aneuploidy_mapq
 */

include { bam_mapq_depth_bin_sample } from '../../../modules/local/germplasm/bam_mapq_depth_bin.nf'

workflow RUN_ANEUPLOIDY_MAPQ {
    main:
    if (params.mod != 'aneuploidy_check') {
        error "RUN_ANEUPLOIDY_MAPQ: params.mod must be aneuploidy_check (got: ${params.mod})."
    }
    if (!params.output_dir || !params.job) {
        error 'RUN_ANEUPLOIDY_MAPQ: params.output_dir and params.job are required.'
    }
    if (!params.taxa_bam_map) {
        error 'RUN_ANEUPLOIDY_MAPQ: params.taxa_bam_map is required.'
    }

    def publish_dir = "${params.output_dir}/${params.job}/stats/${params.mod}/info/mapq_check"
    def skip_existing = params.mapq_skip_existing != false

    def existing_samples = [] as Set
    if (skip_existing) {
        def publish_path = file(publish_dir)
        if (publish_path.exists()) {
            publish_path.eachFileMatch(~/.+\.mapq_depth_bin\.tsv$/) { f ->
                existing_samples << f.name.replace('.mapq_depth_bin.tsv', '')
            }
        }
        log.info "aneuploidy_mapq: skip_existing enabled under ${publish_dir} (${existing_samples.size()} TSVs)"
    }

    def sample_bam = channel
        .fromPath(params.taxa_bam_map, checkIfExists: true)
        .splitCsv(sep: '\t', header: false)
        .filter { row -> row.size() >= 3 && row[0] != 'Taxa' }
        .map { row ->
            def sample = row[0].toString().trim()
            def bam_path = row[2].toString().trim()
            tuple(sample, file(bam_path))
        }
        .filter { _sample, bam -> bam.exists() }
        .filter { sample, _bam -> !skip_existing || !existing_samples.contains(sample) }

    bam_mapq_depth_bin_sample(sample_bam)
}
