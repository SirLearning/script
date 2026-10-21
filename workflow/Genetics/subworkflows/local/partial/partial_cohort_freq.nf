nextflow.enable.dsl=2

/*
 * Cohort-restricted allele frequency for merged subgenome pfiles.
 * Turns the sample_manifest keep files into per-cohort afreq/acount/vmiss so
 * dual_af and variant_sets no longer need hand-made frequency tables.
 *
 * Launch: partial_router.nf --partial_task cohort_freq --mod <pfile mod>
 */

include {
    hasMergedSubgenomeTestPfiles
    listMergedSubgenomeTestPfileTuples
} from '../../../modules/local/infra/infra_plink_reuse.nf'
include { plink2_cohort_freq } from '../../../modules/local/genotype/processor/processor_cohort_freq.nf'

// Cohort name -> keep file written by partial_task sample_manifest
// (assoc.keep / wgs_core.keep).
def cohortKeepFileName(cohort) {
    if (cohort == 'wgs_core') {
        return 'wgs_core.keep'
    }
    if (cohort == 'assoc') {
        return 'assoc.keep'
    }
    return "${cohort}.keep"
}

workflow RUN_COHORT_FREQ {
    main:
    if (!params.output_dir || !params.job) {
        error 'RUN_COHORT_FREQ: params.output_dir and params.job are required.'
    }
    if (!params.user_dir) {
        error 'RUN_COHORT_FREQ: params.user_dir is required (Conda stats env).'
    }
    if (!params.mod) {
        error 'RUN_COHORT_FREQ: params.mod is required (pfile + publish folder).'
    }

    def proc_dir = params.process_dir ?: "${params.output_dir}/${params.job}/process/${params.mod}"
    if (!hasMergedSubgenomeTestPfiles(proc_dir)) {
        error "RUN_COHORT_FREQ: missing merged *_test.plink2 under ${proc_dir}"
    }

    def keep_dir = params.cohort_keep_dir ?: "${params.output_dir}/${params.job}/process/sample_manifest"
    def cohorts = "${params.cohort_freq_cohorts}".split(',').collect { it.trim() }.findAll { it }
    if (!cohorts) {
        error 'RUN_COHORT_FREQ: params.cohort_freq_cohorts is empty.'
    }

    def jobs = []
    cohorts.each { cohort ->
        def keep = file("${keep_dir}/${cohortKeepFileName(cohort)}", checkIfExists: true)
        listMergedSubgenomeTestPfileTuples(proc_dir).each { pfile ->
            // Nextflow tuples are immutable, so rebuild rather than append in place.
            jobs << (pfile.toList() + [cohort, keep])
        }
    }

    log.info "${params.c_green}cohort_freq ${cohorts.join(',')} on ${proc_dir}; keep from ${keep_dir}${params.c_reset}"
    plink2_cohort_freq(channel.from(jobs))
}
