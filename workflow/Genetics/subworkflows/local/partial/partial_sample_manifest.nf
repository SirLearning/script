nextflow.enable.dsl=2

/*
 * Stage 0/1 cohort contract: WAP association panel vs WGS-core, plus RefSeq v1.0
 * FASTA/annotation inventory. Does not model phenotypes.
 *
 * Launch: partial_router.nf --partial_task sample_manifest --mod sample_manifest
 */

include { sample_manifest_build } from '../../../modules/local/genotype/stats/stats_cohort.nf'

workflow RUN_SAMPLE_MANIFEST {
    main:
    if (!params.output_dir || !params.job) {
        error 'RUN_SAMPLE_MANIFEST: params.output_dir and params.job are required.'
    }
    if (!params.user_dir) {
        error 'RUN_SAMPLE_MANIFEST: params.user_dir is required (Conda stats env).'
    }
    if (!params.home_dir) {
        error 'RUN_SAMPLE_MANIFEST: params.home_dir is required (taxaBamMap + RefSeq v1 FASTA).'
    }
    if (!params.mod) {
        error 'RUN_SAMPLE_MANIFEST: params.mod is required (publish folder, e.g. sample_manifest).'
    }

    def group_file = file(
        params.sample_group_file ?: "${params.output_dir}/meta_data/sample_group.txt",
        checkIfExists: true,
    )
    def taxa_dir = params.taxa_bam_map_dir ?: "${params.home_dir}/00data/05taxaBamMap"
    log.info "${params.c_green}sample_manifest from ${group_file}; taxa maps ${taxa_dir}${params.c_reset}"
    sample_manifest_build(group_file, taxa_dir, params.home_dir)
}
