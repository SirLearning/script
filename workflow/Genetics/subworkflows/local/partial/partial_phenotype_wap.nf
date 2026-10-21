nextflow.enable.dsl=2

/*
 * WAP phenotype matrix for association runs: align a raw phenotype table to
 * assoc.keep and publish a plink2 --pheno file plus per-trait coverage QC.
 *
 * Launch: partial_router.nf --partial_task phenotype_wap --mod phenotype_wap \
 *           --phenotype_wap_input ... [--phenotype_wap_traits Trait1,Trait2]
 */

include { wap_phenotype_matrix } from '../../../modules/local/static/phenotype/phenotype_wap.nf'

workflow RUN_PHENOTYPE_WAP {
    main:
    if (!params.output_dir || !params.job) {
        error 'RUN_PHENOTYPE_WAP: params.output_dir and params.job are required.'
    }
    if (!params.user_dir) {
        error 'RUN_PHENOTYPE_WAP: params.user_dir is required (Conda stats env).'
    }
    if (!params.mod) {
        error 'RUN_PHENOTYPE_WAP: params.mod is required (publish folder).'
    }
    if (!params.phenotype_wap_input) {
        error 'RUN_PHENOTYPE_WAP: --phenotype_wap_input is required (wide phenotype table).'
    }

    def keep_dir = params.cohort_keep_dir ?: "${params.output_dir}/${params.job}/process/sample_manifest"
    def keep = file(params.phenotype_wap_keep ?: "${keep_dir}/assoc.keep", checkIfExists: true)
    def phenotype = file(params.phenotype_wap_input, checkIfExists: true)
    def prefix = params.phenotype_wap_prefix ?: "${params.mod}"

    log.info "${params.c_green}phenotype_wap from ${phenotype}; keep ${keep}${params.c_reset}"
    wap_phenotype_matrix(phenotype, keep, prefix)
}
