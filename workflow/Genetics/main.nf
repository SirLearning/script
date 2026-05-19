nextflow.enable.dsl=2

include { header; helpMessage; check_input } from './subworkflows/local/genetics_helpers.nf'
include { RUN_V1_PLINK } from './subworkflows/local/plink_genotype_modes.nf'
include { RUN_TEST_PLINK } from './subworkflows/local/plink_genotype_modes.nf'
include { RUN_TEST_PLINK_CAMP } from './subworkflows/local/plink_genotype_modes.nf'
include { RUN_TEST_COMMON_THIN } from './subworkflows/local/plink_genotype_modes.nf'
include { RUN_WHEAT_INTEGRATED } from './subworkflows/local/wheat_integrated_study.nf'

workflow {
    log.info header()
    if (params.help) {
        helpMessage()
        exit 0
    }
    if (params.mod != null && params.mod.toString().startsWith('wheat_')) {
        log.info "${params.c_purple}Wheat integrated analytics (${params.mod}).${params.c_reset}"
        RUN_WHEAT_INTEGRATED()
    } else {
        def input_out = check_input()
        def ch_vcf = input_out.vcf

        log.info "${params.c_purple}Process all vcf with plink / plink2 tools.${params.c_reset}"
        if (params.mod == 'v1_plink') {
            log.info "${params.c_yellow}Using PLINK PROCESS module.${params.c_reset}"
            RUN_V1_PLINK(ch_vcf)
        } else if (params.mod == 'test_thin') {
            log.info "${params.c_yellow}Using PLINK TEST module.${params.c_reset}"
            RUN_TEST_PLINK(ch_vcf)
        } else if (params.mod == 'test_camp') {
            log.info "${params.c_yellow}Using PLINK TEST CAMP module.${params.c_reset}"
            camp_vmap4_map_tsv = file(params.camp)
            RUN_TEST_PLINK_CAMP(ch_vcf, camp_vmap4_map_tsv)
        } else if (params.mod == 'test_common_thin') {
            log.info "${params.c_yellow}Using PLINK COMMON-THIN TEST module.${params.c_reset}"
            RUN_TEST_COMMON_THIN(ch_vcf)
        } else {
            log.info "${params.c_yellow}No workflow module was chosen.${params.c_reset}"
        }
    }
}
