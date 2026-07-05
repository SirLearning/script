nextflow.enable.dsl=2

include { header; helpMessage; check_input } from './subworkflows/local/entry/genetics_helpers.nf'
include {
    RUN_V1_PLINK
    RUN_TEST_PLINK
    RUN_TEST_PLINK_CAMP
    RUN_TEST_COMMON_THIN
    RUN_TEST_COMMON_ONLY
} from './subworkflows/local/plink/plink_genotype_modes.nf'
include { RUN_WHEAT_INTEGRATED } from './subworkflows/local/wheat/wheat_integrated_study.nf'

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
        } else if (params.mod == 'test_common_only') {
            log.info "${params.c_yellow}Using PLINK MAF-ONLY COMMON TEST module (test_common_only).${params.c_reset}"
            RUN_TEST_COMMON_ONLY(ch_vcf)
        } else {
            log.error "${params.c_red}Unknown or unsupported params.mod: '${params.mod}'${params.c_reset}"
            log.info """
            Supported VCF/PLINK mods:
              - v1_plink
              - test_thin
              - test_camp
              - test_common_thin
              - test_common_only

            Wheat integrated mods (no VCF channel):
              - params.mod starting with wheat_ (e.g. wheat_pca_tsne, wheat_gwas)

            See workflow/Genetics/docs/GENETICS_WORKFLOW.md and doc/PROGRESS_README.md for run examples and replay log.
            """.stripIndent()
            exit 1
        }
    }
}
