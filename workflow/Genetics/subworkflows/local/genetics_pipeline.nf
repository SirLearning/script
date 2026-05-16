nextflow.enable.dsl=2

include { RUN_V1_PLINK } from '../../modules/local/genotype_modules.nf'
include { RUN_TEST_PLINK } from '../../modules/local/genotype_modules.nf'
include { RUN_TEST_PLINK_CAMP } from '../../modules/local/genotype_modules.nf'
include { RUN_TEST_COMMON_THIN } from '../../modules/local/genotype_modules.nf'
include { RUN_WHEAT_INTEGRATED } from '../../modules/local/wheat_modules.nf'

include { getJobConfig } from '../../genotype/utils.nf'
include { getVcfIdFromPath } from '../../genotype/utils.nf'
include { getRefV1SubChr } from '../../genotype/utils.nf'

def header() {
    return """
    ${params.c_blue}=======================================================${params.c_reset}
    ${params.c_green}      GENETIC PIPELINE (Static/Dynamic/Genotype)      ${params.c_reset}
    ${params.c_blue}=======================================================${params.c_reset}
    """.stripIndent()
}

def helpMessage() {
    log.info header()
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --home_dir <dir> --src_dir <dir> --mod <string> --job <string> [options]

    Required params:
        --home_dir <dir>      Home directory for predefined modules
        --src_dir <dir>       Source directory for scripts and resources
        --mod <string>        Predefined module name for input data (v1_plink, test_plink, test_thin,
                              test_plink_camp, test_camp, test_common_thin, or wheat_* integrated modes)
        --job <string>        Job name/ID (default: genotype_<mod>)
        --tool <string>       Tool to use (plink, hail). Default: plink

    Wheat integrated (`wheat_*` mods): `--home_dir` / `--src_dir` are unused; require `--output_dir`,
        `--job`, `--user_dir` (Conda stats env), and the task-specific `wheat_*` input params from
        `nextflow.config` / CLI (see `workflow/Genetics/README.md`, Wheat integrated section).

    Options:
        --output_dir <dir>    Directory to output results
        --help                This usage statement

    Examples:
        nextflow run main.nf --home_dir /path/to/home --src_dir /path/to/src --mod test_plink --job test

    Chronological Nextflow command log (append each new run at the end):
        See repo file doc/NF_CMD.md.
    """.stripIndent()
}

workflow check_input {
    main:
    if (params.help) {
        helpMessage()
        exit 0
    }

    log.info header()
    def summary = [:]
    summary['Pipeline Name']  = 'Genetics Pipeline'
    summary['Run Name']       = workflow.runName
    summary['Launch Dir']     = workflow.launchDir
    summary['Work Dir']       = workflow.workDir
    summary['Script Dir']     = workflow.projectDir
    summary['User']           = workflow.userName
    summary['Config Profile'] = workflow.profile
    summary['Home Dir']       = params.home_dir
    summary['Src Dir']        = params.src_dir
    summary['Run Mode']       = params.mod
    summary['Job ID']         = params.job
    summary['Tool']           = params.tool ?: 'plink'
    if(params.output_dir) summary['Output Dir'] = params.output_dir
    log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
    log.info "-\033[2m--------------------------------------------------\033[0m-"

    def job_config = getJobConfig(params.job, params.home_dir)
    if (job_config.vcf_file) {
        def f = file(job_config.vcf_file)
        log.info "${params.c_green}Using VCF file:${params.c_reset} ${f}"
        if (!f.exists()) {
            log.error "Mod VCF file not found: ${f}"
            exit 1
        }
        def (id, chr) = getVcfIdFromPath(f)
        ch_vcf = channel.of([ id, chr, f ])
    } else if (job_config.vcf_dir) {
        def pattern = "${job_config.vcf_dir}/*.{vcf,vcf.gz,vcf.bgz}"
        log.info "${params.c_green}Using VCF dir:${params.c_reset} ${pattern}"
        ch_vcf = channel.fromPath(pattern, checkIfExists: true)
            .map { vcf ->
            def (id, chr) = getVcfIdFromPath(vcf)
            return [ id, chr, vcf ]
            }
    } else {
        helpMessage()
        log.error "No valid input found for job: ${params.job}"
        exit 1
    }

    def man_chr_filter = getRefV1SubChr("ALL")
    ch_vcf = ch_vcf.filter { vcf_tuple ->
        vcf_tuple[1] in man_chr_filter
    }

    if (params.chr) {
        def chr_filter_list = params.chr.toString().split(',').collect { chr -> chr.trim() }
        log.info "${params.c_green}Filtering for chromosome(s):${params.c_reset} ${chr_filter_list}"
        ch_vcf = ch_vcf.filter { vcf_tuple -> vcf_tuple[1].toString() in chr_filter_list }
    }

    emit:
    vcf = ch_vcf
}

workflow GENETICS_PIPELINE {
    main:
    if (params.mod != null && params.mod.toString().startsWith('wheat_')) {
        log.info header()
        if (params.help) {
            helpMessage()
            exit 0
        }
        log.info "${params.c_purple}Wheat integrated analytics (${params.mod}).${params.c_reset}"
        RUN_WHEAT_INTEGRATED()
    } else {
        def input_out = check_input()
        def ch_vcf = input_out.vcf

        log.info "${params.c_purple}Process all vcf with plink / plink2 tools.${params.c_reset}"
        if (params.mod == 'v1_plink') {
            log.info "${params.c_yellow}Using PLINK PROCESS module.${params.c_reset}"
            RUN_V1_PLINK(ch_vcf)
        } else if (params.mod == 'test_plink' || params.mod == 'test_thin') {
            log.info "${params.c_yellow}Using PLINK TEST module.${params.c_reset}"
            RUN_TEST_PLINK(ch_vcf)
        } else if (params.mod == 'test_plink_camp' || params.mod == 'test_camp') {
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
