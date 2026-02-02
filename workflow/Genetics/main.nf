nextflow.enable.dsl=2

// --- Include workflow modules ---
include { plink_processor as PLINK_PROCESSOR } from './genotype/processor.nf'
include { plink_stats as PLINK_STATS } from './genotype/stats.nf'
include { test_plink as TEST_PLINK } from './genotype/processor.nf'
include { database as DATABASE } from './genotype/database.nf'
include { kinship as KINSHIP } from './dynamic/kinship.nf'
include { population_structure as POPULATION_STRUCTURE } from './dynamic/ps.nf'
include { GWAS } from './static/gwas/gwas.nf'
include { HAIL } from './genotype/hail.nf'
// --- Include util functions ---
include { getJobConfig } from './genotype/utils.nf'
include { getVcfIdFromPath } from './genotype/utils.nf'
include { getRefV1SubChr } from './genotype/utils.nf'

// --- Header ---
def header() {
    return """
    ${params.c_blue}=======================================================${params.c_reset}
    ${params.c_green}      GENETIC PIPELINE (Static/Dynamic/Genotype)      ${params.c_reset}
    ${params.c_blue}=======================================================${params.c_reset}
    """.stripIndent()
}

// --- Help / usage ---
def helpMessage() {
    log.info header()
    log.info """
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf --home_dir <dir> --src_dir <dir> --mod <string> --job <string> [options]

    Required params:
        --home_dir <dir>      Home directory for predefined modules
        --src_dir <dir>       Source directory for scripts and resources
        --mod <string>        Predefined module name for input data (all, gwas)
        --job <string>        Job name/ID (default: genotype_<mod>)
        --tool <string>       Tool to use (plink, hail). Default: plink

    Options:
        --output_dir <dir>    Directory to output results
        --help                This usage statement

    Examples:
        nextflow run main.nf --home_dir /path/to/home --src_dir /path/to/src --mod all --job test
        // 20260130.1
        nextflow run /data/home/tusr1/git/script/workflow/Genetics/main.nf \
            --home_dir /data/home/tusr1/01projects/vmap4 \
            --user_dir /data/home/tusr1 \
            --src_dir /data/home/tusr1/git/script/src \
            --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
            --mod plink_test \
            --job test_plink
        // 20260130.2 - run MAO=1 chr002 only
        nextflow run /data/home/tusr1/git/script/workflow/Genetics/main.nf \
            --home_dir /data/home/tusr1/01projects/vmap4 \
            --user_dir /data/home/tusr1 \
            --src_dir /data/home/tusr1/git/script/src \
            --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
            --mod plink_process \
            --job rebuild
    
    screen prefix commands:
        // test_plink
        screen -dmS test_plink bash -c "\
        cd /data1/dazheng_tusr1/01work && \
        source ~/.bashrc && conda activate run && \
        "
        // rebuild
        screen -dmS rebuild bash -c "\
        cd /data/home/tusr1/01projects/vmap4/05reliable.lib/02rebuild.chr002 && \
        source ~/.bashrc && conda activate run && \
        "
    """.stripIndent()
}

// --- Workflow ---
workflow {
    def input_out = check_input()
    def ch_vcf = input_out.vcf
    // --- Main workflow execution ---
    log.info "${params.c_purple}Process all vcf with plink / plink2 tools.${params.c_reset}"
    if (params.mod == 'plink_process') {
        log.info "${params.c_yellow}Using PLINK PROCESS module.${params.c_reset}"
        PLINK_PROCESSOR(ch_vcf)
    } else if (params.mod == 'plink_test') {
        log.info "${params.c_yellow}Using PLINK TEST module.${params.c_reset}"
        TEST_PLINK(ch_vcf)
    } else {
        log.info "${params.c_yellow}Using default VMap4 v1 plink workflow module.${params.c_reset}"
        def vmap4_v1_plink_out = vmap4_v1_plink(ch_vcf)
    }
    // --- Completion Handler ---
    workflow.onComplete {
        log.info "-\033[2m--------------------------------------------------\033[0m-"
        log.info "Pipeline completed at: $workflow.complete"
        log.info "Duration             : $workflow.duration"
        log.info "Success              : $workflow.success"
        log.info "Exit Status          : $workflow.exitStatus"
        log.info "-\033[2m--------------------------------------------------\033[0m-"
    }
}

workflow check_input {
    main:
    // Show help message
    if (params.help) {
        helpMessage()
        exit 0
    }
    // Print header and summary
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
    // --- Input Handling ---
    def job_config = getJobConfig(params.job, params.home_dir)
    // Build input channel of VCF tuples: [ val(id), path(vcf) ]
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
    // chr filter for manual edit
    def man_chr_filter = getRefV1SubChr("ALL")
    // man_chr_filter -= ["0","43","44"] // Other chromosomes
    // man_chr_filter -= ["1","2","7","8","13","14","19","20","25","26","31","32","37","38"] // A genome
    // man_chr_filter -= ["3","4","9","10","15","16","21","22","27","28","33","34","39","40"] // B genome
    // man_chr_filter -= ["5","6","11","12","17","18","23","24","29","30","35","36","41","42"] // D genome
    ch_vcf = ch_vcf.filter { vcf_tuple ->
        vcf_tuple[1] in man_chr_filter
    }
    // Optional chr param filtering, can use multiple chromosomes separated by comma, like --chr 1,2,3
    if (params.chr) {
        def chr_filter_list = params.chr.toString().split(',').collect { chr -> chr.trim() }
        log.info "${params.c_green}Filtering for chromosome(s):${params.c_reset} ${chr_filter_list}"
        ch_vcf = ch_vcf.filter { vcf_tuple -> vcf_tuple[1].toString() in chr_filter_list }
    }

    emit:
    vcf = ch_vcf
}

workflow vmap4_v1_plink {
    take:
    ch_vcf

    main:
    def processor_out = PLINK_PROCESSOR(ch_vcf)
    def stats_out = PLINK_STATS(
        processor_out.smiss, 
        processor_out.vmiss,
        processor_out.scount,
        processor_out.gcount,
        processor_out.afreq,
        processor_out.hardy,
        processor_out.popdep)
    

    emit:
    out = processor_out
}

