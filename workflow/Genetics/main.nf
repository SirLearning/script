nextflow.enable.dsl=2

// --- Include workflow modules ---
include { processor as PROCESSOR } from './genotype/processor.nf'
include { plink_processor as PLINK_PROCESSOR } from './genotype/processor.nf'
include { database as DATABASE } from './genotype/database.nf'
include { stats as STATS } from './genotype/stats.nf'
include { assess as ASSESS } from './genotype/assess.nf'
include { annotate as ANNOTATE } from './genotype/annotate.nf'
include { kinship as KINSHIP } from './dynamic/kinship.nf'
include { population_structure as POPULATION_STRUCTURE } from './dynamic/ps.nf'
include { GWAS } from './static/gwas/gwas.nf'
include { HAIL } from './genotype/hail.nf'
// --- Include util functions ---
include { getJobConfig } from './genotype/utils.nf'

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
        nextflow run /data/dazheng/git/script/workflow/Genetics/main.nf --home_dir /data/dazheng/01projects/vmap4 --src_dir /data/dazheng/git/script/src --output_dir /data/dazheng/01projects/vmap4/05ana.geno/01chr1.test --mod all --job test
    """.stripIndent()
}


// --- Workflow ---
workflow {
    def input_out = check_input()
    def ch_vcf = input_out.vcf

    // --- Main workflow execution ---
    log.info "${params.c_purple}Assessing genotype process tools${params.c_reset}"
    assess_geno_tools(ch_vcf)

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
        def id = f.baseName.replaceAll(/\.vcf(\.gz|\.bgz)?$/, '')
        ch_vcf = channel.of([ id, f ])
    } else if (job_config.vcf_dir) {
        def pattern = "${job_config.vcf_dir}/*.{vcf,vcf.gz,vcf.bgz}"
        log.info "${params.c_green}Using VCF dir:${params.c_reset} ${pattern}"
        ch_vcf = channel.fromPath(pattern, checkIfExists: true)
            .map { vcf -> [ vcf.baseName.replaceAll(/\.vcf(\.gz|\.bgz)?$/, ''), vcf ] }
    } else {
        helpMessage()
        log.error "No valid input found for job: ${params.job}"
        exit 1
    }

    emit:
    vcf = ch_vcf
}

workflow assess_geno_tools {
    take:
    ch_vcf

    main:
    def assess_out = ASSESS(ch_vcf)

    emit:
    qc_reports = assess_out.qc_reports
}

workflow build_geno_db {
    take:
    ch_vcf

    main:
    // Common processing (Filtering)
    PROCESSOR(ch_vcf)
    def ch_filtered_vcf = PROCESSOR.out.vcf.map{ id, prefix, vcf, tbi -> tuple(id, vcf) }
    
    emit:
    vcf = ch_filtered_vcf
}

workflow build_pheno_db {
    // Placeholder
}

workflow association_study {
    // Placeholder
}


