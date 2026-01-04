nextflow.enable.dsl=2

// --- Colors ---
c_reset  = "\033[0m";
c_dim    = "\033[2m";
c_black  = "\033[0;30m";
c_green  = "\033[0;32m";
c_yellow = "\033[0;33m";
c_blue   = "\033[0;34m";
c_purple = "\033[0;35m";
c_cyan   = "\033[0;36m";
c_white  = "\033[0;37m";

// --- Include modules ---
include { processor as PROCESSOR } from './genotype/processor.nf'
include { database as DATABASE } from './genotype/database.nf'
include { stats as STATS } from './genotype/stats.nf'
include { access as ASSESS } from './genotype/assess.nf'
include { annotate as ANNOTATE } from './genotype/annotate.nf'
include { kinship as KINSHIP } from './dynamic/kinship.nf'
include { population_structure as POPULATION_STRUCTURE } from './dynamic/ps.nf'
include { GWAS } from './static/gwas/gwas.nf'
include { HAIL } from './genotype/hail.nf'

// --- Header ---
def header() {
    return """
    ${c_blue}=======================================================${c_reset}
    ${c_green}      GENETIC PIPELINE (Static/Dynamic/Genotype)      ${c_reset}
    ${c_blue}=======================================================${c_reset}
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

def getJobConfig(job, home_dir) {
    def jobConfigs = [
        "chr1": [
            vcf_dir: "${params.home_dir}/00data/06vcf/01chr1"
        ],
        "test": [
            vcf_dir: "${params.home_dir}/00data/06vcf/02test"
        ],
        "test_first": [
            vcf_file: "${params.home_dir}/00data/06vcf/02test/chr001.f1M.vcf"
        ],
        "test_middle": [
            vcf_file: "${params.home_dir}/00data/06vcf/02test/chr001.m1M.vcf"
        ],
        "test_last": [
            vcf_file: "${params.home_dir}/00data/06vcf/02test/chr001.l1M.vcf"
        ],
        "vmap4": [
            vcf_dir: "${params.home_dir}/00data/06vcf/03vmap4"  // pending
        ]
    ]

    if (!jobConfigs.containsKey(job)) {
        log.error "Unknown job specified: ${job}"
        System.exit(1)
    }

    return jobConfigs[job]
}

// --- Summary ---
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
summary['Module']         = params.mod
summary['Job ID']         = params.job
summary['Tool']           = params.tool ?: 'plink'
if(params.output_dir) summary['Output Dir'] = params.output_dir

// --- Workflow ---
workflow {
    def ch_vcf = check_input.out.vcf

    // --- Main workflow execution ---
    if (params.genotype_tool == 'hail') {
        log.info "${c_purple}Using Hail platform for genotype processing${c_reset}"
        hail_platform(ch_vcf)
    } else {
        log.info "${c_purple}Using default genotype processing pipeline${c_reset}"
        build_genotype_database(ch_vcf)
    }
}

workflow check_input {
    take:
    ch_vcf

    main:
    // Show help message
    if (params.help) {
        helpMessage()
        exit 0
    }
    // Print header and summary
    log.info header()
    log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
    log.info "-\033[2m--------------------------------------------------\033[0m-"
    // --- Input Handling ---
    def job_config = getJobConfig(params.job, params.home_dir)
    // Build input channel of VCF tuples: [ val(id), path(vcf) ]
    def ch_vcf
    if (job_config.vcf_file) {
        def f = file(job_config.vcf_file)
        log.info "${c_green}Using VCF file:${c_reset} ${f}"
        if (!f.exists()) {
            log.error "Mod VCF file not found: ${f}"
            exit 1
        }
        def id = f.baseName.replaceAll(/\.vcf(\.gz|\.bgz)?$/, '')
        ch_vcf = Channel.of([ id, f ])
    } else if (job_config.vcf_dir) {
        def pattern = "${job_config.vcf_dir}/*.{vcf,vcf.gz,vcf.bgz}"
        log.info "${c_green}Using VCF dir:${c_reset} ${pattern}"
        ch_vcf = Channel.fromPath(pattern, checkIfExists: true)
            .map { vcf -> [ vcf.baseName.replaceAll(/\.vcf(\.gz|\.bgz)?$/, ''), vcf ] }
    } else {
        helpMessage()
        log.error "No valid input found for job: ${params.job}"
        exit 1
    }

    emit:
    vcf = ch_vcf
}

workflow build_genotype_database {
    take:
    ch_vcf

    main:
    // Common processing (Filtering)
    PROCESSOR(ch_vcf)
    def ch_filtered_vcf = PROCESSOR.out.vcf.map{ meta, id, vcf -> tuple(meta, vcf) }
    ASSESS(ch_filtered_vcf)
    STATS(ch_filtered_vcf)
    KINSHIP(ch_filtered_vcf)
    POPULATION_STRUCTURE(ch_filtered_vcf)
    
    emit:
    vcf = ch_filtered_vcf
}

workflow association_study {
    // Placeholder
}

workflow hail_platform {
    take:
    ch_vcf

    main:
    // Hail Workflow
    // Ensure these processes are properly imported or defined
    HAIL(ch_vcf)
    
    if (params.pheno && params.trait) {
        ch_pheno = Channel.fromPath(params.pheno)
        ch_covar = params.covar ? Channel.fromPath(params.covar) : Channel.fromPath("NO_FILE").map{it -> file("NO_FILE")}
        
        HAIL_GWAS(ch_vcf, ch_pheno, ch_covar, params.trait)
    }
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


