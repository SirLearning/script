nextflow.enable.dsl=2

// --- Include modules ---
// Avoid conflict with `process` keyword by aliasing the imported workflow
include { process as PROCESS } from './process/process.nf'
include { stats as STATS } from './process/stats.nf'
include { access as ASSESS } from './process/assess.nf'
include { annotate as ANNOTATE } from './process/annotate.nf'
include { kinship as KINSHIP } from './analysis/kinship.nf'
include { population_structure as POPULATION_STRUCTURE } from './analysis/ps.nf'

// --- Required Parameters ---
params.mod = null
params.job = null

// --- Optional Parameters ---
params.process_dir = null
params.filter_dir = null
params.assess_dir = null
params.stats_dir = null
params.kinship_dir = null
params.ps_dir = null

workflow {
    if (params.help) { usage(); System.exit(0) }

    // --- Input Handling ---
    def jobConfig = getJobConfig(params.job, params.home_dir)
    // Resolve output directory
    def process_dir_resolved = params.process_dir ?: "${params.output_dir}/${params.job}/process"
    def filter_dir_resolved = params.filter_dir ?: "${params.output_dir}/${params.job}/filter"
    def assess_dir_resolved = params.assess_dir ?: "${params.output_dir}/${params.job}/assess"
    def stats_dir_resolved = params.stats_dir ?: "${params.output_dir}/${params.job}/stats"
    def kinship_dir_resolved = params.kinship_dir ?: "${params.output_dir}/${params.job}/kinship"
    def ps_dir_resolved = params.ps_dir ?: "${params.output_dir}/${params.job}/ps"

    // Build input channel of tuples: [ val(meta), path(vcf) ]
    def ch_vcf

    if (jobConfig.vcf_file) {
        def f = file(jobConfig.vcf_file)
        log.info "Using VCF file: ${f}"
        if (!f.exists()) {
            log.error "Mod VCF file not found: ${f}"
            System.exit(1)
        }
        def id = f.baseName.replaceAll(/\.vcf(\.gz)?$/, '')
        def meta = [id: id]
        // Emit one tuple [meta, vcf] as a proper channel item
        ch_vcf = Channel.of([ meta, f ])
        // Debug view to confirm tuple structure
        ch_vcf.view { item -> "DEBUG ch_vcf single-file -> ${item}" }
    } else if (jobConfig.vcf_dir) {
        def pattern = "${jobConfig.vcf_dir}/*.vcf.gz"
        ch_vcf = Channel.fromPath(pattern, checkIfExists: true)
            .map { vcf -> [ [id: vcf.baseName.replaceAll(/\.vcf(\.gz)?$/, '')], vcf ] }
        ch_vcf.view { item -> "DEBUG ch_vcf multi-file -> ${item}" }
    } else {
        usage()
        log.error "No valid input found for job: ${params.job}"
        System.exit(1)
    }

    def job_config = getJobConfig(params.mod ?: "all")
    def combined_ch = ch_vcf.combine(Channel.value(job_config))

    if (params.mod == "all") {
        PROCESS(combined_ch)
        ASSESS(PROCESS.out.vcf, job_config)
        STATS(PROCESS.out.vcf, job_config)
        KINSHIP(PROCESS.out.vcf, job_config)
        POPULATION_STRUCTURE(PROCESS.out.vcf, job_config)
    }
}

// --- Help / usage ---
def usage() {
    log.info """
    ================================================
    Genotype pipeline entry (Association/genotype)
    ================================================
    Required params:
        --home_dir <dir>      Home directory for predefined modules
        --src_dir <dir>       Source directory for scripts and resources
        --mod <string>        Predefined module name for input data
        --job <string>        Job name/ID (default: genotype_<mod>)

    Common params (see nextflow.config for more):

    Example: test on 107
        nextflow run /data/dazheng/git/script/Association/genotype/main.nf \
            --home_dir /data/dazheng/01projects/vmap4 \
            --src_dir /data/dazheng/git/script/src \
            --output_dir /data/dazheng/01projects/vmap4/05ana.geno/01chr1.test \
            --mod all --job test

    Examples using screen:
        screen -dmS genotype_pipe bash -c "cd /data/home/dazheng/01projects/vmap4/04chr1Geno && source ~/.bashrc && conda activate stats && nextflow ..."
        screen -dmS test bash -c "cd /data/dazheng/01projects/vmap4/05chr1Geno/02test_wf && source ~/.bashrc && conda activate stats && nextflow ..."
    """
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
        ]
    ]

    if (!jobConfigs.containsKey(job)) {
        log.error "Unknown job specified: ${job}"
        System.exit(1)
    }

    return jobConfigs[job]
}





