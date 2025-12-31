nextflow.enable.dsl=2

// --- Include modules ---
// Avoid conflict with `process` keyword by aliasing the imported workflow
include { processor as PROCESSOR } from './genotype/processor.nf'
include { database as DATABASE } from './genotype/database.nf'
include { stats as STATS } from './genotype/stats.nf'
include { access as ASSESS } from './genotype/assess.nf'
include { annotate as ANNOTATE } from './genotype/annotate.nf'
include { kinship as KINSHIP } from './dynamic/kinship.nf'
include { population_structure as POPULATION_STRUCTURE } from './dynamic/ps.nf'
include { GWAS } from './static/gwas/gwas.nf'

// Hail specific modules
include { HAIL } from './genotype/hail.nf'

workflow {
    if (params.help) { usage(); System.exit(0) }

    // --- Input Handling ---
    def job_config = getJobConfig(params.job, params.home_dir)

    // Build input channel of VCF tuples: [ val(id), path(vcf) ]
    def ch_vcf

    if (job_config.vcf_file) {
        def f = file(job_config.vcf_file)
        log.info "Using VCF file: ${f}"
        if (!f.exists()) {
            log.error "Mod VCF file not found: ${f}"
            System.exit(1)
        }
        def id = f.baseName.replaceAll(/\.vcf(\.gz|\.bgz)?$/, '')
        // Emit one tuple [id, vcf] as a proper channel item
        ch_vcf = Channel.of([ id, f ])
        // Debug view to confirm tuple structure
        ch_vcf.view { item -> "DEBUG ch_vcf single-file -> ${item}" }
    } else if (job_config.vcf_dir) {
        def pattern = "${job_config.vcf_dir}/*.{vcf,vcf.gz,vcf.bgz}"
        ch_vcf = Channel.fromPath(pattern, checkIfExists: true)
            .map { vcf -> [ vcf.baseName.replaceAll(/\.vcf(\.gz|\.bgz)?$/, ''), vcf ] }
        ch_vcf.view { item -> "DEBUG ch_vcf multi-file -> ${item}" }
    } else {
        usage()
        log.error "No valid input found for job: ${params.job}"
        System.exit(1)
    }
}

workflow build_genotype_database {
    // Common processing (Filtering)
    // PROCESSOR handles tool selection internally for filtering
    PROCESSOR(ch_vcf)
    def ch_filtered_vcf = PROCESSOR.out.vcf.map{ meta, id, vcf -> tuple(meta, vcf) }
    ASSESS(ch_filtered_vcf)
    STATS(ch_filtered_vcf)
    KINSHIP(ch_filtered_vcf)
    POPULATION_STRUCTURE(ch_filtered_vcf)
}

workflow association_study {
    if (!params.pheno) { log.error "Missing --pheno"; System.exit(1) }
    if (!params.trait) { log.error "Missing --trait"; System.exit(1) }
    ch_pheno = Channel.fromPath(params.pheno)
    
    if (params.covar) {
        ch_covar = Channel.fromPath(params.covar)
    } else {
        def no_covar = file("NO_FILE")
        if (!no_covar.exists()) no_covar.text = ""
        ch_covar = Channel.of(no_covar)
    }
    
    // If using Hail, we might use raw VCF. If using PLINK, we might want processed VCF.
    // For now, pass raw VCF.
    GWAS(ch_vcf, ch_pheno, ch_covar)
}

workflow hail_platform {
    // Hail Workflow
    HAIL_QC(ch_filtered_vcf)
    HAIL_PCA(ch_filtered_vcf)
    HAIL_KINSHIP(ch_filtered_vcf)
    
    // If phenotype is provided, run GWAS
    if (params.pheno && params.trait) {
        ch_pheno = Channel.fromPath(params.pheno)
        ch_covar = params.covar ? Channel.fromPath(params.covar) : Channel.fromPath("NO_FILE").map{it -> file("NO_FILE")}
        
        HAIL_GWAS(ch_filtered_vcf, ch_pheno, ch_covar, params.trait)
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
        --mod <string>        Predefined module name for input data (all, gwas)
        --job <string>        Job name/ID (default: genotype_<mod>)
        --tool <string>       Tool to use (plink, hail). Default: plink

    Common params (see nextflow.config for more):

    Example: test on 107
        nextflow run /data/dazheng/git/script/workflow/Association/genotype/main.nf \
            --home_dir /data/dazheng/01projects/vmap4 \
            --src_dir /data/dazheng/git/script/src \
            --output_dir /data/dazheng/01projects/vmap4/05ana.geno/01chr1.test \
            --mod all --job test
        nextflow run /data/dazheng/git/script/workflow/Association/genotype/main.nf \
            --home_dir /data/dazheng/01projects/vmap4 \
            --src_dir /data/dazheng/git/script/src \
            --output_dir /data/dazheng/01projects/vmap4/05ana.geno/01chr1.test \
            --mod assess --job test

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





