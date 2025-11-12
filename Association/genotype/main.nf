#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --- Include modules ---
include { FILTER_VCF }           from './process/filter.nf'
include { GENOTYPE_STATS }       from './process/stats.nf'
include { KINSHIP_ANALYSIS }     from './analysis/kinship.nf'
include { POPULATION_STRUCTURE } from './analysis/ps.nf'

// --- Help / usage ---
def usage() {
    log.info """
    ================================================
    Genotype pipeline entry (Association/genotype)
    ================================================
    Required params:
        --home_dir <dir>      Home directory for predefined modules
        --mod <string>        Predefined module name for input data
        --job <string>        Job name/ID (default: genotype_<mod>)

    Common params (see nextflow.config for more):
        --outdir <dir>         Output base directory (default: results/genotype)
        --filter_tool          bcftools|gatk (default: bcftools)
        --pca_k <int>          Number of PCs for PCA (default: 5)
        --enable_genotype_stats  true|false (default: false)

    Example:
      nextflow run Association/genotype/main.nf --mod test_first --home_dir /data/home/user/project --outdir results/genotype --enable_genotype_stats true

    Examples using screen:
        screen -dmS genotype_pipe bash -c "cd /data/home/dazheng/01projects/vmap4/04chr1Geno && source ~/.bashrc && conda activate stats && nextflow ..."
    """
}

params.home_dir = null
params.mod = null
params.job = null

def getModConfig(mod, home_dir) {
    def modConfigs = [
        "chr1": [
            vcf_file: "${home_dir}/00data/06vcf/01chr1/chr001.vcf.gz"
        ],
        "test_first": [
            vcf_file: "${home_dir}/00data/06vcf/02test/chr001.first.1M.vcf.gz"
        ],
        "test_mid": [
            vcf_file: "${home_dir}/00data/06vcf/02test/chr001.mid.1M.vcf.gz"
        ],
        "test_last": [
            vcf_file: "${home_dir}/00data/06vcf/02test/chr001.last.1M.vcf.gz"
        ],
        "vmap4": [
            vcf_path: "${home_dir}/00data/06vcf/03run/"
        ]
    ]

    if (!modConfigs.containsKey(mod)) {
        log.error "Unknown mod specified: ${mod}"
        System.exit(1)
    }

    return modConfigs[mod]
}

def getJobConfig(job) {
    def jobConfigs = [
        "filter": [
            
        ]
    ]

    if (!jobConfigs.containsKey(job)) {
        log.error "Unknown job specified: ${job}"
        System.exit(1)
    }

    return jobConfigs[job]
}

workflow {
    if (params.help) { usage(); System.exit(0) }

    // Build input channel of tuples: [ val(meta), path(vcf) ]
    Channel ch_vcf

    def modConfig = getModConfig(params.mod, params.home_dir)
    if (modConfig.vcf_file) {
        def f = file(modConfig.vcf_file)
        if (!f.exists()) {
            log.error "Mod VCF file not found: ${f}"
            System.exit(1)
        }
        def id = params.mod
        ch_vcf = Channel.of( [ [id: id], f ] )
    } else if (modConfig.vcf_path) {
        def pattern = "${modConfig.vcf_path}/*.vcf.gz"
        ch_vcf = Channel.fromPath(pattern, checkIfExists: true)
            .map { vcf -> [ [id: vcf.baseName.replaceAll(/\.vcf(\.gz)?$/, '')], vcf ] }
    } else {
        usage()
        log.error "No valid input found for mod: ${params.mod}"
        System.exit(1)
    }

    GENOTYPE_PIPELINE(ch_vcf)
}

// --- Genotype Processing Workflow ---
workflow GENOTYPE_PIPELINE {
    take:
        ch_vcf // channel: [ val(meta), path(vcf) ]

    main:
        // 1. Filter/assess VCF file
        FILTER_VCF(ch_vcf)

        // 2. Kinship analysis
        KINSHIP_ANALYSIS(FILTER_VCF.out.vcf)

        // 3. Population structure analysis
        POPULATION_STRUCTURE(FILTER_VCF.out.vcf)

        // 4. Optional genotype summary stats (guarded by params.enable_genotype_stats)
        GENOTYPE_STATS(FILTER_VCF.out.vcf)

    emit:
        filtered_vcf   = FILTER_VCF.out.vcf
        kinship_matrix = KINSHIP_ANALYSIS.out.kinship_matrix
        pca_results    = POPULATION_STRUCTURE.out.pca_results
        geno_stats     = GENOTYPE_STATS.out.geno_stats
}


