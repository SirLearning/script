nextflow.enable.dsl = 2

// Association pipeline: VCF filter -> GWAS -> Post-analysis
// Inputs
//  - params.vcf: path to input VCF[.gz]
//  - params.pheno: phenotype table (CSV/TSV). Must contain a sample ID column (default: IID) and the target trait column.
//  - params.trait: phenotype column name to analyze (required)
//  - params.covariates: optional covariates table (CSV/TSV) with IID and covariate columns
//  - params.filter_tool: 'bcftools' (default) or 'gatk' for variant filtration
//  - params.models: list of GWAS models to run [GLM, MLM, FarmCPU, Blink]
//  - params.gwas_engine: 'rmvp' (default) or 'gapit'
//  - params.outdir: output directory (default: 'Association/output')
//  - params.gff: optional GFF/GTF for annotation in post step

// Filter options (when using bcftools)
params.qual       = params.qual       ?: 30
params.maf        = params.maf        ?: 0.05
params.max_miss   = params.max_miss   ?: 0.2
params.region     = params.region     ?: null
params.filter_tool= params.filter_tool?: 'bcftools'

// GWAS options
params.models       = params.models       ?: ['GLM','MLM','FarmCPU','Blink']
params.gwas_engine  = params.gwas_engine  ?: 'rmvp'  // 'rmvp' or 'gapit'
params.trait        = params.trait        ?: null
params.pheno_id_col = params.pheno_id_col ?: 'IID'   // sample id column in pheno/covar
params.covar_names  = params.covar_names  ?: []      // subset/ordering of covariates; default use all columns except IID and trait
params.plink_threads= params.plink_threads?: 4

// Post options
params.gff         = params.gff         ?: null
params.outdir      = params.outdir      ?: 'Association/output'

// Validate required params
if (!params.vcf)   exit 1, 'Missing required parameter: --vcf PATH_TO_VCF'
if (!params.pheno) exit 1, 'Missing required parameter: --pheno PATH_TO_PHENOTYPE_TABLE'
if (!params.trait) exit 1, 'Missing required parameter: --trait PHENOTYPE_COLUMN_NAME'

// Include subworkflows
include { FILTER_VCF } from './filter/filter.nf'
include { PHENOTYPE } from './pre/phenotype.nf'
include { GWAS       } from './run/gwas.nf'
include { POST       } from './post/post.nf'

workflow {
    main:
    // Create input channels
    ch_vcf   = Channel.fromPath(params.vcf, checkIfExists: true)
    ch_pheno = Channel.fromPath(params.pheno, checkIfExists: true)
    ch_covar = params.covariates ? Channel.value(params.covariates) : Channel.value('')

    // Optional: preprocess phenotype (BLUE/BLUP/summary)
    ch_pheno_for_gwas = (params.pheno_process ?: false) ? PHENOTYPE(ch_pheno) : ch_pheno

    // Step 1: Filter VCF
    filtered_vcf_ch = FILTER_VCF(ch_vcf)

    // Step 2: GWAS (supports PLINK baseline + rMVP/GAPIT models)
    gwas_out_ch = GWAS(
        filtered_vcf_ch,
        ch_pheno_for_gwas,
        ch_covar
    )

    // Step 3: Post-analysis (plots, summary, optional annotation)
    POST(
        gwas_out_ch
    )
}
