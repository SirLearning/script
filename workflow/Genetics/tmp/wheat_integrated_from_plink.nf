nextflow.enable.dsl=2

/*
 * Partial entry: wheat integrated study from existing merged *_test.plink2.
 * Computation: genotype/processor.nf (PLINK2, awk); plots: genotype/stats.nf; GWAS: static/gwas/gwas.nf.
 */
include { WHEAT_STUDY_FROM_PLINK } from '../subworkflows/local/wheat_integrated_study.nf'

workflow {
    if (!(params.mod in ['test_thin', 'test_common_thin'])) {
        error "wheat_integrated_from_plink.nf: params.mod must be test_thin or test_common_thin (got: ${params.mod})."
    }
    if (!params.output_dir || !params.job) {
        error "params.output_dir and params.job are required."
    }
    if (!params.user_dir) {
        error "params.user_dir is required (Conda stats env)."
    }

    def source_mod = params.mod
    def wheat_mod = params.wheat_integrated_mod ?: 'wheat_pca_tsne'
    if (!(wheat_mod in ['wheat_pca_tsne', 'wheat_tagsnp', 'wheat_snp_qc'])) {
        error "params.wheat_integrated_mod must be wheat_pca_tsne, wheat_tagsnp, or wheat_snp_qc (got: ${wheat_mod})."
    }

    WHEAT_STUDY_FROM_PLINK(source_mod, wheat_mod)
}
