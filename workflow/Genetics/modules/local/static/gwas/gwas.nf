nextflow.enable.dsl=2

include { format_vcf_plink } from '../../genotype/processor/processor_vcf.nf'

include {
    compute_pca
    prepare_pheno_and_covar
    run_plink_glm
    run_rmvp_or_gapit
    run_gwas_gapit
    run_gwas_plink
} from './gwas_legacy.nf'

include {
    plink2_gwas_glm
    gcta_gwas
} from './gwas_plink2.nf'

include {
    plot_gwas_association
    plot_gwas
} from './gwas_plot.nf'

workflow GWAS {
    take:
    ch_vcf
    ch_pheno
    ch_covar

    main:
    // TODO: Ensure vcf_to_plink is included or defined
    // ch_plink = vcf_to_plink(ch_vcf).plink
    // For now, assuming VCF_TO_PLINK is defined elsewhere or will be renamed
    ch_plink = format_vcf_plink(ch_vcf).plink
    ch_eigen = compute_pca(ch_plink).eigen

    // fam path for pheno/covar prep
    ch_fam = ch_plink.map { bed,bim,fam -> fam }
    ch_pcv = prepare_pheno_and_covar(ch_pheno, ch_covar, ch_fam, ch_eigen).pcv

    // Pair genotype tuple with pheno/covar tuple (single combination)
    ch_pair = ch_plink.combine(ch_pcv).map { pl, pcv ->
        // pl is [bed,bim,fam], pcv is [pheno,covar]
        tuple(pl[0], pl[1], pl[2], pcv[0], pcv[1])
    }

    ch_run1  = run_plink_glm(ch_pair)
    ch_run2  = run_rmvp_or_gapit(ch_pair)
    results = ch_run1.mix(ch_run2)

    emit:
    results
}
