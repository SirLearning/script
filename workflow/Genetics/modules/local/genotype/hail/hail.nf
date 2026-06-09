nextflow.enable.dsl=2

include {
    vcf_to_mt_hail
    filter_hail
} from './hail_io.nf'

include {
    hail_qc
    hail_kinship
    hail_pca
} from './hail_stats.nf'

include { hail_normal_gwas } from './hail_gwas.nf'

workflow HAIL {
    take:
    // Expect a channel: [ val(id), path(vcf), val(job_config) ]
    vcf_in

    main:
    // 1. Convert VCF to Hail MatrixTable
    mt = vcf_to_mt_hail(vcf_in.map{ id, vcf, job_config -> tuple(id, id, vcf) })

    // 2. Filter VCF using Hail (single channel: tuple id, vcf, job_config)
    filtered_vcf = filter_hail(
        mt.mt.join(vcf_in).map { id, mtx, vcf, job_config -> tuple(id, vcf, job_config) }
    )

    // 3. Hail QC
    qc_stats = hail_qc(filtered_vcf.vcf)

    // 4. Hail Kinship
    kinship = hail_kinship(filtered_vcf.vcf)

    // 5. Hail PCA
    pca_results = hail_pca(filtered_vcf.vcf)

    emit:
    vcf = filtered_vcf.vcf
    qc_stats = qc_stats.qc_stats
    kinship = kinship.kinship
    pca_scores = pca_results.scores
    pca_loadings = pca_results.loadings
    pca_eigenvalues = pca_results.eigenvalues
}
