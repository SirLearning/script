nextflow.enable.dsl=2

include { HAIL_PCA } from './hail_pca.nf'

workflow population_structure {
    take:
    vcf_in
    config

    main:
    if (params.tool == 'hail') {
        HAIL_PCA(vcf_in)
        pca_results = HAIL_PCA.out.scores
    } else {
        // Run PCA for Population Structure Analysis
        pca_results = run_pca(vcf_in).pca_results
    }

    emit:
    pca = pca_results
}

process run_pca {
    tag "PCA: ${meta.id}"
    publishDir "${params.outdir}/genotype/population_structure", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.eigenvec"), emit: pca_results
    path "*.log"

    script:
    def prefix = meta.id
    """
    echo "Running PCA for population structure analysis on ${vcf}" > ${prefix}.log

    # Use plink2 for PCA
    plink2 \\
        --vcf ${vcf} \\
        --pca ${params.pca_k} \\
        --out ${prefix}

    echo "PCA complete. Results are in ${prefix}.eigenvec" >> ${prefix}.log
    """
}
