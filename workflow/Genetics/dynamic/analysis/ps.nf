nextflow.enable.dsl=2

workflow population_structure {
    take:
    // Expect a combined channel: [ val(meta), path(vcf), val(job_config) ]
    vcf_in

    main:
    // Run PCA for Population Structure Analysis
    pca_results = run_pca(vcf_in)

    emit:
    pca = pca_results.pca_results
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
