nextflow.enable.dsl=2

process POPULATION_STRUCTURE {
    tag "PCA: ${meta.id}"
    publishDir "${params.outdir}/genotype/ps", mode: 'copy'

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
