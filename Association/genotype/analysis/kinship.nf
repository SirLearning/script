nextflow.enable.dsl=2

process KINSHIP_ANALYSIS {
    tag "Kinship: ${meta.id}"
    publishDir "${params.outdir}/genotype/kinship", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.kinship.txt"), emit: kinship_matrix
    path "*.log"

    script:
    def prefix = meta.id
    """
    echo "Calculating kinship matrix for ${vcf}" > ${prefix}.log

    # Use plink2 to calculate kinship
    plink2 \\
        --vcf ${vcf} \\
        --make-king-table \\
        --out ${prefix}

    # The output is ${prefix}.kin0, we can rename it
    mv ${prefix}.kin0 ${prefix}.kinship.txt

    echo "Kinship analysis complete." >> ${prefix}.log
    """
}