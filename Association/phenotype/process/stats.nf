nextflow.enable.dsl=2

process PHENOTYPE_STATS {
    tag "Stats Pheno: ${meta.id}"
    publishDir "${params.outdir}/phenotype/stats", mode: 'copy'

    input:
    tuple val(meta), path(pheno_file)

    output:
    tuple val(meta), path("*.blup.txt"), emit: stats_pheno
    path "*.log"

    script:
    def prefix = meta.id
    """
    echo "Running phenotype statistical analysis for ${pheno_file}" > ${prefix}.log
    
    Rscript ${projectDir}/src/r/phenotype_stats.R ${pheno_file} ${prefix} >> ${prefix}.log 2>&1
    """
}