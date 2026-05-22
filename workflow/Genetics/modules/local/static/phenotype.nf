nextflow.enable.dsl=2

process PHENOTYPE_PROCESS {
    tag "Process Pheno: ${meta.id}"
    publishDir "${params.outdir}/phenotype/process", mode: 'copy'

    input:
    tuple val(meta), path(pheno_file)

    output:
    tuple val(meta), path("*.processed.txt"), emit: processed_pheno
    path "*.log"

    script:
    def prefix = meta.id
    """
    echo "Processing phenotype file: ${pheno_file}" > ${prefix}.log

    # This is a placeholder process.
    # In a real scenario, you might connect to a database using a Java/Python script.
    # For now, we just copy the file to demonstrate the workflow.
    cp ${pheno_file} ${prefix}.processed.txt

    echo "Phenotype processing complete." >> ${prefix}.log
    """
}

process multi_traits {
    tag "Multi Traits: ${meta.id}"
    publishDir "${params.outdir}/phenotype/multi_traits", mode: 'copy'

    input:
    tuple val(meta), path(pheno_file)

    output:
    tuple val(meta), path("*.multi_traits.txt"), emit: multi_traits_pheno
    path "*.log"

    script:
    def prefix = meta.id
    """
    echo "Running multi-trait analysis for ${pheno_file}" > ${prefix}.log
    
    Rscript ${projectDir}/src/r/phenotype_multi_traits.R ${pheno_file} ${prefix} >> ${prefix}.log 2>&1
    """
}

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