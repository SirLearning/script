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
