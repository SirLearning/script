nextflow.enable.dsl = 2

params.input = null
params.output_dir = './results/hapmap'
params.window_size = 100000

workflow {
    if (!params.input) error "Missing required param: --input"
    RUN_HAPMAP(file(params.input))
}

process RUN_HAPMAP {
    publishDir params.output_dir, mode: 'copy'

    input:
    path geno_table

    output:
    path '*.tsv'

    script:
    """
    python -m genetics.wheat.hapmap \
      --input ${geno_table} \
      --output-prefix hapmap \
      --window-size ${params.window_size}
    """
}
