nextflow.enable.dsl = 2

params.input = null
params.output_dir = './results/genetic_map'

workflow {
    if (!params.input) error "Missing required param: --input"
    RUN_GENETIC_MAP(file(params.input))
}

process RUN_GENETIC_MAP {
    publishDir params.output_dir, mode: 'copy'

    input:
    path marker_table

    output:
    path '*.tsv'
    path '*.png'

    script:
    """
    python -m genetics.wheat.genetic_map \
      --input ${marker_table} \
      --output-prefix genetic_map
    """
}
