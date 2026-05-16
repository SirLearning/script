nextflow.enable.dsl = 2

params.input = null
params.output_dir = './results/tagsnp'
params.max_tags = 1000
params.ld_threshold = 0.8

workflow {
    if (!params.input) error "Missing required param: --input"
    RUN_TAGSNP(file(params.input))
}

process RUN_TAGSNP {
    publishDir params.output_dir, mode: 'copy'

    input:
    path genotype

    output:
    path '*.tsv'

    script:
    """
    python -m genetics.wheat.tagsnp \
      --input ${genotype} \
      --output-prefix tagsnp \
      --max-tags ${params.max_tags} \
      --ld-threshold ${params.ld_threshold}
    """
}
