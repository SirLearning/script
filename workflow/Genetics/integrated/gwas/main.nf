nextflow.enable.dsl = 2

params.genotype = null
params.phenotype = null
params.output_dir = './results/gwas'
params.trait = 'Trait'

workflow {
    if (!params.genotype) error "Missing required param: --genotype"
    if (!params.phenotype) error "Missing required param: --phenotype"
    RUN_WHEAT_GWAS(file(params.genotype), file(params.phenotype))
}

process RUN_WHEAT_GWAS {
    publishDir params.output_dir, mode: 'copy'

    input:
    path genotype
    path phenotype

    output:
    path '*.tsv'
    path '*.png'

    script:
    """
    python -m genetics.wheat.wheat_gwas \
      --genotype ${genotype} \
      --phenotype ${phenotype} \
      --output-prefix gwas \
      --trait ${params.trait}
    """
}
