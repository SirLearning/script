nextflow.enable.dsl = 2

params.kmer_matrix = null
params.phenotype = null
params.output_dir = './results/kgwas'
params.trait = 'Trait'

workflow {
    if (!params.kmer_matrix) error "Missing required param: --kmer_matrix"
    if (!params.phenotype) error "Missing required param: --phenotype"
    RUN_KGWAS(file(params.kmer_matrix), file(params.phenotype))
}

process RUN_KGWAS {
    publishDir params.output_dir, mode: 'copy'

    input:
    path kmer_matrix
    path phenotype

    output:
    path '*.tsv'
    path '*.png'

    script:
    """
    python -m genetics.wheat.kgwas \
      --kmer-matrix ${kmer_matrix} \
      --phenotype ${phenotype} \
      --output-prefix kgwas \
      --trait ${params.trait}
    """
}
