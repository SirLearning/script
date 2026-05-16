nextflow.enable.dsl = 2

params.input = null
params.output_dir = './results/pca_tsne'
params.n_pcs = 10
params.tsne_perplexity = 30

workflow {
    if (!params.input) error "Missing required param: --input"
    RUN_PCA_TSNE(file(params.input))
}

process RUN_PCA_TSNE {
    publishDir params.output_dir, mode: 'copy'

    input:
    path geno_matrix

    output:
    path '*.tsv'
    path '*.png'

    script:
    """
    python -m genetics.wheat.population_structure \
      --input ${geno_matrix} \
      --output-prefix population_structure \
      --n-pcs ${params.n_pcs} \
      --tsne-perplexity ${params.tsne_perplexity}
    """
}
