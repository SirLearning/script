nextflow.enable.dsl = 2

params.input = null
params.output_dir = './results/cnv'
params.del_z = -2
params.dup_z = 2

workflow {
    if (!params.input) error "Missing required param: --input"
    RUN_CNV(file(params.input))
}

process RUN_CNV {
    publishDir params.output_dir, mode: 'copy'

    input:
    path depth_matrix

    output:
    path '*.tsv'
    path '*.png'

    script:
    """
    python -m genetics.wheat.cnv \
      --input ${depth_matrix} \
      --output-prefix cnv \
      --del-z ${params.del_z} \
      --dup-z ${params.dup_z}
    """
}
