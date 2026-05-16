nextflow.enable.dsl = 2

params.input = null
params.output_dir = './results/snp_qc'
params.maf = 0.05
params.max_missing = 0.1
params.min_qual = 30

workflow {
    if (!params.input) error "Missing required param: --input"
    RUN_SNP_QC(file(params.input))
}

process RUN_SNP_QC {
    publishDir params.output_dir, mode: 'copy'

    input:
    path snp_table

    output:
    path '*.tsv'
    path '*.png'

    script:
    """
    python -m genetics.wheat.snp_qc \
      --input ${snp_table} \
      --output-prefix snp_qc \
      --maf ${params.maf} \
      --max-missing ${params.max_missing} \
      --min-qual ${params.min_qual}
    """
}
