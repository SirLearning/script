nextflow.enable.dsl = 2

process PHENO_PROCESS {
    tag "phenotype process ${params.trait}"
    publishDir params.outdir + "/phenotype", mode: 'copy'

    input:
    path pheno

    output:
    path 'pheno_processed.tsv', emit: pheno
    path 'blue.tsv', optional: true
    path 'blup.tsv', optional: true
    path 'pheno_summary.tsv'
    path 'qc_hist.png', optional: true
    path 'qc_box.png', optional: true

    script:
    def idcol = params.pheno_id_col
    def trait = params.trait
    def envcol = params.envcol ?: ''
    def repcol = params.repcol ?: ''
    def fixed  = params.fixed_effects ?: ''
    def random = params.random_effects ?: ''
    def effect = (params.pheno_effect ?: 'raw').toString().toUpperCase()
    def use_blue = (params.use_blue ?: (effect=='BLUE'))
    def use_blup = (params.use_blup ?: (effect=='BLUP'))
    def scriptPath = "${projectDir}/src/r/gwas/phenotype_process.R"
    """
    set -euo pipefail
    Rscript ${scriptPath} \
      --pheno ${pheno} \
      --idcol ${idcol} \
      --trait ${trait} \
      --envcol "${envcol}" \
      --repcol "${repcol}" \
      --fixed "${fixed}" \
      --random "${random}" \
      $( [ "${use_blue}" = "true" ] && echo "--blue" ) \
      $( [ "${use_blup}" = "true" ] && echo "--blup" ) \
      --outdir .

    # choose processed phenotype for GWAS
    case "${effect}" in
      BLUE)
        if [ -s blue.tsv ]; then cp blue.tsv pheno_processed.tsv; else cp ${pheno} pheno_processed.tsv; fi
        ;;
      BLUP)
        if [ -s blup.tsv ]; then cp blup.tsv pheno_processed.tsv; else cp ${pheno} pheno_processed.tsv; fi
        ;;
      *)
        cp ${pheno} pheno_processed.tsv
        ;;
    esac
    """
}

workflow PHENOTYPE {
    take:
    ch_pheno

    main:
    out = PHENO_PROCESS(ch_pheno)

    emit:
    out.pheno
}
