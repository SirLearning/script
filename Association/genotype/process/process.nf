nextflow.enable.dsl=2

process VCF_TO_PLINK {
    tag "plink from vcf"
    publishDir params.outdir + '/run', mode: 'copy'
    cpus params.plink_threads as int

    input:
    path vcf

    output:
    tuple path('geno.bed'), path('geno.bim'), path('geno.fam'), emit: plink

    script:
    """
    set -euo pipefail
    plink --vcf ${vcf} --double-id --allow-extra-chr --make-bed --out geno --threads ${task.cpus}
    """
}