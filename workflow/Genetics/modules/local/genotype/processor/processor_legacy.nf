nextflow.enable.dsl=2

process mk_vcftools_basic_info {
    tag "vcftools stats: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/sample", mode: 'copy', pattern: "*.{imiss, idepth}"
    publishDir "${params.output_dir}/${params.job}/process/variant", mode: 'copy', pattern: "*.{frq,hwe,lmiss,ldepth.mean,lqual}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), path(vcf)

    output:
    tuple val(id), val(chr), path("${id}.frq"), emit: frq
    tuple val(id), val(chr), path("${id}.hwe"), emit: hwe
    tuple val(id), val(chr), path("${id}.lmiss"), emit: lmiss
    tuple val(id), val(chr), path("${id}.imiss"), emit: imiss
    tuple val(id), val(chr), path("${id}.ldepth.mean"), emit: ldepth_mean
    tuple val(id), val(chr), path("${id}.idepth"), emit: idepth
    tuple val(id), val(chr), path("${id}.lqual"), emit: lqual
    path "${chr}.vcftools.log", emit: log

    script:
    """
    set -euo pipefail
    exec > ${chr}.vcftools.log 2>&1

    vcftools --gzvcf ${vcf} \\
        --freq \\
        --hardy \\
        --missing-site \\
        --missing-indv \\
        --site-mean-depth \\
        --depth \\
        --site-quality \\
        --out ${id}.vcftools
    """
}
