nextflow.enable.dsl=2

/*
 * Cohort-restricted allele frequency and missingness from merged test pfiles.
 * One PLINK2 pass per (subgenome, cohort keep file) pair; the keep files come
 * from partial_task sample_manifest (assoc.keep / wgs_core.keep).
 *
 * Downstream: dual_af_from_afreq.nf and partial_task variant_sets.
 */

process plink2_cohort_freq {
    tag "cohort freq ${cohort}: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/variant", mode: 'copy', pattern: "*.{afreq,acount,vmiss}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_16'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar), val(cohort), path(keep)

    output:
    tuple val(cohort), val(id), path("${id}.${cohort}.afreq"), emit: afreq
    tuple val(cohort), val(id), path("${id}.${cohort}.acount"), emit: acount
    tuple val(cohort), val(id), path("${id}.${cohort}.vmiss"), emit: vmiss
    path "*.log", emit: logs

    script:
    """
    set -euo pipefail
    exec > ${id}.${cohort}.cohort_freq.log 2>&1

    plink2 --pfile ${prefix} \\
        --keep ${keep} \\
        --allow-extra-chr \\
        --freq \\
        --missing variant-only \\
        --threads ${task.cpus} \\
        --out ${id}.${cohort}

    plink2 --pfile ${prefix} \\
        --keep ${keep} \\
        --allow-extra-chr \\
        --freq counts \\
        --threads ${task.cpus} \\
        --out ${id}.${cohort}.counts

    mv ${id}.${cohort}.counts.acount ${id}.${cohort}.acount
    """
}
