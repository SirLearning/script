nextflow.enable.dsl=2

/*
 * WAP association compute on the genotype publish tree.
 *
 * Cohort contract (doc/project_knowledge/domain/sample_cohorts.yaml): only WAP
 * accessions carry phenotypes, so association runs restrict genotypes to
 * assoc.keep. Variant restriction comes from partial_task variant_sets
 * (*.extract lists), not from ad-hoc MAF filters here.
 *
 * Compute stays in native tools (plink2 / gcta). Plotting lives in gwas_plot.nf.
 */

/*
 * Structure covariates for the WAP panel.
 *
 * Population structure is a property of the individual, not of a subgenome, so
 * PCs are computed once from one pfile (params.assoc_wap_pca_source) and reused
 * for every subgenome's GLM. Per-subgenome PCA also fails outright on the
 * unanchored `Others` pfile, where almost nothing survives the MAF and LD filters.
 */
process plink2_wap_pca_covar {
    tag "wap pca: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/assoc/pca", mode: 'copy', pattern: "*.{eigenvec,eigenval}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)
    path keep

    output:
    path "${id}.wap.pca.eigenvec", emit: covar
    path "*.eigenval", emit: eigenval
    path "*.log", emit: logs

    script:
    """
    set -euo pipefail
    exec > ${id}.wap_pca.log 2>&1

    plink2 --pfile ${prefix} \\
        --keep ${keep} \\
        --allow-extra-chr \\
        --maf ${params.assoc_wap_pca_maf} \\
        --indep-pairwise ${params.pca_ld_window} ${params.pca_ld_step} ${params.pca_ld_r2} \\
        --threads ${task.cpus} \\
        --out ${id}.wap

    plink2 --pfile ${prefix} \\
        --keep ${keep} \\
        --allow-extra-chr \\
        --extract ${id}.wap.prune.in \\
        --pca approx ${params.pc_num} \\
        --threads ${task.cpus} \\
        --out ${id}.wap.pca
    """
}

process plink2_wap_glm {
    tag "wap glm ${trait}: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/assoc", mode: 'copy', pattern: "*.glm.*"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_16'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)
    path keep
    path phenotype
    val trait
    path covar
    path extract

    output:
    tuple val(id), val(trait), path("${id}.${trait}.glm.*"), emit: glm
    path "*.log", emit: logs

    script:
    def covar_arg = covar.name == 'NO_COVAR' ? '' : "--covar ${covar} --covar-variance-standardize"
    def extract_arg = extract.name == 'NO_EXTRACT' ? '' : "--extract ${extract}"
    """
    set -euo pipefail
    exec > ${id}.${trait}.wap_glm.log 2>&1

    plink2 --pfile ${prefix} \\
        --keep ${keep} \\
        --allow-extra-chr \\
        --pheno ${phenotype} \\
        --pheno-name ${trait} \\
        ${covar_arg} ${extract_arg} \\
        --glm omit-ref no-x-sex hide-covar \\
        --threads ${task.cpus} \\
        --out ${id}
    """
}

/*
 * Rare-variant gene burden for one subgenome and one variant_sets set.
 *
 * plink2 has no gene-collapsing test, so the regression runs in Python. plink2
 * still does the genotype work: it applies the cohort keep and the set's extract
 * list and writes an additive matrix, which the burden module then folds into one
 * score per gene. Gene membership comes from {id}.variant_sets.tsv (VEP Gene) —
 * the *.extract list carries variant IDs only.
 */
process wap_burden_test {
    tag "wap burden ${trait}: ${id} ${set_name}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/assoc", mode: 'copy', pattern: "*.burden.tsv"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_16'

    input:
    // pfile and set travel as one tuple: two parallel channels would be zipped by
    // arrival order, which risks pairing a subgenome with another's extract list.
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar),
        val(set_name), path(extract), path(variant_sets)
    path keep
    path phenotype
    val trait
    path covar

    output:
    tuple val(id), val(trait), val(set_name), path("*.burden.tsv"), emit: burden
    path "*.log", emit: logs

    script:
    def covar_arg = covar.name == 'NO_COVAR' ? 'None' : "\"${covar}\""
    """
    set -euo pipefail
    exec > ${id}.${trait}.${set_name}.burden.log 2>&1

    plink2 --pfile ${prefix} \\
        --keep ${keep} \\
        --extract ${extract} \\
        --allow-extra-chr \\
        --export Av \\
        --threads ${task.cpus} \\
        --out ${id}.${set_name}

    python - <<'PY'
import sys

if "${params.src_dir}" and "${params.src_dir}" != "null":
    sys.path.insert(0, "${params.src_dir}/python")

from genetics.gwas.burden import run_burden_test

run_burden_test(
    "${id}.${set_name}.traw",
    "${variant_sets}",
    "${set_name}",
    "${phenotype}",
    "${trait}",
    "${id}.${trait}.${set_name}",
    covar_table=${covar_arg},
    collapse="${params.assoc_wap_burden_collapse}",
)
PY
    """
}
