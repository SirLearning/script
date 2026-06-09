nextflow.enable.dsl=2

process subsampling_pfile_for_test {
    tag "subsample pfile for test: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.{pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)
    val thin_rate

    output:
    tuple val(id), val(chr), val("${id}.thin"), path("${id}.thin.pgen"), path("${id}.thin.psam"), path("${id}.thin.pvar"), emit: pfile
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > subsample_pfile_${id}.log 2>&1

    plink2 --pfile ${prefix} \\
        --thin ${thin_rate} \\
        --seed \$(date +%s) \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${id}.thin
    """
}

process subsampling_common_variant_pfile_for_test {
    tag "common-thin pfile for test: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.{pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), val("${id}.common.thin"), path("${id}.common.thin.pgen"), path("${id}.common.thin.psam"), path("${id}.common.thin.pvar"), emit: pfile
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > common_thin_pfile_${id}.log 2>&1

    plink2 --pfile ${prefix} \\
        --geno ${params.hf_geno} \\
        --maf ${params.hf_maf} \\
        --thin ${params.hf_thin_rate} \\
        --seed ${params.ld_decay_random_seed} \\
        --allow-extra-chr \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${id}.common.thin
    """
}

process merge_subgenome_test_pfile {
    tag "merge plink2 test subgenome pfile: ${subgenome}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.{bed,bim,fam,pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    // Receives lists of files/values grouped by subgenome
    tuple val(subgenome), val(ids), val(chrs), val(prefixes), path(pgens), path(psams), path(pvars)

    output:
    tuple val(subgenome), val("sub_${subgenome}"), val("${subgenome}_test.plink"), path("${subgenome}_test.plink.bed"), path("${subgenome}_test.plink.bim"), path("${subgenome}_test.plink.fam"), emit: bfile
    tuple val(subgenome), val("sub_${subgenome}"), val("${subgenome}_test.plink2"), path("${subgenome}_test.plink2.pgen"), path("${subgenome}_test.plink2.psam"), path("${subgenome}_test.plink2.pvar"), emit: pfile
    tuple val(subgenome), path("*.log"), emit: logs

    script:
    def has_multiple_prefixes = prefixes.size() > 1
    def merge_list_block = has_multiple_prefixes
        ? """echo "${prefixes[1..-1].join('\n')}" > ${subgenome}.merge_list.txt

    plink2 --pfile ${prefixes[0]} \\
        --pmerge-list ${subgenome}.merge_list.txt \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${subgenome}_test.plink2"""
        : """plink2 --pfile ${prefixes[0]} \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${subgenome}_test.plink2"""
    """
    set -euo pipefail
    exec > merge_subgenome_test.${subgenome}.log 2>&1
    
    ${merge_list_block}

    plink2 --pfile ${subgenome}_test.plink2 \\
        --make-bed \\
        --threads ${task.cpus} \\
        --out ${subgenome}_test.plink
    """
}

