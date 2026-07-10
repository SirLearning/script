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

process subsampling_maf_only_pfile_for_test {
    tag "maf-only common pfile for test: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.{pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), val("${id}.common.only"), path("${id}.common.only.pgen"), path("${id}.common.only.psam"), path("${id}.common.only.pvar"), emit: pfile
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > common_only_pfile_${id}.log 2>&1

    plink2 --pfile ${prefix} \\
        --maf ${params.hf_maf} \\
        --thin ${params.hf_thin_rate} \\
        --seed ${params.ld_decay_random_seed} \\
        --allow-extra-chr \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${id}.common.only
    """
}

process build_chr_rare_inter_extract {
    tag "rare-inter extract: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/extract", mode: 'copy', pattern: "*.txt"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/extract/info", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/extract/afreq/${id}", mode: 'copy', pattern: "${id}.*.afreq"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'rare_inter_extract'

    input:
    path group_file
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar), path("${id}.rare_inter.txt"), emit: pfile_with_extract

    script:
    def published_afreq_dir = "${params.output_dir}/${params.job}/process/${params.mod}/extract/afreq"
    def published_extract = "${params.output_dir}/${params.job}/process/${params.mod}/extract/${id}.rare_inter.txt"
    def published_summary = "${params.output_dir}/${params.job}/process/${params.mod}/extract/info/${id}.rare_inter.summary.tsv"
    """
    #!/usr/bin/env python
    import sys

    sys.stdout = open("build_rare_inter_${id}.log", "w")
    sys.stderr = sys.stdout

    from genetics.genomics.variant.group_maf_intersection import write_rare_inter_extract_for_chromosome

    print("Building rare-inter extract for ${id} (all sample_group MAF < ${params.hf_maf}) ...")
    summary = write_rare_inter_extract_for_chromosome(
        "${prefix}",
        "${group_file}",
        "${id}.rare_inter.txt",
        chr_id="${id}",
        maf_threshold=${params.hf_maf},
        summary_path="${id}.rare_inter.summary.tsv",
        published_afreq_dir="${published_afreq_dir}",
        published_extract_path="${published_extract}",
        published_summary_path="${published_summary}",
        afreq_emit_dir=".",
        threads=${task.cpus},
        group_workers=${params.rare_inter_group_workers},
        plink_threads=${params.rare_inter_plink_threads},
    )
    print(summary)
    """
}

process subsampling_common_inter_pfile_for_test {
    tag "common-inter thin pfile for test: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.{pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar), path(extract)

    output:
    tuple val(id), val(chr), val("${id}.common.inter"), path("${id}.common.inter.pgen"), path("${id}.common.inter.psam"), path("${id}.common.inter.pvar"), emit: pfile
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > common_inter_pfile_${id}.log 2>&1

    if [ ! -s "${extract}" ]; then
        echo "SKIP: empty rare-inter extract for ${id}; no all-group-rare sites to thin."
        exit 0
    fi

    plink2 --pfile ${prefix} \\
        --extract ${extract} \\
        --thin ${params.hf_common_inter_thin_rate} \\
        --seed ${params.ld_decay_random_seed} \\
        --allow-extra-chr \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${id}.common.inter
    """
}

process subsampling_rare_only_pfile_for_test {
    tag "rare-only pfile for test: ${id}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}", mode: 'copy', pattern: "*.{pgen,psam,pvar}"
    publishDir "${params.output_dir}/${params.job}/process/${params.mod}/logs", mode: 'copy', pattern: "*.log"
    conda "${params.user_dir}/miniconda3/envs/stats"
    label 'cpus_32'

    input:
    tuple val(id), val(chr), val(prefix), path(pgen), path(psam), path(pvar)

    output:
    tuple val(id), val(chr), val("${id}.rare.only"), path("${id}.rare.only.pgen"), path("${id}.rare.only.psam"), path("${id}.rare.only.pvar"), emit: pfile
    tuple val(id), val(chr), path("*.log"), emit: logs

    script:
    """
    set -euo pipefail
    exec > rare_only_pfile_${id}.log 2>&1

    plink2 --pfile ${prefix} \\
        --max-maf ${params.hf_maf} \\
        --thin ${params.thin_rate} \\
        --seed ${params.ld_decay_random_seed} \\
        --allow-extra-chr \\
        --make-pgen \\
        --threads ${task.cpus} \\
        --out ${id}.rare.only
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

