nextflow.enable.dsl=2

process ct_md5_wtk_fq {
    tag "ct_fq_md5_${id}"
    publishDir "results/fastq_md5sums", mode: 'copy'

    input:
    tuple val(id), path(dir)

    output:
    path "${id}.md5"

    script:
    f1_file = file("${dir}/${id}_f1.fastq.gz").get()
    r2_file = file("${dir}/${id}_r2.fastq.gz").get()
    """
    set -euo pipefail
    exec > ct_fq_md5_${id}.log 2>&1

    echo "Checking MD5 sums for ${dir}..."
    md5sum ${f1_file} ${r2_file} > ${id}.md5
    echo "MD5 sums for ${dir} saved to ${id}.md5"
    """
}

process ct_md5_wtk_bam {
    tag "ct_bam_md5_${id}"
    publishDir "results/bam_md5sums", mode: 'copy'

    input:
    tuple val(id), path(bam_file)

    output:
    tuple val(id), path("${bam_file}.md5"), path(bam_file)

    script:
    """
    set -euo pipefail
    exec > ct_bam_md5_${id}.log 2>&1

    echo "Checking MD5 sum for BAM file ${bam_file}..."
    md5sum ${bam_file} > ${bam_file}.md5
    echo "MD5 sum for BAM file ${bam_file} saved to ${bam_file}.md5"
    """
}

process ck_md5_wtk {
    tag "ck_md5_${bam_id}"

    input:
    tuple val(server), val(bam_id), val(bam_file), val(bai_file), val(md5_file)

    output:
    tuple val(server), path("${bam_id}_md5_check.txt")

    script:
    """
    set -euo pipefail
    exec > ck_md5_${bam_id}.log 2>&1
    
    # 【重点】这里使用 val 而不是 path 获取文件路径，为了防止 Nextflow 将 USB 中几十G的 BAM 文件“多此一举”复制回本地工作目录
    if [[ ! -f "${bam_file}" ]]; then
        echo -e "${bam_file}\tMISSING_BAM\tFAIL" > "${bam_id}_md5_check.txt"
        exit 0
    fi
    if [[ ! -f "${md5_file}" ]]; then
        echo -e "${bam_file}\tMISSING_MD5\tFAIL" > "${bam_id}_md5_check.txt"
        exit 0
    fi

    expected_md5=\$(awk 'NR==1 {print \$1}' "${md5_file}")
    actual_md5=\$(md5sum "${bam_file}" | awk '{print \$1}')

    if [[ "\$expected_md5" == "\$actual_md5" ]]; then
        echo -e "${bam_file}\t${md5_file}\tPASS" > "${bam_id}_md5_check.txt"
    else
        echo -e "${bam_file}\t${md5_file}\tFAIL" > "${bam_id}_md5_check.txt"
    fi
    """
}

process collect_md5_results {
    tag "collect_${server}"
    publishDir "${params.output_dir}/logs", mode: 'copy'

    input:
    tuple val(server), path(check_files)

    output:
    path "${server}_md5_check.txt"

    script:
    """
    cat ${check_files} > ${server}_md5_check.txt
    echo "MD5 check results collected and saved to ${server}_md5_check.txt"
    """
}
