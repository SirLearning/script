#! /usr/bin/env nextflow

nextflow.enable.dsl=2

params.bams_dir = "data/bams"
params.reference_genome = "data/reference/genome.fa"
params.usb_mnt_dir = "/mnt"

def helpMessage() {
    log.info """
    Usage: nextflow run align.nf [options]

    Options:

    Example:
    screen -dmS cp_115 bash -c "\
        cd /data/home/tuser1/run && \
        source ~/.bashrc && conda activate run && \
        nextflow run /data/home/tuser1/git/script/workflow/Genetics/genotype/align.nf \
            --user_dir /data/home/tuser1 \
            --output_dir /data/home/tuser1 \
            --bams_dir /data/home/tuser1/01bam \
            --usb_mnt_dir /mnt "
    """
}

workflow {
    files_ch = channel.fromPath(params.bams_dir + "/CRR*.rmdup.bam", type: 'file').map { bam -> 
        def id = bam.baseName.replace(".rmdup.bam", "")
        def bai = file("${bam}.bai")
        def md5 = file("${bam}.md5")
        [id, bam, bai, md5]
    }

    // Transpose list of lists to tuple of lists to match process input
    files = files_ch.collect(flat: false).map { list ->
        if (!list) return tuple([], [], [], [])
        def t = list.transpose()
        tuple(t[0], t[1], t[2], t[3])
    }

    usb_dirs_list = ["usb", "usb-2", "usb-3"]
    dirs = channel.from(usb_dirs_list).map { dir -> 
        def dir_path = file("${params.usb_mnt_dir}/${dir}")
        [dir, dir_path]
    }.collect(flat: false).map { list ->
        if (!list) return tuple([], [])
        def t = list.transpose()
        tuple(t[0], t[1])
    }

    _cp_out = cp_bams_based_on_usb_size(files, dirs)
}

process ck_md5sum_wtk_fq {
    tag "ck_fq_md5_${id}"
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
    exec > ck_fq_md5_${id}.log 2>&1

    echo "Checking MD5 sums for ${dir}..."
    md5sum ${f1_file} ${r2_file} > ${id}.md5
    echo "MD5 sums for ${dir} saved to ${id}.md5"
    """
}

process ck_md5sum_wtk_bam {
    tag "ck_bam_md5_${id}"
    publishDir "results/bam_md5sums", mode: 'copy'

    input:
    tuple val(id), path(bam_file)

    output:
    tuple val(id), path("${bam_file}.md5"), path(bam_file)

    script:
    """
    set -euo pipefail
    exec > ck_bam_md5_${id}.log 2>&1

    echo "Checking MD5 sum for BAM file ${bam_file}..."
    md5sum ${bam_file} > ${bam_file}.md5
    echo "MD5 sum for BAM file ${bam_file} saved to ${bam_file}.md5"
    """
}

process prepare_reference_genome {
    input:
    tuple val(ref_id), path(reference_genome)

    output:
    tuple val(ref_id), path(reference_genome), path("${ref_id}_bwa_mem2_idx")

    script:
    """
    set -euo pipefail
    exec > prepare_reference_${ref_id}.log 2>&1

    echo "Preparing reference genome ${ref_id} for BWA-MEM2..."
    bwa-mem2 index -p ${ref_id}_bwa_mem2_idx ${reference_genome}
    echo "Reference genome ${ref_id} indexed for BWA-MEM2 and saved as ${ref_id}_bwa_mem2_idx.*"
    """
}

process align_with_bwa_mem2 {
    tag "align_${fastq_id}"

    input:
    tuple val(fastq_id), val(sample_id), path(f1_file), path(r2_file)
    tuple val(ref_id), path(reference_genome), path(bwa_mem2_idx)

    output:
    tuple val(fastq_id), path("${fastq_id}.bam")

    script:
    """
    set -euo pipefail
    exec > align_${fastq_id}.log 2>&1

    echo "Aligning ${fastq_id} with BWA-MEM2..."
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R "@RG\\tID:${fastq_id}\\tSM:${sample_id}\\tPL:illumina" \\
        ${reference_genome} \\
        ${f1_file} ${r2_file} | samtools view -bS - > ${fastq_id}.bam
    echo "Alignment for ${fastq_id} completed and saved to ${fastq_id}.bam"
    """
}

process filter_bam_with_samtools {
    tag "filter_${bam_id}"

    input:
    tuple val(bam_id), path(input_bam)

    output:
    tuple val(bam_id), path("${bam_id}.rmdup.bam")

    script:
    """
    set -euo pipefail
    exec > filter_${bam_id}.log 2>&1

    echo "Filtering BAM file ${input_bam}..."
    samtools sort \\
        -n -m 4G \\
        -@ ${task.cpus} ${input_bam} | \\
    samtools fixmate \\
        -m -u \\
        -@ ${task.cpus} - - | \\
    samtools sort \\
        -m 4G \\
        -@ ${task.cpus} - - | \\
    samtools markdup \\
        -r -@ ${task.cpus} - ${bam_id}.rmdup.bam
    echo "Filtered BAM file saved as ${bam_id}.rmdup.bam"
    """
}

process filter_bam_with_samtools_no_rmdup {
    tag "filter_no_rmdup_${bam_id}"

    input:
    tuple val(bam_id), path(input_bam)

    output:
    tuple val(bam_id), path("${bam_id}.markdup.bam")

    script:
    """
    set -euo pipefail
    exec > filter_no_rmdup_${bam_id}.log 2>&1

    echo "Filtering BAM file ${input_bam} without removing duplicates..."
    samtools sort \\
        -n -m 4G \\
        -@ ${task.cpus} ${input_bam} | \\
    samtools fixmate \\
        -m -u \\
        -@ ${task.cpus} - - | \\
    samtools sort \\
        -m 4G \\
        -@ ${task.cpus} - - | \\
    samtools markdup \\
        -@ ${task.cpus} - ${bam_id}.markdup.bam
    echo "Filtered BAM file (no rmdup) saved as ${bam_id}.markdup.bam"
    """
}

process prepare_bam_index_with_samtools {
    tag "index_${bam_id}"

    input:
    tuple val(bam_id), path(input_bam)

    output:
    tuple val(bam_id), path("${input_bam}.bai"), path(input_bam)

    script:
    """
    set -euo pipefail
    exec > index_${bam_id}.log 2>&1

    echo "Indexing BAM file ${input_bam} with samtools..."
    samtools index -@ ${task.cpus} ${input_bam}
    echo "BAM file ${input_bam} indexed and saved as ${input_bam}.bai"
    """
}

process ct_depth_with_mosdepth {
    tag "ct_depth_${prefix}"

    input:
    tuple val(prefix), path(input_bam)

    output:
    tuple val(prefix), path("${prefix}.mosdepth")

    script:
    """
    set -euo pipefail
    exec > ct_depth_${prefix}.log 2>&1

    echo "Calculating coverage depth for ${input_bam} with mosdepth..."
    mosdepth \\
        -t ${task.cpus} \\
        -n ${prefix} \\
        ${input_bam}
    echo "Coverage depth calculated and saved as ${prefix}.mosdepth.*"
    """
}

process cp_bams_based_on_usb_size {
    tag "cp_bams"
    publishDir "${params.output_dir}/logs", mode: 'copy'
    conda "${params.user_dir}/miniconda3/envs/stats"

    input:
    tuple val(bam_ids), path(bams), path(bai_files), path(md5_files)
    tuple val(dir_names), val(dir_paths)

    output:
    path "usb_transfer_log.txt"

    script:
    """
    #!/usr/bin/env python3
    import sys
    from infra.server.cp import run_copy_process

    bams_str = "${bams}"
    bams = bams_str.split()

    bais_str = "${bai_files}"
    bais = bais_str.split()

    md5s_str = "${md5_files}"
    md5s = md5s_str.split()
    
    files = bams + bais + md5s

    usb_dirs_str = "${dir_paths}"
    usb_dirs = usb_dirs_str.split()

    log_file = "usb_transfer_log.txt"
    threads = int("${task.cpus}")

    print(f"Running copy process for bams with {threads} threads...")
    run_copy_process(files, usb_dirs, log_file, threads)
    """
}
