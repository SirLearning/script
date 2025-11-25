version 1.0

workflow wgs {
    input {
        String sample_id
        String read_group

        File fastq1
        File fastq2
        Array[File] genomes

        String docker = "registry.cn-hangzhou.aliyuncs.com/damo-gene/dna-variant-calling-cpu@sha256:64c131e1c08fed5c7399ddef7b9aeb87eea0eaacaa8f5dde0ab93019b5a7230e"
    }

    call FastqToBam {
        input:
            fastq1=fastq1,
            fastq2=fastq2,
            genomes=genomes,
            sample_id=sample_id,
            read_group=read_group,
            docker=docker
    }

    output {
        File sorted_bam = FastqToBam.bam_out
    }
}

task FastqToBam {
    input {
        String sample_id
        String read_group

        File fastq1
        File fastq2
        Array[File] genomes

        Int cpu = 32
        String memory = "64G"
        String disk = "local-disk 100 cloud_essd_pl1"
        String docker
    }

    String sorted_bam = sample_id + ".sort.bam"

    command <<<
        set -euxo pipefail
        threads=$(nproc)
        
        damo-bwa-mem2 mem -R ~{read_group} \
        ~{genomes[0]} \
        ~{fastq1} ~{fastq2} \
        -t $threads | damo_samtools sort -o ~{sorted_bam} -T tmp -l3 -@8 --write-index
    >>>

    runtime {
        cpu: cpu
        memory: memory
        disks: disk
        docker: docker
    }

    output {
        File bam_out = sorted_bam
    }
}
