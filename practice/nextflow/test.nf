nextflow.enable.dsl = 2

// 默认参数设置
params.reads = "$projectDir/data/*_{1,2}.fastq.gz"
params.reference = "$projectDir/reference/genome.fasta"
params.outdir = "$projectDir/results"
params.threads = 8
params.help = false

// 定义打印帮助信息的函数
def printHelp() {
    log.info """
    ==================================================
    DNA序列比对和变异检测流程
    ==================================================
    
    用法: nextflow run test.nf [选项]
    
    选项:
      --reads         输入测序数据路径 (默认: $params.reads)
      --reference     参考基因组路径 (默认: $params.reference)
      --outdir        输出目录 (默认: $params.outdir)
      --threads       CPU线程数 (默认: $params.threads)
      --help          显示帮助信息
    """
}

// 定义流程

// 1. 质量控制
process fastqc {
    tag "FastQC: ${sample_id}"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "${sample_id}*_fastqc.{zip,html}"
    
    script:
    """
    fastqc -t ${params.threads} ${reads}
    """
}

// 2. 比对到参考基因组
process alignReads {
    tag "BWA: ${sample_id}"
    publishDir "${params.outdir}/aligned", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    path reference
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")
    
    script:
    """
    # 如果索引不存在则创建索引
    if [ ! -f ${reference}.bwt ]; then
        bwa index ${reference}
    fi
    
    # 进行序列比对
    bwa mem -t ${params.threads} ${reference} ${reads[0]} ${reads[1]} | \
    samtools sort -@ ${params.threads} -o ${sample_id}.sorted.bam -
    
    # 为BAM文件创建索引
    samtools index ${sample_id}.sorted.bam
    """
}

// 3. 标记重复序列
process markDuplicates {
    tag "MarkDup: ${sample_id}"
    publishDir "${params.outdir}/dedup", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bai")
    
    script:
    """
    gatk MarkDuplicates \
        -I ${bam} \
        -O ${sample_id}.dedup.bam \
        -M ${sample_id}.metrics.txt \
        --CREATE_INDEX true
    """
}

// 4. 变异检测
process callVariants {
    tag "GATK: ${sample_id}"
    publishDir "${params.outdir}/vcf", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    
    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi")
    
    script:
    """
    # 创建参考基因组字典
    if [ ! -f ${reference.baseName}.dict ]; then
        gatk CreateSequenceDictionary -R ${reference}
    fi
    
    # 调用变异
    gatk HaplotypeCaller \
        -R ${reference} \
        -I ${bam} \
        -O ${sample_id}.vcf.gz
    """
}

// 5. 变异过滤
process filterVariants {
    tag "Filter: ${sample_id}"
    publishDir "${params.outdir}/filtered_vcf", mode: 'copy'
    
    input:
    tuple val(sample_id), path(vcf), path(tbi)
    path reference
    
    output:
    path "${sample_id}.filtered.vcf.gz"
    
    script:
    """
    gatk VariantFiltration \
        -R ${reference} \
        -V ${vcf} \
        -O ${sample_id}.filtered.vcf.gz \
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0" \
        --filter-name "basic_filter"
    """
}

// 定义工作流
workflow {
    // 检查参考基因组文件是否存在
    reference_file = file(params.reference)
    if (!reference_file.exists()) {
        error "参考基因组文件不存在: ${params.reference}"
    }
    
    // 创建输入通道
    read_pairs_ch = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .ifEmpty { error "未找到匹配的测序数据: ${params.reads}" }

    // 显示帮助信息
    if (params.help) {
        printHelp()
        exit 0
    }
    
    // 打印参数信息
    log.info """
    ==================================================
    DNA序列比对和变异检测流程
    ==================================================
    测序数据: ${params.reads}
    参考基因组: ${params.reference}
    输出目录: ${params.outdir}
    线程数: ${params.threads}
    """
    
    // 检查参考基因组文件是否存在
    reference_file = file(params.reference)
    // 进行序列比对
    alignReads(read_pairs_ch, reference_file)
    
    // 标记重复序列
    markDuplicates(alignReads.out)
    
    // 变异检测
    callVariants(markDuplicates.out, reference_file)
    
    // 变异过滤
    filterVariants(callVariants.out, reference_file)

    // 工作流完成通知
    workflow.onComplete = {
        log.info "流程完成时间: $workflow.complete"
        log.info "执行状态: ${ workflow.success ? '成功' : '失败' }"
        log.info "执行时长: $workflow.duration"
    }
}

