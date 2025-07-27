nextflow.enable.dsl = 2

// 默认参数设置
params.jar = "$projectDir/TIGER_20250526.jar"
params.reference = "$projectDir/data/input/chr1_10M.fa"
params.bam_maps = "$projectDir/data/input/taxaBamMap.txt"
params.samtools = "$projectDir/bin/samtools"
params.outdir = "$projectDir/results"
params.threads = 32
params.memory = "100g"
params.help = false

// FastCall2 参数
params.min_depth = 30
params.min_qual = 20
params.ploidy = 2
params.min_allele_freq = 0.2
params.max_missing = 3
params.min_maf = 0.8
params.min_call_rate = 0.35
params.min_het_freq = 0.2
params.chunk_size = 1

// 性能分析参数
params.perf_events = "cycles,instructions,cache-misses,branches,branch-misses"
params.sample_sizes = "1,2,3,4,5,6"
params.generate_flamegraph = true

// 定义打印帮助信息的函数
def printHelp() {
    log.info """
    ==================================================
    FastCall2 性能分析流程
    ==================================================
    
    用法: nextflow run fastcall.nf [选项]
    
    必要文件:
      --jar           TIGER jar文件路径 (默认: $params.jar)
      --reference     参考基因组路径 (默认: $params.reference)
      --bam_maps      BAM文件映射表 (默认: $params.bam_maps)
      --samtools      samtools路径 (默认: $params.samtools)
    
    输出设置:
      --outdir        输出目录 (默认: $params.outdir)
      --threads       CPU线程数 (默认: $params.threads)
      --memory        内存设置 (默认: $params.memory)
    
    性能分析:
      --perf_events   perf监控事件 (默认: $params.perf_events)
      --sample_sizes  样本大小列表 (默认: $params.sample_sizes)
      --generate_flamegraph  生成火焰图 (默认: $params.generate_flamegraph)
      
    FastCall2参数:
      --min_depth     最小深度 (默认: $params.min_depth)
      --min_qual      最小质量 (默认: $params.min_qual)
      --ploidy        倍性 (默认: $params.ploidy)
      --help          显示帮助信息
    """
}

// 1. 准备不同样本大小的BAM映射文件
process prepareBamMaps {
    tag "Prepare BAM maps for ${sample_size} samples"
    publishDir "${params.outdir}/bam_maps", mode: 'copy'
    
    input:
    val sample_size
    path original_bam_map
    
    output:
    tuple val(sample_size), path("${sample_size}.tbm.txt")
    
    script:
    """
    head -n ${sample_size} ${original_bam_map} > ${sample_size}.tbm.txt
    """
}

// 2. 执行FastCall2 disc模式的性能分析
process perfAnalysisDisc {
    tag "Perf analysis disc: ${sample_size} samples"
    publishDir "${params.outdir}/perf_logs", mode: 'copy'
    
    input:
    tuple val(sample_size), path(bam_map)
    path reference
    path jar
    
    output:
    tuple val(sample_size), val("disc"), path("perf.${sample_size}.disc.log")
    path "${sample_size}disc/"
    
    script:
    """
    mkdir -p ${sample_size}disc
    
    perf stat -e ${params.perf_events} \\
        java -Xmx${params.memory} -jar ${jar} \\
        -app FastCall2 -mod disc \\
        -a ${reference} \\
        -b ${bam_map} \\
        -c 0 -d ${params.min_depth} -e ${params.min_qual} \\
        -f ${params.ploidy} -g ${params.min_allele_freq} \\
        -h ${params.max_missing} -i ${params.min_maf} \\
        -j ${params.min_call_rate} -k ${params.min_het_freq} \\
        -l ${params.chunk_size} -m ${params.threads} \\
        -n ${sample_size}disc/ \\
        -o ${params.samtools} \\
        > perf.${sample_size}.disc.log 2>&1
    """
}

// 3. 执行FastCall2 blib模式的性能分析
process perfAnalysisBlib {
    tag "Perf analysis blib: ${sample_size} samples"
    publishDir "${params.outdir}/perf_logs", mode: 'copy'
    
    input:
    tuple val(sample_size), path(disc_dir)
    path reference
    path jar
    
    output:
    tuple val(sample_size), val("blib"), path("perf.${sample_size}.blib.log")
    path "${sample_size}blib/"
    
    script:
    """
    mkdir -p ${sample_size}blib
    
    perf stat -e ${params.perf_events} \\
        java -Xmx${params.memory} -jar ${jar} \\
        -app FastCall2 -mod blib \\
        -a ${reference} \\
        -b 1 -c 2 -d ${params.threads} \\
        -e ${disc_dir} \\
        -f ${sample_size}blib/ \\
        > perf.${sample_size}.blib.log 2>&1
    """
}

// 4. 执行FastCall2 scan模式的性能分析
process perfAnalysisScan {
    tag "Perf analysis scan: ${sample_size} samples"
    publishDir "${params.outdir}/perf_logs", mode: 'copy'
    
    input:
    tuple val(sample_size), path(bam_map), path(blib_dir)
    path reference
    path jar
    
    output:
    tuple val(sample_size), val("scan"), path("perf.${sample_size}.scan.log")
    path "${sample_size}scan/"
    
    script:
    """
    mkdir -p ${sample_size}scan
    
    # 查找lib文件
    lib_file=\$(find ${blib_dir} -name "*.lib.gz" | head -1)
    
    perf stat -e ${params.perf_events} \\
        java -Xmx${params.memory} -jar ${jar} \\
        -app FastCall2 -mod scan \\
        -a ${reference} \\
        -b ${bam_map} \\
        -c \$lib_file \\
        -d 1 -e 0 -f ${params.min_depth} -g ${params.min_qual} \\
        -h 0.05 -i ${params.samtools} \\
        -j ${params.threads} \\
        -k ${sample_size}scan/ \\
        > perf.${sample_size}.scan.log 2>&1
    """
}

// 5. 生成详细的perf record分析
process perfRecord {
    tag "Perf record: ${sample_size} samples"
    publishDir "${params.outdir}/perf_record", mode: 'copy'
    
    input:
    tuple val(sample_size), path(bam_map)
    path reference
    path jar
    
    output:
    tuple val(sample_size), path("record.${sample_size}.perf.data"), path("record.${sample_size}.log")
    
    when:
    params.generate_flamegraph
    
    script:
    """
    mkdir -p ${sample_size}disc_record
    
    perf record -o record.${sample_size}.perf.data \\
        java -Xmx${params.memory} -jar ${jar} \\
        -app FastCall2 -mod disc \\
        -a ${reference} \\
        -b ${bam_map} \\
        -c 0 -d ${params.min_depth} -e ${params.min_qual} \\
        -f ${params.ploidy} -g ${params.min_allele_freq} \\
        -h ${params.max_missing} -i ${params.min_maf} \\
        -j ${params.min_call_rate} -k ${params.min_het_freq} \\
        -l ${params.chunk_size} -m ${params.threads} \\
        -n ${sample_size}disc_record/ \\
        -o ${params.samtools} \\
        > record.${sample_size}.log 2>&1
    """
}

// 6. 生成火焰图
process generateFlameGraph {
    tag "Flame graph: ${sample_size} samples"
    publishDir "${params.outdir}/flamegraphs", mode: 'copy'
    
    input:
    tuple val(sample_size), path(perf_data), path(log_file)
    
    output:
    path "flame.${sample_size}.svg"
    
    when:
    params.generate_flamegraph
    
    script:
    """
    # 检查是否安装了FlameGraph工具
    if ! command -v stackcollapse-perf.pl &> /dev/null; then
        echo "Warning: FlameGraph tools not found. Please install from:"
        echo "https://github.com/brendangregg/FlameGraph"
        touch flame.${sample_size}.svg
    else
        perf script -i ${perf_data} | \\
        stackcollapse-perf.pl | \\
        flamegraph.pl > flame.${sample_size}.svg
    fi
    """
}

// 7. 分析性能指标
process analyzePerformance {
    tag "Analyze performance metrics"
    publishDir "${params.outdir}/analysis", mode: 'copy'
    
    input:
    path "perf_logs/*"
    
    output:
    path "ipc_analysis.txt"
    path "branch_analysis.txt"
    path "timing_analysis.txt"
    path "performance_summary.md"
    
    script:
    """
    # 提取IPC信息
    grep "insn per cycle" perf_logs/*.disc.log | awk -F':' '{print \$1, \$2}' > ipc_disc.txt
    grep "insn per cycle" perf_logs/*.blib.log | awk -F':' '{print \$1, \$2}' > ipc_blib.txt
    grep "insn per cycle" perf_logs/*.scan.log | awk -F':' '{print \$1, \$2}' > ipc_scan.txt
    cat ipc_*.txt > ipc_analysis.txt
    
    # 提取分支预测信息
    grep "of all branches" perf_logs/*.disc.log | awk -F':' '{print \$1, \$2}' > branch_disc.txt
    grep "of all branches" perf_logs/*.blib.log | awk -F':' '{print \$1, \$2}' > branch_blib.txt
    grep "of all branches" perf_logs/*.scan.log | awk -F':' '{print \$1, \$2}' > branch_scan.txt
    cat branch_*.txt > branch_analysis.txt
    
    # 提取时间信息
    grep "seconds time elapsed" perf_logs/*.log | awk -F':' '{print \$1, \$2}' > timing_analysis.txt
    
    # 生成性能总结报告
    cat > performance_summary.md << 'EOF'
# FastCall2 性能分析报告

## 分析概述
本报告包含了FastCall2三个模式(disc/blib/scan)在不同样本大小下的性能表现。

## 关键指标
- **IPC (Instructions Per Cycle)**: 指令执行效率，通常>1.0为好
- **Branch Miss Rate**: 分支预测失败率，低于1%为佳
- **Cache Miss**: 缓存未命中次数
- **执行时间**: 各模式的实际运行时间

## 优化建议
1. 如果IPC较低，考虑优化算法逻辑
2. 如果Cache Miss较高，考虑优化数据结构和访问模式
3. 如果Branch Miss较高，考虑减少条件分支

## 详细数据
请查看对应的分析文件：
- ipc_analysis.txt: IPC详细数据
- branch_analysis.txt: 分支预测详细数据
- timing_analysis.txt: 时间详细数据
EOF
    """
}

// 定义工作流
workflow {
    // 显示帮助信息
    if (params.help) {
        printHelp()
        exit 0
    }
    
    // 检查必要文件
    jar_file = file(params.jar)
    reference_file = file(params.reference)
    bam_map_file = file(params.bam_maps)
    
    if (!jar_file.exists()) {
        error "TIGER jar文件不存在: ${params.jar}"
    }
    if (!reference_file.exists()) {
        error "参考基因组文件不存在: ${params.reference}"
    }
    if (!bam_map_file.exists()) {
        error "BAM映射文件不存在: ${params.bam_maps}"
    }
    if (!file(params.samtools).exists()) {
        error "samtools文件不存在: ${params.samtools}"
    }
    
    // 打印参数信息
    log.info """
    ==================================================
    FastCall2 性能分析流程
    ==================================================
    JAR文件: ${params.jar}
    参考基因组: ${params.reference}
    BAM映射文件: ${params.bam_maps}
    输出目录: ${params.outdir}
    线程数: ${params.threads}
    内存: ${params.memory}
    样本大小: ${params.sample_sizes}
    性能事件: ${params.perf_events}
    生成火焰图: ${params.generate_flamegraph}
    """
    
    // 创建样本大小通道
    sample_sizes_ch = Channel
        .from(params.sample_sizes.split(','))
        .map { it.toInteger() }
    
    // 1. 准备不同样本大小的BAM映射文件
    prepareBamMaps(sample_sizes_ch, bam_map_file)
    
    // 2. 执行disc模式性能分析
    perfAnalysisDisc(prepareBamMaps.out, reference_file, jar_file)
    
    // 3. 执行blib模式性能分析
    blib_input = perfAnalysisDisc.out
        .map { sample_size, _mode, _log -> 
            tuple(sample_size, file("${params.outdir}/perf_logs/${sample_size}disc"))
        }
    perfAnalysisBlib(blib_input, reference_file, jar_file)
    
    // 4. 执行scan模式性能分析
    scan_input = prepareBamMaps.out
        .join(perfAnalysisBlib.out.map { sample_size, _mode, _log -> 
            tuple(sample_size, file("${params.outdir}/perf_logs/${sample_size}blib"))
        })
    perfAnalysisScan(scan_input, reference_file, jar_file)
    
    // 5. 生成perf record数据（可选）
    if (params.generate_flamegraph) {
        perfRecord(prepareBamMaps.out, reference_file, jar_file)
        generateFlameGraph(perfRecord.out)
    }
    
    // 6. 收集所有性能日志并分析
    all_perf_logs = perfAnalysisDisc.out
        .mix(perfAnalysisBlib.out)
        .mix(perfAnalysisScan.out)
        .map { _sample_size, _mode, log -> log }
        .collect()
    
    analyzePerformance(all_perf_logs)
    
    // 工作流完成通知
    workflow.onComplete = {
        log.info "性能分析流程完成时间: $workflow.complete"
        log.info "执行状态: ${ workflow.success ? '成功' : '失败' }"
        log.info "执行时长: $workflow.duration"
        log.info "结果目录: ${params.outdir}"
        
        if (workflow.success) {
            log.info """
            
            分析结果已生成:
            - 性能日志: ${params.outdir}/perf_logs/
            - 性能分析: ${params.outdir}/analysis/
            ${params.generate_flamegraph ? "- 火焰图: ${params.outdir}/flamegraphs/" : ""}
            ${params.generate_flamegraph ? "- perf记录: ${params.outdir}/perf_record/" : ""}
            """
        }
    }
}
