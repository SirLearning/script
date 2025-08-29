nextflow.enable.dsl = 2

// 输入文件夹参数
params.input_dir = "$projectDir/input"

// 默认参数设置 - 基于输入文件夹
params.jar = "${params.input_dir}/TIGER_20250526.jar"
params.reference = "${params.input_dir}/chr1_10M.fa"
params.bam_maps = "${params.input_dir}/taxaBamMap.txt"
params.samtools = System.getenv("CONDA_PREFIX") ? "${System.getenv("CONDA_PREFIX")}/bin/samtools" : "samtools"
params.outdir = "${params.input_dir}/results_callpath"
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

// 调用路径分析参数
params.call_graph_mode = "dwarf"  // dwarf, fp, lbr
params.sampling_frequency = 4000  // 采样频率
params.call_graph_depth = 10      // 调用图深度
params.cache_events = "L1-dcache-loads,L1-dcache-load-misses,L1-icache-loads,L1-icache-load-misses,LLC-loads,LLC-load-misses"
params.branch_events = "branches,branch-misses,branch-loads,branch-load-misses"
params.cpu_events = "cycles,instructions,cache-references,cache-misses"

// 可视化参数
params.generate_callgraph = true
params.generate_dotgraph = true
params.generate_heatmap = true
params.min_percentage = 0.5       // 最小显示百分比
params.sample_count = 3           // 分析的样本数量

// 定义打印帮助信息的函数
def printHelp() {
    log.info """
    ==================================================
    FastCall2 调用路径和缓存分析流程
    ==================================================
    
    用法: nextflow run perf_callpath_analysis.nf [选项]
    
    输入文件夹:
      --input_dir     输入文件夹路径 (默认: $params.input_dir)
                      该文件夹应包含以下文件:
                      - TIGER_20250526.jar (或其他版本的TIGER jar文件)
                      - chr1_10M.fa (或其他参考基因组文件)
                      - taxaBamMap.txt (BAM文件映射表)
    
    单独文件设置 (可选，会覆盖input_dir中的默认设置):
      --jar           TIGER jar文件路径 (默认: ${params.input_dir}/TIGER_20250526.jar)
      --reference     参考基因组路径 (默认: ${params.input_dir}/chr1_10M.fa)
      --bam_maps      BAM文件映射表 (默认: ${params.input_dir}/taxaBamMap.txt)
      --samtools      samtools路径 (默认: $params.samtools)
    
    输出设置:
      --outdir        输出目录 (默认: ${params.input_dir}/results_callpath)
      --threads       CPU线程数 (默认: $params.threads)
      --memory        内存设置 (默认: $params.memory)
    
    调用路径分析:
      --call_graph_mode     调用图模式 (默认: $params.call_graph_mode)
      --sampling_frequency  采样频率 (默认: $params.sampling_frequency)
      --call_graph_depth    调用图深度 (默认: $params.call_graph_depth)
      --sample_count        分析样本数 (默认: $params.sample_count)
      
    性能事件:
      --cache_events  缓存事件 (默认: $params.cache_events)
      --branch_events 分支事件 (默认: $params.branch_events)
      --cpu_events    CPU事件 (默认: $params.cpu_events)
      
    可视化选项:
      --generate_callgraph  生成调用图 (默认: $params.generate_callgraph)
      --generate_dotgraph   生成DOT图 (默认: $params.generate_dotgraph)
      --generate_heatmap    生成热力图 (默认: $params.generate_heatmap)
      --min_percentage      最小显示百分比 (默认: $params.min_percentage)
      
    FastCall2参数:
      --min_depth     最小深度 (默认: $params.min_depth)
      --min_qual      最小质量 (默认: $params.min_qual)
      --ploidy        倍性 (默认: $params.ploidy)
      --help          显示帮助信息
    """
}

// 1. 准备样本BAM映射文件
process prepareSamples {
    tag "Prepare samples: ${sample_count}"
    publishDir "${params.outdir}/samples", mode: 'copy'
    
    input:
    val sample_count
    path original_bam_map
    
    output:
    tuple val(sample_count), path("${sample_count}.tbm.txt")
    
    script:
    """
    head -n ${sample_count} ${original_bam_map} > ${sample_count}.tbm.txt
    """
}

// 2. 调用图性能记录 - 详细调用路径
process perfCallGraphRecord {
    tag "Call graph record: ${sample_count} samples"
    publishDir "${params.outdir}/callgraph_data", mode: 'copy'
    
    input:
    tuple val(sample_count), path(bam_map)
    path reference
    path jar
    
    output:
    tuple val(sample_count), path("callgraph.${sample_count}.perf.data"), path("callgraph.${sample_count}.log")
    
    script:
    """
    mkdir -p ${sample_count}_callgraph
    
    # 记录详细的调用图信息
    perf record \\
        -F ${params.sampling_frequency} \\
        -g --call-graph ${params.call_graph_mode},${params.call_graph_depth} \\
        --branch-filter any,u \\
        -o callgraph.${sample_count}.perf.data \\
        java -Xmx${params.memory} -jar ${jar} \\
        -app FastCall2 -mod disc \\
        -a ${reference} \\
        -b ${bam_map} \\
        -c 0 -d ${params.min_depth} -e ${params.min_qual} \\
        -f ${params.ploidy} -g ${params.min_allele_freq} \\
        -h ${params.max_missing} -i ${params.min_maf} \\
        -j ${params.min_call_rate} -k ${params.min_het_freq} \\
        -l ${params.chunk_size} -m ${params.threads} \\
        -n ${sample_count}_callgraph/ \\
        -o ${params.samtools} \\
        > callgraph.${sample_count}.log 2>&1
    """
}

// 3. 缓存性能记录 - CPU缓存分析
process perfCacheRecord {
    tag "Cache record: ${sample_count} samples"
    publishDir "${params.outdir}/cache_data", mode: 'copy'
    
    input:
    tuple val(sample_count), path(bam_map)
    path reference
    path jar
    
    output:
    tuple val(sample_count), path("cache.${sample_count}.perf.data"), path("cache.${sample_count}.log")
    
    script:
    """
    mkdir -p ${sample_count}_cache
    
    # 记录缓存性能数据
    perf record \\
        -e ${params.cache_events} \\
        -g --call-graph ${params.call_graph_mode} \\
        -o cache.${sample_count}.perf.data \\
        java -Xmx${params.memory} -jar ${jar} \\
        -app FastCall2 -mod disc \\
        -a ${reference} \\
        -b ${bam_map} \\
        -c 0 -d ${params.min_depth} -e ${params.min_qual} \\
        -f ${params.ploidy} -g ${params.min_allele_freq} \\
        -h ${params.max_missing} -i ${params.min_maf} \\
        -j ${params.min_call_rate} -k ${params.min_het_freq} \\
        -l ${params.chunk_size} -m ${params.threads} \\
        -n ${sample_count}_cache/ \\
        -o ${params.samtools} \\
        > cache.${sample_count}.log 2>&1
    """
}

// 4. 分支预测记录 - 执行路径分析
process perfBranchRecord {
    tag "Branch record: ${sample_count} samples"
    publishDir "${params.outdir}/branch_data", mode: 'copy'
    
    input:
    tuple val(sample_count), path(bam_map)
    path reference
    path jar
    
    output:
    tuple val(sample_count), path("branch.${sample_count}.perf.data"), path("branch.${sample_count}.log")
    
    script:
    """
    mkdir -p ${sample_count}_branch
    
    # 记录分支预测和执行路径
    perf record \\
        -e ${params.branch_events} \\
        -j any,u \\
        -g --call-graph ${params.call_graph_mode} \\
        -o branch.${sample_count}.perf.data \\
        java -Xmx${params.memory} -jar ${jar} \\
        -app FastCall2 -mod disc \\
        -a ${reference} \\
        -b ${bam_map} \\
        -c 0 -d ${params.min_depth} -e ${params.min_qual} \\
        -f ${params.ploidy} -g ${params.min_allele_freq} \\
        -h ${params.max_missing} -i ${params.min_maf} \\
        -j ${params.min_call_rate} -k ${params.min_het_freq} \\
        -l ${params.chunk_size} -m ${params.threads} \\
        -n ${sample_count}_branch/ \\
        -o ${params.samtools} \\
        > branch.${sample_count}.log 2>&1
    """
}

// 5. 生成调用图报告
process generateCallGraphReport {
    tag "Call graph report: ${sample_count} samples"
    publishDir "${params.outdir}/callgraph_reports", mode: 'copy'
    
    input:
    tuple val(sample_count), path(perf_data), path(log_file)
    
    output:
    path "callgraph_report_${sample_count}.txt"
    path "callgraph_tree_${sample_count}.txt"
    path "callgraph_folded_${sample_count}.txt"
    
    script:
    """
    # 生成调用图报告
    perf report -i ${perf_data} \\
        --stdio \\
        --sort=overhead,symbol \\
        --call-graph=graph,${params.min_percentage} \\
        > callgraph_report_${sample_count}.txt
    
    # 生成调用树
    perf report -i ${perf_data} \\
        --stdio \\
        --call-graph=tree,${params.min_percentage} \\
        > callgraph_tree_${sample_count}.txt
    
    # 生成折叠格式（用于后续可视化）
    perf script -i ${perf_data} | \\
    awk '
    /^[[:space:]]*[0-9a-f]+/ {
        if (NF >= 2) {
            symbol = \$2
            gsub(/^.*\\//, "", symbol)  # 移除路径前缀
            if (symbol != "" && symbol !~ /^[0-9a-f]+\$/) {
                stack[depth] = symbol
                depth++
            }
        }
    }
    /^[[:space:]]*\$/ {
        if (depth > 0) {
            call_stack = ""
            for (i = depth-1; i >= 0; i--) {
                if (call_stack == "") {
                    call_stack = stack[i]
                } else {
                    call_stack = stack[i] ";" call_stack
                }
            }
            if (call_stack != "") {
                print call_stack " 1"
            }
        }
        depth = 0
        delete stack
    }
    ' > callgraph_folded_${sample_count}.txt
    """
}

// 6. 生成DOT图
process generateDotGraph {
    tag "DOT graph: ${sample_count} samples"
    publishDir "${params.outdir}/dot_graphs", mode: 'copy'
    
    input:
    tuple val(sample_count), path(perf_data), path(log_file)
    path python_script
    
    output:
    path "callgraph_${sample_count}.dot"
    path "callgraph_${sample_count}.png"
    path "callgraph_${sample_count}.svg"
    
    when:
    params.generate_dotgraph
    
    script:
    """
    # 使用独立的Python脚本生成DOT格式的调用图
    perf script -i ${perf_data} | \\
    python3 ${python_script} callgraph > callgraph_${sample_count}.dot
    
    # 生成PNG和SVG图像
    if command -v dot &> /dev/null; then
        dot -Tpng callgraph_${sample_count}.dot -o callgraph_${sample_count}.png
        dot -Tsvg callgraph_${sample_count}.dot -o callgraph_${sample_count}.svg
    else
        echo "Graphviz not found. Only DOT file generated."
        touch callgraph_${sample_count}.png callgraph_${sample_count}.svg
    fi
    """
}

// 7. 缓存热力图分析
process generateCacheHeatmap {
    tag "Cache heatmap: ${sample_count} samples"
    publishDir "${params.outdir}/cache_analysis", mode: 'copy'
    
    input:
    tuple val(sample_count), path(cache_data), path(log_file)
    path python_script
    
    output:
    path "cache_report_${sample_count}.txt"
    path "cache_heatmap_${sample_count}.html"
    path "cache_analysis_${sample_count}.json"
    
    when:
    params.generate_heatmap
    
    script:
    """
    # 生成缓存报告
    perf report -i ${cache_data} \\
        --stdio \\
        --sort=overhead,symbol \\
        > cache_report_${sample_count}.txt
    
    # 使用独立的Python脚本分析缓存数据并生成热力图
    perf script -i ${cache_data} | \\
    python3 ${python_script} cache ${sample_count}
    """
}

// 8. 分支预测分析
process analyzeBranchPrediction {
    tag "Branch analysis: ${sample_count} samples"
    publishDir "${params.outdir}/branch_analysis", mode: 'copy'
    
    input:
    tuple val(sample_count), path(branch_data), path(log_file)
    path python_script
    
    output:
    path "branch_report_${sample_count}.txt"
    path "branch_flow_${sample_count}.dot"
    path "branch_flow_${sample_count}.svg"
    
    script:
    """
    # 生成分支预测报告
    perf report -i ${branch_data} \\
        --stdio \\
        --sort=overhead,symbol \\
        > branch_report_${sample_count}.txt
    
    # 使用独立的Python脚本分析分支流
    perf script -i ${branch_data} | \\
    python3 ${python_script} branch > branch_flow_${sample_count}.dot
    
    # 生成SVG
    if command -v dot &> /dev/null; then
        dot -Tsvg branch_flow_${sample_count}.dot -o branch_flow_${sample_count}.svg
    else
        touch branch_flow_${sample_count}.svg
    fi
    """
}

// 9. 综合分析报告
process generateComprehensiveReport {
    tag "Comprehensive analysis report"
    publishDir "${params.outdir}/comprehensive_analysis", mode: 'copy'
    
    input:
    path "callgraph_reports/*"
    path "cache_analysis/*"
    path "branch_analysis/*"
    path "dot_graphs/*"
    path python_script
    
    output:
    path "comprehensive_report.html"
    path "performance_insights.md"
    path "optimization_guide.md"
    
    script:
    """
    # 使用独立的Python脚本生成综合报告
    python3 ${python_script} reports
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
    python_script = file("$projectDir/perf_analysis_generator.py")
    
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
    if (!python_script.exists()) {
        error "Python分析脚本不存在: $projectDir/perf_analysis_generator.py"
    }
    
    // 打印参数信息
    log.info """
    ==================================================
    FastCall2 调用路径和缓存分析流程
    ==================================================
    JAR文件: ${params.jar}
    参考基因组: ${params.reference}
    BAM映射文件: ${params.bam_maps}
    输出目录: ${params.outdir}
    样本数量: ${params.sample_count}
    调用图模式: ${params.call_graph_mode}
    采样频率: ${params.sampling_frequency}
    生成DOT图: ${params.generate_dotgraph}
    生成热力图: ${params.generate_heatmap}
    """
    
    // 1. 准备样本
    prepareSamples(params.sample_count, bam_map_file)
    
    // 2. 并行执行三种性能记录
    perfCallGraphRecord(prepareSamples.out, reference_file, jar_file)
    perfCacheRecord(prepareSamples.out, reference_file, jar_file)
    perfBranchRecord(prepareSamples.out, reference_file, jar_file)
    
    // 3. 生成调用图报告
    generateCallGraphReport(perfCallGraphRecord.out)
    
    // 4. 生成DOT图（可选）
    if (params.generate_dotgraph) {
        generateDotGraph(perfCallGraphRecord.out, python_script)
    }
    
    // 5. 生成缓存热力图（可选）
    if (params.generate_heatmap) {
        generateCacheHeatmap(perfCacheRecord.out, python_script)
    }
    
    // 6. 分析分支预测
    analyzeBranchPrediction(perfBranchRecord.out, python_script)
    
    // 7. 收集所有分析结果并生成综合报告
    callgraph_reports = generateCallGraphReport.out.collect()
    cache_analysis = params.generate_heatmap ? generateCacheHeatmap.out.collect() : Channel.of([]).collect()
    branch_analysis = analyzeBranchPrediction.out.collect()
    dot_graphs = params.generate_dotgraph ? generateDotGraph.out.collect() : Channel.of([]).collect()
    
    generateComprehensiveReport(
        callgraph_reports,
        cache_analysis,
        branch_analysis,
        dot_graphs,
        python_script
    )
    
    // 工作流完成通知
    workflow.onComplete = {
        log.info "调用路径分析流程完成时间: $workflow.complete"
        log.info "执行状态: ${ workflow.success ? '成功' : '失败' }"
        log.info "执行时长: $workflow.duration"
        log.info "结果目录: ${params.outdir}"
        
        if (workflow.success) {
            log.info """
            
            分析结果已生成:
            - 调用图数据: ${params.outdir}/callgraph_data/
            - 缓存分析数据: ${params.outdir}/cache_data/
            - 分支分析数据: ${params.outdir}/branch_data/
            - 调用图报告: ${params.outdir}/callgraph_reports/
            ${params.generate_heatmap ? "- 缓存热力图: ${params.outdir}/cache_analysis/" : ""}
            - 分支分析: ${params.outdir}/branch_analysis/
            ${params.generate_dotgraph ? "- DOT图表: ${params.outdir}/dot_graphs/" : ""}
            - 综合报告: ${params.outdir}/comprehensive_analysis/
            
            主要报告文件:
            - ${params.outdir}/comprehensive_analysis/comprehensive_report.html (主报告)
            - ${params.outdir}/comprehensive_analysis/performance_insights.md (性能洞察)
            - ${params.outdir}/comprehensive_analysis/optimization_guide.md (优化指南)
            ${params.generate_dotgraph ? "- ${params.outdir}/dot_graphs/*.svg (调用图可视化)" : ""}
            ${params.generate_heatmap ? "- ${params.outdir}/cache_analysis/*.html (缓存热力图)" : ""}
            """
        }
    }
}
