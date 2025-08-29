nextflow.enable.dsl = 2

// 输入文件夹参数
params.input_dir = "$projectDir/input"

// 默认参数设置 - 基于输入文件夹
params.jar = "${params.input_dir}/TIGER_20250526.jar"
params.reference = "${params.input_dir}/chr1_10M.fa"
params.bam_maps = "${params.input_dir}/taxaBamMap.txt"
params.samtools = System.getenv("CONDA_PREFIX") ? "${System.getenv("CONDA_PREFIX")}/bin/samtools" : "samtools"
params.outdir = "${params.input_dir}/results"
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
params.basic_events = "cycles,instructions,cache-references,cache-misses,branches,branch-misses"
params.cpu_events = "cpu-cycles,ref-cycles,instructions,stalled-cycles-frontend,stalled-cycles-backend"
params.memory_events = "cache-references,cache-misses,LLC-loads,LLC-load-misses,LLC-stores,LLC-store-misses,dTLB-loads,dTLB-load-misses,iTLB-loads,iTLB-load-misses"
params.branch_events = "branches,branch-misses,branch-loads,branch-load-misses"
params.advanced_events = "context-switches,cpu-migrations,page-faults,minor-faults,major-faults"
params.sample_sizes = "1,2,3,4,5,6"
params.generate_flamegraph = true
params.comprehensive_analysis = true

// 定义打印帮助信息的函数
def printHelp() {
    log.info """
    ==================================================
    FastCall2 性能分析流程
    ==================================================
    
    用法: nextflow run basic.nf [选项]
    
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
      --outdir        输出目录 (默认: ${params.input_dir}/results)
      --threads       CPU线程数 (默认: $params.threads)
      --memory        内存设置 (默认: $params.memory)
    
    性能分析:
      --perf_events   perf监控事件 (默认: $params.perf_events)
      --sample_sizes  样本大小列表 (默认: $params.sample_sizes)
      --generate_flamegraph  生成火焰图 (默认: $params.generate_flamegraph)
      --comprehensive_analysis  启用全面分析 (默认: $params.comprehensive_analysis)
      
    详细perf事件:
      --basic_events    基础性能事件 (默认: $params.basic_events)
      --cpu_events      CPU性能事件 (默认: $params.cpu_events) 
      --memory_events   内存性能事件 (默认: $params.memory_events)
      --branch_events   分支预测事件 (默认: $params.branch_events)
      --advanced_events 高级系统事件 (默认: $params.advanced_events)
      
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
    def events = params.comprehensive_analysis ? "${params.cpu_events},${params.memory_events},${params.advanced_events}" : params.perf_events
    """
    mkdir -p ${sample_size}disc
    
    perf stat -e ${events} \\
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
    def events = params.comprehensive_analysis ? "${params.cpu_events},${params.memory_events},${params.advanced_events}" : params.perf_events
    """
    mkdir -p ${sample_size}blib
    
    perf stat -e ${events} \\
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
    def events = params.comprehensive_analysis ? "${params.cpu_events},${params.memory_events},${params.advanced_events}" : params.perf_events
    """
    mkdir -p ${sample_size}scan
    
    # 查找lib文件
    lib_file=\$(find ${blib_dir} -name "*.lib.gz" | head -1)
    
    perf stat -e ${events} \\
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

// 5.1 增强的perf record进程 - 多维度分析
process enhancedPerfRecord {
    tag "Enhanced perf record: ${sample_size} samples"
    publishDir "${params.outdir}/perf_record_enhanced", mode: 'copy'
    
    input:
    tuple val(sample_size), path(bam_map)
    path reference
    path jar
    
    output:
    tuple val(sample_size), path("cpu.${sample_size}.perf.data"), path("memory.${sample_size}.perf.data"), path("branch.${sample_size}.perf.data"), path("enhanced.${sample_size}.log")
    
    when:
    params.comprehensive_analysis
    
    script:
    """
    mkdir -p ${sample_size}disc_enhanced
    
    # CPU性能记录 - 高频采样
    perf record -F 997 -g --call-graph dwarf \\
        -o cpu.${sample_size}.perf.data \\
        java -Xmx${params.memory} -jar ${jar} \\
        -app FastCall2 -mod disc \\
        -a ${reference} -b ${bam_map} \\
        -c 0 -d ${params.min_depth} -e ${params.min_qual} \\
        -f ${params.ploidy} -g ${params.min_allele_freq} \\
        -h ${params.max_missing} -i ${params.min_maf} \\
        -j ${params.min_call_rate} -k ${params.min_het_freq} \\
        -l ${params.chunk_size} -m ${params.threads} \\
        -n ${sample_size}disc_enhanced/ \\
        -o ${params.samtools} \\
        > cpu_record.${sample_size}.log 2>&1
    
    # 内存访问记录
    perf record -e cache-misses -g \\
        -o memory.${sample_size}.perf.data \\
        java -Xmx${params.memory} -jar ${jar} \\
        -app FastCall2 -mod disc \\
        -a ${reference} -b ${bam_map} \\
        -c 0 -d ${params.min_depth} -e ${params.min_qual} \\
        -f ${params.ploidy} -g ${params.min_allele_freq} \\
        -h ${params.max_missing} -i ${params.min_maf} \\
        -j ${params.min_call_rate} -k ${params.min_het_freq} \\
        -l ${params.chunk_size} -m ${params.threads} \\
        -n ${sample_size}disc_enhanced/ \\
        -o ${params.samtools} \\
        > memory_record.${sample_size}.log 2>&1
    
    # 分支预测记录
    perf record -e branch-misses -g \\
        -o branch.${sample_size}.perf.data \\
        java -Xmx${params.memory} -jar ${jar} \\
        -app FastCall2 -mod disc \\
        -a ${reference} -b ${bam_map} \\
        -c 0 -d ${params.min_depth} -e ${params.min_qual} \\
        -f ${params.ploidy} -g ${params.min_allele_freq} \\
        -h ${params.max_missing} -i ${params.min_maf} \\
        -j ${params.min_call_rate} -k ${params.min_het_freq} \\
        -l ${params.chunk_size} -m ${params.threads} \\
        -n ${sample_size}disc_enhanced/ \\
        -o ${params.samtools} \\
        > branch_record.${sample_size}.log 2>&1
    
    # 合并日志
    cat cpu_record.${sample_size}.log memory_record.${sample_size}.log branch_record.${sample_size}.log > enhanced.${sample_size}.log
    """
}

// 5.2 全面的性能统计分析
process comprehensivePerfStat {
    tag "Comprehensive perf stat: ${sample_size} samples, ${mode} mode"
    publishDir "${params.outdir}/perf_comprehensive", mode: 'copy'
    
    input:
    tuple val(sample_size), val(mode), path(bam_map)
    path reference
    path jar
    
    output:
    tuple val(sample_size), val(mode), path("comprehensive.${sample_size}.${mode}.json")
    
    when:
    params.comprehensive_analysis
    
    script:
    def cmd_args = ""
    if (mode == "disc") {
        cmd_args = "-app FastCall2 -mod disc -a ${reference} -b ${bam_map} -c 0 -d ${params.min_depth} -e ${params.min_qual} -f ${params.ploidy} -g ${params.min_allele_freq} -h ${params.max_missing} -i ${params.min_maf} -j ${params.min_call_rate} -k ${params.min_het_freq} -l ${params.chunk_size} -m ${params.threads} -n ${sample_size}${mode}/ -o ${params.samtools}"
    } else if (mode == "blib") {
        cmd_args = "-app FastCall2 -mod blib -a ${reference} -b 1 -c 2 -d ${params.threads} -e ${sample_size}disc/ -f ${sample_size}${mode}/"
    } else if (mode == "scan") {
        cmd_args = "-app FastCall2 -mod scan -a ${reference} -b ${bam_map} -c lib_file -d 1 -e 0 -f ${params.min_depth} -g ${params.min_qual} -h 0.05 -i ${params.samtools} -j ${params.threads} -k ${sample_size}${mode}/"
    }
    
    """
    mkdir -p ${sample_size}${mode}
    
    # 全面的性能统计 - JSON格式输出
    perf stat -j \\
        -e cycles,instructions,cache-references,cache-misses \\
        -e branches,branch-misses,stalled-cycles-frontend,stalled-cycles-backend \\
        -e L1-dcache-loads,L1-dcache-load-misses,L1-dcache-stores \\
        -e L1-icache-loads,L1-icache-load-misses \\
        -e LLC-loads,LLC-load-misses,LLC-stores,LLC-store-misses \\
        -e dTLB-loads,dTLB-load-misses,iTLB-loads,iTLB-load-misses \\
        -e context-switches,cpu-migrations,page-faults \\
        -e task-clock,cpu-clock \\
        java -Xmx${params.memory} -jar ${jar} ${cmd_args} \\
        > comprehensive.${sample_size}.${mode}.json 2>&1
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

// 6.1 增强的火焰图生成
process enhancedFlameGraph {
    tag "Enhanced flame graph: ${sample_size} samples"
    publishDir "${params.outdir}/visualizations", mode: 'copy'
    
    input:
    tuple val(sample_size), path(cpu_data), path(memory_data), path(branch_data), path(log_file)
    
    output:
    path "flame_*.svg"
    path "icicle_*.svg"
    path "diff_flame_*.svg"
    
    when:
    params.comprehensive_analysis
    
    script:
    """
    # 生成CPU火焰图
    if command -v stackcollapse-perf.pl &> /dev/null; then
        perf script -i ${cpu_data} | \\
        stackcollapse-perf.pl | \\
        flamegraph.pl --title="FastCall2 CPU Profile (${sample_size} samples)" > flame_cpu_${sample_size}.svg
        
        # 生成倒置火焰图（icicle图）
        perf script -i ${cpu_data} | \\
        stackcollapse-perf.pl | \\
        flamegraph.pl --inverted --title="FastCall2 CPU Icicle (${sample_size} samples)" > icicle_cpu_${sample_size}.svg
        
        # 生成内存访问火焰图
        perf script -i ${memory_data} | \\
        stackcollapse-perf.pl | \\
        flamegraph.pl --title="FastCall2 Memory Profile (${sample_size} samples)" > flame_memory_${sample_size}.svg
        
        # 生成分支预测火焰图
        perf script -i ${branch_data} | \\
        stackcollapse-perf.pl | \\
        flamegraph.pl --title="FastCall2 Branch Profile (${sample_size} samples)" > flame_branch_${sample_size}.svg
    else
        echo "Warning: FlameGraph tools not found"
        touch flame_cpu_${sample_size}.svg icicle_cpu_${sample_size}.svg flame_memory_${sample_size}.svg flame_branch_${sample_size}.svg diff_flame_placeholder.svg
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

// 7.1 高级性能分析和报告生成
process advancedPerformanceAnalysis {
    tag "Advanced performance analysis"
    publishDir "${params.outdir}/advanced_analysis", mode: 'copy'
    
    input:
    path "perf_data/*"
    path "perf_logs/*"
    path "perf_comprehensive/*"
    
    output:
    path "performance_report.html"
    path "hotspots_analysis.txt"
    path "memory_analysis.txt"
    path "cpu_utilization.txt"
    path "optimization_recommendations.md"
    
    when:
    params.comprehensive_analysis
    
    script:
    """
    # 生成热点分析
    echo "=== FastCall2 热点函数分析 ===" > hotspots_analysis.txt
    if ls perf_data/*.perf.data 1> /dev/null 2>&1; then
        for data_file in perf_data/*.perf.data; do
            echo "=== 分析文件: \$data_file ===" >> hotspots_analysis.txt
            perf report -i \$data_file --stdio --sort=overhead | head -30 >> hotspots_analysis.txt
            echo "" >> hotspots_analysis.txt
        done
    else
        echo "No perf data files found" >> hotspots_analysis.txt
    fi
    
    # 内存访问模式分析
    echo "=== 内存访问模式分析 ===" > memory_analysis.txt
    if ls perf_data/memory.*.perf.data 1> /dev/null 2>&1; then
        for data_file in perf_data/memory.*.perf.data; do
            echo "--- 文件: \$data_file ---" >> memory_analysis.txt
            perf report -i \$data_file --stdio --sort=overhead | head -20 >> memory_analysis.txt
            echo "" >> memory_analysis.txt
        done
    else
        echo "No memory perf data files found" >> memory_analysis.txt
    fi
    
    # CPU利用率分析
    echo "=== CPU利用率分析 ===" > cpu_utilization.txt
    if ls perf_comprehensive/*.json 1> /dev/null 2>&1; then
        grep -E "(task-clock|cpu-clock|cycles)" perf_comprehensive/*.json | \\
            awk '{print \$1, \$2}' >> cpu_utilization.txt
    else
        grep -E "(task-clock|cpu-clock|cycles)" perf_logs/*.log | \\
            awk '{print \$1, \$2}' >> cpu_utilization.txt
    fi
    
    # 生成优化建议
    cat > optimization_recommendations.md << 'EOF'
    # FastCall2 性能优化建议

    ## 1. CPU性能优化
    - **IPC分析**: 目标IPC > 1.0，当前需查看分析结果
    - **前端停顿**: 如果stalled-cycles-frontend过高，考虑优化指令获取
    - **后端停顿**: 如果stalled-cycles-backend过高，考虑优化执行单元使用

    ## 2. 内存性能优化
    - **L1缓存命中率**: L1-dcache-load-misses/L1-dcache-loads应 < 5%
    - **LLC命中率**: LLC-load-misses/LLC-loads应 < 20%
    - **TLB命中率**: dTLB-load-misses/dTLB-loads应 < 1%

    ## 3. 分支预测优化
    - **分支预测失败率**: branch-misses/branches应低于1%
    - **热点函数**: 重点优化调用频繁的函数

    ## 4. 并行性优化
    - **线程利用率**: 监控多线程效率
    - **上下文切换**: context-switches应保持较低水平
    - **CPU迁移**: cpu-migrations过多表明负载不均

    ## 5. I/O性能优化
    - **页错误**: page-faults中major-faults应尽量少
    - **内存使用**: 避免内存不足导致的swap
    - **文件I/O**: 优化BAM文件读取模式

    ## 6. Java特定优化
    - **垃圾回收**: 监控GC频率和停顿时间
    - **堆内存**: 调整-Xmx参数避免OOM
    - **JIT编译**: 关注热点方法的编译优化

    EOF
        
        # 生成HTML报告
        cat > performance_report.html << 'EOF'
    <!DOCTYPE html>
    <html>
    <head>
        <title>FastCall2 Performance Analysis Report</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 20px; line-height: 1.6; }
            .header { background: #f4f4f4; padding: 20px; border-radius: 5px; }
            .metric { background: #f9f9f9; padding: 15px; margin: 10px 0; border-left: 4px solid #007acc; }
            .warning { color: #ff6600; font-weight: bold; }
            .good { color: #00aa00; font-weight: bold; }
            .critical { color: #ff0000; font-weight: bold; }
            table { border-collapse: collapse; width: 100%; margin: 20px 0; }
            th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }
            th { background-color: #007acc; color: white; }
            .section { margin: 30px 0; }
            .code { background: #f5f5f5; padding: 10px; border-radius: 3px; font-family: monospace; }
        </style>
    </head>
    <body>
        <div class="header">
            <h1>FastCall2 Performance Analysis Report</h1>
            <p>Generated on: $(date)</p>
            <p>Analysis Type: Comprehensive Performance Profiling</p>
        </div>
        
        <div class="section">
            <h2>Executive Summary</h2>
            <div class="metric">
                <h3>分析范围</h3>
                <p>本报告包含FastCall2三个执行模式的全面性能分析：</p>
                <ul>
                    <li><strong>disc模式</strong>: 基因型发现和初步处理</li>
                    <li><strong>blib模式</strong>: 库文件构建</li>
                    <li><strong>scan模式</strong>: 最终扫描和输出</li>
                </ul>
            </div>
            
            <div class="metric">
                <h3>关键性能指标</h3>
                <ul>
                    <li>CPU效率 (IPC - Instructions Per Cycle)</li>
                    <li>内存访问模式 (Cache Miss Rates)</li>
                    <li>分支预测准确性 (Branch Miss Rate)</li>
                    <li>并行执行效率 (Thread Utilization)</li>
                    <li>I/O性能 (Page Faults, Context Switches)</li>
                </ul>
            </div>
        </div>
        
        <div class="section">
            <h2>Performance Metrics Overview</h2>
            <table>
                <tr><th>指标类别</th><th>测量项目</th><th>目标值</th><th>评估标准</th></tr>
                <tr><td rowspan="3">CPU性能</td><td>IPC</td><td>&gt; 1.0</td><td>指令执行效率</td></tr>
                <tr><td>前端停顿率</td><td>&lt; 10%</td><td>指令获取效率</td></tr>
                <tr><td>后端停顿率</td><td>&lt; 20%</td><td>执行单元利用率</td></tr>
                <tr><td rowspan="3">内存性能</td><td>L1缓存命中率</td><td>&gt; 95%</td><td>数据访问局部性</td></tr>
                <tr><td>LLC缓存命中率</td><td>&gt; 80%</td><td>内存层次效率</td></tr>
                <tr><td>TLB命中率</td><td>&gt; 99%</td><td>地址转换效率</td></tr>
                <tr><td rowspan="2">分支预测</td><td>分支预测命中率</td><td>&gt; 99%</td><td>控制流预测准确性</td></tr>
                <tr><td>间接分支命中率</td><td>&gt; 95%</td><td>虚函数调用效率</td></tr>
                <tr><td rowspan="2">并发性能</td><td>上下文切换率</td><td>适中</td><td>线程调度开销</td></tr>
                <tr><td>CPU迁移次数</td><td>较低</td><td>NUMA感知调度</td></tr>
            </table>
        </div>
        
        <div class="section">
            <h2>Analysis Results</h2>
            <div class="metric">
                <h3>数据文件位置</h3>
                <p>详细分析数据请查看以下文件：</p>
                <ul>
                    <li><strong><a href="hotspots_analysis.txt">hotspots_analysis.txt</a></strong> - 热点函数分析</li>
                    <li><strong><a href="memory_analysis.txt">memory_analysis.txt</a></strong> - 内存访问模式分析</li>
                    <li><strong><a href="cpu_utilization.txt">cpu_utilization.txt</a></strong> - CPU利用率统计</li>
                    <li><strong><a href="optimization_recommendations.md">optimization_recommendations.md</a></strong> - 优化建议</li>
                </ul>
            </div>
            
            <div class="metric">
                <h3>可视化结果</h3>
                <p>如果启用了可视化选项，请查看：</p>
                <ul>
                    <li>火焰图 (FlameGraphs) - CPU热点可视化</li>
                    <li>冰柱图 (Icicle Graphs) - 调用栈倒置视图</li>
                    <li>差异图 (Diff Graphs) - 不同样本大小的性能对比</li>
                </ul>
            </div>
        </div>
        
        <div class="section">
            <h2>Quick Actions</h2>
            <div class="metric">
                <h3>立即检查项目</h3>
                <ol>
                    <li>查看 <span class="code">hotspots_analysis.txt</span> 中的Top 10热点函数</li>
                    <li>检查 <span class="code">memory_analysis.txt</span> 中的缓存未命中模式</li>
                    <li>参考 <span class="code">optimization_recommendations.md</span> 中的具体优化建议</li>
                    <li>如有火焰图，重点关注宽度最大的函数调用</li>
                </ol>
            </div>
        </div>
        
        <div class="section">
            <h2>Troubleshooting</h2>
            <div class="metric">
                <h3>常见性能问题及解决方案</h3>
                <table>
                    <tr><th>问题症状</th><th>可能原因</th><th>解决方案</th></tr>
                    <tr><td>IPC &lt; 0.5</td><td>大量内存等待</td><td>优化数据结构，增加内存</td></tr>
                    <tr><td>高L1缓存未命中</td><td>数据访问模式差</td><td>重组数据结构，提高局部性</td></tr>
                    <tr><td>高分支预测失败</td><td>复杂条件逻辑</td><td>简化分支，使用查找表</td></tr>
                    <tr><td>高上下文切换</td><td>线程竞争激烈</td><td>减少锁竞争，优化线程数</td></tr>
                </table>
            </div>
        </div>
    </body>
    </html>
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
    
    // 5.1 增强的perf record分析（如果启用全面分析）
    if (params.comprehensive_analysis) {
        enhancedPerfRecord(prepareBamMaps.out, reference_file, jar_file)
        enhancedFlameGraph(enhancedPerfRecord.out)
        
        // 全面的性能统计分析
        disc_input_for_comprehensive = prepareBamMaps.out.map { sample_size, bam_map -> 
            tuple(sample_size, "disc", bam_map) 
        }
        comprehensivePerfStat(disc_input_for_comprehensive, reference_file, jar_file)
    }
    
    // 6. 收集所有性能日志并分析
    all_perf_logs = perfAnalysisDisc.out
        .mix(perfAnalysisBlib.out)
        .mix(perfAnalysisScan.out)
        .map { _sample_size, _mode, log -> log }
        .collect()
    
    analyzePerformance(all_perf_logs)
    
    // 6.1 高级性能分析（如果启用全面分析）
    if (params.comprehensive_analysis) {
        // 收集所有perf数据文件
        all_perf_data = Channel.empty()
        all_comprehensive_data = Channel.empty()
        
        if (params.generate_flamegraph) {
            all_perf_data = enhancedPerfRecord.out
                .map { _sample_size, cpu_data, memory_data, branch_data, _log -> 
                    [cpu_data, memory_data, branch_data] 
                }
                .flatten()
                .collect()
        } else {
            all_perf_data = Channel.of([]).collect()
        }
        
        all_comprehensive_data = comprehensivePerfStat.out
            .map { _sample_size, _mode, json -> json }
            .collect()
        
        advancedPerformanceAnalysis(all_perf_data, all_perf_logs, all_comprehensive_data)
    }
    
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
            ${params.comprehensive_analysis ? "- 增强记录: ${params.outdir}/perf_record_enhanced/" : ""}
            ${params.comprehensive_analysis ? "- 全面统计: ${params.outdir}/perf_comprehensive/" : ""}
            ${params.comprehensive_analysis ? "- 可视化图表: ${params.outdir}/visualizations/" : ""}
            ${params.comprehensive_analysis ? "- 高级分析: ${params.outdir}/advanced_analysis/" : ""}
            
            主要报告文件:
            ${params.comprehensive_analysis ? "- ${params.outdir}/advanced_analysis/performance_report.html (主报告)" : ""}
            - ${params.outdir}/analysis/performance_summary.md (基础总结)
            ${params.comprehensive_analysis ? "- ${params.outdir}/advanced_analysis/optimization_recommendations.md (优化建议)" : ""}
            """
        }
    }
}
