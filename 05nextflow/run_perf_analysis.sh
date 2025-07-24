#!/bin/bash

# FastCall2 性能分析启动脚本
# 使用方法: ./run_perf_analysis.sh [选项]

set -e

# 默认参数
JAR_FILE="./TIGER_20250526.jar"
REFERENCE="./data/input/chr1_10M.fa"
BAM_MAPS="./data/input/taxaBamMap.txt"
SAMTOOLS="./bin/samtools"
OUTDIR="./results"
THREADS=32
MEMORY="100g"
SAMPLE_SIZES="1,2,3,4,5,6"
GENERATE_FLAMEGRAPH=true
PROFILE="standard"

# 显示帮助信息
show_help() {
    cat << EOF
FastCall2 性能分析启动脚本

用法: $0 [选项]

选项:
    -j, --jar FILE          TIGER jar文件路径 (默认: $JAR_FILE)
    -r, --reference FILE    参考基因组文件路径 (默认: $REFERENCE)  
    -b, --bam-maps FILE     BAM映射文件路径 (默认: $BAM_MAPS)
    -s, --samtools FILE     samtools可执行文件路径 (默认: $SAMTOOLS)
    -o, --outdir DIR        输出目录 (默认: $OUTDIR)
    -t, --threads NUM       线程数 (默认: $THREADS)
    -m, --memory SIZE       内存大小 (默认: $MEMORY)
    -n, --samples LIST      样本大小列表 (默认: $SAMPLE_SIZES)
    -f, --no-flamegraph     不生成火焰图
    -p, --profile NAME      使用的配置profile (默认: $PROFILE)
    -h, --help              显示此帮助信息

示例:
    # 基本运行
    $0
    
    # 自定义参数运行
    $0 -j /path/to/TIGER.jar -t 64 -m 200g -n "1,2,4,8"
    
    # 不生成火焰图（快速运行）
    $0 -f
    
    # 使用集群配置
    $0 -p cluster

EOF
}

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        -j|--jar)
            JAR_FILE="$2"
            shift 2
            ;;
        -r|--reference)
            REFERENCE="$2"
            shift 2
            ;;
        -b|--bam-maps)
            BAM_MAPS="$2"
            shift 2
            ;;
        -s|--samtools)
            SAMTOOLS="$2"
            shift 2
            ;;
        -o|--outdir)
            OUTDIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -m|--memory)
            MEMORY="$2"
            shift 2
            ;;
        -n|--samples)
            SAMPLE_SIZES="$2"
            shift 2
            ;;
        -f|--no-flamegraph)
            GENERATE_FLAMEGRAPH=false
            shift
            ;;
        -p|--profile)
            PROFILE="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo "未知选项: $1"
            show_help
            exit 1
            ;;
    esac
done

# 检查必要文件
echo "检查必要文件..."

if [[ ! -f "$JAR_FILE" ]]; then
    echo "错误: TIGER jar文件不存在: $JAR_FILE"
    exit 1
fi

if [[ ! -f "$REFERENCE" ]]; then
    echo "错误: 参考基因组文件不存在: $REFERENCE"
    exit 1
fi

if [[ ! -f "$BAM_MAPS" ]]; then
    echo "错误: BAM映射文件不存在: $BAM_MAPS"
    exit 1
fi

if [[ ! -f "$SAMTOOLS" ]]; then
    echo "错误: samtools可执行文件不存在: $SAMTOOLS"
    exit 1
fi

# 检查perf工具
if ! command -v perf &> /dev/null; then
    echo "错误: 未找到perf工具，请先安装perf"
    echo "在CentOS/RHEL上: sudo yum install perf"
    echo "在Ubuntu/Debian上: sudo apt install linux-tools-generic"
    exit 1
fi

# 检查Nextflow
if ! command -v nextflow &> /dev/null; then
    echo "错误: 未找到nextflow，请先安装Nextflow"
    echo "安装方法: curl -s https://get.nextflow.io | bash"
    exit 1
fi

# 检查Java
if ! command -v java &> /dev/null; then
    echo "错误: 未找到Java，请先安装Java"
    exit 1
fi

# 创建输出目录
mkdir -p "$OUTDIR"

# 显示配置信息
cat << EOF

==================================================
FastCall2 性能分析配置
==================================================
JAR文件:      $JAR_FILE
参考基因组:    $REFERENCE
BAM映射文件:  $BAM_MAPS
samtools:     $SAMTOOLS
输出目录:     $OUTDIR
线程数:       $THREADS
内存:         $MEMORY
样本大小:     $SAMPLE_SIZES
生成火焰图:   $GENERATE_FLAMEGRAPH
配置profile:  $PROFILE
==================================================

EOF

# 构建nextflow命令
NEXTFLOW_CMD="nextflow run fastcall.nf"
NEXTFLOW_CMD="$NEXTFLOW_CMD --jar $JAR_FILE"
NEXTFLOW_CMD="$NEXTFLOW_CMD --reference $REFERENCE"
NEXTFLOW_CMD="$NEXTFLOW_CMD --bam_maps $BAM_MAPS"
NEXTFLOW_CMD="$NEXTFLOW_CMD --samtools $SAMTOOLS"
NEXTFLOW_CMD="$NEXTFLOW_CMD --outdir $OUTDIR"
NEXTFLOW_CMD="$NEXTFLOW_CMD --threads $THREADS"
NEXTFLOW_CMD="$NEXTFLOW_CMD --memory $MEMORY"
NEXTFLOW_CMD="$NEXTFLOW_CMD --sample_sizes $SAMPLE_SIZES"
NEXTFLOW_CMD="$NEXTFLOW_CMD --generate_flamegraph $GENERATE_FLAMEGRAPH"
NEXTFLOW_CMD="$NEXTFLOW_CMD -profile $PROFILE"
NEXTFLOW_CMD="$NEXTFLOW_CMD -with-report $OUTDIR/reports/execution_report.html"
NEXTFLOW_CMD="$NEXTFLOW_CMD -with-timeline $OUTDIR/reports/timeline.html"
NEXTFLOW_CMD="$NEXTFLOW_CMD -with-dag $OUTDIR/reports/dag.html"

echo "执行命令: $NEXTFLOW_CMD"
echo ""

# 记录开始时间
START_TIME=$(date)
echo "开始时间: $START_TIME"

# 执行nextflow
eval $NEXTFLOW_CMD

# 记录结束时间
END_TIME=$(date)
echo ""
echo "结束时间: $END_TIME"

echo ""
echo "性能分析完成！"
echo "结果保存在: $OUTDIR"
echo ""
echo "主要输出文件:"
echo "- 性能日志: $OUTDIR/perf_logs/"
echo "- 分析结果: $OUTDIR/analysis/"
echo "- 执行报告: $OUTDIR/reports/"
if [[ "$GENERATE_FLAMEGRAPH" == "true" ]]; then
    echo "- 火焰图: $OUTDIR/flamegraphs/"
fi
