#!/bin/bash

# Linux/MacOS脚本用于运行转座子分析流程

set -e  # 遇到错误时退出

echo "========================================"
echo "转座子分析流程 - Unix运行脚本"
echo "========================================"
echo

# 函数定义
show_help() {
    cat << EOF

用法: $0 [选项]

选项:
  --input FILE     指定输入数据文件 (默认: newplot.txt)
  --output DIR     指定输出目录 (默认: results)
  --profile PROF   指定执行配置 (默认: standard)
                   可选: standard, docker, conda, hpc, test
  --resume         从上次中断的地方恢复运行
  --help, -h       显示此帮助信息

示例:
  $0
  $0 --input mydata.txt --output my_results
  $0 --profile docker
  $0 --resume

EOF
}

# 检查Nextflow是否已安装
check_nextflow() {
    if ! command -v nextflow &> /dev/null; then
        echo "错误: 未找到Nextflow。请先安装Nextflow。"
        echo "安装指南: https://www.nextflow.io/docs/latest/getstarted.html"
        echo
        echo "快速安装:"
        echo "curl -s https://get.nextflow.io | bash"
        exit 1
    fi
    
    echo "Nextflow版本: $(nextflow -version | head -n1)"
}

# 检查必需文件
check_files() {
    echo "当前工作目录: $(pwd)"
    echo
    echo "检查必需文件..."
    
    if [[ ! -f "transposon_analysis_simple.nf" ]]; then
        echo "错误: 未找到主流程文件 transposon_analysis_simple.nf"
        exit 1
    fi
    
    if [[ ! -f "01compo_count_v1.py" ]]; then
        echo "错误: 未找到脚本文件 01compo_count_v1.py"
        exit 1
    fi
    
    if [[ ! -f "newplot.txt" ]] && [[ "$INPUT_DATA" == "newplot.txt" ]]; then
        echo "警告: 未找到默认输入文件 newplot.txt"
        echo "请确保数据文件在当前目录中或使用 --input 指定文件"
    fi
    
    echo "文件检查完成。"
}

# 默认参数
INPUT_DATA="newplot.txt"
OUTPUT_DIR="results"
PROFILE="standard"
RESUME=""

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case $1 in
        --input)
            INPUT_DATA="$2"
            shift 2
            ;;
        --output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --profile)
            PROFILE="$2"
            shift 2
            ;;
        --resume)
            RESUME="-resume"
            shift
            ;;
        --help|-h)
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

# 主程序
main() {
    check_nextflow
    check_files
    
    echo
    echo "运行参数:"
    echo "  输入文件: $INPUT_DATA"
    echo "  输出目录: $OUTPUT_DIR"
    echo "  执行配置: $PROFILE"
    if [[ -n "$RESUME" ]]; then
        echo "  恢复模式: 是"
    fi
    echo
    
    # 确认运行
    read -p "是否继续运行分析? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "分析已取消。"
        exit 0
    fi
    
    echo
    echo "开始运行转座子分析流程..."
    echo "========================================"
    
    # 记录开始时间
    start_time=$(date)
    
    # 运行Nextflow流程
    nextflow run transposon_analysis_simple.nf \
        --input_data "$INPUT_DATA" \
        --output_dir "$OUTPUT_DIR" \
        -profile "$PROFILE" \
        $RESUME
    
    exit_code=$?
    
    # 记录结束时间
    end_time=$(date)
    
    echo
    echo "========================================"
    if [[ $exit_code -eq 0 ]]; then
        echo "分析完成！"
        echo "开始时间: $start_time"
        echo "结束时间: $end_time"
        echo "结果保存在: $OUTPUT_DIR"
        echo
        if [[ -d "$OUTPUT_DIR" ]]; then
            echo "生成的文件:"
            find "$OUTPUT_DIR" -type f | head -20
            total_files=$(find "$OUTPUT_DIR" -type f | wc -l)
            if [[ $total_files -gt 20 ]]; then
                echo "... 还有 $((total_files - 20)) 个文件"
            fi
        fi
    else
        echo "分析失败！错误代码: $exit_code"
        echo "请检查日志文件以获取详细信息："
        echo "  nextflow log last"
    fi
    echo "========================================"
    
    exit $exit_code
}

# 运行主程序
main "$@"
