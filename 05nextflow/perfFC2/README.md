# FastCall2 性能分析 Nextflow 流程

这个 Nextflow 工作流程专门用于分析 FastCall2 程序在不同样本大小和运行模式下的性能表现。

## 功能特性

- **多模式分析**: 支持 FastCall2 的三个运行模式（disc/blib/scan）
- **可扩展性测试**: 自动测试不同样本大小下的性能表现
- **全面的性能指标**: 收集 CPU 周期、指令数、缓存命中率、分支预测等指标
- **火焰图生成**: 可选生成详细的性能火焰图
- **自动化报告**: 生成性能分析总结报告

## 目录结构

```
项目目录/
├── fastcall.nf           # 主要工作流文件
├── nextflow.config       # 配置文件
├── data/
│   └── input/
│       ├── chr1_10M.fa   # 参考基因组
│       └── taxaBamMap.txt # BAM文件映射表
├── TIGER_20250526.jar    # FastCall2 程序
├── bin/
│   └── samtools          # samtools工具
└── results/              # 输出目录
    ├── bam_maps/         # 不同样本大小的BAM映射文件
    ├── perf_logs/        # 性能分析日志
    ├── perf_record/      # perf record数据
    ├── flamegraphs/      # 火焰图
    ├── analysis/         # 性能分析结果
    └── reports/          # 执行报告
```

## 使用方法

### 1. 基本使用

```bash
# 查看帮助信息
nextflow run fastcall.nf --help

# 使用默认参数运行
nextflow run fastcall.nf

# 生成执行报告
nextflow run fastcall.nf -with-report -with-timeline -with-dag
```

### 2. 自定义参数

```bash
# 指定自定义路径和参数
nextflow run fastcall.nf \
  --jar /path/to/TIGER.jar \
  --reference /path/to/reference.fa \
  --bam_maps /path/to/taxaBamMap.txt \
  --outdir /path/to/results \
  --threads 64 \
  --memory "200g" \
  --sample_sizes "1,2,4,8,16"

# 不生成火焰图（加快执行速度）
nextflow run fastcall.nf --generate_flamegraph false

# 使用不同的性能事件
nextflow run fastcall.nf \
  --perf_events "cycles,instructions,cache-references,cache-misses,L1-dcache-loads,L1-dcache-load-misses"
```

### 3. 使用配置文件

```bash
# 使用特定profile
nextflow run fastcall.nf -profile cluster
nextflow run fastcall.nf -profile debug

# 恢复中断的运行
nextflow run fastcall.nf -resume
```

## 输出文件说明

### 性能日志文件
- `perf.{N}.disc.log`: N个样本的disc模式性能数据
- `perf.{N}.blib.log`: N个样本的blib模式性能数据  
- `perf.{N}.scan.log`: N个样本的scan模式性能数据

### 分析结果文件
- `ipc_analysis.txt`: IPC（每周期指令数）分析
- `branch_analysis.txt`: 分支预测分析
- `timing_analysis.txt`: 执行时间分析
- `performance_summary.md`: 性能分析总结报告

### 火焰图文件
- `flame.{N}.svg`: N个样本的性能火焰图

## 参数说明

### 必要参数
- `--jar`: TIGER jar文件路径
- `--reference`: 参考基因组FASTA文件路径
- `--bam_maps`: BAM文件映射表路径
- `--samtools`: samtools可执行文件路径

### 可选参数
- `--outdir`: 输出目录（默认: `$projectDir/results`）
- `--threads`: CPU线程数（默认: 32）
- `--memory`: 内存大小（默认: "100g"）
- `--sample_sizes`: 测试的样本大小列表（默认: "1,2,3,4,5,6"）
- `--perf_events`: perf监控的事件（默认: "cycles,instructions,cache-misses,branches,branch-misses"）
- `--generate_flamegraph`: 是否生成火焰图（默认: true）

### FastCall2参数
- `--min_depth`: 最小深度（默认: 30）
- `--min_qual`: 最小质量（默认: 20）
- `--ploidy`: 倍性（默认: 2）
- `--min_allele_freq`: 最小等位基因频率（默认: 0.2）

## 前置要求

### 软件依赖
1. **Nextflow** (>= 21.04.0)
2. **perf** 工具（Linux性能分析工具）
3. **Java** (>= 8) 用于运行FastCall2
4. **FlameGraph** 工具（可选，用于生成火焰图）

### 安装 FlameGraph（可选）
```bash
git clone https://github.com/brendangregg/FlameGraph.git
export PATH=$PATH:/path/to/FlameGraph
```

### 系统配置
```bash
# 配置perf权限（非root用户）
echo "kernel.perf_event_paranoid = -1" | sudo tee /etc/sysctl.d/perf.conf
sudo sysctl -p /etc/sysctl.d/perf.conf
```

## 性能指标解释

### IPC (Instructions Per Cycle)
- **含义**: 每个CPU周期执行的指令数
- **理想值**: > 1.0，通常1.0-2.0为良好表现
- **优化**: 如果过低，考虑算法优化或减少分支

### Cache Miss Rate
- **含义**: 缓存未命中率
- **理想值**: 越低越好，通常 < 5%
- **优化**: 优化数据结构和访问模式

### Branch Miss Rate
- **含义**: 分支预测失败率
- **理想值**: < 1%
- **优化**: 减少条件分支，使用更可预测的代码结构

## 故障排除

### 常见问题

1. **perf权限问题**
   ```bash
   # 解决方案：配置perf权限
   sudo sysctl kernel.perf_event_paranoid=-1
   ```

2. **内存不足**
   ```bash
   # 解决方案：减少内存使用或增加系统内存
   nextflow run fastcall.nf --memory "50g"
   ```

3. **FlameGraph工具未找到**
   ```bash
   # 解决方案：安装FlameGraph或禁用火焰图生成
   nextflow run fastcall.nf --generate_flamegraph false
   ```

## 示例运行

```bash
# 完整的性能分析运行示例
nextflow run fastcall.nf \
  --jar ./TIGER_20250526.jar \
  --reference ./data/input/chr1_10M.fa \
  --bam_maps ./data/input/taxaBamMap.txt \
  --samtools ./bin/samtools \
  --outdir ./results \
  --threads 32 \
  --memory "100g" \
  --sample_sizes "1,2,4,8" \
  --generate_flamegraph true \
  -with-report results/reports/report.html \
  -with-timeline results/reports/timeline.html \
  -profile standard
```

这个流程将帮助您全面分析 FastCall2 的性能特征，识别性能瓶颈，并为优化提供数据支持。
