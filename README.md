# Vmap4 Data Processing Pipeline

基于 `00seq_data_process.ipynb` 的方法，这是一个完整的小麦基因组数据处理 Nextflow 工作流系统。

## 项目结构

```
01vmap4/
├── main.nf                          # 主入口文件
├── nextflow.config                  # 主配置文件
├── run_vmap4_pipeline.sh           # 运行脚本
├── DataProcess/
│   ├── alignment/
│   │   ├── sequence_alignment.nf   # 序列比对工作流
│   │   └── nextflow.config         # 比对工作流配置
│   ├── calling/
│   │   ├── perf/
│   │   │   ├── fastcall.nf         # 性能分析工作流
│   │   │   └── nextflow.config     # 性能分析配置
│   │   └── run/
│   │       ├── runFastCall2.nf     # FastCall2 变异检测
│   │       └── nextflow.config     # 变异检测配置
│   └── overall/
│       ├── data_processing.nf      # 数据下载处理工作流
│       ├── depth_qc.nf            # 深度计算和质控工作流
│       └── nextflow.config        # 数据处理配置
├── output/                         # 输出目录
└── seq/                           # 序列处理模块
```

## 工作流类型

### 1. 数据处理工作流 (`data_processing`)
基于 `00seq_data_process.ipynb` 中的数据下载和预处理方法：
- SRA 数据下载（使用 `prefetch` 和 `fasterq-dump`）
- FASTQ 文件提取和压缩
- 质量控制（FastQC）
- 子采样（可选）
- MD5 校验和生成

### 2. 序列比对工作流 (`alignment`)
基于 notebook 中的比对方法：
- 支持 BWA 和 BWA-MEM2
- 自动化的 BAM 处理流程（排序、去重、索引）
- 比对统计报告生成

### 3. 深度计算和质控工作流 (`depth_qc`)
基于深度计算部分的方法：
- BAM 文件索引
- 使用 mosdepth 计算深度
- samtools flagstat 质量控制
- 深度统计汇总

### 4. FastCall2 变异检测工作流 (`fastcall2`)
基于 FastCall2 运行方法：
- disc、blib、scan 三步骤流程
- 支持多染色体并行处理
- 自动化参数优化

### 5. 性能分析工作流 (`performance`)
新增的性能监控功能：
- CPU、内存、I/O 使用监控
- 执行时间分析
- 性能瓶颈识别
- 详细的性能报告

## 使用方法

### 快速开始

1. **使用运行脚本（推荐）**
```bash
# 序列比对工作流
./run_vmap4_pipeline.sh -w alignment \
    --fastq-dir ./fastq \
    --reference genome.fa \
    --bam-list samples.txt

# 深度计算工作流
./run_vmap4_pipeline.sh -w depth_qc \
    --bam-dir ./bam \
    --bam-list bam_files.txt

# FastCall2 变异检测
./run_vmap4_pipeline.sh -w fastcall2 \
    --reference genome.fa \
    --taxa-bam-map taxa.txt \
    --tiger-jar TIGER.jar \
    --samtools-path samtools
```

2. **直接使用 Nextflow**
```bash
# 序列比对
nextflow run DataProcess/alignment/sequence_alignment.nf \
    --fastq_dir ./fastq \
    --reference genome.fa \
    --sample_list samples.txt \
    --fq_list fqlist.txt

# 深度计算
nextflow run DataProcess/overall/depth_qc.nf \
    --bam_dir ./bam \
    --bam_list bam_files.txt
```

### 配置文件

每个工作流都有对应的配置文件，可以根据需要修改：

- **资源配置**：CPU 核数、内存大小、执行时间
- **软件路径**：各种生物信息学工具的路径
- **参数优化**：质量阈值、深度阈值等
- **执行环境**：本地、集群、容器等

### 输入文件格式

1. **SRA 列表文件** (`sra_accessions.txt`)
```
SRR7164604
SRR7164620
SRR7164628
```

2. **样本列表文件** (`samples.txt`)
```
sample1
sample2
sample3
```

3. **FASTQ 列表文件** (`fqlist.txt`)
```
SRR7164604
SRR7164620
SRR7164628
```

4. **BAM 列表文件** (`bam_files.txt`)
```
sample1.rmdup.bam
sample2.rmdup.bam
sample3.rmdup.bam
```

5. **Taxa-BAM 映射文件** (`taxa.txt`)
```
sample1	/path/to/sample1.bam	5.2
sample2	/path/to/sample2.bam	4.8
sample3	/path/to/sample3.bam	6.1
```

## 输出结果

每个工作流会生成以下类型的输出：

### 通用输出
- `reports/`：执行报告、时间线、流程图
- `logs/`：详细的执行日志
- `work/`：临时工作目录

### 工作流特定输出

**数据处理工作流**：
- `fastq/`：提取的 FASTQ 文件
- `quality_control/`：FastQC 报告
- `checksums/`：MD5 校验文件

**序列比对工作流**：
- `bam/`：比对后的 BAM 文件
- `stats/`：比对统计报告

**深度计算工作流**：
- `depth/`：深度计算结果
- `qc/`：质量控制报告
- `plots/`：覆盖度图表（可选）

**FastCall2 工作流**：
- `disc/`：disc 步骤输出
- `blib/`：blib 步骤输出
- `scan/`：最终 VCF 文件
- `final/`：合并的变异数据

**性能分析工作流**：
- `performance/`：性能监控数据
- `reports/`：性能分析报告

## 性能优化

### 资源配置建议

1. **小型数据集**（< 10 个样本）
   - CPU: 8-16 核
   - 内存: 16-32 GB
   - 存储: SSD 推荐

2. **中型数据集**（10-100 个样本）
   - CPU: 32-64 核
   - 内存: 64-128 GB
   - 存储: 高速 SSD

3. **大型数据集**（> 100 个样本）
   - 建议使用集群环境
   - 每个节点: 32+ 核, 128+ GB 内存
   - 分布式存储系统

### 集群配置

支持多种集群管理器：
- SLURM
- PBS/Torque
- SGE
- AWS Batch

使用 `-profile cluster` 启用集群配置。

## 故障排除

### 常见问题

1. **内存不足错误**
   - 增加内存分配：`--memory 32G`
   - 使用 `high_performance` profile

2. **磁盘空间不足**
   - 清理临时文件：`nextflow clean -f`
   - 使用更大的工作目录

3. **网络连接问题**（SRA 下载）
   - 设置重试次数：配置文件中的 `maxRetries`
   - 使用本地数据源

4. **软件路径错误**
   - 检查配置文件中的软件路径
   - 使用绝对路径

### 调试模式

使用调试配置运行：
```bash
nextflow run workflow.nf -profile debug
```

## 依赖软件

确保以下软件已安装并在 PATH 中：

- **必需软件**：
  - Nextflow (>= 21.04.0)
  - samtools (>= 1.18)
  - bwa 或 bwa-mem2

- **数据处理工作流**：
  - sra-tools (prefetch, fasterq-dump)
  - fastqc
  - pigz
  - seqtk

- **深度计算工作流**：
  - mosdepth

- **FastCall2 工作流**：
  - TIGER jar 文件
  - bcftools（可选，用于 VCF 合并）

- **性能分析**：
  - R (可选，用于图表生成)
  - ggplot2, dplyr 等 R 包

## 版本信息

- 版本: 1.0.0
- 基于: 00seq_data_process.ipynb 方法
- Nextflow 最低版本: 21.04.0
- 兼容性: Linux, macOS

## 许可证

本项目遵循开源许可证，具体请查看 LICENSE 文件。

## 贡献

欢迎提交 Issue 和 Pull Request 来改进这个工作流系统。

## 支持

如有问题，请联系开发团队或查看 GitHub Issues。


# old script

There are 4 main script dictions:

- **00basic**: basic genomic analysis scripts
- **01vmap4**: jobs I did for vmap4 project
- **02transposon**: jobs I did in transposon field
- **03JM44**: reference jobs in JM44 project

You can also find other dictions like **practice**, which contains scripts I did when I was learning python or R.