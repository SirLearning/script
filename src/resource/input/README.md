# 输入文件夹示例文件

## 说明

这个文件夹是输入文件的存放位置。请将以下文件放入此文件夹：

### 必需文件

1. **TIGER_20250526.jar** (或其他版本的TIGER jar文件)
   - FastCall2的主程序
   - 文件大小通常在几MB到几十MB

2. **chr1_10M.fa** (或其他参考基因组文件)
   - FASTA格式的参考基因组
   - 可以是全基因组或部分染色体
   - 确保同时有对应的.fai索引文件

3. **taxaBamMap.txt** (BAM文件映射表)
   - 制表符分隔的文本文件
   - 第一列：样本名称
   - 第二列：BAM文件的完整路径

### 示例文件结构

```
input/
├── TIGER_20250526.jar
├── chr1_10M.fa
├── chr1_10M.fa.fai      # 参考基因组索引文件
├── taxaBamMap.txt
├── README.md            # 本文件
└── results/             # 输出目录（自动创建）
    ├── results/         # 基础性能分析输出
    └── results_callpath/ # 调用路径分析输出
```

### taxaBamMap.txt 格式示例

```
sample01	/path/to/bam/sample01.sorted.bam
sample02	/path/to/bam/sample02.sorted.bam
sample03	/path/to/bam/sample03.sorted.bam
sample04	/path/to/bam/sample04.sorted.bam
```

### 使用方法

当所有文件就位后，您可以运行：

```bash
# 基础性能分析
nextflow run ../perf_basic_analysis.nf --input_dir input

# 调用路径分析
nextflow run ../perf_callpath_analysis.nf --input_dir input
```

或者使用相对路径：

```bash
# 从项目根目录运行
nextflow run perf_basic_analysis.nf --input_dir input
nextflow run perf_callpath_analysis.nf --input_dir input
```
