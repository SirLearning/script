# FastCall2 高通量 Nextflow 流程文件说明

## 文件结构

```
01vmap4/runFC2/
├── runFastCall2.nf      # 主要的Nextflow流程文件
├── nextflow.config      # 配置文件
├── README.md           # 详细文档
├── run_example.sh      # Linux/Unix运行示例
├── run_example.bat     # Windows运行示例
└── SUMMARY.md          # 本文件
```

## 核心特性

### 1. 高度并行化
- 按染色体并行处理，每个染色体独立运行disc、blib、scan三个步骤
- 支持21个小麦染色体同时处理
- 可根据系统资源调整并行度

### 2. 灵活配置
- 所有FastCall2参数都可通过命令行自定义
- 支持多种执行环境配置（本地、集群、测试）
- 内存和线程数可动态调整

### 3. 流程设计
```
输入文件 → disc (发现变异位点) → blib (构建变异库) → scan (基因型检测) → 结果合并
```

### 4. 错误处理
- 自动重试失败的任务
- 详细的日志记录
- 参数验证和错误提示

## 快速开始

### 准备工作
1. 安装Nextflow
2. 准备参考基因组文件
3. 使用Workshop.jar创建taxaBamMap.txt文件
4. 确保samtools可用

### 运行命令
```bash
nextflow run runFastCall2.nf \
  --reference /path/to/genome.fa \
  --taxaBamMap /path/to/taxaBamMap.txt \
  --tiger_jar /path/to/TIGER.jar \
  --samtools_path /path/to/samtools
```

## 输出结构
```
fastcall2_output/
├── disc/          # 每个染色体的.ing文件
├── blib/          # 每个染色体的.lib.gz文件  
├── scan/          # 每个染色体的VCF文件
├── final/         # 合并后的结果和统计
└── pipeline_info/ # 流程执行报告
```

## 高级用法

### 自定义参数
```bash
nextflow run runFastCall2.nf \
  --reference genome.fa \
  --taxaBamMap map.txt \
  --tiger_jar TIGER.jar \
  --samtools_path samtools \
  --disc_min_depth 20 \
  --scan_p_value 0.01 \
  --threads 48 \
  --memory "200g"
```

### 集群运行
```bash
nextflow run runFastCall2.nf \
  --reference genome.fa \
  --taxaBamMap map.txt \
  --tiger_jar TIGER.jar \
  --samtools_path samtools \
  -profile hpc
```

### 测试运行
```bash
nextflow run runFastCall2.nf \
  --reference genome.fa \
  --taxaBamMap map.txt \
  --tiger_jar TIGER.jar \
  --samtools_path samtools \
  -profile test
```

## 优势

1. **效率提升**: 相比串行处理，并行化可将处理时间减少到原来的1/21
2. **资源优化**: 可根据系统资源动态调整参数
3. **易于管理**: 统一的配置和监控
4. **可重复性**: 完整的执行记录和版本控制
5. **扩展性**: 容易适配不同的计算环境

## 注意事项

1. 确保有足够的磁盘空间存储中间文件
2. 根据样本数量调整内存分配
3. 验证BAM文件的索引完整性
4. 检查染色体命名一致性

## 技术细节

- 基于Nextflow DSL2语法
- 支持容器化部署
- 完整的错误处理和重试机制
- 自动生成执行报告和统计信息

这个流程已经在你的工作环境中配置完成，可以直接使用！
