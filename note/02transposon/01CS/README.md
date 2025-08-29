# 转座子分析流程 (统一版本)

## 概述
这是一个用于分析转座子数据的统一流程，将原来的多个Python脚本合并为一个功能强大的分析工具。包含以下功能：
- 转座子组成分析
- GFF3文件类别统计
- 库文件统计分析
- 数据质量控制
- 生成综合分析报告（HTML和文本格式）

## 文件结构
```
02transposon/01CS/
├── transposon_analyzer.py            # 统一的Python分析脚本 ⭐
├── transposon_analysis_simple.nf     # Nextflow流程文件（可选）
├── nextflow.config                   # 配置文件
├── run_analysis.bat                  # Windows运行脚本
├── run_analysis.sh                   # Linux/Mac运行脚本
├── test_environment.bat              # 环境测试脚本
├── newplot.txt                       # 输入数据文件
└── README.md                         # 本文件
```

## 新版本特性

### ✨ 主要改进
- **统一脚本**: 将原来的3个Python脚本(`01compo_count_v1.py`, `01compo_count.py`, `02lib_stats.py`)合并为一个功能完整的`transposon_analyzer.py`
- **模块化分析**: 支持单独运行不同类型的分析或运行完整分析
- **更好的错误处理**: 完善的异常处理和日志记录
- **丰富的输出**: 生成JSON、HTML、文本多种格式的报告
- **参数化**: 完全支持命令行参数配置

### 🔄 兼容性
- 保持与原有Nextflow流程的兼容性
- 可以选择直接使用Python脚本或通过Nextflow运行
- 支持可选的输入文件（GFF和库文件）

## 安装要求

### 必需组件
- **Python 3.7+** 
  ```bash
  python --version  # 检查版本
  ```

### 可选组件
- **Nextflow** (如果要使用流程管理功能)
  ```bash
  # 通过Java安装
  curl -s https://get.nextflow.io | bash
  
  # 或通过conda安装 (Linux/Mac)
  conda install -c bioconda nextflow
  ```

## 使用方法

### 🚀 方式1: 直接使用Python脚本 (推荐)

#### 基本用法
```bash
# 运行完整分析
python transposon_analyzer.py --mode all --input newplot.txt --output results

# 仅分析组成
python transposon_analyzer.py --mode composition --input newplot.txt --output results

# 包含GFF和库文件的完整分析
python transposon_analyzer.py --mode all --input newplot.txt --gff chr1A.gff3 --lib library.lib --output results
```

#### 参数说明
| 参数 | 必需 | 描述 | 默认值 |
|------|------|------|--------|
| `--mode` | 否 | 分析模式: composition, gff, library, qc, all | all |
| `--input` | 是 | 主输入数据文件 | - |
| `--gff` | 否 | GFF3注释文件 | - |
| `--lib` | 否 | 库文件 | - |
| `--output` | 否 | 输出目录 | results |

#### 分析模式详解
- **composition**: 转座子组成分析（基于newplot.txt格式）
- **gff**: GFF3文件类别统计
- **library**: 库文件统计分析
- **qc**: 数据质量控制
- **all**: 运行所有可用的分析

### 🚀 方式2: 使用批处理脚本

#### Windows用户
```cmd
# 环境测试
test_environment.bat

# 运行分析（交互式）
run_analysis.bat

# 带参数运行
run_analysis.bat --input mydata.txt --output my_results --mode all
```

#### Linux/Mac用户
```bash
# 给脚本执行权限
chmod +x run_analysis.sh

# 运行分析
./run_analysis.sh --input newplot.txt --output results
```

### 🚀 方式3: 使用Nextflow流程

```bash
# 基本运行
nextflow run transposon_analysis_simple.nf

# 自定义参数
nextflow run transposon_analysis_simple.nf \
    --input_data "mydata.txt" \
    --gff_file "annotation.gff3" \
    --output_dir "my_results"

# 使用不同配置
nextflow run transposon_analysis_simple.nf -profile docker
```

## 输出结果

### 📁 输出目录结构
```
results/
├── composition_stats.txt              # 组成分析详细结果
├── composition_summary.txt            # 组成分析摘要
├── gff_category_counts.txt            # GFF类别计数（如果有GFF文件）
├── gff_analysis_report.txt            # GFF分析报告
├── lib_stats_*.txt                    # 库统计文件（如果有库文件）
├── qc_report.txt                      # 质量控制报告
├── data_validation.txt                # 数据验证结果
├── transposon_analysis_report.html    # 📊 HTML综合报告
├── transposon_analysis_report.txt     # 📄 文本综合报告
└── analysis_results.json              # 🔧 JSON格式结果（机器可读）
```

### 📊 报告内容
- **组成分析**: 转座子类别统计和百分比分析
- **GFF分析**: 基因组注释中的转座子分类
- **库统计**: 库文件中各染色体的统计信息
- **质量控制**: 数据格式验证和基本统计

## 故障排除

### 常见问题

1. **Python脚本语法错误**
   ```bash
   # 检查脚本语法
   python -m py_compile transposon_analyzer.py
   ```

2. **文件未找到**
   ```bash
   # 确保所有文件在正确位置
   ls -la transposon_analyzer.py newplot.txt
   ```

3. **权限错误 (Linux/Mac)**
   ```bash
   # 给脚本执行权限
   chmod +x *.sh *.py
   ```

4. **模块导入错误**
   ```bash
   # 确保Python版本正确
   python --version
   # 应该是 3.7 或更高版本
   ```

### 🔍 调试技巧

1. **查看详细错误信息**
   ```bash
   python transposon_analyzer.py --mode all --input newplot.txt --output results
   ```

2. **测试单个模式**
   ```bash
   # 先测试组成分析
   python transposon_analyzer.py --mode composition --input newplot.txt --output test_results
   ```

3. **检查输入文件格式**
   ```bash
   head -5 newplot.txt  # 查看前5行
   wc -l newplot.txt    # 检查行数
   ```

## 高级功能

### 🔧 自定义分析

如果需要修改分析逻辑，可以直接编辑`transposon_analyzer.py`文件中的相应方法：

- `analyze_composition()`: 组成分析逻辑
- `analyze_gff()`: GFF分析逻辑
- `analyze_library()`: 库分析逻辑
- `quality_control()`: 质量控制逻辑

### 📈 批量处理

```bash
# 处理多个文件
for file in data/*.txt; do
    python transposon_analyzer.py --mode all --input "$file" --output "results_$(basename $file .txt)"
done
```

### 🔄 集成到其他工具

由于输出包含JSON格式，可以轻松集成到其他分析流程中：

```python
import json
with open('results/analysis_results.json', 'r') as f:
    results = json.load(f)
    # 使用结果数据
```

## 版本历史

- **v2.0** (当前版本): 统一Python脚本，改进用户体验
- **v1.0**: 原始版本，多个独立Python脚本

## 贡献指南

1. 备份原始文件
2. 在测试数据上验证修改
3. 更新文档
4. 提交更改

## 联系信息

如有问题或建议，请通过以下方式联系：
- 创建GitHub Issue
- 邮箱咨询

---

**🎉 快速开始**: 运行 `test_environment.bat` 检查环境，然后执行 `run_analysis.bat` 开始分析！
