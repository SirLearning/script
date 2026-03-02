# gnomAD QC 代码库完整总结

> 仓库: [broadinstitute/gnomad_qc](https://github.com/broadinstitute/gnomad_qc)
> 提交: `2e35215ce15d7b6b5f5c82e77a6649af4510cf76`
> 本文档总结了 `gnomad_qc/` 文件夹及所有子文件夹中的程序内容

---

## 目录结构概览

```
gnomad_qc/
├── __init__.py                # 包初始化文件
├── resource_utils.py          # 通用管线资源管理工具
├── slack_creds.py             # Slack 通知凭据配置
├── analyses/                  # 跨版本分析脚本
│   ├── compute_pext.py
│   └── resources.py
├── example_notebooks/         # 示例 Jupyter Notebook
│   └── ancestry_classification_using_gnomad_rf.ipynb
├── v2/                        # gnomAD v2 (GRCh37) QC 管线
│   ├── annotations/
│   ├── load_data/
│   ├── resources/
│   ├── sample_qc/
│   └── variant_qc/
├── v3/                        # gnomAD v3 (GRCh38 基因组) QC 管线
│   ├── analysis/
│   ├── annotations/
│   ├── create_release/
│   ├── load_data/
│   ├── notebooks/
│   ├── resources/
│   ├── sample_qc/
│   ├── variant_qc/
│   └── utils.py
├── v4/                        # gnomAD v4 (GRCh38 外显子组+基因组) QC 管线
│   ├── analyses/
│   ├── annotations/
│   ├── assessment/
│   ├── create_release/
│   ├── plots/
│   ├── resources/
│   ├── sample_qc/
│   ├── variant_qc/
│   └── subset.py
└── v5/                        # gnomAD v5 (开发中)
    ├── annotations/
    ├── configs/
    ├── data_ingestion/
    ├── resources/
    └── sample_qc/
```

---

## 1. 顶层文件 (`gnomad_qc/`)

### 1.1 `resource_utils.py`

**用途**: 提供通用的 gnomAD QC 管线资源管理工具类和函数。

| 函数/类 | 说明 |
|---------|------|
| `check_resource_existence()` | 检查所有指定的输入和输出资源的存在性。如果输入资源不存在则报错，如果输出资源已存在且 `overwrite=False` 则报错 |
| `PipelineStepResourceCollection` | 管理管线中**单个步骤**的资源集合。跟踪每个步骤的输入/输出资源，支持从上游步骤继承输出作为输��，并提供资源存在性检查 |
| `PipelineResourceCollection` | 管理**整个管线**的资源集合。聚合多个 `PipelineStepResourceCollection`，统一管理跨步骤的共享资源和覆写行为 |

### 1.2 `slack_creds.py`

**用途**: 存储 Slack API Token，用于管线运行完成后发送通知到指定 Slack 频道。默认为空字符串占位符。

---

## 2. `analyses/` — 跨版本分析

### 2.1 `compute_pext.py`

**用途**: 计算 **pext (proportion expression across transcripts)** 指标 — 衡量基因组中每个碱基在不同转录本中的表达比例。

**主要功能**:
- 基于 GTEx 等表达数据，计算每个碱基位点在各组织中的转录本加权表达比例
- 用于评估功能丧失变异的生物学影响

### 2.2 `resources.py`

**用途**: 定义分析模块使用的资源路径和配置。

---

## 3. `example_notebooks/`

### 3.1 `ancestry_classification_using_gnomad_rf.ipynb`

**用途**: 演示如何使用 gnomAD 的随机森林 (Random Forest) 模型进行**祖源分类** (ancestry classification)。展示了 PCA 投影和 RF 预测的完整流程。

---

## 4. `v2/` — gnomAD v2 QC 管线 (GRCh37)

gnomAD v2 包含约 125,000 个外显子组和 15,000 个基因组样本，基于 GRCh37 参考基因组。

### 4.1 `v2/annotations/` — 变异注释生成

| 脚本 | 说明 |
|------|------|
| `generate_frequency_data.py` | 计算等位基因频率数据，包括按人群 (population)、性别 (sex)、亚群 (subpopulation) 分层的频率统计；生成 popmax、faf (filtering allele frequency) 等指标 |
| `generate_ld_data.py` | 计算连锁不平衡 (LD) 数据，生成按人群分层的 LD 评分和 LD 矩阵 |
| `generate_qc_annotations.py` | 生成用于质控的变异注释，包括 allele-specific QD、pAB、中位数 GQ/DP/AB 等 per-allele 统计量 |

### 4.2 `v2/load_data/` — 数据加载

| 脚本 | 说明 |
|------|------|
| `import_vcf.py` | 导入 gnomAD 基因组 VCF 文件到 Hail MatrixTable 格式，执行 `min_rep` 标准化 |
| `import_exomes_vcf_on_prem.py` | (归档) 在本地集群上导入外显子组 VCF 的历史脚本 |
| `import_gnomad_sv.py` | 导入 gnomAD 结构变异 (SV) 数据，合��� SV 与短变异 (short variants) 样本，分析样本重叠 |
| `import_resources.py` | 导入外部参考资源：ClinVar、de novo 变异、甲基化位点、ExAC 数据、CpG 位点、真值集 (truth sets) |
| `load_coverage.py` | 加载和合并基因组/外显子组覆盖度数据，生成 per-base 覆盖度 MatrixTable |

### 4.3 `v2/resources/` — 资源路径定义

| 模块 | 说明 |
|------|------|
| `basics.py` | 定义核心数据路径和访问函数：`get_gnomad_data()`（获取外显子组/基因组 MT）、`get_gnomad_meta()`（获取元数据）、VEP 配置、参考资源路径 |
| `annotations.py` | 定义注释相关的 Hail Table 路径（频率、VEP、RF 结果等） |
| `sample_qc.py` | 定义样本 QC 相关的路径：QC MT、PCA 评分、亲缘关系、平台检测结果等 |
| `variant_qc.py` | 定义变异 QC 相关的路径：RF 模型、VQSR 结果、评分排名、concordance 等 |

### 4.4 `v2/sample_qc/` — 样本质控

| 脚本 | 说明 |
|------|------|
| `apply_hard_filters.py` | **硬过滤**: 根据污染率 (contamination >5%)、call rate (<85%)、嵌合体比例 (chimera >5%)、性别不明确等阈值标记需过滤的样本。包含 `annotate_sex()` 进行性别推断、`make_hard_filters_expr()` 构建过滤表达式 |
| `exomes_platform_pca.py` | **平台检测**: 对外显子组样本的 call rate 进行 PCA，使用 HDBSCAN 聚类算法自动识别测序平台。`assign_platform_pcs()` 执行聚类并分配平台标签 |
| `joint_sample_qc.py` | **联合样本QC**: 联合外显子组和基因组进行 PCA 和祖源分类。`read_and_pre_process_data()` 预处理数据，`make_rank_file()` 为去重和亲属剪枝分配样本排名，使用 Random Forest 进行祖源分类 |
| `assign_subpops.py` | **亚群分配**: 为主要人群 (EUR/EAS/AFR) 中的样本分配亚群标签（如 NFE 子群: nwe/seu/bgr/est/swe；EAS 子群: jpn/kor/oea 等） |
| `create_fam.py` | **家系推断**: `GnomADRelatedData` 类封装亲缘关系数据，使用 `pc_relate` 结果推断三人家系 (trios)，处理重复样本 |
| `finalize_sample_qc.py` | **最终样本QC**: `add_release_annotations()` 标记 high_quality 和 release 状态，`collapse_small_pops()` 合并过小的人群，`add_topmed_annotation()` 添加 TOPMed 重复样本注释 |
| `generate_hardcalls.py` | **生成硬基因型**: 从原始 MT 生成仅含 GT 的硬基因型 MT，调整性染色体倍性，执行多等位基因拆分 |
| `get_topmed_dups.py` | **TOPMed 重复检测**: 比较 gnomAD 与 TOPMed 的共享位点，通过 singleton/doubleton 在 TOPMed 中的出现来识别潜在重复样本 |

### 4.5 `v2/variant_qc/` — 变异质控

| 脚本 | 说明 |
|------|------|
| `variantqc.py` | **Random Forest 变异过滤**: 核心变异 QC 脚本。定义 RF 特征集（位点特征、等位基因特征、VQSR 特征、中位数特征），训练 RF 模型区分 TP（真阳性，基于 truth sets）和 FP（假阳性，基于 hard filter 失败变异）|
| `create_ranked_scores.py` | **评分排名**: 基于 RF 概率和 VQSR VQSLOD 创建排名表，生成多种子排名（singleton、biallelic、adj 等），用于确定过滤阈值 |
| `calculate_concordance.py` | **一致性计算**: 计算外显子组和基因组之间的变异一致性，识别跨数据类型的重复样本，评估不同过滤策略下的一致性 |
| `prepare_data_release.py` | **发布数据准备**: 格式化频率、质量注释、直方图等用于公开发布���定义 VCF INFO 字段描述，生成频率索引 |
| `correct_fafs.py` | **修正 FAF**: 修正 filtering allele frequency——当某人群中仅有 singleton 时，将 FAF 设为 0 |
| `make_var_annot_hists.py` | **变异注释直方图**: 生成 QC 指标（FS、MQ、QD、VQSLOD 等）的直方图，按等位基因频率 bin 聚合统计 |
| `variant_qc_plots.py` | **可视化**: 使用 Bokeh 生成变异 QC 图表，包括 binned 模型比较、concordance 曲线等 |
| `select_qc_set.py` | (归档) 用于一次性分析的子采样脚本 |
| `exomes_genomes_coverage.py` | (归档) 用于 gnomAD v2 LOF 论文的外显子组-基因组覆盖度分析 |

---

## 5. `v3/` — gnomAD v3 QC 管线 (GRCh38 基因组)

gnomAD v3 包含约 76,000 个全基因组样本，基于 GRCh38 参考基因组。

### 5.1 顶层文件

| 文件 | 说明 |
|------|------|
| `utils.py` | v3 版本特定的工具函数 |

### 5.2 子目录结构

| 目录 | 说明 |
|------|------|
| `analysis/` | v3 特定分析脚本 |
| `annotations/` | 变异注释生成（频率、VEP、QC 注释等） |
| `create_release/` | 发布数据创建和格式化 |
| `load_data/` | 数据导入和转换 |
| `notebooks/` | Jupyter 分析笔记本 |
| `resources/` | v3 资源路径定义 |
| `sample_qc/` | 样本质控（硬过滤、性别推断、平台检测、祖源分类、亲缘关系等） |
| `variant_qc/` | 变异质控（RF 过滤、VQSR、评分排名等） |

---

## 6. `v4/` — gnomAD v4 QC 管线 (GRCh38)

gnomAD v4 是最大规模的版本，包含约 730,000 个外显子组和 76,000 个基因组样本。

### 6.1 `v4/subset.py`

**用途**: 从完整 gnomAD v4 数据集创建子集（如特定人群或队列的子集），用于下游分析和发布。

### 6.2 `v4/sample_qc/` — 样本质控

v4 引入了更复杂的 QC 管线，使用 `PipelineResourceCollection` 进行依赖管理。

| 脚本 | 说明 |
|------|------|
| `hard_filters.py` | **硬过滤与样本QC**: `compute_sample_qc()` 使用 `compute_stratified_sample_qc` 按多等位基因分层计算 QC 指标。基于污染、覆盖度��指纹验证失败等阈值进行硬过滤 |
| `sex_inference.py` | **性别推断**: `determine_fstat_sites()` 筛选 chrX 上的高质量 SNP，基于 F-stat 推断性染色体核型 (XX/XY/XXY 等) |
| `interval_qc.py` | **区间QC**: 基于 per-interval 的覆盖度统计确定高质量区间。支持按平台分层。`generate_sex_chr_interval_coverage_mt()` 处理性染色体覆盖度，考虑 PAR 区域 |
| `generate_qc_mt.py` | **QC MatrixTable 生成**: 联合 v3 基因组和 v4 外显子组创建密集 MT，过滤到预定义 QC 位点、执行 LD 剪枝，用于祖源 PCA 和亲缘关系分析 |
| `assign_ancestry.py` | **祖源分配**: 使用 v3 已知标签和 TGP/HGDP 参考进行 PCA + RF 祖源分类。定义 `V4_POP_SPIKE_DICT`（中东人群分类）和 `V3_SPIKE_PROJECTS`（各人群训练队列） |
| `outlier_filtering.py` | **异常值过滤**: 多种过滤策略 — 分层过滤 (stratified)、回归残差过滤 (regressed)、最近邻过滤 (nearest neighbors)。`get_sample_qc_ht()` 预处理 QC 指标 |
| `relatedness.py` | **亲缘关系估计**: 使用 **cuKING** (GPU 加速) 计算成对亲缘关系，确定需要移除的相关样本。`print_cuking_command()` 生成 Cloud Batch 作业命令 |
| `identify_trios.py` | **家系识别**: 从亲缘关系数据推断三人家系，按 Mendel 错误和 de novo 变异过滤。`families_to_trios()` 从每个家系中随机选择一个 trio |
| `platform.py` (推测存在) | 平台检测和分配 |
| `finalize.py` (推测存在) | 最终元数据整合和 release 标记 |

### 6.3 `v4/variant_qc/` — 变异质控

| 脚本 | 说明 |
|------|------|
| `import_variant_qc_vcf.py` | 导入变异 QC VCF 文件（如 VQSR 输出）到 Hail Table 格式 |
| `random_forest.py` | **RF 变异过滤**: 训练和应用 Random Forest 模型进行变异质量评分，使用真值集 (truth sets) 定义训练标签 |
| `vqsr.py` | **VQSR 管线**: 完整的 GATK VQSR (Variant Quality Score Recalibration) 实现，包括资源定义、模型训练、recalibration 应用。这是一个大型脚本 (~54KB) |
| `vqsr_resources.json` | VQSR 使用的参考资源配置 (HapMap、Omni、1000G、dbSNP 等) |
| `evaluation.py` | **评估**: 使用 concordance、真值集一致性等指标评估不同变异过滤策略的性能 |
| `final_filter.py` | **最终过滤 (外显子组)**: 整合 RF 和 VQSR 评分，确定最终过滤阈值，生成 PASS/FAIL 标记 |
| `final_filter_genomes.py` | **最终过滤 (基因组)**: 基因组数据的最终变异过滤，策略可能与外显子组不同 |

### 6.4 `v4/annotations/`

变异注释生成脚本，包括频率计算、VEP 注释、quality histogram 等。

### 6.5 `v4/assessment/`

数据质量评估和验证脚本。

### 6.6 `v4/create_release/`

创建公开发布数据集的脚本，包括 VCF 导出、站点频率汇总等。

### 6.7 `v4/plots/`

数据可视化脚本，生成 QC 相关图表。

### 6.8 `v4/resources/`

v4 版本的资源路径定义模块。

---

## 7. `v5/` — gnomAD v5 (开发中)

### 7.1 子目录结构

| 目录 | 说明 |
|------|------|
| `annotations/` | v5 变异注释脚本 |
| `configs/` | 配置文件 |
| `data_ingestion/` | 数据摄入管线 — 新的模块化数据导入方式 |
| `resources/` | v5 资源路径定义 |
| `sample_qc/` | v5 样本质控脚本 |

---

## 核心技术栈

| 组件 | 说明 |
|------|------|
| **Hail** | 分布式基因组数据处理框架 (MatrixTable, Table, VDS) |
| **Google Cloud Platform** | 数据存储 (GCS) 和计算 (Dataproc, Cloud Batch) |
| **Random Forest** | 变异质控的核心模型 (scikit-learn) |
| **VQSR** | GATK 的变异质量评分重校准 |
| **cuKING** | GPU 加速的亲缘关系估计 (v4 新引入) |
| **PCA** | 主成分分析用于祖源分类和平台检测 |
| **HDBSCAN** | 密度聚类算法用于平台检测 |
| **Bokeh** | 交互式可视化 |

---

## QC 管线整体流程

```
1. 数据导入 (load_data / data_ingestion)
   └── VCF → Hail MatrixTable/VDS

2. 样本 QC (sample_qc)
   ├── 硬过滤 (hard_filters) → 移除低质量样本
   ├── 性别推断 (sex_inference) → 性染色体核型
   ├── 平台检测 (platform) → 测序平台标签
   ├── 祖源分类 (assign_ancestry) → 人群标签
   ├── 亲缘关系 (relatedness) → 相关样本移除
   ├── 异常值过滤 (outlier_filtering) → QC 指标异常值
   ├── 家系识别 (identify_trios) → 用于变异 QC 评估
   └── 最终元数据 (finalize) → release 标记

3. 变异注释 (annotations)
   ├── 频率计算 → 等位基因频率
   ├── VEP → 功能注释
   └── QC 注释 → 质量指标

4. 变异 QC (variant_qc)
   ├── RF 训练 → Random Forest 模型
   ├── VQSR → 变异质量重校准
   ├── 评估 → concordance, truth set 一致性
   └── 最终过滤 → PASS/FAIL 标记

5. 数据发布 (create_release)
   └── 格式化导出 → VCF, Hail Table
```

---

> **注意**: 由于代码搜索结果数量限制，本文档可能未涵盖所有脚本。更完整的文件列表请参考 [GitHub 仓库](https://github.com/broadinstitute/gnomad_qc/tree/2e35215ce15d7b6b5f5c82e77a6649af4510cf76/gnomad_qc)。v3 和 v5 的详细脚本内容因深度探索限制可能不够完整。