# gnomAD v4 QC 管线 — 全部函数完整参考手册

> 仓库: [broadinstitute/gnomad_qc](https://github.com/broadinstitute/gnomad_qc)
> 提交: `2e35215ce15d7b6b5f5c82e77a6649af4510cf76`
> 路径: `gnomad_qc/v4/`

---

## 目录

- [1. 顶层文件](#1-顶层文件)
  - [1.1 subset.py](#11-subsetpy)
- [2. sample_qc/ — 样本质控](#2-sample_qc--样本质控)
  - [2.1 hard_filters.py](#21-hard_filterspy)
  - [2.2 sex_inference.py](#22-sex_inferencepy)
  - [2.3 interval_qc.py](#23-interval_qcpy)
  - [2.4 generate_qc_mt.py](#24-generate_qc_mtpy)
  - [2.5 assign_ancestry.py](#25-assign_ancestrypy)
  - [2.6 outlier_filtering.py](#26-outlier_filteringpy)
  - [2.7 relatedness.py](#27-relatednesspy)
  - [2.8 identify_trios.py](#28-identify_triospy)
  - [2.9 platform_inference.py](#29-platform_inferencepy)
  - [2.10 create_sample_qc_metadata_ht.py](#210-create_sample_qc_metadata_htpy)
  - [2.11 create_sample_qc_metadata_ht_genomes.py](#211-create_sample_qc_metadata_ht_genomespy)
- [3. variant_qc/ — 变异质控](#3-variant_qc--变异质控)
  - [3.1 import_variant_qc_vcf.py](#31-import_variant_qc_vcfpy)
  - [3.2 random_forest.py](#32-random_forestpy)
  - [3.3 vqsr.py](#33-vqsrpy)
  - [3.4 evaluation.py](#34-evaluationpy)
  - [3.5 final_filter.py](#35-final_filterpy)
  - [3.6 final_filter_genomes.py](#36-final_filter_genomespy)
- [4. annotations/ — 变异注释](#4-annotations--变异注释)
  - [4.1 generate_freq.py](#41-generate_freqpy)
  - [4.2 generate_freq_genomes.py](#42-generate_freq_genomespy)
  - [4.3 generate_variant_qc_annotations.py](#43-generate_variant_qc_annotationspy)
  - [4.4 compute_coverage.py](#44-compute_coveragepy)
  - [4.5 fix_freq_an.py](#45-fix_freq_anpy)
  - [4.6 insilico_predictors.py](#46-insilico_predictorspy)
  - [4.7 vep_context_ht.py](#47-vep_context_htpy)
  - [4.8 recover_and_complete_vep115.py](#48-recover_and_complete_vep115py)
  - [4.9 vrs_annotation_batch.py](#49-vrs_annotation_batchpy)
- [5. create_release/ — 发布数据创建](#5-create_release--发布数据创建)
  - [5.1 create_release_utils.py](#51-create_release_utilspy)
  - [5.2 create_release_sites_ht.py](#52-create_release_sites_htpy)
  - [5.3 create_combined_faf_release_ht.py](#53-create_combined_faf_release_htpy)
  - [5.4 create_de_novo_release.py](#54-create_de_novo_releasepy)
  - [5.5 create_false_dup_liftover.py](#55-create_false_dup_liftoverpy)
  - [5.6 validate_and_export_vcf.py](#56-validate_and_export_vcfpy)
  - [5.7 make_var_annot_hists.py](#57-make_var_annot_histspy)
- [6. assessment/ — 数据评估](#6-assessment--数据评估)
  - [6.1 calculate_per_sample_stats.py](#61-calculate_per_sample_statspy)
  - [6.2 summary_stats.py](#62-summary_statspy)
- [7. analyses/ — 特定分析](#7-analyses--特定分析)
  - [7.1 grpmax_comps.py](#71-grpmax_compspy)
- [8. resources/ — 资源路径定义](#8-resources--资源路径定义)
  - [8.1 basics.py](#81-basicspy)
  - [8.2 constants.py](#82-constantspy)
  - [8.3 annotations.py](#83-annotationspy)
  - [8.4 sample_qc.py](#84-sample_qcpy)
  - [8.5 variant_qc.py](#85-variant_qcpy)
  - [8.6 meta.py](#86-metapy)
  - [8.7 release.py](#87-releasepy)
  - [8.8 assessment.py](#88-assessmentpy)
- [9. plots/ — 可视化](#9-plots--可视化)
  - [9.1 sample_growth.R](#91-sample_growthr)

---

## 1. 顶层文件

### 1.1 `subset.py`

**文件用途**: 从 gnomAD v4 VariantDataset 中筛选指定样本子集，输出 VDS 或 VCF 格式。要求 Hail 0.2.120 版本以优化性能和成本。

#### 类: `ProcessingConfig`

数据处理配置的 dataclass，封装了所有筛选、输出和数据加载选项。

| 方法/属性 | 说明 |
|-----------|------|
| `__post_init__()` | 验证参数逻辑一致性：`add_variant_qc` 需 `split_multi=True`；`pass_only` 需 `split_multi=True`；VCF 导出需 `split_multi=True` |
| `from_args(cls, args)` | 类方法，从 argparse 参数创建 `ProcessingConfig` 实例 |

#### 函数

| 函数 | 说明 |
|------|------|
| `get_gnomad_datasets(data_type, n_partitions, test)` | 根据数据类型获取 v4 VDS 和元数据 Table。支持 exomes 和 genomes |
| `filter_to_samples(vds, meta_ht, config)` | 根据样本列表或 workspace 信息过滤 VDS 中的样本 |
| `process_and_export(vds, meta_ht, config)` | 执行多等位基因拆分、PASS-only 过滤、变异 QC 注释添加，并导出为 VDS 或 VCF |
| `main(args)` | 命令行入口：加载数据 → 过滤样本 → 处理导出 |

---

## 2. `sample_qc/` — 样本质控

### 2.1 `hard_filters.py`

**文件用途**: 确定未通过硬过滤阈值的样本。基于污染率、覆盖度、指纹验证等指标进行初步样本过滤。

| 函数 | 说明 |
|------|------|
| `compute_sample_qc(n_partitions, test, n_alt_alleles_strata, n_alt_alleles_strata_name)` | 使用 `compute_stratified_sample_qc` 按替代等位基因数量进行分层样本 QC 计算。按 bi-allelic / multi-allelic 或自定义阈值分层。移除着丝粒和端粒区域、过滤至常染色体，计算 QC 指标 |
| `compute_mean_dp(test)` | 计算每个样本在 chr20 上的平均 DP (深度)，用于基线覆盖度评估 |
| `compute_sample_qc_mt_callrate(test)` | 计算每个样本在预定义 QC 位点上的 call rate (调用率) |
| `determine_hard_filter_expr(ht)` | 构建硬过滤表达式：基于污染率 (>5%)、覆盖度 (<15x chr20 mean DP)、指纹验证失败等条件 |
| `main(args)` | 主入口：执行 sample QC → 计算 mean DP → 计算 callrate → 应用硬过滤 → 标记不合格样本 |

---

### 2.2 `sex_inference.py`

**文件用途**: 推断每个样本的染色体性别核型 (XX/XY/XXY/X0 等)。基于 X 染色体 F-statistic 和性染色体倍性。

| 函数 | 说明 |
|------|------|
| `determine_fstat_sites(vds, approx_af_and_no_callrate, min_af, min_callrate)` | 筛选 chrX 上高质量 SNP (双等位基因、常见、高 callrate) 用于 F-stat 计算。如果 `approx_af_and_no_callrate=True` 则使用 AC/(n_samples×2) 近似 AF 并跳过 callrate 过滤 |
| `compute_sex_ploidy(vds, calling_intervals_ht, normalization_contig)` | 使用 `hl.vds.impute_sex_chromosome_ploidy` 计算性染色体倍性和 F-stat |
| `get_ploidy_cutoffs(ht, f_stat_cutoff)` | 基于 F-stat 分布确定 XX/XY 分类的 F-stat 阈值 |
| `infer_sex_karyotype(ht, ploidy_cutoffs)` | 基于倍性和 F-stat 阈值推断完整的性别核型 (XX/XY/XXY/X0/ambiguous) |
| `compute_sex_chr_interval_coverage(vds, calling_intervals_ht)` | 计算性染色体上 per-interval 的覆盖度用于后续分析 |
| `main(args)` | 主入口：确定 F-stat 位点 → 计算倍性 → 确定阈值 → 推断核型 → 输出性别注释 |

---

### 2.3 `interval_qc.py`

**文件用途**: 基于 per-interval 聚合统计定义高质量捕获区间。支持两种方法：mean fraction over DP 0；fraction of samples with coverage over threshold。支持按平台分层。

| 函数 | 说明 |
|------|------|
| `generate_sex_chr_interval_coverage_mt(vds, calling_intervals_ht)` | 创建性染色体上的 interval-by-sample 覆盖度 MatrixTable。在 PAR 区域边界处分割区间，并标注每个区间是否与 PAR 重叠 |
| `filter_to_test(mt, sex_mt, num_partitions)` | 过滤 MT 为测试子集 (少量分区) 用于快速调试 |
| `compute_interval_qc(mt, autosome_par_only, per_platform, platform_ht, sex_ht, by_mean_fraction_over_dp_0, frac_samples_pass_threshold, mean_dp_thresholds, autosome_dp_threshold, sex_chr_dp_threshold)` | 核心函数：计算 per-interval QC 统计指标。支持按平台分层，区分常染色体和性染色体阈值。可选 "mean fraction over DP 0" 或 "fraction of samples passing DP threshold" 两种方法 |
| `get_high_qual_cutoff_dict(cutoff, by_mean_fraction_over_dp_0, per_platform, autosome_par_only)` | 构建高质量区间判定阈值字典 |
| `get_interval_qc_pass(interval_qc_ht, cutoff_dict, per_platform)` | 基于阈值字典标注每个区间是否通过 QC。返回包含 `pass_interval_qc` 注释的 Table |
| `main(args)` | 主入口：计算区间覆盖度 → 计算区间 QC → 标记高质量区间 |

---

### 2.4 `generate_qc_mt.py`

**文件用途**: 创建密集 (dense) MatrixTable，用于亲缘关系估计和祖源 PCA。联合 v3 基因组和 v4 外显子组，过滤到预定义 QC 位点集。

| 函数 | 说明 |
|------|------|
| `create_filtered_dense_mt(mtds, split)` | 将稀疏 MT 或 VDS 过滤到预定义 QC 变异集，然后密集化 (densify)。保留 GT/GQ/DP/AD 条目 |
| `generate_qc_mt(v3_mt, v4_mt, bi_allelic_only, min_af, min_callrate, min_inbreeding_coeff_threshold, ld_r2, filter_lcr, filter_decoy, filter_segdup)` | 核心函数：合并 v3 和 v4 数据，应用 AF/callrate/IC 过滤，执行 LD 剪枝 (pruning)，生成最终 QC MT |
| `create_joint_qc_meta(v3_meta_ht, v4_meta_ht, v4_hard_filter_ht, v4_sex_ht)` | 创建联合 (v3+v4) 元数据 Table，包含样本来源、硬过滤状态、性别信息 |
| `main(args)` | 主入口：密集化数据 → 合并 → QC 过滤 → LD 剪枝 → 输出 QC MT |

---

### 2.5 `assign_ancestry.py`

**文件用途**: 使用已知 v3 人群标签或 TGP/HGDP 标签为样本分配全球祖源标签。

#### 常量

| 常量 | 说明 |
|------|------|
| `V4_POP_SPIKE_DICT` | 字典：定义特殊人群 (Arab/Bedouin/Persian/Qatari) → `mid` (中东) 的映射 |
| `V3_SPIKE_PROJECTS` | 字典：每个人群对应的 v3 训练队列列表 (asj→Jewish_Genome_Project, ami→NHLBI 等) |

#### 函数

| 函数 | 说明 |
|------|------|
| `get_training_samples(joint_meta_ht, pop_ht, v3_spike_projects, v4_spike_dict, hgdp_tgp_outliers_ht, min_pop_prob)` | 确定用于 RF 训练的已知祖源样本集。整合 v3 已知标签、HGDP/TGP 参考、v4 race/ethnicity spike-in；移除 HGDP/TGP 异常值 |
| `run_pca(qc_mt, related_samples_ht, n_pcs, joint_meta_ht)` | 使用 `run_pca_with_relateds` 执行 PCA (移除相关样本后)，返回特征值、评分、loadings |
| `assign_ancestry(pca_scores_ht, training_ht, n_pcs, min_prob)` | 使用 `assign_genetic_ancestry_pcs` (Random Forest) 将 PCA 评分映射到祖源标签。支持自定义每人群最小 RF 概率阈值 |
| `main(args)` | 主入口：获取训练样本 → PCA → RF 分类 → 输出祖源标签 |

---

### 2.6 `outlier_filtering.py`

**文件用途**: 确定样本 QC 指标异常值并标记过滤。实现三种过滤策略：分层过滤、回归残差过滤、最近邻过滤。

| 函数 | 说明 |
|------|------|
| `get_sample_qc_ht(sample_qc_ht, test, seed)` | 预处理 QC Table：添加 `r_snp_indel` 指标、展开 `bases_over_dp_threshold` 数组为独立注释、移除硬过滤样本 |
| `apply_filter(sample_qc_ht, qc_metrics, filtering_method, apply_r_ti_tv_singleton_filter, pop_scores_ht, pop_ht, platform_ht, ...)` | 核心过滤调度器。根据 `filtering_method` 参数调用不同策略 |
| `apply_stratified_filter(sample_qc_ht, qc_metrics, pop_ht, platform_ht, ...)` | **分层过滤**: 按人群 × 平台分层，计算每层的均值±N个标准差，标记超出阈值的样本 |
| `apply_regressed_filter(sample_qc_ht, qc_metrics, pop_scores_ht, platform_ht, ...)` | **回归残差过滤**: 使用 PCA scores 回归 QC 指标，对残差进行分层过滤。消除人群结构对 QC 指标的混杂影响 |
| `apply_nearest_neighbors_filter(sample_qc_ht, qc_metrics, pop_scores_ht, ...)` | **最近邻过滤**: 基于 PCA 空间中的最近邻 (nearest neighbors) 确定局部参考分布，标记显著偏离邻居的样本 |
| `finalize_outlier_filtering(stratified_ht, regressed_ht, nn_ht)` | 合并三种过滤策略的结果，生成最终的异常值过滤标记 |
| `main(args)` | 主入口：预处理 → 分层过滤 → 回归过滤 → NN 过滤 → 合并 → 输出 |

---

### 2.7 `relatedness.py`

**文件用途**: 计算所有样本对之间的亲缘关系估计值。使用 cuKING (GPU 加速) 和/或 Hail pc_relate。

| 函数 | 说明 |
|------|------|
| `print_cuking_command(cuking_input_path, cuking_output_path, min_emission_kinship, cuking_split_factor)` | 生成并打印用于在 Cloud Batch 上运行 cuKING 的命令行。支持自定义最小发射 kinship 阈值和矩阵分割因子 |
| `prepare_cuking_inputs(mt, n_pcs)` | 准备 cuKING 所需的 Parquet 格式输入文件：从 QC MT 中提取基因型和 PC 相关信息 |
| `load_cuking_outputs(cuking_output_path)` | 加载 cuKING 的 Parquet 输出为 Hail Table，包含成对亲缘系数估计 |
| `run_pc_relate(qc_mt, related_samples_ht, n_pcs)` | 使用 Hail `pc_relate` 计算成对 IBD (identity by descent) 估计 |
| `compute_related_samples_to_drop(relatedness_ht, sample_rankings_ht)` | 使用最大独立集 (maximal independent set) 算法确定需移除的相关样本，优先保留高质量样本 |
| `create_sample_rankings(meta_ht)` | 创建样本排名 Table：基于数据类型、硬过滤状态、覆盖度等为每个样本分配优先级分数 |
| `run_ibd(qc_mt, related_samples_ht)` | 运行 IBD 估计 (Z0/Z1/Z2) 用于进一步亲缘关系分类 |
| `main(args)` | 主入口：准备输入 → cuKING/pc_relate → 加载结果 → 排名 → 确定移除集 |

---

### 2.8 `identify_trios.py`

**文件用途**: 从亲缘关系数据推断三人家系 (trios)，并基于 Mendel 错误和 de novo 变异进行过滤验证。

| 函数 | 说明 |
|------|------|
| `families_to_trios(ped, seed)` | 将包含多个 trios 的家系 Pedigree 转换为每个家系仅保留一个随机 trio 的 Pedigree |
| `filter_relatedness_ht(ht, filter_ht)` | 过滤亲缘关系 Table：仅保留两个样本均为外显子组且未被 QC 过滤的配对 |
| `platform_table_to_dict(platform_ht)` | 将平台 Table 转换为 `{sample_id: platform}` 字典格式 |
| `infer_pedigree(relatedness_ht, sex_ht, duplicate_ht)` | 使用 `infer_families` 从亲缘关系数据推断家系结构 |
| `validate_trios_with_mendel(trio_mt, pedigree)` | 计算每个 trio 的 Mendel 错误率和 de novo 变异数，用于验证家系正确性 |
| `filter_trios(ped, mendel_ht, max_mendel_errors, max_de_novos)` | 基于 Mendel 错误和 de novo 阈值过滤不合格的 trios |
| `main(args)` | 主入口：推断家系 → 识别重复 → 推断 trios → Mendel 验证 → 过滤 |

---

### 2.9 `platform_inference.py`

**文件用途**: 基于 per-interval DP 覆盖度的 PCA 结果，使用 HDBSCAN 聚类算法分配测序平台标签。

| 函数 | 说明 |
|------|------|
| `main(args)` | 完整平台推断流程：①加载区间覆盖度 MT → ②移除硬过滤样本 → ③过滤至常染色体 → ④以 "fraction over dp 0" 为 callrate → ⑤使用 `run_platform_pca` 执行 PCA → ⑥使用 `assign_platform_from_pcs` + HDBSCAN 聚类分配平台 |

---

### 2.10 `create_sample_qc_metadata_ht.py`

**文件用途**: 将所有样本 QC 模块的输出合并为单一元数据 Table (外显子组)。

| 函数 | 说明 |
|------|------|
| `get_project_meta()` | 加载项目元数据并添加 GATK 版本注释。标记 UKB 样本和 fixed homalt 模型状态 |
| `create_sex_imputation_struct(sex_ht)` | 从性别推断 Table 创建结构化的性别推断字段 |
| `create_population_inference_struct(pop_ht, pop_pr_ht)` | 从祖源分类结果创建人群推断结构体，包含 RF 概率和 PCA 得分 |
| `create_sample_qc_struct(sample_qc_ht)` | 从样本 QC 结果创建结构化 QC 指标字段 |
| `create_hard_filter_struct(hard_filter_ht, hard_filter_no_sex_ht)` | 合并有/无性别推断的硬过滤结果 |
| `create_outlier_filter_struct(outlier_ht)` | 从异常值过滤结果创建过滤标记结构体 |
| `create_relatedness_struct(relatedness_ht, related_drop_ht, rankings_ht)` | 创建亲缘关系信息结构体：最近关系度、是否需移除、排名等 |
| `merge_all_metadata(...)` | 核心合并函数：将所有 QC 模块输出按样本 key 联合，生成最终元数据 Table |
| `main(args)` | 主入口：加载所有子模块输出 → 创建结构体 → 合并 → 标记 release 状态 → 输出 |

---

### 2.11 `create_sample_qc_metadata_ht_genomes.py`

**文件用途**: 创建基因组样本的 QC 元数据 Table。处理 HGDP/TGP 子集从 v3.1 到 v4.0 的样本更新 (169 新增、110 移除)。

| 常量 | 说明 |
|------|------|
| `N_DIFF_SAMPLES` | 更新差异样本数: 169 - 110 = 59 |

| 函数 | 说明 |
|------|------|
| `import_updated_annotations(ht, subset_ht)` | 从更新的 HGDP/TGP 子集元数据 HT 导入更新注释 (亚群、freemix、硬过滤、release 状态等)，标注更新后的元数据到完整 v3.1 基因组 meta HT |
| `main(args)` | 主入口：加载 v3.1 元数据 → 导入 HGDP/TGP 更新 → 合并验证 → 输出 v4.0 基因组元数据 |

---

## 3. `variant_qc/` — 变异质控

### 3.1 `import_variant_qc_vcf.py`

**文件用途**: 将变异 QC 结果 VCF (VQSR/RF/IF 输出) 导入为 Hail Table。

| 函数 | 说明 |
|------|------|
| `import_variant_qc_vcf(vcf_path, model_id, num_partitions, import_header_path, array_elements_required, is_split, deduplicate_check)` | 导入变异 QC VCF。根据 model_id 前缀 (`rf_`/`vqsr_`/`if_`) 决定处理方式。解析 AS_VQSLOD、AS_QUALapprox、AS_VarDP、AS_SB_TABLE 等字段；对 VQSR 结果执行多等位基因拆分 |
| `main(args)` | 主入口：导入 VCF → 拆分 → 写入 HT |

---

### 3.2 `random_forest.py`

**文件用途**: 在 gnomAD v4 变异 QC 数据上训练和应用 Random Forest 模型。

#### 常量

| 常量 | 说明 |
|------|------|
| `FEATURES` | RF 特征列表: `AS_MQRankSum`, `AS_pab_max`, `AS_QD`, `AS_ReadPosRankSum`, `AS_SOR`, `allele_type`, `has_star`, `n_alt_alleles`, `variant_type` |
| `LABEL_COL` / `PREDICTION_COL` / `PROBABILITY_COL` / `TRAIN_COL` | RF 训练/预测列名 |

| 函数 | 说明 |
|------|------|
| `train_rf(ht, test, features, fp_to_tp, num_trees, max_depth, transmitted_singletons, sibling_singletons, adj, filter_centromere_telomere, test_intervals, interval_qc_pass_ht)` | 训练 RF 模型。使用 `train_rf_model` 函数，支持自定义 FP:TP 比例、树数量、最大深度；可选 transmitted/sibling singletons 作为 TP 训练数据；可过滤着丝粒/端粒 |
| `get_rf_resources(test, overwrite, model_id)` | 获取 RF 管线的 `PipelineResourceCollection`，管理输入/输出依赖 |
| `main(args)` | 主入口：①训练 RF → ②保存模型 → ③应用模型到完整数据集 → ④中位数插补缺失特征 → ⑤输出预测结果 |

---

### 3.3 `vqsr.py`

**文件用途**: 完整的 GATK VQSR (Variant Quality Score Recalibration) 实现 (~54KB 大文件)。管理从资源准备到模型训练再到重校准应用的完整流程。

> 此文件非常大，包含多个函数和类，负责 VQSR 的完整管线。使用参考资源 (HapMap、Omni、1000G、dbSNP)，通过 `vqsr_resources.json` 配置。

---

### 3.4 `evaluation.py`

**文件用途**: 创建用于评估图表的聚合变异统计 Table，通过评分 bin 进行聚合。

| 函数 | 说明 |
|------|------|
| `create_bin_ht(ht, info_ht, rf_annotations_ht, n_bins, model_type)` | 创建 bin 注释 Table。根据模型类型 (vqsr/rf/if) 添加对应的评分和训练位点标记。过滤 lowqual 变异和低置信区域 |
| `get_evaluation_resources(test, overwrite, model_id)` | 获取评估管线的 `PipelineResourceCollection` |
| `main(args)` | 主入口：①创建 bin HT → ②创建聚合 bin HT → ③真值样本 concordance → ④grouped binned HT |

---

### 3.5 `final_filter.py`

**文件用途**: 创建外显子组的最终过滤 Table。整合 RF/VQSR/IF 评分、阈值确定、PASS/FAIL 标记。

#### 常量

| 常量 | 说明 |
|------|------|
| `REGION_INFO_FIELDS` | 区域信息字段: `non_lcr`, `in_calling_intervals`, `pass_interval_qc` |
| `TRUTH_SET_FIELDS` | 真值集字段: `hapmap`, `omni`, `mills`, `kgp_phase1_hc`, `transmitted_singleton_*`, `sibling_singleton_*` |
| `TRAINING_INFO_FIELDS` | 按模型类型的训练信息字段 (RF/AS_VQSR/IF) |
| `VARIANT_QC_RESULT_FIELDS` | 按模型类型的 QC 结果字段 (RF概率/VQSLOD/IF SCORE) |
| `FINAL_FILTER_FIELDS` | 最终过滤 Table 的顶层注释列表 |

| 函数 | 说明 |
|------|------|
| `process_score_cutoffs(score_cutoffs)` | 处理评分阈值配置，支持 SNP 和 indel 不同阈值 |
| `create_final_filter_ht(ht, info_ht, freq_ht, model_id, score_cutoffs, inbreeding_coeff_cutoff)` | 核心函数：整合所有注释、应用评分阈值、设置 InbreedingCoeff 硬过滤、标记 PASS/FAIL。使用 `add_filters_expr` 构建过滤标记 |
| `get_final_variant_qc_resources(test, overwrite, model_id)` | 获取最终过滤管线的 `PipelineResourceCollection` |
| `main(args)` | 主入口：加载 bin HT 和评分 → 确定阈值 → 生成最终过滤 HT |

---

### 3.6 `final_filter_genomes.py`

**文件用途**: 创建基因组的最终过滤 Table (与外显子组策略不同)。

#### 常量

| 常量 | 说明 |
|------|------|
| `ALLELE_TYPE_FIELDS` | 等位基因类型字段 (移除 `original_alleles` 和 `has_star`) |
| `REGION_INFO_FIELDS` | 区域信息字段: 仅 `non_lcr` |
| `TRUTH_SET_FIELDS` | 真值集字段 (比外显子组少 sibling_singleton) |

| 函数 | 说明 |
|------|------|
| `get_final_variant_qc_resources(test, overwrite, model_id)` | 获取基因组最终过滤管线的 `PipelineResourceCollection` |
| `create_final_filter_ht_genomes(ht, vqsr_ht, freq_ht, score_cutoffs, inbreeding_coeff_cutoff)` | 基因组特定的最终过滤函数。直接使用 v3 VQSR 评分和 bin，添加端粒/着丝粒过滤 |
| `main(args)` | 主入口 |

---

## 4. `annotations/` — 变异注释

### 4.1 `generate_freq.py`

**文件用途**: 生成 v4 外显子组频率数据注释。将 VDS 按子集拆分 → 密集化 → 计算频率/直方图 → 合并 → 修正高 AB het 伪影 → 计算 InbreedingCoeff/FAF/grpmax。

#### 常量

| 常量 | 说明 |
|------|------|
| `AGE_HISTS` | 年龄直方图: `age_hist_het`, `age_hist_hom` |
| `QUAL_HISTS` | 质量直方图: `gq_hist_all`, `dp_hist_all`, `gq_hist_alt`, `dp_hist_alt`, `ab_hist_alt` |
| `FREQ_HIGH_AB_HET_ROW_FIELDS` | 高 AB het 修正相关字段 |
| `FREQ_ROW_FIELDS` | 频率 HT 核心行字段: `freq`, `qual_hists`, `raw_qual_hists`, `age_hists` |

| 函数 | 说明 |
|------|------|
| `split_vds_by_strata(vds, meta_ht, strata_expr)` | 按分层表达式 (如 UKB vs non-UKB) 拆分 VDS |
| `compute_freq_by_strata_from_vds(vds, meta_ht, strata_list, ...)` | 核心频率计算：从 VDS 密集化 → 按分层计算 call stats → 生成频率数组和直方图 |
| `correct_for_high_ab_hets(freq_ht, af_threshold)` | 修正高 allele balance het GATK 伪影：当 AF > 阈值时，从 AC/AN 中移除高 AB het 贡献 |
| `compute_inbreeding_coeff(freq_ht)` | 使用原始 (raw) call stats 计算双等位基因位点的 InbreedingCoeff |
| `generate_faf_grpmax(freq_ht)` | 使用 AB 修正后的频率计算 filtering allele frequency (FAF) 和 grpmax |
| `create_final_freq_ht(freq_ht)` | 创建最终频率 HT：合并所有分层频率、修正、FAF/grpmax |
| `main(args)` | 主入口：拆分 VDS → 计算频率 → 合并 → AB 修正 → IC → FAF/grpmax → 输出 |

---

### 4.2 `generate_freq_genomes.py`

**文件用途**: 创建 v4.0 基因组的频率 HT。专为处理 HGDP/TGP 子集的样本增删更新而设计。合并 v3.1 release、更新后的 HGDP/TGP 以及新增样本的频率数据。

#### 常量

| 常量 | 说明 |
|------|------|
| `POP_MAP` | HGDP 人群名称映射 (如 `bantusafrica` → `bantusouthafrica`) |

| 函数 | 说明 |
|------|------|
| `get_hgdp_tgp_meta(v3_meta_ht, updates, ...)` | 获取更新后的 HGDP/TGP 元数据，整合污染标记、亚群 PCA 异常值、关系推断更新 |
| `compute_new_variant_an(vds, meta_ht, ...)` | 为仅存在于更新子集中的新变异计算 allele number |
| `merge_freq_data(v3_freq_ht, subset_freq_ht, new_an_ht)` | 合并 v3.1、HGDP/TGP 子集、新变异的频率数据 |
| `compute_faf_grpmax_genomes(freq_ht)` | 计算基因组的 FAF 和 grpmax |
| `main(args)` | 主入口：加载 v3.1 频率 → 更新子集频率 → 合并 → FAF/grpmax → 输出 |

---

### 4.3 `generate_variant_qc_annotations.py`

**文件用途**: 生成用于变异 QC 的注释 (info 统计量、VEP、真值集、trio/sibling 统计)。

#### 常量

| 常量 | 说明 |
|------|------|
| `INFO_METHODS` | Info 计算方法: `AS` (allele-specific)、`quasi` (准计算)、`set_long_AS_missing` |
| `INFO_FEATURES` | 用于 QC 的 info 特征: `AS_MQRankSum`, `AS_pab_max`, `AS_MQ`, `AS_QD`, `AS_ReadPosRankSum`, `AS_SOR`, `AS_FS` |
| `NON_INFO_FEATURES` | 非 info 特征: `variant_type`, `allele_type`, `n_alt_alleles`, `was_mixed`, `has_star` |
| `TRUTH_DATA` | 真值数据集: `hapmap`, `omni`, `mills`, `kgp_phase1_hc` |

| 函数 | 说明 |
|------|------|
| `extract_as_pls(lpl_expr, allele_idx)` | 从 LPL (Local PhredLikelihoods) 数组中提取特定等位基因的 PL 值 |
| `get_as_info_expr_from_entry(mt)` | 从 MatrixTable 条目中计算 allele-specific info 注释 (AS_QD, AS_FS, AS_MQ 等) |
| `compute_info(vds, method, ...)` | 核心 info 计算函数：支持三种方法。`AS` 方法使用 allele-specific 聚合；`quasi` 方法使用近似计算；`set_long_AS_missing` 处理超长多等位基因 |
| `generate_trio_stats(mt, pedigree)` | 使用 trios 生成传递/非传递变异统计 (transmitted_singleton) |
| `generate_sib_stats(mt, relatedness_ht)` | 使用同胞对生成兄弟姐妹 singleton 统计 |
| `run_vep(ht, vep_version)` | 使用 VEP (v105 或 v115) 注释变异的功能影响 |
| `create_variant_qc_annotation_ht(info_ht, truth_ht, trio_stats_ht, sib_stats_ht)` | 整合所有 QC 注释为单一 Table |
| `main(args)` | 主入口：计算 info → 拆分 → VEP → 真值集 → trio/sib 统计 → 合并 |

---

### 4.4 `compute_coverage.py`

**文件用途**: 计算 gnomAD v4 外显子组和基因组的覆盖度统计。

| 函数 | 说明 |
|------|------|
| `get_exomes_group_membership_ht(meta_ht, ds_ht, non_ukb_ds_ht)` | 获取外显子组分组成员 HT，包含全量和 non-UKB 的频率分层信息 |
| `adjust_interval_padding(calling_intervals_ht, padding)` | 调整捕获区间的 padding (扩展/收缩边界) |
| `compute_coverage_stats_from_vds(vds, group_membership_ht, ...)` | 核心覆盖度计算：从 VDS 计算 per-base 的覆盖度统计 (mean, median, 各 DP 阈值的 fraction over) |
| `compute_all_sites_an(vds, group_membership_ht, ...)` | 计算所有位点的 allele number，用于频率修正 |
| `main(args)` | 主入口：计算覆盖度 → 计算 all-sites AN → 导出 TSV |

---

### 4.5 `fix_freq_an.py`

**文件用途**: 修正 v4.0 频率 HT 中的 AN 错误 (由于 `hail.vds.filter_samples` 的行为导致 UKB/non-UKB 子集特有变异的 AN 不正确)。

| 函数 | 说明 |
|------|------|
| `prep_vds_for_all_sites_stats(vds)` | 预处理 VDS 以计算所有位点的统计：处理性别倍性调整、AD/PL 缺失填充等 |
| `main(args)` | 主入口：生成全位点 AN → 生成频率修正 AN → 应用修正 → 重新计算 InbreedingCoeff/FAF/grpmax |

---

### 4.6 `insilico_predictors.py`

**文件用途**: 生成 in silico 预测器注释 (SIFT、PolyPhen、CADD、REVEL、SpliceAI、PrimateAI)。

| 函数 | 说明 |
|------|------|
| `get_sift_polyphen_from_vep(ht)` | 从 VEP 105 注释中提取 MANE Select (或 canonical) 转录本的最大 SIFT 和 PolyPhen 评分 |
| `create_cadd_grch38_ht()` | 创建 GRCh38 CADD 评分 Table，合并全基因组 SNV、gnomAD v3.0/v3.1 indel 等多个来源 |
| `_load_cadd_raw(cadd_tsv)` | 内部函数：加载单个 CADD TSV 文件 |
| `main(args)` | 主入口：CADD HT 创建 → SIFT/PolyPhen 提取 → REVEL/SpliceAI/PrimateAI 导入 → 合并输出 |

---

### 4.7 `vep_context_ht.py`

**文件用途**: 对 GRCh38 context Table (所有可能的 SNV) 添加 VEP 注释。

| 函数 | 说明 |
|------|------|
| `main(args)` | 加载现有 VEP 版本的 context HT → 删除旧 VEP 注释 → 运行新版本 VEP (v105/v115) → 添加全局 VEP 版本/配置信息 → 写入结果 |
| `get_script_argument_parser()` | 创建命令行参数解析器 |

---

### 4.8 `recover_and_complete_vep115.py`

**文件用途**: 从 VEP 115 运行失败中恢复并完成注释。处理 chr18 着丝粒区域的 VEP JSON 解析错误。

#### 常量

| 常量 | 说明 |
|------|------|
| `PARTIAL_HT_PATH` | 部分完成的 HT 路径 |
| `CHR18_CENTROMERE_INTERVAL` | chr18 着丝粒排除区间: `chr18:15460900-20861207` |

| 函数 | 说明 |
|------|------|
| `_read_schema_metadata(schema_ref_path)` | 读取参考 HT 的 schema 元数据 |
| `_get_context_ht_path(version)` | 获取 context HT 路径 |
| `_read_json_metadata(file_path)` | 读取 JSON 元数据文件 |
| `_get_sorted_partition_indices(index_to_filename)` | 从索引文件获取已完成分区的排序索引 |
| `_drop_vep_proc_id_if_present(ht)` | 如果存在则删除 `vep_proc_id` 字段 |
| `main(args)` | 主入口 7 步流程：拷贝部分 HT → 提取元数据 → 重建 HT → 过滤未注释变异 → VEP 主体 → VEP chr18 着丝粒 → 合并 |

---

### 4.9 `vrs_annotation_batch.py`

**文件用途**: 批处理脚本：通过创建分片 VCF → 对每个分片运行 vrs-annotation → 合并结果回 Hail Table 的方式添加 VRS (Variant Representation Specification) ID。

#### 常量

| 常量 | 说明 |
|------|------|
| `VRS_SCHEMA_VERSION` | VRS schema 版本: `1.3.0` |
| `VRS_PYTHON_VERSION` | VRS Python 版本: `0.8.4` |
| `SEQREPO_VERSION` | SeqRepo 版本: `2018-11-26` |

| 函数 | 说明 |
|------|------|
| `init_job(batch, name, image, cpu, memory, disk_size)` | 初始化 Hail Batch 作业，设置默认参数 (CPU/内存/磁盘) |
| `export_sharded_vcfs(ht, output_path, n_shards)` | 将 HT 分片导出为多个 VCF 文件 |
| `run_vrs_annotation(batch, vcf_path, output_path, image, header_path)` | 在单个 VCF 分片上运行 VRS 注释 |
| `merge_vrs_results(original_ht, vrs_results_paths)` | 合并 VRS 注释结果回原始 HT |
| `main(args)` | 主入口：分片 VCF → VRS 注释 → 合并 |

---

## 5. `create_release/` — 发布数据创建

### 5.1 `create_release_utils.py`

**文件用途**: 创建 v4 发布所需的通用工具函数和版本常量。

#### 常量

| 常量 | 说明 |
|------|------|
| `DBSNP_VERSION` | `b156` |
| `SIFT_VERSION` | `5.2.2` |
| `POLYPHEN_VERSION` | `2.2.2` |
| `VRS_SCHEMA_VERSION` | `2.0.1` |
| `VRS_PYTHON_VERSION` | `2.2.0` |
| `SEQREPO_VERSION` | `2024-12-20` |
| `GENCODE_VERSIONS` | VEP 版本到 GENCODE 版本映射: `105→Release 39`, `115→Release 49` |
| `MANE_SELECT_VERSIONS` | VEP 版本到 MANE Select 版本映射: `105→v0.95`, `115→v1.4` |
| `VEP_VERSIONS_TO_ADD` | 发布版本到新增 VEP 版本映射: `4.1.1→["115"]` |

| 函数 | 说明 |
|------|------|
| `remove_missing_vep_fields(vep_expr)` | 移除 VEP 105 注释中已排除或全行缺失的字段：`colocated_variants`、`context`、`minimised`、`swissprot`、`trembl`、`uniparc` |

---

### 5.2 `create_release_sites_ht.py`

**文件用途**: 创建用于公开发布的 v4.0 外显子组和基因组 sites HT。

#### 常量

| 常量 | 说明 |
|------|------|
| `SUBSETS_TO_DROP` | 需从 v4.0 基因组中移除的 v3 子集 (除 HGDP/TGP 外) |

| 函数 | 说明 |
|------|------|
| `get_tables_for_release(version)` | 根据版本获取需纳入 release 的 Table 列表 |
| `prepare_freq_for_release(freq_ht, data_type)` | 准备频率数据：过滤分层、重命名字段、确保格式一致 |
| `prepare_vep_for_release(vep_ht, vep_version)` | 准备 VEP 注释：移除缺失字段、更新 LOFTEE end_trunc 过滤器 |
| `create_release_ht(freq_ht, info_ht, vep_ht, filter_ht, coverage_ht, ...)` | 核心函数：整合频率、info、VEP、过滤、覆盖度、in silico 预测器、VRS、区域标记等所有注释。添加全局元数据 (版本号、工具版本、频率索引字典) |
| `main(args)` | 主入口：准备各组件 → 合并 → 添加全局信息 → 输出 release sites HT |

---

### 5.3 `create_combined_faf_release_ht.py`

**文件用途**: 创建 v4 外显子组+基因组联合频率和 FAF Table。包含卡方/Fisher 精确检验和 Cochran–Mantel–Haenszel 检验。

| 函数 | 说明 |
|------|------|
| `filter_gene_to_test(ht, pcsk9, zfy)` | 过滤到 PCSK9 (chr1) 和/或 ZFY (chrY) 基因区间用于测试 |
| `merge_exome_genome_freq(exome_ht, genome_ht)` | 合并外显子组和基因组频率数据 |
| `compute_joint_freq(merged_ht)` | 计算联合频率 (merged AC/AN) |
| `compute_freq_comparison_tests(merged_ht)` | 执行频率比较统计检验：① Hail contingency_table_test (卡方/Fisher) ② CMH 检验 |
| `compute_joint_faf_grpmax(joint_freq_ht)` | 计算联合 FAF 和 grpmax |
| `main(args)` | 主入口：合并频率 → 联合频率 → 统计检验 → FAF/grpmax → 输出 |

---

### 5.4 `create_de_novo_release.py`

**文件用途**: 创建 v4 外显子组的 de novo 变异发布 Table。

| 函数 | 说明 |
|------|------|
| `get_releasable_de_novo_calls_ht(mt, priors_ht, ped, test, n_partitions)` | 从可发布 trios 的密集 MT 获取 de novo 调用。执行 split-multi、annotate_adj、近似缺失的 AD/PL 字段、调用 `default_get_de_novo_expr` |
| `annotate_de_novo_with_vep(ht, vep_ht)` | 使用 VEP 注释 de novo 变异的功能影响 (CSQ 类别、LOFTEE 等) |
| `create_release_de_novo_ht(ht)` | 选择发布所需字段，添加公开频率，过滤低置信区域 |
| `main(args)` | 主入口：de novo 调用 → VEP 注释 → 选择字段 → 输出 |

---

### 5.5 `create_false_dup_liftover.py`

**文件用途**: 为三个临床重要的 GRCh38 假重复基因 (KCNE1、CBS、CRYAA) 创建自定义 liftover 文件。

| 常量 | 说明 |
|------|------|
| `FALSE_DUP_GENES` | `["KCNE1", "CBS", "CRYAA"]` |

| 函数 | 说明 |
|------|------|
| `filter_liftover_to_false_dups(data_type)` | 读取 v2 liftover Table，过滤到 chr21 上三个目标基因的变异 |
| `main(args)` | 主入口：分别过滤外显子组和基因组 liftover → 联合 → 合并频率 → 计算 FAF/grpmax → 输出 |

---

### 5.6 `validate_and_export_vcf.py`

**文件用途**: 验证并导出 gnomAD VCF。定义 VCF INFO 字段、执行完整性检查、生成 VCF header 和数据文件。

#### 常量

| 常量 | 说明 |
|------|------|
| `NEW_SITE_FIELDS` | 新增站点字段: `monoallelic`, `only_het`, `transmitted_singleton` |
| `SITE_FIELDS` | 按数据类型的站点字段 (exomes 额外包含 `sibling_singleton`) |
| `REGION_FLAG_FIELDS` | 按数据类型的区域标记字段 (exomes 包含 `fail_interval_qc` 等) |

| 函数 | 说明 |
|------|------|
| `validate_release_ht(ht, data_type, ...)` | 对 release HT 执行全面验证：检查注释完整性、频率数组长度、全局注释等 |
| `build_vcf_info_dict(data_type, ...)` | 构建 VCF INFO 字段描述字典 |
| `export_vcf(ht, output_path, data_type, ...)` | 导出 VCF 文件：调整类型兼容性、添加 VEP CSQ header、处理 VRS 字段 |
| `main(args)` | 主入口：验证 → 构建 header → 导出 VCF |

---

### 5.7 `make_var_annot_hists.py`

**文件用途**: 生成变异注释的质量直方图 JSON 文件。

| 常量 | 说明 |
|------|------|
| `LOG10_ANNOTATIONS` | 需要对数缩放的注释: `AS_VarDP`, `QUALapprox`, `AS_QUALapprox` |

| 函数 | 说明 |
|------|------|
| `create_frequency_bins_expr_inbreeding(AF)` | 为 InbreedingCoeff 创建频率 bin 表达式: `< 0.0005` 和 `>= 0.0005` |
| `main(args)` | 加载 release HT → 加载直方图参数 → 按频率 bin 聚合注释直方图 → 单独处理 InbreedingCoeff → 导出 JSON |

---

## 6. `assessment/` — 数据评估

### 6.1 `calculate_per_sample_stats.py`

**文件用途**: 计算每个样本的变异计数和聚合样本统计。

#### 常量

| 常量 | 说明 |
|------|------|
| `SUM_STAT_FILTERS` | 过滤组合: variant_qc (none/pass)、capture (ukb/broad/intersect/union)、frequency (all/rare)、consequence (all/lof/missense/synonymous) |
| `LOFTEE_LABELS` | LOFTEE 标签 (v4 移除了 `OS`) |

| 函数 | 说明 |
|------|------|
| `get_per_sample_counts(vds, meta_ht, release_ht, filter_group, ...)` | 核心函数：计算每个样本在各过滤组合下的变异计数 (het/hom/non-ref/singleton 等) |
| `compute_aggregate_stats(per_sample_ht, ancestry_ht)` | 聚合样本统计：均值、分位数 (0/25/50/75/100)，可按祖源分层 |
| `main(args)` | 主入口：计算 per-sample 计数 → 聚合统计 → 按祖源分层 → 输出 |

---

### 6.2 `summary_stats.py`

**文件用途**: 运行 gnomAD v4 数据的汇总统计。

| 函数 | 说明 |
|------|------|
| `main(args)` | 加载 release sites HT → 可选过滤 (interval QC pass / calling interval) → ��用 `get_summary_counts` 按变异类别获取汇总计数 → 支持 100k downsampling → 输出 |

---

## 7. `analyses/` — 特定分析

### 7.1 `grpmax_comps.py`

**文件用途**: 比较 gnomAD v4 和 v2 的 grpmax 统计数据。分析欧洲人群与其他多元化人群之间的等位基因频率差异。

#### 常量

| 常量 | 说明 |
|------|------|
| `DIVERSE_GRPS` | `["afr", "amr", "eas", "mid", "sas"]` |
| `EUR_GRPS` | `{"all_eur": ["nfe", "fin", "asj"], "nfe_only": ["nfe"]}` |
| `AF_THRESHOLDS` | `[0.0001, 0.001, 0.01]` |
| `NS_CONSEQ_TERMS` | 非同义结果项 (HIGH + MEDIUM 影响) |

| 函数 | 说明 |
|------|------|
| `get_eur_freq(ht, eur_grps, version)` | 计算欧洲超级人群的 AF。合并指定子人群的 AC/AN |
| `filter_to_threshold(ht, af_threshold, version, eur_filter)` | 过滤到欧洲 AF < 阈值但 grpmax AF > 阈值的变异 |
| `main(args)` | 主入口：加载 v4/v2 release → 计算欧洲 AF → 跨阈值比较 → 按结果分类 → 输出表格 |

---

## 8. `resources/` — 资源路径定义

### 8.1 `basics.py`

**文件用途**: 定义核心数据路径和数据访问函数。

| 函数/资源 | 说明 |
|-----------|------|
| `get_gnomad_v4_vds(split, remove_hard_filtered_samples, high_quality_only, release_only, test, chrom, autosomes_only, ...)` | **核心数据访问函数**: 获取 gnomAD v4 外显子组 VDS，支持丰富的过滤选项 (样本过滤、染色体过滤、质量过滤、分区控制等) |
| `get_gnomad_v4_genomes_vds(...)` | 获取 gnomAD v4 基因组 VDS |
| `calling_intervals(interval_name, padding)` | 获取捕获区间 Table (UKB/Broad/Union/Intersect) |
| `get_checkpoint_path(name, mt, ht)` | 获取 checkpoint 文件路径 |
| `get_logging_path(name)` | 获取日志文件路径 |
| `qc_temp_prefix(version)` | 获取临时文件前缀 |
| `gnomad_v4_testset_meta` | 测试数据集元数据 |
| `ukb_f_stat` | UKB F-statistic 资源 |
| `all_ukb_samples_to_remove` | 需移除的 UKB 样本列表 |

---

### 8.2 `constants.py`

**文件用途**: 定义版本号和全局常量。

| 常量 | 说明 |
|------|------|
| `CURRENT_VERSION` | 当前版本号 |
| `CURRENT_RAW_VERSION` | 原始数据版本 |
| `CURRENT_SAMPLE_QC_VERSION` | 样本 QC 版本 |
| `CURRENT_RELEASE` | 当前发布版本 |
| `CURRENT_ANNOTATION_VERSION` | 注释版本 |
| `CURRENT_FREQ_VERSION` | 频率版本 |
| `DEFAULT_VEP_VERSION` | 默认 VEP 版本 (`105`) |
| `DATA_TYPES` | 数据类型: `exomes`, `genomes` |
| `RELEASE_DATA_TYPES` | 发布数据类型: `exomes`, `genomes`, `joint` |

---

### 8.3 `annotations.py`

**文件用途**: 定义注释相关的资源路径函数。

| 函数 | 说明 |
|------|------|
| `_annotations_root(version, test, data_type)` | 获取注释文件根路径 |
| `get_info(split, test)` | 获取 info HT 资源 |
| `get_vep(test, data_type, vep_version)` | 获取 VEP 注释 HT 资源 |
| `get_freq(data_type, version, test)` | 获取频率 HT 资源 |
| `get_split_vds(strata)` | 获取分片 VDS 资源 |
| `get_downsampling(data_type)` | 获取 downsampling 资源 |
| `get_variant_qc_annotations(test)` | 获取变异 QC 注释 HT |
| `get_trio_stats(test)` | 获取 trio 统计 HT |
| `get_sib_stats(test)` | 获取 sibling 统计 HT |
| `get_insilico_predictors(predictor)` | 获取 in silico 预测器 HT (CADD/REVEL/SpliceAI/PrimateAI/SIFT/PolyPhen) |
| `get_vrs(data_type)` | 获取 VRS 注释 HT |
| `get_all_sites_an_and_qual_hists(data_type)` | 获取全位点 AN 和质量直方图 |
| `get_combined_frequency()` | 获取联合频率 HT |
| `get_freq_comparison()` | 获取频率比较检验 HT |
| `hgdp_tgp_updated_callstats` | HGDP/TGP 更新的 call stats 资源 |

---

### 8.4 `sample_qc.py`

**文件用途**: 定义样本 QC 相关的所有资源路径。

| 函数/资源 | 说明 |
|-----------|------|
| `get_sample_qc_root(version, test, data_type)` | 获取样本 QC 根路径 |
| `get_sample_qc(strat, test, data_type)` | 获取样本 QC HT (bi_allelic/multi_allelic/all) |
| `fingerprinting_failed` | 指纹验证失败样本 |
| `sample_chr20_mean_dp` | 样本 chr20 平均 DP |
| `sample_qc_mt_callrate` | 样本 QC MT callrate |
| `contamination` | 污染估计 |
| `hard_filtered_samples` | 硬过滤后的不合格样本 |
| `hard_filtered_samples_no_sex` | 无性别推断的硬过滤样本 |
| `sex` | 性别推断结果 |
| `f_stat_sites` | F-stat 计算位点 |
| `ploidy` | 倍性估计 |
| `sex_chr_coverage` | 性染色体覆盖度 |
| `interval_coverage` | 区间覆盖度 MT |
| `interval_qc` | 区间 QC 结果 |
| `interval_qc_pass` | 通过区间 QC 的区间 |
| `platform` | 平台分配结果 |
| `platform_pca_*` | 平台 PCA 相关资源 (eigenvalues/loadings/scores) |
| `get_joint_qc(test)` | 联合 QC MT |
| `joint_qc_meta` | 联合 QC 元数据 |
| `predetermined_qc_sites` | 预定义 QC 位点集 |
| `ancestry_pca_*` | 祖源 PCA 相关资源 |
| `get_pop_ht(test)` | 人群分配 HT |
| `pop_rf_path` | RF 模型路径 |
| `relatedness` | 亲缘关系 HT |
| `related_samples_to_drop` | 需移除的相关样本 |
| `sample_rankings` | 样本排名 HT |
| `pedigree` | 家系 Pedigree |
| `trios` | 已验证 trios |
| `duplicates` | 重复样本 |
| `stratified_filtering` / `regressed_filtering` / `nearest_neighbors_filtering` | 三种异常值过滤 HT |
| `finalized_outlier_filtering` | 最终异常值过滤 HT |
| `dense_trio_mt(releasable, split, test)` | 密集 trio MT |
| `get_cuking_input_path()` / `get_cuking_output_path()` | cuKING I/O 路径 |
| HGDP/TGP 相关资源 | `hgdp_tgp_meta_updated`, `hgdp_tgp_populations_updated`, `hgdp_tgp_pop_outliers`, `hgdp_tgp_related_*` 等 |

---

### 8.5 `variant_qc.py`

**文件用途**: 定义变异 QC 相关的资源路径。

| 函数/资源 | 说明 |
|-----------|------|
| `TRUTH_SAMPLES` | 真值样本字典 (NA12878 等)，包含 s、data_type 等信息 |
| `VQSR_FEATURES` | VQSR 使用的特征列表 |
| `get_variant_qc_result(model_id)` | 获取变异 QC 结果 HT |
| `get_rf_model_path(model_id)` | 获取 RF 模型文件路径 |
| `get_rf_run_path()` | 获取 RF 运行信息路径 |
| `get_rf_training(model_id)` | 获取 RF 训练数据 HT |
| `get_score_bins(model_id, aggregated)` | 获取评分 bin HT |
| `get_binned_concordance(model_id, truth_sample)` | 获取 binned concordance HT |
| `get_callset_truth_data(truth_sample)` | 获取 callset 真值数据 HT |
| `final_filter` | 最终过滤 HT |

---

### 8.6 `meta.py`

**文件用途**: 定义元数据相关资源。

| 资源 | 说明 |
|------|------|
| `meta(data_type)` | 样本元数据 HT |
| `project_meta` | 项目元数据 HT |
| `gatk_versions` | GATK 版本信息 HT |
| `meta_tsv_path(data_type)` | 元数据 TSV 路径 |

---

### 8.7 `release.py`

**文件用途**: 定义发布数据相关的资源路径。

| 函数/资源 | 说明 |
|-----------|------|
| `release_sites(data_type, public)` | 发布 sites HT |
| `release_coverage(data_type)` | 覆盖度发布 HT |
| `release_all_sites_an` | 全位点 AN 发布 HT |
| `release_vcf_path(data_type)` | VCF 输出路径 |
| `release_header_path(data_type)` | VCF header 路径 |
| `validated_release_ht(data_type)` | 验证后的 release HT |
| `release_de_novo` | de novo 发布 HT |
| `release_ht_path(data_type, public)` | release HT 路径 |
| `qual_hists_json_path(data_type)` | 质量直方图 JSON 路径 |
| `annotation_hists_params_path(data_type)` | 注释直方图参数路径 |
| `get_false_dup_genes_path()` | 假重复基因 liftover 路径 |
| `get_freq_array_readme()` | 频率数组 README 路径 |
| `get_combined_frequency()` | 联合频率 HT |

---

### 8.8 `assessment.py`

**文件用途**: 定义评估相关的资源路径。

| 函数/资源 | 说明 |
|-----------|------|
| `get_per_sample_counts(data_type, filter_group, test)` | 每样本变异计数 HT |
| `get_summary_stats_filtering_groups(data_type)` | 汇总统计过滤组 |
| `release_summary_stats(test, data_type, filter_name)` | 发布汇总统计 HT |

---

## 9. `plots/` — 可视化

### 9.1 `sample_growth.R`

**文件用途**: R 脚本，使用 ggplot2 绘制 gnomAD 各版本的样本增长曲线图。展示外显子组和基因组样本数量在 v2/v3/v4 各版本间的增长趋势。

---

## 附录: v4 QC 管线完整执行流程

```
╔══════════════════════════════════════════════════════════════════╗
║                      gnomAD v4 QC 管线流程                       ║
╠══════════════════════════════════════════════════════════════════╣
║                                                                  ║
║  ┌─ sample_qc ────────────────────────────────────────────────┐ ║
║  │  1. hard_filters.py        → 硬过滤 (污染/覆盖度/指纹)     │ ║
║  │  2. sex_inference.py       → 性别推断 (F-stat/倍性/核型)    │ ║
║  │  3. interval_qc.py         → 区间 QC (高质量区间标记)       │ ║
║  │  4. platform_inference.py  → 平台检测 (PCA + HDBSCAN)      │ ║
║  │  5. generate_qc_mt.py      → QC MT 生成 (v3+v4密集��)      │ ║
║  │  6. assign_ancestry.py     → 祖源分类 (PCA + RF)           │ ║
║  │  7. relatedness.py         → 亲缘关系 (cuKING/pc_relate)   │ ║
║  │  8. outlier_filtering.py   → 异常值过滤 (3种策略)          │ ║
║  │  9. identify_trios.py      → 家系识别 (Mendel验证)         │ ║
║  │ 10. create_sample_qc_metadata_ht.py → 合并元数据           │ ║
║  └────────────────────────────────────────────────────────────┘ ║
║                            ↓                                     ║
║  ┌─ annotations ──────────────────────────────────────────────┐ ║
║  │  1. generate_variant_qc_annotations.py → info/VEP/真值集   │ ║
║  │  2. generate_freq.py        → 外显子组频率                  │ ║
║  │  3. generate_freq_genomes.py → 基因组频率                   │ ║
║  │  4. fix_freq_an.py          → AN 修正 (v4.1)               │ ║
║  │  5. compute_coverage.py     → 覆盖度统计                    │ ║
║  │  6. insilico_predictors.py  → CADD/REVEL/SpliceAI          │ ║
║  │  7. vep_context_ht.py       → VEP context                  │ ║
║  │  8. vrs_annotation_batch.py → VRS ID                       │ ║
║  └─────────────────────────────────────���──────────────────────┘ ║
║                            ↓                                     ║
║  ┌─ variant_qc ──────────────────────────────────────────────┐ ║
║  │  1. import_variant_qc_vcf.py → 导入 VQSR/RF VCF           │ ║
║  │  2. random_forest.py        → RF 训练和应用                 │ ║
║  │  3. vqsr.py                 → VQSR 管线                     │ ║
║  │  4. evaluation.py           → 评估 (binned concordance)    │ ║
║  │  5. final_filter.py         → 外显子组最终过滤              │ ║
║  │  6. final_filter_genomes.py → 基因组最终过滤                │ ║
║  └────────────────────────────────────────────────────────────┘ ║
║                            ↓                                     ║
║  ┌─ create_release ──────────────────────────────────────────┐ ║
║  │  1. create_release_sites_ht.py → 整合 release HT           │ ║
║  │  2. create_combined_faf_release_ht.py → 联合 FAF            │ ║
║  │  3. validate_and_export_vcf.py → 验证并导出 VCF            │ ║
║  │  4. make_var_annot_hists.py → 质量直方图                    │ ║
║  │  5. create_de_novo_release.py → de novo 发布               │ ║
║  │  6. create_false_dup_liftover.py → 假重复基因 liftover      │ ║
║  └────────────────────────────────────────────────────────────┘ ║
║                            ↓                                     ║
║  ┌─ assessment ──────────────────────────────────────────────┐ ║
║  │  1. calculate_per_sample_stats.py → 样本级统计              │ ║
║  │  2. summary_stats.py         → 汇总统计                     │ ║
║  └────────────────────────────────────────────────────────────┘ ║
║                                                                  ║
╚══════════════════════════════════════════════════════════════════╝
```

---

> **注意**: 本文档基于对仓库源代码的深度分析编写。由于部分大文件 (如 `vqsr.py` 53KB、`outlier_filtering.py` 55KB、`generate_freq_genomes.py` 80KB) 内容截断，部分内部辅助函数可能未被完全列出。完整源代码请参考 [GitHub 仓库](https://github.com/broadinstitute/gnomad_qc/tree/2e35215ce15d7b6b5f5c82e77a6649af4510cf76/gnomad_qc/v4)。