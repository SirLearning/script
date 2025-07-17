# 小麦GWAS分析报告

## 分析概述
本报告基于GWAS Catalog数据库中的小麦相关遗传变异数据进行分析。

## 数据概况
- **总SNP数量**: 45
- **显著性SNP数量** (P < 5×10⁻⁸): 43
- **涉及染色体数量**: 5
- **研究数量**: 4
- **基因组膨胀因子** (λ): 88.950

## 染色体分布
CHR_ID
6     34
3      2
8      1
16     1
4      1

## 研究性状分布
DISEASE/TRAIT
seropositivity for triticum aestivum (common wheat) peptide (twist_14325)               10
seropositivity for triticum urartu (wild einkorn wheat) peptide (twist_11627)           10
seropositivity for triticum urartu (wild einkorn wheat) peptide (twist_12614)            4
never eat sugar vs no eggs, dairy, wheat or sugar restrictions (ukb data field 6144)     4
seropositivity for triticum aestivum wheat peptide (twist_41103)                         3
wheat-dependent exercise-induced anaphylaxis                                             3
never eat wheat vs no wheat restrictions (ukb data field 6144)                           3
seropositivity for triticum monococcum (einkorn wheat) peptide (twist_2100)              2
hydrolysed wheat protein allergy                                                         2
seropositivity for triticum monococcum (einkorn wheat) peptide (twist_25044)             1
never eat wheat vs no eggs, dairy, wheat or sugar restrictions (ukb data field 6144)     1

## 统计分析
- **P值范围**: 1e-151 - 9e-06
- **显著性阈值**: P < 5×10⁻⁸

## 输出文件
- 完整GWAS数据: `wheat_gwas_data.csv`
- 显著性SNP: `significant_snps.csv`
- 曼哈顿图: `wheat_gwas_manhattan_plot.png`
- QQ图: `wheat_gwas_qq_plot.png`

## 数据来源
- **数据库**: GWAS Catalog
- **分析工具**: Python pandas, matplotlib
- **数据文件**: gwas_catalog.pkl from biomni data lake

---
*报告生成时间: 2025-07-13 15:54:18*
