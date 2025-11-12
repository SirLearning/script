# GWAS Benchmark (Python)

Tools to benchmark different GWAS software outputs, perform model selection, and evaluate results using scikit-learn metrics.

## Install (recommended env)

```bash
pip install pandas numpy scikit-learn matplotlib seaborn pyyaml scipy
```

## Inputs
- GWAS result tables from different tools (PLINK, GEMMA, FarmCPU, etc.). Columns auto-detected.
- Ground truth file (causal SNPs or regions):
  - Point variants: columns `CHR POS`
  - Regions: columns `CHR START END`

## Manifest
YAML example:
```yaml
methods:
  - name: GEMMA_LMM
    file: /path/to/gemma.assoc.txt
  - name: PLINK_MLM
    file: /path/to/plink.mlm.txt
truth: /path/to/causal_truth.tsv
truth_window: 50000
phenotype: height
```

CSV example:
```csv
name,file
GEMMA_LMM,/path/to/gemma.assoc.txt
PLINK_MLM,/path/to/plink.mlm.txt
```

## Run
```bash
python -m gwas.benchmark --manifest bench.yaml --outdir out/bench
# or
python -m gwas.benchmark --manifest bench.csv --truth truth.tsv --outdir out/bench --topk 10,20,50,100
```

## Outputs
- `metrics_summary.tsv`: per-method ROC-AUC, PR-AUC, recall@K, enrichment@K, lambda_GC
- `<method>.labeled.tsv`: harmonized table with SCORE and HIT columns
- Plots: `roc_auc_bar.png`, `pr_auc_bar.png`, `qq_<method>.png`
- `best_model.txt`: recommendation by PR-AUC

## Notes
- Matching truth to GWAS hits uses ±window around point truths and interval overlap otherwise.
- If your GWAS tables lack CHR/POS but have SNP like `1:12345`, they will be parsed automatically.
- For very large files, consider pre-filtering by P-value to speed up plotting.
