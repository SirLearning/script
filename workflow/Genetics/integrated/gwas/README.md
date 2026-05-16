# GWAS Workflow

- **Nextflow entry**: `workflow/Genetics/integrated/gwas/main.nf`
- **Python module entry**: `python -m genetics.wheat.wheat_gwas`

## Usage

```bash
nextflow run /home/runner/work/script/script/workflow/Genetics/integrated/gwas/main.nf \
  --genotype <genotype_matrix.tsv> \
  --phenotype <phenotype.tsv> \
  --trait Trait \
  --output_dir <out_dir>
```

## Parameters

- `--genotype`: genotype matrix with `Sample` and marker columns
- `--phenotype`: phenotype table with `Sample` and trait column
- `--trait`: trait column name
- `--output_dir`: output directory
