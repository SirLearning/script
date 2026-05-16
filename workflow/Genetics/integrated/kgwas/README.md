# kGWAS Workflow

- **Nextflow entry**: `workflow/Genetics/integrated/kgwas/main.nf`
- **Python module entry**: `python -m genetics.wheat.kgwas`

## Usage

```bash
nextflow run /home/runner/work/script/script/workflow/Genetics/integrated/kgwas/main.nf \
  --kmer_matrix <kmer_matrix.tsv> \
  --phenotype <phenotype.tsv> \
  --trait Trait \
  --output_dir <out_dir>
```

## Parameters

- `--kmer_matrix`: sample-by-k-mer matrix with `Sample` column
- `--phenotype`: phenotype table with `Sample` and trait column
- `--trait`: phenotype trait column
- `--output_dir`: output directory
