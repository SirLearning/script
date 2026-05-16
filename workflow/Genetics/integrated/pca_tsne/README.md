# PCA and t-SNE Workflow

- **Nextflow entry**: `workflow/Genetics/integrated/pca_tsne/main.nf`
- **Python module entry**: `python -m genetics.wheat.population_structure`

## Usage

```bash
nextflow run /home/runner/work/script/script/workflow/Genetics/integrated/pca_tsne/main.nf \
  --input <genotype_matrix.tsv> \
  --output_dir <out_dir> \
  --n_pcs 10 \
  --tsne_perplexity 30
```

## Parameters

- `--input`: sample-by-marker numeric matrix with `Sample` column
- `--output_dir`: output directory
- `--n_pcs`: PCA component number
- `--tsne_perplexity`: t-SNE perplexity
