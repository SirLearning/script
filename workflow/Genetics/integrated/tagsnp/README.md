# tagSNP Workflow

- **Nextflow entry**: `workflow/Genetics/integrated/tagsnp/main.nf`
- **Python module entry**: `python -m genetics.wheat.tagsnp`

## Usage

```bash
nextflow run /home/runner/work/script/script/workflow/Genetics/integrated/tagsnp/main.nf \
  --input <genotype_matrix.tsv> \
  --max_tags 1000 \
  --ld_threshold 0.8 \
  --output_dir <out_dir>
```

## Parameters

- `--input`: genotype matrix with `Sample` and marker columns
- `--max_tags`: upper limit of selected tagSNP count
- `--ld_threshold`: LD correlation threshold for pruning
- `--output_dir`: output directory
