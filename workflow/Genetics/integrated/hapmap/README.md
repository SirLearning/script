# HAPMAP Construction Workflow

- **Nextflow entry**: `workflow/Genetics/integrated/hapmap/main.nf`
- **Python module entry**: `python -m genetics.wheat.hapmap`

## Usage

```bash
nextflow run /home/runner/work/script/script/workflow/Genetics/integrated/hapmap/main.nf \
  --input <genotype_table.tsv> \
  --output_dir <out_dir> \
  --window_size 100000
```

## Parameters

- `--input`: genotype table containing `CHR`, `POS`, and sample genotype columns
- `--output_dir`: output directory
- `--window_size`: block window size for haplotype block stitching
