# Genetic Map Workflow

- **Nextflow entry**: `workflow/Genetics/integrated/genetic_map/main.nf`
- **Python module entry**: `python -m genetics.wheat.genetic_map`

## Usage

```bash
nextflow run /home/runner/work/script/script/workflow/Genetics/integrated/genetic_map/main.nf \
  --input <marker_genotype.tsv> \
  --output_dir <out_dir>
```

## Parameters

- `--input`: marker table with `Marker`, `CHR`, `POS`, and sample genotype columns
- `--output_dir`: output directory
