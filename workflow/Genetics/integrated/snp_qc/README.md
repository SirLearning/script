# SNP Calling QC Workflow

- **Nextflow entry**: `workflow/Genetics/integrated/snp_qc/main.nf`
- **Python module entry**: `python -m genetics.wheat.snp_qc`

## Usage

```bash
nextflow run /home/runner/work/script/script/workflow/Genetics/integrated/snp_qc/main.nf \
  --input <snp_table.tsv> \
  --output_dir <out_dir> \
  --maf 0.05 \
  --max_missing 0.1 \
  --min_qual 30
```

## Parameters

- `--input`: SNP table with `MAF`, `MISSING_RATE`, `QUAL`
- `--output_dir`: output directory
- `--maf`: minimum MAF threshold
- `--max_missing`: maximum missing rate
- `--min_qual`: minimum quality threshold
