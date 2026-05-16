# CNV Workflow

- **Nextflow entry**: `workflow/Genetics/integrated/cnv/main.nf`
- **Python module entry**: `python -m genetics.wheat.cnv`

## Usage

```bash
nextflow run /home/runner/work/script/script/workflow/Genetics/integrated/cnv/main.nf \
  --input <normalized_depth.tsv> \
  --del_z -2 \
  --dup_z 2 \
  --output_dir <out_dir>
```

## Parameters

- `--input`: normalized depth table with `CHR`, `START`, `END`, and sample depth columns
- `--del_z`: deletion z-score cutoff
- `--dup_z`: duplication z-score cutoff
- `--output_dir`: output directory
