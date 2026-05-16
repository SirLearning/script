# Genetics Nextflow Pipeline (nf-core style layout)

## Pipeline purpose

`workflow/Genetics` provides genotype and wheat integrated workflows for population genetics and downstream statistics.

## Entry point

```bash
nextflow run /absolute/path/to/workflow/Genetics/main.nf \
  -c /absolute/path/to/workflow/Genetics/nextflow.config \
  --output_dir <out> --job <job> --mod <mod> [other params]
```

## Key parameters

- `mod`: execution mode (`v1_plink`, `test_plink`, `test_thin`, `test_plink_camp`, `test_camp`, `test_common_thin`, `wheat_*`)
- `job`: run identifier
- `output_dir`: output root
- `home_dir`, `src_dir`: required for VCF-based (non-`wheat_*`) modes
- `user_dir`: conda stats environment base, required by wheat integrated modes
- `camp`: required for `test_plink_camp` / `test_camp`

## nf-core style directories

- `main.nf`: top-level entry workflow
- `subworkflows/local/genetics_pipeline.nf`: router + input checks
- `modules/local/*.nf`: reusable workflow modules
- `conf/*.config`: split config files (`base`, `resources`, `modules`)
- `nextflow_schema.json`: parameter schema
- `tests/nextflow.config`: test-only config placeholder

## Notes

- Functionality is preserved from the previous layout.
- Existing domain workflows remain in `genotype/`, `dynamic/`, `static/`, `integrated/`, and `tmp/`.
