# Auxiliary Nextflow entry scripts

Small workflows in this folder reuse the same **`workflow/Genetics/nextflow.config`** as `main.nf`. Pass the config with an **absolute** path so it works no matter what the launch directory is:

```bash
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/ld_plots_redraw.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  ...
```

Entry scripts use `include { ... } from '../genotype/...'` paths relative to this directory.

- `ld_plots_redraw.nf` — LD decay + cross-chr plots from existing `.vcor`.
- `assess_plink_debug.nf` — §9 tier-1 assess on `test_thin` / `test_common_thin`: export narrow VCF from `*_test.plink2`, then `quick_count` + `bcftools_qc_assess` + MAF-bin TSV under `assess/<mod>/`.
