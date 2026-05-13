# Nextflow run log (`NF_CMD`)

Chronological log of **Nextflow** command lines for `workflow/Genetics`. Append **new runs at the end** with a `### YYYY-MM-DD` (or dated tag) heading, optional one-line working directory / outcome, then a fenced `bash` block. See `.cursor/rules/workstation-nextflow.mdc`.

**Note:** Stats processes often use `publishDir` with `mode: 'copy'`. Changing `output_prefix` produces **new** filenames under `stats/.../plots` (and related dirs) but **does not delete** older PNG/TSV names. If plots “look unchanged,” confirm you are opening files whose basenames match the **current** prefix (e.g. `*.variant.ld_decay.*`), not leftover files from an older prefix.

---

### 2026-01-30 — `test_plink` (example 20260130.1)

```bash
nextflow run /data/home/tusr1/git/script/workflow/Genetics/main.nf \
    --home_dir /data/home/tusr1/01projects/vmap4 \
    --user_dir /data/home/tusr1 \
    --src_dir /data/home/tusr1/git/script/src \
    --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
    --mod test_plink \
    --job test_plink
```

---

### 2026-01-30 — `v1_plink` / `rebuild` (example 20260130.2)

```bash
nextflow run /data/home/tusr1/git/script/workflow/Genetics/main.nf \
    --home_dir /data/home/tusr1/01projects/vmap4 \
    --user_dir /data/home/tusr1 \
    --src_dir /data/home/tusr1/git/script/src \
    --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
    --mod v1_plink \
    --job rebuild
```

---

### 2026-02-05 — `test_plink` + `process_dir` (example 20260205.1)

```bash
nextflow run /data/home/tusr1/git/script/workflow/Genetics/main.nf \
    --home_dir /data/home/tusr1/01projects/vmap4 \
    --user_dir /data/home/tusr1 \
    --src_dir /data/home/tusr1/git/script/src \
    --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
    --process_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process \
    --mod test_plink \
    --job test_plink
```

---

### 2026-02-09 — `v1_plink` / `rebuild` + `process_dir` (example 20260209.1)

```bash
nextflow run /data/home/tusr1/git/script/workflow/Genetics/main.nf \
    --home_dir /data/home/tusr1/01projects/vmap4 \
    --user_dir /data/home/tusr1 \
    --src_dir /data/home/tusr1/git/script/src \
    --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
    --process_dir /data1/dazheng_tusr1/vmap4.VCF.v1/rebuild/process \
    --mod v1_plink \
    --job rebuild
```

---

### 2026-04-01 — `test_plink_camp` + `camp` (example 20260401.1)

```bash
nextflow run /data/home/tusr1/git/script/workflow/Genetics/main.nf \
    --home_dir /data/home/tusr1/01projects/vmap4 \
    --user_dir /data/home/tusr1 \
    --src_dir /data/home/tusr1/git/script/src \
    --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
    --process_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process \
    --camp /data1/dazheng_tusr1/vmap4.VCF.v1/camp_vmap4_map.tsv \
    --mod test_plink_camp \
    --job test_plink
```

---

### 2026-05-13 — `test_thin` restats (reuse merged `*_test.plink2`)

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/11run_test_thin_restats`. Outcome: exit 0, 52 tasks, ~21m 50s.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/11run_test_thin_restats && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/main.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --process_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/test_thin \
  --mod test_thin \
  --job test_plink
```

---

### 2026-05-13 — `test_common_thin` restats (reuse merged `*_test.plink2`)

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/12run_test_common_thin_restats`. Outcome: exit 0, 52 tasks, ~1m 38s.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/12run_test_common_thin_restats && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/main.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --process_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/test_common_thin \
  --mod test_common_thin \
  --job test_plink
```

---

### 2026-05-13 — LD plots only (`tmp/ld_plots_redraw.nf`), `test_thin`

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/14run_ld_plots_redraw`. Outcome: exit 0, 8 tasks, ~15m 22s.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/14run_ld_plots_redraw && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/ld_plots_redraw.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_thin
```

---

### 2026-05-13 — LD plots only (`tmp/ld_plots_redraw.nf`), `test_common_thin` (planned)

Use a **new** numbered run folder (e.g. `15run_ld_plots_redraw`).

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/15run_ld_plots_redraw && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/ld_plots_redraw.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_common_thin
```
