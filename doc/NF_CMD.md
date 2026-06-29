# Nextflow run log (`NF_CMD`)

Chronological log of **Nextflow** command lines for `workflow/Genetics`. Append **new runs at the end** with a `### YYYY-MM-DD` (or dated tag) heading, optional one-line working directory / outcome, then a fenced `bash` block. See `.cursor/rules/workstation-nextflow.mdc`.

**Router (VCF/PLINK branch, `main.nf`):** As of **2026-05-16**, only `--mod v1_plink`, `test_thin`, `test_camp`, and `test_common_thin` are dispatched (plus any `wheat_*` integrated mod). Older blocks below may still use `--mod test_plink` or `test_plink_camp`; for **new** full-pipeline runs, substitute **`test_thin`** and **`test_camp`** respectively (same processor workflows; `job` / `output_dir` paths are unchanged).

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

### 2026-05-13 — LD plots only (`tmp/ld_plots_redraw.nf`), `test_thin` (re-run, `*.variant.ld_decay` prefix)

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/16run_ld_plots_redraw_test_thin`. Outcome: exit 0, 8 tasks succeeded, ~15m 57s.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/16run_ld_plots_redraw_test_thin && \
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

### 2026-05-13 — LD plots only (`tmp/ld_plots_redraw.nf`), `test_common_thin` (re-run, `*.variant.ld_decay` prefix)

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/17run_ld_plots_redraw_test_common_thin`. Outcome: exit 0, 8 tasks succeeded (same session immediately after `test_thin`).

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/17run_ld_plots_redraw_test_common_thin && \
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

---

### 2026-05-13 — Assess debug (`tmp/assess_plink_debug.nf`), `test_thin`

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/21run_assess_debug_test_thin`. Outcome: exit 0, 16 tasks succeeded, ~2m 42s. Outputs under `.../test_plink/assess/test_thin/` (`export/`, `info/*.mac_site_histogram.tsv`, per-subgenome `*.maf_missing.tsv`, `*.counts.tsv`).

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/21run_assess_debug_test_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/assess_plink_debug.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_thin
```

---

### 2026-05-13 — Assess debug (`tmp/assess_plink_debug.nf`), `test_common_thin`

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/22run_assess_debug_test_common_thin`. Outcome: exit 0, 16 tasks succeeded, ~43s wall clock. Outputs under `.../test_plink/assess/test_common_thin/`.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/22run_assess_debug_test_common_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/assess_plink_debug.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_common_thin
```

---

### 2026-05-14 — Assess debug (`tmp/assess_plink_debug.nf`) — PLINK2 slice + `assess_slice.py` plots

Working directory (example): `/data/home/tusr1/01projects/vmap4/08stats.genome/23run_assess_plink2_debug_stub`. Agent session: `nextflow ... -preview` exit 0 against real `--output_dir` / `--job` (no Graphviz). Re-run without `-preview` after `test_plink` processor outputs exist. Published layout: `assess/<mod>/*.afreq`, `*.vmiss`, `*.tsv`, `plots/*.png`, `info/*.tsv`, `logs/*.log`.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/23run_assess_plink2_debug_stub && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/assess_plink_debug.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_thin
```

---

### 2026-05-15 — Assess debug (`tmp/assess_plink_debug.nf`) — full run `test_thin` (PLINK2 + plots)

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/24run_assess_plink2_full_test_thin`. Outcome: exit 0, **8** succeeded (4× slice + 4× Python plots); wall ~65 s; Nextflow trace aggregate succeedDuration ~7m 25s. Outputs under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/assess/test_thin/`. Graphviz warning only.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/24run_assess_plink2_full_test_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/assess_plink_debug.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_thin
```

---

### 2026-05-15 — Assess debug (`tmp/assess_plink_debug.nf`) — full run `test_common_thin` (PLINK2 + plots)

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/25run_assess_plink2_full_test_common_thin`. Outcome: exit 0, **8** succeeded, succeedDuration ~1m 25s. Outputs under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/assess/test_common_thin/`.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/25run_assess_plink2_full_test_common_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/assess_plink_debug.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_common_thin
```

---

### 2026-05-16 — Wheat PCA/t-SNE from merged PLINK2 (`tmp/wheat_integrated_from_plink.nf`), `test_thin`

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/26run_wheat_pca_test_thin`. Outcome: **8** succeeded (~12m 30s). Genotype matrices: `test_plink/process/test_thin/integrated/geno_matrix/`; PCA/t-SNE: `test_plink/integrated/test_thin/{info,plots}/`.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/26run_wheat_pca_test_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/wheat_integrated_from_plink.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --user_dir /data/home/tusr1 \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_thin
```

---

### 2026-05-16 — Wheat PCA/t-SNE from merged PLINK2 (`tmp/wheat_integrated_from_plink.nf`), `test_common_thin`

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/27run_wheat_pca_test_common_thin`. Outcome: **8** succeeded (~11m 1s). Genotype matrices: `test_plink/process/test_common_thin/integrated/geno_matrix/`; PCA/t-SNE: `test_plink/integrated/test_common_thin/{info,plots}/`.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/27run_wheat_pca_test_common_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/wheat_integrated_from_plink.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --user_dir /data/home/tusr1 \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_common_thin
```

---

### 2026-05-17 — Assess debug (`tmp/assess_plink_debug.nf`) — revalidation `test_thin`

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/28run_assess_revalidate_test_thin`. Outcome: exit 0, **8** succeeded; wall ~1m 14s; Nextflow trace `succeedDuration` ~9m 19s (aggregate). Published refreshed artefacts under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/assess/test_thin/`. Graphviz warning only.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/28run_assess_revalidate_test_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/assess_plink_debug.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_thin
```

---

### 2026-05-17 — Assess debug (`tmp/assess_plink_debug.nf`) — revalidation `test_common_thin`

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/29run_assess_revalidate_test_common_thin`. Outcome: exit 0, **8** succeeded; wall ~23s; Nextflow trace `succeedDuration` ~1m 9s (aggregate). Published refreshed artefacts under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/assess/test_common_thin/`. Graphviz warning only.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/29run_assess_revalidate_test_common_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/assess_plink_debug.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_common_thin
```

---

### 2026-05-17 — Assess debug (`tmp/assess_plink_debug.nf`) — PLINK2 `--freq counts` + MAC / singleton summaries (`test_thin`)

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/30run_assess_singleton_mac_test_thin`. Outcome: exit 0, **8** succeeded; aggregate `succeedDuration` ~1m 32s. Publishes `.acount`, `*.singleton_mac.summary.tsv`, `*.mac_category_counts.tsv`, `*.mac_category.bar.png`, `*.mac.dist.png`, plus existing assess plot/table names under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/assess/test_thin/`. Graphviz warning only.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/30run_assess_singleton_mac_test_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/assess_plink_debug.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_thin
```

---

### 2026-05-17 — Assess debug (`tmp/assess_plink_debug.nf`) — PLINK2 `--freq counts` + MAC / singleton summaries (`test_common_thin`)

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/31run_assess_singleton_mac_test_common_thin`. Outcome: exit 0, **8** succeeded. Same artefact layout under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/assess/test_common_thin/`. Graphviz warning only.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/31run_assess_singleton_mac_test_common_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/assess_plink_debug.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_common_thin
```

---

### 2026-05-17 — Output naming refactor: assess (subgenome-only) + wheat-from-PLINK (`integrated/<plink_mod>/<wheat_mod>/`)

**Intent:** Drop redundant `job` / PLINK-`mod` / `test_*__` prefixes from artefact **basenames**; **paths** carry `job` and assess PLINK `mod`. Wheat-from-PLINK **plots/tables** publish under **`…/integrated/{plink_source_mod}/{wheat_task_mod}/`** so `test_thin` and `test_common_thin` do not collide. **`plink2_pca`** uses **`--pca approx`** (valid PLINK 2.0 syntax; avoids GRM NaN with extreme missingness). Prior `test_plink/integrated/*` and assess outputs for the two test mods were removed before regenerate.

**Assess** — cwds `…/32run_assess_naming_test_thin`, `…/33run_assess_naming_test_common_thin`; example `A.assess.maf.dist.png` under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/assess/{test_thin,test_common_thin}/`.

**Wheat PCA** — cwds `…/39run_wheat_layout_test_thin`, `…/40run_wheat_layout_test_common_thin`; example `…/integrated/test_thin/wheat_pca_tsne/plots/A.pca.png`.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/32run_assess_naming_test_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/assess_plink_debug.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_thin
```

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/33run_assess_naming_test_common_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/assess_plink_debug.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_common_thin
```

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/39run_wheat_layout_test_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/wheat_integrated_from_plink.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --user_dir /data/home/tusr1 \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_thin
```

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/40run_wheat_layout_test_common_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/wheat_integrated_from_plink.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --user_dir /data/home/tusr1 \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_common_thin
```

---

### 2026-05-20 — Wheat PCA/t-SNE from merged PLINK2 (`tmp/wheat_integrated_from_plink.nf`), `test_thin`

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/41run_wheat_pca_tsne_test_thin`. Outcome: exit 0, **8** succeeded (~12m 36s). Outputs under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/integrated/test_thin/wheat_pca_tsne/{info,plots}/` (includes `*.tsne.png` / `*.tsne.tsv` per subgenome). Graphviz warning only.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/41run_wheat_pca_tsne_test_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/wheat_integrated_from_plink.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --user_dir /data/home/tusr1 \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_thin
```

---

### 2026-05-20 — Wheat PCA/t-SNE from merged PLINK2 (`tmp/wheat_integrated_from_plink.nf`), `test_common_thin`

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/42run_wheat_pca_tsne_test_common_thin`. Outcome: exit 0, **8** succeeded (~2m 42s). Outputs under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/integrated/test_common_thin/wheat_pca_tsne/{info,plots}/`. Graphviz warning only.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/42run_wheat_pca_tsne_test_common_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/wheat_integrated_from_plink.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --user_dir /data/home/tusr1 \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_common_thin
```

---

### 2026-05-20 — Wheat PCA/t-SNE with sample `Group` coloring (`tmp/wheat_integrated_from_plink.nf`), `test_common_thin`

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/44run_wheat_pca_grp_test_common_thin`. Outcome: exit 0, **8** succeeded (~13m 43s). Tables include `Sample` + `Group`; plots use group hue. Outputs under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/integrated/test_common_thin/wheat_pca_tsne/{info,plots}/`. Graphviz warning only.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/44run_wheat_pca_grp_test_common_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/wheat_integrated_from_plink.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --user_dir /data/home/tusr1 \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_common_thin
```

---

### 2026-05-20 — Wheat PCA/t-SNE with sample `Group` coloring (`tmp/wheat_integrated_from_plink.nf`), `test_thin`

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/43run_wheat_pca_grp_test_thin`. Launched in `screen` session `wheat43` after Cursor interrupt (prior launch never created run dirs). **In progress** at log time — expect ~12–15m total; same publish layout under `…/integrated/test_thin/wheat_pca_tsne/`.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/43run_wheat_pca_grp_test_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/wheat_integrated_from_plink.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --user_dir /data/home/tusr1 \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_thin
```

---

### 2026-05-25 — GWH genome FASTA upload (`tmp/gwh_ftp_upload.nf`), Jm229

Working directory: `/data/dazheng/01projects/gwh_upload/Jm229/01run`. Outcome: exit 0, **2** succeeded (~8m 28s). Resumed partial remote copy (~7.5 GB) and uploaded remaining ~8.3 GB of `Jm229.final.fasta` (15782270833 bytes) to `ftp://submit.big.ac.cn/GWH/Batch0092978/`. Summary: `gwh_upload_summary.tsv` status **PASS**; logs under `logs/`.

```bash
cd /data/dazheng/01projects/gwh_upload/Jm229/01run && \
source ~/.bashrc && conda activate run && \
nextflow run /data/dazheng/git/script/workflow/Genetics/tmp/gwh_ftp_upload.nf \
  -c /data/dazheng/git/script/workflow/Genetics/nextflow.config \
  --output_dir /data/dazheng/01projects/gwh_upload/Jm229/01run \
  -resume
```

---

### 2026-05-25 — GWH genome MD5 verify (`tmp/gwh_ftp_upload.nf --verify_only`), Jm229

Working directory: `/data/dazheng/01projects/gwh_upload/Jm229/02verify`. Outcome: exit 0, **2** succeeded (~42m). Local and remote sizes match (15782270833 bytes); MD5 **f1bdc85c49789c21f94888045bcca818** on both sides. Summary: `gwh_md5_summary.tsv` status **PASS**; logs under `logs/`.

```bash
cd /data/dazheng/01projects/gwh_upload/Jm229/02verify && \
source ~/.bashrc && conda activate run && \
nextflow run /data/dazheng/git/script/workflow/Genetics/tmp/gwh_ftp_upload.nf \
  -c /data/dazheng/git/script/workflow/Genetics/nextflow.config \
  --output_dir /data/dazheng/01projects/gwh_upload/Jm229/02verify \
  --verify_only
```

---

### 2026-06-09 — Per-chromosome variant counts (`tmp/chr_variant_counts.nf`), `test_common_thin`

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/46run_chr_variant_counts_test_common_thin`. Outcome: exit 0, **1** succeeded (~16s). Published `test_common_thin.chr_variant_counts.tsv` and `.by_ref.tsv` under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/test_common_thin/info/` (total **1,914,487** variants across PLINK chr 0–44).

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/46run_chr_variant_counts_test_common_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/chr_variant_counts.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_common_thin \
  --process_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/test_common_thin
```

---

### 2026-06-09 — Per-chromosome variant counts (`tmp/chr_variant_counts.nf`), `test_thin` + `test_common_thin` (parallel)

Working directories: `…/47run_chr_variant_counts_test_thin`, `…/48run_chr_variant_counts_test_common_thin`. Outcome: both **exit 0**, **1** succeeded each (~19s wall, launched in parallel). Published under `…/test_plink/process/{test_thin,test_common_thin}/info/` and `…/logs/`:

| `mod` | Total variants |
| --- | --- |
| `test_thin` | **2,820,532** |
| `test_common_thin` | **1,914,487** |

```bash
# test_thin
cd /data/home/tusr1/01projects/vmap4/08stats.genome/47run_chr_variant_counts_test_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/chr_variant_counts.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_thin \
  --process_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/test_thin

# test_common_thin (parallel in a second shell)
cd /data/home/tusr1/01projects/vmap4/08stats.genome/48run_chr_variant_counts_test_common_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/chr_variant_counts.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_common_thin \
  --process_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/test_common_thin
```

---

### 2026-06-09 — `test_thin` vs `test_common_thin` chr variant compare plots (`tmp/chr_variant_compare_plots.nf`)

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/51run_chr_variant_compare_plots`. Outcome: exit 0, **1** succeeded (~16s). Line plot (variant count by ref chromosome) and bar plot (`test_common_thin / test_thin` fraction 0–1) under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/stats/thin_common_compare/{plots,info,logs}/`.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/51run_chr_variant_compare_plots && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/chr_variant_compare_plots.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink
```

---

### 2026-06-09 — `test_thin` vs `test_common_thin` genome 1 Mb density plots (refresh)

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/52run_chr_variant_compare_mb_density`. Outcome: exit 0, **1** succeeded (~22s). Line plot uses **1 Mb bins** along PLINK 0–44 vmap4 sequence order (~14,548 bins); bar plot unchanged. Outputs under `…/test_plink/stats/thin_common_compare/` (`*.variant.genome_mb_density.line.png`, `*.variant.genome_mb_density.info.tsv`, `*.variant.common_fraction.bar.png`).

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/52run_chr_variant_compare_mb_density && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/chr_variant_compare_plots.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink
```

---

### 2026-06-09 — `test_rebulld_lib` stats from migrated chr002 process tables

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/56run_test_rebulld_lib_stats`. Outcome: exit 0, **8** succeeded (~49m). Migrated legacy `rebuild/` artefacts into `test_plink/process/test_rebulld_lib/` and `stats/test_rebulld_lib/`; reran `TEST_PLINK_STATS` sample + variant QC (LD plots skipped — single chr002 ~46M variants). Outputs: `stats/test_rebulld_lib/{plots,info,logs,thresholds}/chr002.*`.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/56run_test_rebulld_lib_stats && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/test_rebulld_lib_stats.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_rebulld_lib
```

---

### 2026-06-09 — Variant MAC stats (`tmp/mac_stats_from_gcount.nf`), `test_thin` / `test_common_thin` / `test_rebulld_lib`

Working directories: `…/57run_mac_stats_test_thin`, `…/58run_mac_stats_test_common_thin`, `…/59run_mac_stats_test_rebulld_lib`. Outcome: **exit 0** — `test_thin` **4/4** (~37s), `test_common_thin` **4/4** (~19s retry), `test_rebulld_lib` **1/1** (~20m). Publishes per dataset under `stats/<mod>/{info,plots,thresholds,logs}/`: `*.variant.mac.info.tsv` (all MAC→site counts), `*.variant.mac.dist.1_100.png`, `*.variant.mac.heatmap_hobs.1_100.png`, `*.variant.mac.th.tsv`. Note: `test_common_thin` has no sites with MAC 1–100 after common filtering (min MAC ≈615); plots are zero-filled / placeholder heatmaps.

```bash
# example: test_thin
cd /data/home/tusr1/01projects/vmap4/08stats.genome/57run_mac_stats_test_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/mac_stats_from_gcount.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_thin
```

---

### 2026-06-09 — Partial rerun path migration (`partial_router.nf`)

**Note:** Genetics partial reruns no longer use `workflow/Genetics/tmp/*.nf` (removed). Use **`subworkflows/local/entry/partial_router.nf`** with **`--partial_task`** (`assess_plink`, `assess_vcf`, `ld_redraw`, `mac_stats`, `mac_dist_redraw`, `chr_counts`, `chr_compare`, `rebuild_lib_stats`, `wheat_from_plink`, `abstract_mq_50_bams`). FTP-only scripts moved to **`subworkflows/tmp/ops/`**. Historic command blocks above still reference old `tmp/` paths for audit.

**Example** (`assess_plink`, `test_thin`):

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/23run_assess_plink2_debug_stub && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/local/entry/partial_router.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --partial_task assess_plink \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --mod test_thin
```

---

### 2026-06-10 — `abstract_mq_50_bams` (chr002 validation, 50 BAMs)

Working directory: `/data/home/tusr1/01projects/vmap4/05reliable.lib/03run_abstract_mq_50_bams` (screen `mq50_chr2`).

```bash
cd /data/home/tusr1/01projects/vmap4/05reliable.lib/03run_abstract_mq_50_bams && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/local/entry/partial_router.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --partial_task abstract_mq_50_bams \
  --mod abstract_mq_50_bams \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data/home/tusr1/01projects/vmap4/05reliable.lib/03run_abstract_mq_50_bams \
  --job abstract_mq_50_bams \
  --chr 2
```

All chromosomes (omit `--chr`; use a new numbered run folder):

```bash
cd /data/home/tusr1/01projects/vmap4/05reliable.lib/<NNrun_abstract_mq_50_bams_allchr> && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/local/entry/partial_router.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --partial_task abstract_mq_50_bams \
  --mod abstract_mq_50_bams \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data/home/tusr1/01projects/vmap4/05reliable.lib/<NNrun_abstract_mq_50_bams_allchr> \
  --job abstract_mq_50_bams
```

Output: `process/abstract_mq_50_bams/reference/chrNNN.site_mq.ref.txt.gz` (full-chr grid: chrom, pos, MQ or `.`) plus `*.site_mq.ref.info.tsv`.

---

### 2026-06-10 — SUPP — `abstract_mq_50_bams` full-chromosome MQ reference (not variant targets)

**Change:** `calc_site_mq_bcftools` scans `-r CHR:1-LEN` (free mpileup), pads to PopDep-like length, publishes **`process/abstract_mq_50_bams/reference/chrNNN.site_mq.ref.txt.gz`**. Does **not** use legacy variant list `2_1_122798052.pos.txt`. Params: `mq_pad_all_positions`, `mq_ref_gzip`, `mq_keep_sparse`.

**Note:** Re-run chr002+ after stopping the prior `03run_abstract_mq_50_bams` screen job (old output name `variant/*.site_mq.txt`).

---

### 2026-06-10 — `abstract_mq_50_bams` all chromosomes (50 BAMs, full-chr reference)

Working directory: `/data/home/tusr1/01projects/vmap4/05reliable.lib/04run_abstract_mq_50_bams_ref` (screen `mq50_allchr`, 45 tasks, `maxForks=4`).

```bash
cd /data/home/tusr1/01projects/vmap4/05reliable.lib/04run_abstract_mq_50_bams_ref && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/local/entry/partial_router.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --partial_task abstract_mq_50_bams \
  --mod abstract_mq_50_bams \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data/home/tusr1/01projects/vmap4/05reliable.lib/04run_abstract_mq_50_bams_ref \
  --job abstract_mq_50_bams \
  -process.withLabel:site_mq_bcftools.maxForks=4
```

Log: `run_logs/nextflow.mq50_allchr.log`. Outputs: `process/abstract_mq_50_bams/reference/chr{000..044}.site_mq.ref.txt.gz`.

---

### 2026-06-10 — `abstract_mq_50_bams` all chromosomes (rerun after code update)

Working directory: `/data/home/tusr1/01projects/vmap4/05reliable.lib/06run_abstract_mq_50_bams_ref` (screen `mq50_allchr`, 45 tasks, `maxForks=4` in `resources.config`).

Stopped prior attempt: `04run_abstract_mq_50_bams_ref` (old code: single `*.site_mq.txt`, `sort -u`).

```bash
cd /data/home/tusr1/01projects/vmap4/05reliable.lib/06run_abstract_mq_50_bams_ref && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/local/entry/partial_router.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --partial_task abstract_mq_50_bams \
  --mod abstract_mq_50_bams \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data/home/tusr1/01projects/vmap4/05reliable.lib/06run_abstract_mq_50_bams_ref \
  --job abstract_mq_50_bams
```

Log: `run_logs/nextflow.mq50_allchr.log`. Outputs: `process/abstract_mq_50_bams/reference/chrNNN.site_mq.calls.tsv.gz` (CHROM POS REF ALT MQ) and `chrNNN.site_mq.ref.txt.gz` (full-chr grid).

---

### 2026-06-10 — `abstract_mq_50_bams` all chromosomes (I16 mpileup-only, float MQ)

Stopped prior attempt: `06run_abstract_mq_50_bams_ref` (screen `mq50_allchr`, bcftools call + giant `*.calls.vcf`).

Working directory: `/data/home/tusr1/01projects/vmap4/05reliable.lib/08run_abstract_mq_50_bams_i16` (screen `mq50_i16_allchr`, 45 tasks, `maxForks=4`).

```bash
cd /data/home/tusr1/01projects/vmap4/05reliable.lib/08run_abstract_mq_50_bams_i16 && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/local/entry/partial_router.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --partial_task abstract_mq_50_bams \
  --mod abstract_mq_50_bams \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data/home/tusr1/01projects/vmap4/05reliable.lib/08run_abstract_mq_50_bams_i16 \
  --job abstract_mq_50_bams
```

Log: `run_logs/nextflow.mq50_i16_allchr.log`. MQ from mpileup INFO/I16 mean MAPQ (float, no bcftools call). Outputs: `chrNNN.site_mq.calls.tsv.gz` + `chrNNN.site_mq.ref.txt.gz`.

---

### 2026-06-10 — `abstract_mq_50_bams` resume (I16, 16-way I/O-bound tuning)

Stop `08run` after the four large in-flight chromosomes finish (chr000, chr029, chr033, chr041); chr002 and chr018 are already done. Do **not** delete `08run/work/`.

Resource change (`conf/resources.config`, label `site_mq_bcftools`): `cpus=4`, `memory=64.GB`, `maxForks=16` (was 32 / 256.GB / 4).

**Resume note:** Launch cwd must be **`08run`** (not `09run`) for `-resume`; `-work-dir` alone does not restore cache. Also skip already-published chr via `--mq_chr_exclude` because `cpus` 32→4 changes task hash and would rerun completed chromosomes.

Working directory (launch cwd): `/data/home/tusr1/01projects/vmap4/05reliable.lib/08run_abstract_mq_50_bams_i16` (screen `mq50_i16_resume`). Log tee: `09run_.../run_logs/nextflow.mq50_i16_resume.log`. Publish tree: `08run/abstract_mq_50_bams/process/abstract_mq_50_bams/reference/`.

```bash
source ~/.bashrc && conda activate run && \
screen -dmS mq50_i16_resume bash -c 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate run && /data/home/tusr1/01projects/vmap4/05reliable.lib/09run_abstract_mq_50_bams_i16/resume_from_08.sh'
```

Skip list: `--mq_chr_exclude '0,2,18,24,29,41,44'` (38 remaining tasks). Log: `09run_.../run_logs/nextflow.mq50_i16_resume.log`.

---

### 2026-06-14 — `abstract_mq_50_bams` frozen publish under test_plink/process

All 45 chr published to **`/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/abstract_mq_50_bams/{reference,logs}/`**. Workflow now uses **`params.mq_dir`** (not `output_dir/job/process/mod`). Re-runs skip existing `*.site_mq.ref.*` when **`mq_force_rerun=false`**. Use **`--job test_plink`** for operator consistency; main pipeline does not invoke this partial task.

```bash
nextflow run .../partial_router.nf -c .../nextflow.config \
  --partial_task abstract_mq_50_bams --mod abstract_mq_50_bams --job test_plink \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 --src_dir /data/home/tusr1/git/script/src \
  -preview
```

Python: `from genetics.genomics.variant.mq import site_mq_ref_path, site_mq_calls_path, load_mq_data`

---

### 2026-06-11 — MAC vs missing regression (`mac_stats`), `test_rebulld_lib` (background)

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/67run_mac_miss_reg_test_rebulld_lib`. Screen session: `mac_miss_rebulld`. Reruns `variant_mac_stats` + `variant_mac_missing_reg` on chr002 gcount/vmiss (includes full-range and MAC 0–100 / 0–500 / 0–1000 regression plots). Publishes under `test_plink/stats/test_rebulld_lib/{info,plots,thresholds,logs}/`.

```bash
source ~/.bashrc && conda activate run && \
screen -dmS mac_miss_rebulld bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/08stats.genome/67run_mac_miss_reg_test_rebulld_lib && nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/local/entry/partial_router.nf -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config --partial_task mac_stats --home_dir /data/home/tusr1/01projects/vmap4 --user_dir /data/home/tusr1 --src_dir /data/home/tusr1/git/script/src --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 --job test_plink --mod test_rebulld_lib 2>&1 | tee run_logs/nextflow.mac_miss_reg_test_rebulld_lib.log'
```

Log: `67run_mac_miss_reg_test_rebulld_lib/run_logs/nextflow.mac_miss_reg_test_rebulld_lib.log`. Attach: `source ~/.bashrc && conda activate run && screen -r mac_miss_rebulld`.

---

### 2026-06-11 — MAC vs MAF regression (`mac_stats` + `variant_mac_maf_reg`)

Working directories: `70run_mac_maf_reg_test_thin` (**12/12**, ~2m), `71run_mac_maf_reg_test_common_thin` (**12/12**, ~1m), `72run_mac_maf_reg_test_rebulld_lib` (screen `mac_maf_rebulld`, background). Publishes `{id}.variant.mac_maf.reg.png` plus `.reg.mac0_{100,500,1000}.png` under `stats/<mod>/plots/` (MAF from `.gcount` only).

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/70run_mac_maf_reg_test_thin && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/local/entry/partial_router.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --partial_task mac_stats --mod test_thin --job test_plink \
  --home_dir /data/home/tusr1/01projects/vmap4 --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1
```

`test_rebulld_lib` (background):

```bash
source ~/.bashrc && conda activate run && \
screen -dmS mac_maf_rebulld bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/08stats.genome/72run_mac_maf_reg_test_rebulld_lib && nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/local/entry/partial_router.nf -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config --partial_task mac_stats --home_dir /data/home/tusr1/01projects/vmap4 --user_dir /data/home/tusr1 --src_dir /data/home/tusr1/git/script/src --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 --job test_plink --mod test_rebulld_lib 2>&1 | tee run_logs/nextflow.log'
```

---

### 2026-06-11 — MAC miss bin50 refactor (full MAC only + mac_an full in bin50 task)

Working directory: `81run_mac_miss_bin50s_v2_test_thin` (**4/4**, ~1m35s). Bin50 outputs use full MAC range only (no `mac0_100/500/1000` subsets).

**Note (later same day):** Removed duplicate `mac_miss_full` partial task and `reg.full.*` MAC vs F_MISS plots; canonical full-data regression is `variant_mac_missing_reg` → `reg.png` / `reg.mac0_{100,500,1000}.png` via `mac_stats` or main `TEST_PLINK_STATS`.

---

### 2026-06-11 — MAC miss naming cleanup + R_mac full-data plots

Working directories: `82run_mac_miss_reg_rmac_test_thin` (`mac_stats`, **12/12**, ~3m15s), `83run_mac_miss_bin50s_v3_test_thin` (`mac_miss_bin50_sample`, **4/4**, ~2m15s).

**Full-data** (`variant_mac_missing_reg` / `mac_stats`): `reg.png`, `reg.mac0_{100,500,1000}.png`, `reg.mac_an.png`, `reg.mac_an.mac0_{100,500,1000}.png`, `R_mac.png`, `R_mac.mac0_{100,500,1000}.png`.

**Bin50 sample** (`mac_miss_bin50_sample`): `reg.bin50s.png`, `reg.mac_an.bin50s.png`, `R_mac.bin50s.png` only (no `full`/`all` suffixes).

---

### 2026-06-15 — `main_raw_popdepth` chr36 smoke (`10stats.genome`)

Working directory: `/data/home/tusr1/01projects/vmap4/10stats.genome/01run_main_raw_popdepth` (screen `popdep_chr36`). New stats module root: **`10stats.genome/`** (numbered runs from `01`; `08stats.genome` kept for history). Process label **`popdep_tiger`**: 16 CPUs, 64 GB, `maxForks=1`. Publish: `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/main_raw/variant/chr036.popdep.txt`.

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdep_chr36 bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/01run_main_raw_popdepth && nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/local/entry/partial_router.nf -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config --partial_task main_raw_popdepth --mod main_raw_popdepth --job test_plink --chr 36 --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink --home_dir /data/home/tusr1/01projects/vmap4 --user_dir /data/home/tusr1 --src_dir /data/home/tusr1/git/script/src 2>&1 | tee run_logs/nextflow.popdep_chr36.log'
```

**Follow-up (full 44 chr, after chr36 OK):** same cwd, drop `--chr`, add `-resume`, new screen e.g. `popdep_all`:

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdep_all bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/01run_main_raw_popdepth && nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/local/entry/partial_router.nf -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config --partial_task main_raw_popdepth --mod main_raw_popdepth --job test_plink --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink --home_dir /data/home/tusr1/01projects/vmap4 --user_dir /data/home/tusr1 --src_dir /data/home/tusr1/git/script/src -resume 2>&1 | tee run_logs/nextflow.popdep_all.log'
```

Existing `*.popdep.txt` under `main_raw/variant/` are skipped (`popdep_force_rerun=false`).

---

### 2026-06-15 — `PopDepFull` chr36 benchmark (`10stats.genome`)

Working directory: `/data/home/tusr1/01projects/vmap4/10stats.genome/02run_popdepfull_chr36_benchmark` (screen `popdepfull_chr36`). Jar `TIGER_PD_20260615.jar`, app `PopDepFull`, `-e 32`, `Xmx128G`, `maxForks=1`. Benchmark publish: `02run_popdepfull_chr36_benchmark/publish/variant/` (production `main_raw/variant/` unchanged).

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdepfull_chr36 bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/02run_popdepfull_chr36_benchmark && bash run_logs/monitor_io.sh run_logs/io_monitor.tsv & nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/local/entry/partial_router.nf -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config -c /data/home/tusr1/01projects/vmap4/10stats.genome/02run_popdepfull_chr36_benchmark/popdepfull_benchmark.config --partial_task main_raw_popdepth --mod main_raw_popdepth --job test_plink --chr 36 --popdep_tiger_jar TIGER_PD_20260615.jar --popdep_tiger_app PopDepFull --popdep_force_rerun true --popdep_publish_dir /data/home/tusr1/01projects/vmap4/10stats.genome/02run_popdepfull_chr36_benchmark/publish --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink --home_dir /data/home/tusr1/01projects/vmap4 --user_dir /data/home/tusr1 --src_dir /data/home/tusr1/git/script/src 2>&1 | tee run_logs/nextflow.popdepfull_chr36.log'
```

**Production full sweep (after benchmark OK):** `05run_main_raw_popdepth_full`, drop `--chr` and `--popdep_publish_dir`, keep PopDepFull params:

```bash
--popdep_tiger_jar TIGER_PD_20260615.jar
--popdep_tiger_app PopDepFull
--popdep_force_rerun false
```

### 2026-06-16 — `main_raw_popdepth` full PopDepFull sweep (`10stats.genome/06run`)

Working directory: `/data/home/tusr1/01projects/vmap4/10stats.genome/06run_main_raw_popdepth_full` (screen `popdepfull_all`). PopDepFull, `maxForks=2`, 32 CPU / 128 GB per task. Publish BGZF+tabix grids to `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/main_raw/variant/chrNNN.popdep.txt.bgz` (+ `.tbi`); logs under `.../main_raw/logs/`. Skips chr32 (no `chr032.vcf.gz`).

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdepfull_all bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/06run_main_raw_popdepth_full && date -Iseconds | tee run_logs/start.txt && nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/local/entry/partial_router.nf -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config --partial_task main_raw_popdepth --mod main_raw_popdepth --job test_plink --popdep_tiger_jar TIGER_PD_20260615.jar --popdep_tiger_app PopDepFull --popdep_tiger_max_forks 2 --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink --home_dir /data/home/tusr1/01projects/vmap4 --user_dir /data/home/tusr1 --src_dir /data/home/tusr1/git/script/src 2>&1 | tee run_logs/nextflow.log; echo $? | tee run_logs/exit_code.txt; date -Iseconds | tee run_logs/end.txt'
```

### 2026-06-16 — `main_raw_popdepth` PopDepCrossChr production (`10stats.genome/07run`)

Working directory: `/data/home/tusr1/01projects/vmap4/10stats.genome/07run_main_raw_popdepth_crosschr` (screen `popdep_crosschr`). Single JVM task: `TIGER_PD_20260616.jar`, app `PopDepCrossChr`, 64 CPU / 256 GB (`popdep_tiger_crosschr` label). One BAM pass per taxon for all chr; publish BGZF+tabix to `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/main_raw/variant/chrNNN.popdep.txt.bgz` (+ `.tbi`); logs under `.../main_raw/logs/`. Skips chr32 (no `chr032.vcf.gz`). Replaces stopped 06run PopDepFull sweep.

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdep_crosschr bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/07run_main_raw_popdepth_crosschr && date -Iseconds | tee run_logs/start.txt && nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/local/entry/partial_router.nf -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config --partial_task main_raw_popdepth --mod main_raw_popdepth --job test_plink --popdep_tiger_jar TIGER_PD_20260616.jar --popdep_tiger_app PopDepCrossChr --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink --home_dir /data/home/tusr1/01projects/vmap4 --user_dir /data/home/tusr1 --src_dir /data/home/tusr1/git/script/src 2>&1 | tee run_logs/nextflow.log; echo $? | tee run_logs/exit_code.txt; date -Iseconds | tee run_logs/end.txt'
```

### 2026-06-16 — `main_raw_popdepth` PopDepCrossChr v2 (`10stats.genome/08run`)

Working directory: `/data/home/tusr1/01projects/vmap4/10stats.genome/08run_main_raw_popdepth_crosschr` (screen `popdep_crosschr_v2`). Updated `TIGER_PD_20260616.jar` (segmentIndex long fix + length file `Chr\\tLength\\tnTaxa`). Single task: `tb.ALL.txt`, 3-column length (nTaxa 7675/7594/7738/8285 by A/B/D/Others), **16 threads / 384 GB**. Publish to `main_raw/variant/chrNNN.popdep.txt.bgz` (+ `.tbi`); skip chr32.

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdep_crosschr_v2 bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/08run_main_raw_popdepth_crosschr && date -Iseconds | tee run_logs/start.txt && nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/local/entry/partial_router.nf -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config --partial_task main_raw_popdepth --mod main_raw_popdepth --job test_plink --popdep_tiger_jar TIGER_PD_20260616.jar --popdep_tiger_app PopDepCrossChr --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink --home_dir /data/home/tusr1/01projects/vmap4 --user_dir /data/home/tusr1 --src_dir /data/home/tusr1/git/script/src 2>&1 | tee run_logs/nextflow.log; echo $? | tee run_logs/exit_code.txt; date -Iseconds | tee run_logs/end.txt'
```

### 2026-06-16 — PopDepCrossChr thread benchmark (`10stats.genome/09run`)

Working directory: `/data/home/tusr1/01projects/vmap4/10stats.genome/09run_popdep_crosschr_threads_bench` (screen `popdep_bench9`). Sequential **60 min** runs for threads **8, 16, 24, 32, 48, 64, 96, 128** (`384G`, `TIGER_PD_20260616.jar`, skip chr32). Bench publish under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/_bench_popdep_crosschr_threads/e{threads}/`. Summary: `run_logs/bench_summary.tsv`.

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdep_bench9 bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/09run_popdep_crosschr_threads_bench && bash run_all_bench.sh 2>&1 | tee run_logs/run_all.log'
```

### 2026-06-17 — PopDepFull vs PopDepCrossChr head-to-head (`10stats.genome/10run`)

Working directory: `/data/home/tusr1/01projects/vmap4/10stats.genome/10run_popdep_head2head`. Production **ABD** BAMs via taxa panels `bench_assets/tb.N3_micro.txt` (3 taxa) and `tb.N20.txt` (20 taxa); tiers T1 (seek, CrossChr `-e 1` vs Full `maxForks=1`), T2 (throughput, `-e 8` vs `maxForks=2`), optional T3 chr36-only (biased). Publish under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/_bench_popdep_head2head/`. Summary: `run_logs/head2head_summary.tsv`. New params: `--popdep_taxa_bam_file`, `--popdep_chr_list`.

```bash
source ~/.bashrc && conda activate run && \
cd /data/home/tusr1/01projects/vmap4/10stats.genome/10run_popdep_head2head && \
bash run_all_head2head.sh 2>&1 | tee run_logs/run_all.log
```

Optional: `HEAD2HEAD_TIERS=1,2,3` to include chr36-only biased tier; `HEAD2HEAD_JAR=TIGER_PD_YYYYMMDD.jar` when testing optimized builds.

### 2026-06-17 — PopDep truncated extract + head2head (`10stats.genome/10run`)

Extract cwd: `.../10run_popdep_head2head/bench_assets` — `N_TAXA=20 REGION_BP=10000000 bash extract_popdep_bench_bams.sh` (~112 min, 13G under `00data/popdep_bench/chr10000000/`).

Bench cwd: `.../10run_popdep_head2head` (screen `popdep_h2h_trunc`). Truncated T1+T2: `N_TAXA=20 HEAD2HEAD_TIERS=1,2 bash run_truncated_head2head.sh`.

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdep_h2h_trunc bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/10run_popdep_head2head && N_TAXA=20 HEAD2HEAD_TIERS=1,2 bash run_truncated_head2head.sh 2>&1 | tee run_logs/run_truncated_all.log'
```

### 2026-06-17 — PopDep truncated head2head 11run (bgzip/tabix split)

Working directory: `/data/home/tusr1/01projects/vmap4/10stats.genome/11run_popdep_head2head_bgz_split`. Re-run of truncated **10 Mb × 44 chr × 20 taxa** T1+T2 after `popdep_tiger_gz_to_bgzip_tabix` was split from TIGER compute. **10run results retained** for comparison (`10run_popdep_head2head/run_logs/head2head_summary.tsv`). Publish under `.../_bench_popdep_head2head/11run_bgz_split/trunc_10000000_N20/` (distinct from 10run `.../trunc_10000000_N20/`). Shared BAMs: `00data/popdep_bench/chr10000000/`.

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdep_h2h_11run bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/11run_popdep_head2head_bgz_split && N_TAXA=20 HEAD2HEAD_TIERS=1,2 bash run_truncated_head2head.sh 2>&1 | tee run_logs/run_truncated_all.log'
```

### 2026-06-17 — PopDep CrossChr flatMap fix + 11run bench retry

**Fix:** `partial_main_raw_popdepth.nf` — `calc_population_depth_crosschr.out.tiger_gz` is a file list; fan-out to `popdep_tiger_gz_to_bgzip_tabix` via `.flatMap` (was `.map` + `getName` on `ArrayList`).

**Smoke (pass):** `-resume` on `11run.../t1_crosschr_micro` — TIGER cached, 44× `popdep_bgz_tabix` exit 0; 44 `*.popdep.txt.bgz` under publish `.../t1_crosschr_micro/variant/`.

screen -dmS popdep_h2h_11run bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/11run_popdep_head2head_bgz_split && N_TAXA=20 HEAD2HEAD_TIERS=1,2 bash run_truncated_head2head.sh 2>&1 | tee run_logs/run_truncated_all.log'
```

### 2026-06-17 — 11run T2-only monitor补跑

Cwd: `10stats.genome/11run_popdep_head2head_bgz_split`. New `monitor_head2head.sh` (wait for TIGER + full NF lifecycle). **T2 only** (`HEAD2HEAD_TIERS=2`, N20 truncated); prior timing rows in `run_logs/head2head_summary_timing_only.tsv`. Screen `popdep_h2h_t2mon`. Use `HEAD2HEAD_BUST_CACHE=1` so Nextflow re-executes TIGER (not resume-only).

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdep_h2h_t2mon bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/11run_popdep_head2head_bgz_split && N_TAXA=20 HEAD2HEAD_TIERS=2 HEAD2HEAD_BUST_CACHE=1 bash run_truncated_head2head.sh 2>&1 | tee run_logs/run_t2_monitor.log'
```

### 2026-06-18 — PopDep bench N200 truncated BAM extract (32 parallel)

Extract cwd / log: `10stats.genome/12run_popdep_extract_n200`. Data: `00data/popdep_bench/chr10000000_N200/` (44 chr × 10 Mb, 200 ABD taxa). Screen `popdep_extract_n200`. Script: `10run_popdep_head2head/bench_assets/extract_popdep_bench_bams.sh`. **Retry:** `SAMTOOLS_THREADS=32 EXTRACT_JOBS=4` (4 parallel taxa × samtools `view`/`index -@32` ≈ 128 threads).

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdep_extract_n200 bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/12run_popdep_extract_n200 && N_TAXA=200 REGION_BP=10000000 EXTRACT_JOBS=4 SAMTOOLS_THREADS=32 bash /data/home/tusr1/01projects/vmap4/10stats.genome/10run_popdep_head2head/bench_assets/extract_popdep_bench_bams.sh 2>&1 | tee -a run_logs/extract_n200.log'
```

### 2026-06-18 — N200 PopDep bench harness + `popdep_tiger_threads` smoke

Repo: `params.popdep_tiger_threads`, `params.popdep_tiger_memory_gb`; `resources.config` wires them to `popdep_tiger` cpus/memory; `partial_main_raw_popdepth.nf` logs resolved PopDepFull settings.

**Param smoke (pass):** cwd `10stats.genome/13run_popdep_n200_bench` — `bash verify_tiger_threads_smoke.sh` (chr36, N200 `tb.N200.txt`). Confirmed Nextflow log `PopDepFull: threads=32|64 … maxForks=2|1` and `work/…/.command.sh` TIGER `-e 32` / `-e 64`.

**Full bench (sequential, not yet launched):** same cwd, screen recommended.

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdep_n200_bench bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/13run_popdep_n200_bench && bash run_n200_bench.sh 2>&1 | tee run_logs/run_all.log'
```

Phases: T1 PopDepFull `fork2×32` then `fork1×64`; T2 PopDepCrossChr threads 8→128 (128G, run to completion). Publish `_bench_popdep_n200/trunc_10000000_N200/`. Summaries: `run_logs/full_fork_summary.tsv`, `run_logs/crosschr_threads_summary.tsv`.

### 2026-06-18 — `test_thin` vs `test_common_thin` chr compare plots (centromere shading)

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/94run_chr_compare_centromere`. Outcome: exit 0, **1** succeeded (~1m 7s). Refreshed genome density line plot with vmap4 centromere bands under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/stats/thin_common_compare/plots/test_thin_vs_test_common_thin.variant.genome_mb_density.line.png`.

```bash
cd /data/home/tusr1/01projects/vmap4/08stats.genome/94run_chr_compare_centromere && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/local/entry/partial_router.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --partial_task chr_compare
```

  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --partial_task chr_compare
```

### 2026-06-19 — N200 PopDep bench preflight (14run, checkpoint ci=10)

Working directory: `/data/home/tusr1/01projects/vmap4/10stats.genome/14run_popdep_n200_bench`. Outcome: **skipped remaining smoke** after `checkpoint.bin` appeared under `_smoke_preflight/crosschr_ckpt/checkpoint` with `-ci 10` (TIGER log: `Checkpoint written after 10 / 200 taxa`). Jar `TIGER_PD_20260619.jar`; taxa map `tb.N200.taxaBamMap.txt`.

### 2026-06-19 — N200 PopDep bench full run (14run, TIGER_PD_20260619)

Working directory: `/data/home/tusr1/01projects/vmap4/10stats.genome/14run_popdep_n200_bench`. Screen session `popdep_n200_14run`. Phases: PopDepFull `fork2×32` / `fork1×64`; CrossChr e8–e128; extras `crosschr_e64_shuffle`, `crosschr_e64_ckpt_ci100`. Publish `_bench_popdep_n200_14run/trunc_10000000_N200/`.

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdep_n200_14run bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/14run_popdep_n200_bench && bash run_n200_bench.sh 2>&1 | tee -a run_logs/run_all.log'
```

### 2026-06-19 — `main_raw_popdepth` PopDepCrossChr production (`10stats.genome/15run`)

Working directory: `/data/home/tusr1/01projects/vmap4/10stats.genome/15run_main_raw_popdepth_crosschr`. Screen `popdep_crosschr_15run`. Full **main_raw** (~8285 taxa, skip chr32): **`TIGER_PD_20260619.jar`**, **96 threads / 640 GB**, **`-o shuffle`**, checkpoint **`-ci 500`** under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/main_raw/checkpoint/popdep_crosschr`. Monitor `run_logs/monitor.tsv` every **600 s** (sdb). Publish `main_raw/variant/chrNNN.popdep.txt.bgz`.

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdep_crosschr_15run bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/15run_main_raw_popdepth_crosschr && bash run_popdep_crosschr.sh 2>&1 | tee run_logs/run_all.log'
```

---

### 2026-06-21 — CS_mp_2018_8X FastCall3 `scan_only` debug (chr2)

Working directory: `/data/home/tusr1/01projects/vmap4/04runScreens/07run_cs_mp2018_scan_debug`. Outcome: **pass** (~4 min, exit 0). Single sample **CS_mp_2018_8X**; taxaBamMap line extracted from canonical `ABD.taxaBamMap.txt` → `chr1.A.taxaBamMap.txt`. vLib read-only from `04runScreens/all/vLib`. chr2 VCF **99.97%** `./.` (`DP=0`) — root cause traced to mpileup excluding non–proper-pair reads (see progress log 2026-06-22).

```bash
RUN=/data/home/tusr1/01projects/vmap4/04runScreens/07run_cs_mp2018_scan_debug
cd "$RUN"
source ~/.bashrc && conda activate run
nextflow run /data/home/tusr1/git/script/workflow/Genetics/modules/local/genotype/calling/caller.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir "$RUN" \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir "$RUN" \
  --job chr1 \
  --mod scan_only \
  --chr 2 \
  --server s243 \
  --vlib_dir /data/home/tusr1/01projects/vmap4/04runScreens/all/vLib \
  --gen_dir "$RUN/chr1/gen" \
  --tiger_jar TIGER_F3_20251121.jar
```

---

### 2026-06-22 — CS_mp_2018_8X BAM proper-pair relabel (scheme B, non-Nextflow)

Working directory: `/data/home/tusr1/01projects/vmap4/04runScreens/08run_cs_mp2018_properflag`. Input: `00data/02bam/bam1/ABD/CS_mp_2018_8X.bam` (~72 GB). Output: `CS_mp_2018_8X.properflag.bam` (~66 GB). **64 threads**; launched in `screen -S cs_mp_properflag`. Outcome: **pass** (~1 h 25 min). Properly paired **78.49%** (was ~0.0003%). Logs: `logs/properflag_*.log`, `logs/flagstat_properflag_*.txt`.

```bash
source ~/.bashrc && conda activate run
screen -dmS cs_mp_properflag bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/04runScreens/08run_cs_mp2018_properflag && THREADS=64 bash run_properflag.sh 2>&1 | tee logs/full_properflag_launch.log'
```

chr2 pilot (subset) before full BAM:

```bash
cd /data/home/tusr1/01projects/vmap4/04runScreens/08run_cs_mp2018_properflag
source ~/.bashrc && conda activate run
samtools view -b -@ 64 /data/home/tusr1/01projects/vmap4/00data/02bam/bam1/ABD/CS_mp_2018_8X.bam 2 -o pilot_chr2_input.bam
samtools index -@ 64 pilot_chr2_input.bam
THREADS=64 bash run_properflag.sh pilot_chr2_input.bam pilot_chr2_properflag.bam
```

---

### 2026-06-22 — CS_mp_2018_8X same-chr TLEN QC (non-Nextflow, `stats` env)

Working directory: `/data/home/tusr1/01projects/vmap4/04runScreens/08run_cs_mp2018_properflag/tlen_qc`. BAM: `CS_mp_2018_8X.properflag.bam`. **64 threads** for `samtools view`. Key outputs: `same_chr_tlen_0-5kb_bin20bp.tsv/png`, stacked proper/improper plot, GMM deconv plot, FR vs RF overlay.

```bash
source ~/.bashrc && conda activate stats
RUN=/data/home/tusr1/01projects/vmap4/04runScreens/08run_cs_mp2018_properflag
python3 "$RUN/plot_same_chr_tlen_fine.py" "$RUN/CS_mp_2018_8X.properflag.bam" \
  --bin-size 20 --max-tlen 5000 --threads 64 \
  --out-tsv "$RUN/tlen_qc/same_chr_tlen_0-5kb_bin20bp.tsv" \
  --out-png "$RUN/tlen_qc/same_chr_tlen_0-5kb_bin20bp.png"
python3 "$RUN/plot_same_chr_tlen_stacked.py" "$RUN/CS_mp_2018_8X.properflag.bam" \
  --bin-size 20 --max-tlen 5000 --threads 64 \
  --out-tsv "$RUN/tlen_qc/same_chr_tlen_0-5kb_bin20bp_stacked.tsv" \
  --out-png "$RUN/tlen_qc/same_chr_tlen_0-5kb_bin20bp_stacked.png"
python3 "$RUN/deconv_tlen_gmm.py" "$RUN/tlen_qc/same_chr_tlen_0-5kb_bin20bp.tsv" \
  --out-dir "$RUN/tlen_qc"
```

---

### 2026-06-23 — Sample cov hybrid heatmaps and plot regeneration (Python, `stats` env)

Working directory: repo root not required; reads published info TSVs under data1. Outcome: **pass** — hybrid log10(Depth) heatmaps for A/B/D/Others; standard heatmaps and GAM panels regenerated for `test_thin` and `test_common_thin`. Registry: `doc/project_knowledge/domain/sample_cov_stats.yaml`.

```bash
source ~/.bashrc && conda activate stats
export PYTHONPATH=/data/home/tusr1/git/script/src/python
BASE=/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/stats

# Hybrid: F_MISS from test_thin, IBS_Ref from test_common_thin, shared Coverage
python <<'PY'
from pathlib import Path
from genetics.genomics.sample.cov import heatmap_thin_miss_common_ibs
base = Path("/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/stats")
for mod in ["A", "B", "D", "Others"]:
    thin = base / "test_thin/info" / f"{mod}.sample.cov.ibs_depth_miss.info.tsv"
    common = base / "test_common_thin/info" / f"{mod}.sample.cov.ibs_depth_miss.info.tsv"
    prefix = base / "test_thin/plots" / f"{mod}.sample.cov"
    heatmap_thin_miss_common_ibs(str(thin), str(common), output_prefix=str(prefix))
PY

# Full regen: standard heatmaps + LOESS + GAM + ABD panels (both datasets)
python <<'PY'
from pathlib import Path
import shutil
from genetics.genomics.sample.cov import heatmap_ibs_depth_missing, compare_subgenome_ibs_depth_gam
base_root = Path("/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/stats")
for job in ["test_thin", "test_common_thin"]:
    base = base_root / job
    for mod in ["A", "B", "D", "Others"]:
        info = base / "info" / f"{mod}.sample.cov.ibs_depth_miss.info.tsv"
        heatmap_ibs_depth_missing(tsv_file=str(info), output_prefix=str(base / "plots" / f"{mod}.sample.cov"), save_info=False)
    abd = [str(base / "info" / f"{m}.sample.cov.ibs_depth_miss.info.tsv") for m in ["A", "B", "D"]]
    compare_subgenome_ibs_depth_gam(abd, output_prefix=str(base / "plots" / "ABD.sample.cov"))
    summary = base / "plots" / "ABD.sample.cov.gam_subgenome_summary.info.tsv"
    if summary.exists():
        shutil.move(summary, base / "info" / summary.name)
PY
```

Note: the 2026-06-23 entry above referenced `sample_cov_stats.yaml`; that file was removed — publish layout is **`doc/project_knowledge/domain/data_publish_tree.yaml`** only.

---

### 2026-06-23 — `sg_cov_replot.nf` (Watkins mosdepth → subgenome depth refresh)

Republish **`*.sample.sg_cov.*`** stats after copying Watkins mosdepth into `00data/04depth/09Watkins` and rebuilding `sg_mean_depth` taxaBamMap files. Entry: **`workflow/Genetics/subworkflows/tmp/sg_cov_replot.nf`**.

**test_thin** — working directory: `/data/home/tusr1/01projects/vmap4/10stats.genome/17run_sg_cov_replot_watkins`. Outcome: exit 0, **5** succeeded (4× `sample_sg_coverage_stats` + 1× `plot_subgenome_gam_ibs_depth_compare_sg`), ~6m 47s.

```bash
cd /data/home/tusr1/01projects/vmap4/10stats.genome/17run_sg_cov_replot_watkins && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/tmp/sg_cov_replot.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --process_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/test_thin \
  --mod test_thin \
  --job test_plink
```

**test_common_thin** — working directory: `/data/home/tusr1/01projects/vmap4/10stats.genome/18run_sg_cov_replot_watkins_common`. Outcome: exit 0, **5** succeeded, ~8m 57s.

```bash
cd /data/home/tusr1/01projects/vmap4/10stats.genome/18run_sg_cov_replot_watkins_common && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/tmp/sg_cov_replot.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --process_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/test_common_thin \
  --mod test_common_thin \
  --job test_plink
```

Published under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/stats/{test_thin,test_common_thin}/{info,plots,thresholds,logs}/` with prefix `{A,B,D,Others}.sample.sg_cov` and `ABD.sample.sg_cov.gam_*`.

---

### 2026-06-23 — `sg_cov_hybrid_thinmiss_commonibs.nf` (hybrid thinmiss + common IBS heatmaps)

Hybrid heatmap for **subgenome mosdepth depth** (`sg_cov`): F_MISS from `test_thin`, IBS from `test_common_thin`, Y-axis log10(subgenome depth). Entry: **`workflow/Genetics/subworkflows/tmp/sg_cov_hybrid_thinmiss_commonibs.nf`**.

Working directory: `/data/home/tusr1/01projects/vmap4/10stats.genome/20run_sg_cov_hybrid_thinmiss`. Outcome: exit 0, **4** succeeded (A/B/D/Others), ~29s.

```bash
cd /data/home/tusr1/01projects/vmap4/10stats.genome/20run_sg_cov_hybrid_thinmiss && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/tmp/sg_cov_hybrid_thinmiss_commonibs.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --thin_stats_mod test_thin \
  --common_stats_mod test_common_thin
```

Published under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/stats/test_thin/`:

- `{A,B,D,Others}.sample.sg_cov.heatmap_thinmiss_commonibs.logdepth.png`
- `{A,B,D,Others}.sample.sg_cov.hybrid_thinmiss_commonibs.info.tsv`
- matching logs under `logs/`

---

### 2026-06-23 — `sg_cov_gam_residual.nf` (GAM residual outlier diagnostics)

GAM te residual analysis for **sg_cov**: residual = observed F_MISS − predicted F_MISS; flag top `gam_residual_outlier_frac` (default 2%) by |residual|; group bar charts + highlighted scatter plots. Entry: **`workflow/Genetics/subworkflows/tmp/sg_cov_gam_residual.nf`**.

Working directory: `/data/home/tusr1/01projects/vmap4/10stats.genome/21run_sg_cov_gam_residual`. Outcome: exit 0, **4** succeeded, ~8m 2s.

```bash
cd /data/home/tusr1/01projects/vmap4/10stats.genome/21run_sg_cov_gam_residual && \
source ~/.bashrc && conda activate run && \
nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/tmp/sg_cov_gam_residual.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir /data/home/tusr1/01projects/vmap4 \
  --user_dir /data/home/tusr1 \
  --src_dir /data/home/tusr1/git/script/src \
  --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 \
  --job test_plink \
  --stats_mod test_thin \
  --gam_residual_outlier_frac 0.02
```

Published under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/stats/test_thin/` for `{A,B,D,Others}`:

- `info/{mod}.sample.sg_cov.gam_residual.info.tsv`
- `plots/{mod}.sample.sg_cov.gam_residual_outlier_group_frac.bar.png`
- `plots/{mod}.sample.sg_cov.gam_residual_outlier_rate_by_group.bar.png`
- `plots/{mod}.sample.sg_cov.gam_residual_outliers_vs_logdepth.png`
- `plots/{mod}.sample.sg_cov.gam_residual_outliers_vs_ibs.png`

---

### 2026-06-25 — PopDepCrossChr N1000 GAM-credible deep panel (full genome)

Working directory: `/data/home/tusr1/01projects/vmap4/10stats.genome/30run_popdep_n1000_gam_credible`. Screen `popdep_n1000_gam`. **1000** taxa from GAM credible band (`tb.N1000_gam_credible_deep.taxaBamMap.txt`); **44** chr (skip chr32); **`TIGER_PD_20260619.jar`**, PopDepCrossChr **96 threads / 640 GB**, taxa order **shuffle**, checkpoint **`-ci 500`**. Publish: `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/_popdep_n1000_gam_credible/variant/chrNNN.popdep.txt.bgz` (+ `.tbi`); logs under `.../_popdep_n1000_gam_credible/logs/`. Monitor: `run_logs/monitor.tsv` every 600 s.

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdep_n1000_gam bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/30run_popdep_n1000_gam_credible && bash run_popdep_n1000_gam_credible.sh 2>&1 | tee run_logs/run_all.log'
```

Attach / progress:

```bash
source ~/.bashrc && conda activate run && screen -r popdep_n1000_gam
tail -f /data/home/tusr1/01projects/vmap4/10stats.genome/30run_popdep_n1000_gam_credible/run_logs/monitor.tsv
grep -i 'Finished taxa' /data/home/tusr1/01projects/vmap4/10stats.genome/30run_popdep_n1000_gam_credible/work/*/popdep_crosschr.log 2>/dev/null | tail -3
```

---

### 2026-06-25 — PopDepCrossChr N500 gam-credible deep (64 threads, ci=100)

Working directory: `/data/home/tusr1/01projects/vmap4/10stats.genome/31run_popdep_n500_gam_credible`. Screen `popdep_n500_gam`. **500** taxa (`tb.N500_gam_credible_deep.taxaBamMap.txt`, WGS depth **9.26–48.08×**); **44** chr (skip chr32); **`TIGER_PD_20260619.jar`**, PopDepCrossChr **64 threads / 640 GB**, shuffle, checkpoint **`-ci 100`** under **`/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/_popdep_n500_gam_credible/checkpoint/popdep_crosschr`** (separate from `main_raw/checkpoint/popdep_crosschr` ~8285 taxa and `_popdep_n1000_gam_credible/`). Publish: `.../_popdep_n500_gam_credible/variant/`.

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdep_n500_gam bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/31run_popdep_n500_gam_credible && bash run_popdep_n500_gam_credible.sh 2>&1 | tee run_logs/run_all.log'
```

Monitor:

```bash
tail -f /data/home/tusr1/01projects/vmap4/10stats.genome/31run_popdep_n500_gam_credible/run_logs/monitor.tsv
grep -i 'Finished taxa\|Checkpoint written' /data/home/tusr1/01projects/vmap4/10stats.genome/31run_popdep_n500_gam_credible/work/*/popdep_crosschr.log 2>/dev/null | tail -5
```

---

### 2026-06-28 — PopDepFull N500 chr32 only (PopDepFull supplement)

Working directory: `/data/home/tusr1/01projects/vmap4/10stats.genome/32run_popdep_n500_chr32_popdepfull`. Screen `popdep_n500_chr32`. **500** taxa (`tb.N500_gam_credible_deep.taxaBamMap.txt`); **chr32 only**; **`TIGER_PD_20260619.jar`**, **PopDepFull** (not CrossChr) **32 threads / 640 GB**, `maxForks=2`. Input VCF gate: `main_raw/chr032.vcf.gz` (symlink to job-root). Publish merges into `.../_popdep_n500_gam_credible/variant/chr032.popdep.txt.bgz` (+ `.tbi`).

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdep_n500_chr32 bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/32run_popdep_n500_chr32_popdepfull && bash run_popdep_n500_chr32_popdepfull.sh 2>&1 | tee run_logs/run_all.log'
```

Monitor:

```bash
tail -f /data/home/tusr1/01projects/vmap4/10stats.genome/32run_popdep_n500_chr32_popdepfull/run_logs/monitor.tsv
grep -i 'Finished taxa' /data/home/tusr1/01projects/vmap4/10stats.genome/32run_popdep_n500_chr32_popdepfull/work/*/popdep_32.log 2>/dev/null | tail -5
```

---

### 2026-06-28 — Single-sample mosdepth D_0132 (10 threads)

Working directory: `/data/home/tusr1/01projects/vmap4/12depth/01run_D_0132_mosdepth`. Screen `mosdepth_D_0132`. **1** BAM (`D_0132.bam` + `.bai` staged); **10** threads via `conf/mosdepth_single_bam.config`. Publish: `00data/04depth/04D/D_0132.bam.mosdepth.*`. Also created `final.{A,B,D,ALL}.taxaBamMap.txt` under `00data/05taxaBamMap/` (copy of `all.*` with **D_0132** inserted at coverage **9.44×**).

```bash
source ~/.bashrc && conda activate run && \
screen -dmS mosdepth_D_0132 bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/12depth/01run_D_0132_mosdepth && nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/tmp/mosdepth_single_bam.nf -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config -c /data/home/tusr1/git/script/workflow/Genetics/conf/mosdepth_single_bam.config --user_dir /data/home/tusr1 --home_dir /data/home/tusr1/01projects/vmap4 --depth_publish_dir /data/home/tusr1/01projects/vmap4/00data/04depth/04D --bam_list /data/home/tusr1/01projects/vmap4/12depth/01run_D_0132_mosdepth/bam_list.tsv -process.cache=false 2>&1 | tee run_logs/nextflow_run4.log'
```

Outcome: **1 succeeded**, ~3m46s; `D_0132.bam.mosdepth.summary.txt` total mean **9.44**.

---

### 2026-06-28 — Single-sample mosdepth ABD_0915 + ABD_0938 (10 threads)

Working directory: `/data/home/tusr1/01projects/vmap4/12depth/02run_ABD_0915_0938_mosdepth`. Screen `mosdepth_ABD0915_0938`. **2** BAMs; publish `00data/04depth/03ABD/`. Rebuilt `final.{A,B,D,ALL}.taxaBamMap.txt` from `all.*` with **ABD_0915** / **ABD_0938** in all four maps (D_0132 retained in `final.D` / `final.ALL`).

```bash
source ~/.bashrc && conda activate run && \
screen -dmS mosdepth_ABD0915_0938 bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/12depth/02run_ABD_0915_0938_mosdepth && nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/tmp/mosdepth_single_bam.nf -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config -c /data/home/tusr1/git/script/workflow/Genetics/conf/mosdepth_single_bam.config --user_dir /data/home/tusr1 --home_dir /data/home/tusr1/01projects/vmap4 --depth_publish_dir /data/home/tusr1/01projects/vmap4/00data/04depth/03ABD --bam_list /data/home/tusr1/01projects/vmap4/12depth/02run_ABD_0915_0938_mosdepth/bam_list.tsv -process.cache=false 2>&1 | tee run_logs/nextflow_run1.log'
```

Outcome: **2 succeeded**, ~8m43s; `ABD_0915` total mean **6.41**, `ABD_0938` total mean **6.23** (summary files now 47 lines incl. `total`).

---

### 2026-06-28 — N500 popdep variant annotate + stats (test_thin / test_common_thin)

Four runs under `10stats.genome/` using tmp `popdep_annotate.nf` + `popdep_stats.nf`. Reference: `--popdep_dir .../_popdep_n500_gam_credible`. Publish: `--job test_plink`.

| Run folder | Screen | Mod | Entry | Outcome |
|------------|--------|-----|-------|---------|
| `33run_popdep_n500_annotate_thin` | `popdep_n500_ann_thin` | `test_thin` | `popdep_annotate.nf` | **4/4** ~7m20s → `process/test_thin/variant/*.popdep.info.tsv` |
| `34run_popdep_n500_stats_thin` | `popdep_n500_stats_thin` | `test_thin` | `popdep_stats.nf` | Mahalanobis **4/4**; missing_reg **3/4** (`Others` vmiss `Position` KeyError) |
| `35run_popdep_n500_annotate_common_thin` | `popdep_n500_ann_common` | `test_common_thin` | `popdep_annotate.nf` | **4/4** ~2m24s |
| `36run_popdep_n500_stats_common_thin` | `popdep_n500_stats_common` | `test_common_thin` | `popdep_stats.nf` | Mahalanobis **4/4**; missing_reg **3/4** (`Others` same vmiss issue) |

Run 1 (annotate thin):

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdep_n500_ann_thin bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/33run_popdep_n500_annotate_thin && nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/tmp/popdep_annotate.nf -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config --home_dir /data/home/tusr1/01projects/vmap4 --user_dir /data/home/tusr1 --src_dir /data/home/tusr1/git/script/src --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 --job test_plink --mod test_thin --process_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/test_thin --popdep_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/_popdep_n500_gam_credible 2>&1 | tee run_logs/nextflow.popdep_n500_annotate_thin.log'
```

Run 2 (stats thin):

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdep_n500_stats_thin bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/34run_popdep_n500_stats_thin && nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/tmp/popdep_stats.nf -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config --home_dir /data/home/tusr1/01projects/vmap4 --user_dir /data/home/tusr1 --src_dir /data/home/tusr1/git/script/src --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 --job test_plink --mod test_thin 2>&1 | tee run_logs/nextflow.popdep_n500_stats_thin.log'
```

Run 3 (annotate common_thin):

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdep_n500_ann_common bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/35run_popdep_n500_annotate_common_thin && nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/tmp/popdep_annotate.nf -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config --home_dir /data/home/tusr1/01projects/vmap4 --user_dir /data/home/tusr1 --src_dir /data/home/tusr1/git/script/src --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 --job test_plink --mod test_common_thin --process_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/test_common_thin --popdep_dir /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/_popdep_n500_gam_credible 2>&1 | tee run_logs/nextflow.popdep_n500_annotate_common_thin.log'
```

Run 4 (stats common_thin):

```bash
source ~/.bashrc && conda activate run && \
screen -dmS popdep_n500_stats_common bash -lc 'source ~/.bashrc && conda activate run && cd /data/home/tusr1/01projects/vmap4/10stats.genome/36run_popdep_n500_stats_common_thin && nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/tmp/popdep_stats.nf -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config --home_dir /data/home/tusr1/01projects/vmap4 --user_dir /data/home/tusr1 --src_dir /data/home/tusr1/git/script/src --output_dir /data1/dazheng_tusr1/vmap4.VCF.v1 --job test_plink --mod test_common_thin 2>&1 | tee run_logs/nextflow.popdep_n500_stats_common_thin.log'
```

Publish layout: `test_plink/process/{mod}/variant/*.popdep.info.tsv`; `test_plink/stats/{mod}/plots/*.variant.popdep_miss.reg_*.png` (A/B/D only); `test_plink/stats/{mod}/info/*.variant.mahalanobis.info.tsv` + `plots/*.mahalanobis.png` (all subgenomes).
