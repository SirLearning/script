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

Working directory: `/data/home/tusr1/01projects/vmap4/05reliable.lib/09run_abstract_mq_50_bams_i16` (screen e.g. `mq50_i16_resume`). Reuses `08run` work dir and publish tree (`--output_dir` → 08run).

```bash
screen -S mq50_i16_resume
/data/home/tusr1/01projects/vmap4/05reliable.lib/09run_abstract_mq_50_bams_i16/resume_from_08.sh
```

Log: `09run_.../run_logs/nextflow.mq50_i16_resume.log`. Completed chromosomes from 08 are skipped via `-resume`; remaining 39 tasks run at `maxForks=16`.
