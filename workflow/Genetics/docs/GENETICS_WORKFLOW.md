# Genetics Nextflow workflow (`workflow/Genetics/`)

**Documentation:** This file (`docs/GENETICS_WORKFLOW.md`) is the **only** narrative operator guide for this tree. Previous theme trees may live under `old/` for reference; **active** process libraries are **`modules/local/{genotype,dynamic,static,integrated}/`**.

This directory holds the **Genotype** branch of the pipeline: PLINK/PLINK2 preprocessing, per-mod processors, stats (Python/R via `src/`), **composed subworkflows** under `subworkflows/local/`, ops scripts under `subworkflows/tmp/ops/`, and **wheat** table/matrix analytics.

---

## Entry point and config

| Item | Role |
| --- | --- |
| `main.nf` | Anonymous top-level `workflow { }`: for non-`wheat_*` mods, builds `ch_vcf` via `check_input` from `subworkflows/local/entry/genetics_helpers.nf`, then dispatches on `params.mod`. For `params.mod` starting with `wheat_`, skips VCF input and runs `RUN_WHEAT_INTEGRATED` (`subworkflows/local/wheat/wheat_integrated_study.nf`). |
| `subworkflows/local/entry/genetics_helpers.nf` | `header()`, `helpMessage()`, and `workflow check_input` (VCF channel + summaries). |
| `subworkflows/local/plink/plink_genotype_modes.nf` | Composed workflows `RUN_V1_PLINK`, `RUN_TEST_PLINK`, `RUN_TEST_PLINK_CAMP`, `RUN_TEST_COMMON_THIN` (processor + stats only). Router-gap modules: `plink/analysis_extensions.nf`. |
| `subworkflows/local/entry/partial_router.nf` | Partial reruns via `--partial_task` (assess, stats redraw, wheat-from-plink). |
| `nextflow.config` | Default `params`, resource labels, conda hints. Pass with **absolute** `-c` when launching auxiliary scripts from arbitrary cwd. |

Typical invocation shape (full examples and chronological runs: **`doc/NF_CMD.md`**):

```text
nextflow run <path-to>/workflow/Genetics/main.nf \
  -c <absolute-path-to>/workflow/Genetics/nextflow.config \
  --home_dir … --user_dir … --src_dir … --output_dir … --mod <mod> --job <job> [optional flags]
```

---

## Active `params.mod` values (`main.nf` router)

These branches are wired in the **anonymous** top-level workflow in `main.nf`. PLINK tracks below run end-to-end genotype processors plus stats; `wheat_*` modes run the integrated table/matrix workflow only (no VCF channel). Details: **Wheat integrated workflows** (section below).

| `params.mod` | Processor workflow | Stats |
| --- | --- | --- |
| `v1_plink` | `PLINK_PROCESSOR` (`plink_processor`) | `PLINK_STATS` (`plink_stats`) |
| `test_thin` | `TEST_PLINK_PROCESSOR` (`test_plink_processor`) | `TEST_PLINK_STATS` (`test_plink_stats`) |
| `test_camp` | `TEST_PLINK_CAMP_PROCESSOR` + **`params.camp`** TSV | `TEST_PLINK_STATS` |
| `test_common_thin` | `TEST_COMMON_THIN_PROCESSOR` (hard-filter + thin from `hf_*` params) | `TEST_PLINK_STATS` |
| `wheat_*` | *(none — table/matrix integrated track)* | `integrated_wheat` (`modules/local/integrated/integrated_wheat.nf`) |

Single-table wheat modes use `wheat_table_input`; `wheat_gwas` / `wheat_kgwas` use dedicated genotype or k-mer + phenotype params (see **Wheat integrated workflows** below).

Any other `params.mod` causes `main.nf` to **exit 1** with a list of supported mods (no silent no-op).

### Included but not routed (known gap)

`analysis_extensions.nf` collects router-gap modules (`database`, `kinship`, `ps`, `GWAS`, `HAIL`) for future `main.nf` branches — they are **not** loaded by the active PLINK router. See **`doc/TODO.md`** §2 (“Router gap”).

### Process libraries (`modules/local/`)

**Active** Nextflow process definitions live under **`workflow/Genetics/modules/local/genotype/`**, **`.../dynamic/`**, **`.../static/`**, and **`.../integrated/`**. Shared DSL helpers live under **`.../infra/`**. Earlier copies or experiments may remain under **`workflow/Genetics/old/`**; new work should `include` from **`modules/local/...`** only. Cross-theme imports (e.g. GWAS pulling `format_vcf_plink` from genotype) use paths relative to the importing file under `modules/local/`.

**Infra library** (`modules/local/infra/`) — cross-cutting Groovy helpers (no `process` blocks):

| File | Contents |
| --- | --- |
| `infra/utils.nf` | Index / compose-only barrel (do not `include { def }` from here; use `infra_*.nf` directly) |
| `infra/infra_tools.nf` | Java setup, software/TIGER/samtools paths, `validatePaths` |
| `infra/infra_job_config.nf` | `getJobConfig`, calling/pop/chr config, `getVcfIdFromPath` |
| `infra/infra_tiger.nf` | TIGER JAR config, FastCall server population lists |
| `infra/infra_ref_v1.nf` | vmap4 IWGSC v1 chr maps, offsets, lengths, taxaBam paths |
| `infra/infra_plink_reuse.nf` | Merged test pfile reuse tuples (`hasMergedSubgenomeTestPfiles`, …) |

**Genotype processor library** (`modules/local/genotype/processor/`) — heavy `plink2` / TIGER / bcftools compute; workflows only in the barrel file:

| File | Contents |
| --- | --- |
| `processor/processor.nf` | Composed workflows `test_plink_processor`, `test_plink_camp`, `test_common_thin_processor`, `plink_processor`, `plink_preprocess`, `vcf_arrange_merge`, `vcf_bcftools_filter` |
| `processor/processor_vcf.nf` | VCF format, arrange, merge, `format_vcf_plink` |
| `processor/processor_test.nf` | Test subsampling and subgenome merge |
| `processor/processor_plink2.nf` | PLINK2 basic info, PCA, tag-SNP prune, CNV awk, LD matrices |
| `processor/processor_legacy.nf` | Legacy vcftools basic info |
| `processor/processor_depth.nf` | Population depth (TIGER) |
| `processor/processor_filter.nf` | Sample/variant PLINK filters; VCF v0 filters |
| `processor/processor_assess.nf` | Assess compute (`plink2_assess_debug_slice`, `quick_count`, `bcftools_qc_assess`) |

**Genotype stats library** (`modules/local/genotype/stats/`) — plotting and summarization only; no heavy `plink2` compute:

| File | Contents |
| --- | --- |
| `stats/stats.nf` | Composed workflows `test_plink_stats`, `plink_stats` (includes all sub-libraries below) |
| `stats/stats_variant.nf` | Variant QC, MAC, popdep, LD plot processes |
| `stats/stats_sample.nf` | Sample QC, kinship/IBS, PCA helper processes |
| `stats/stats_assess.nf` | Assess debug plot processes (`assess_plink_debug_plots`, `dumpnice_vcf_qc_assess`) |
| `stats/stats_integrated.nf` | Wheat integrated PLINK2 plot/report processes |
| `stats/stats_chr_report.nf` | Chr variant count tables and thin vs common-thin compare plots |

**Upstream genotype libraries** (not routed from `main.nf` yet — standalone entry scripts; composed in `subworkflows/local/upstream/raw_data_upstream.nf`):

| Path | Entry | Contents |
| --- | --- | --- |
| `genotype/calling/` | `calling/caller.nf` | FastCall3 workflows (`run_FastCall3`, `load_lib_files`) + `caller_prep.nf` / `caller_fastcall.nf` |
| `genotype/align/` | `align/align.nf` | `RUN_ALIGN_USB_TRANSFER` workflow; BAM USB transfer + MD5 (`align_transfer.nf`, `align_md5.nf`); BWA-MEM2 (`align_bwa.nf`) |

**Hail library** (`genotype/hail/`) — included as `HAIL` workflow (router gap; see §2 TODO):

| File | Contents |
| --- | --- |
| `hail/hail.nf` | Composed workflow `HAIL` |
| `hail/hail_io.nf` | `vcf_to_mt_hail`, `filter_hail` |
| `hail/hail_stats.nf` | QC, kinship, PCA processes |
| `hail/hail_gwas.nf` | `hail_normal_gwas` |

**Genotype database library** (`genotype/database/`) — VCF → Hail MT + bcftools annotation helpers:

| File | Contents |
| --- | --- |
| `database/database.nf` | Composed workflow `database` |
| `database/database_build.nf` | `build_hail_from_vcf` |
| `database/database_annotate.nf` | bcftools annotate / tab / BED helpers |

**Static GWAS library** (`modules/local/static/gwas/`):

| File | Contents |
| --- | --- |
| `gwas/gwas.nf` | Composed workflow `GWAS` (legacy PLINK1 track) |
| `gwas/gwas_legacy.nf` | PCA, pheno/covar prep, PLINK1 / GAPIT / rMVP processes |
| `gwas/gwas_plink2.nf` | `plink2_gwas_glm`, `gcta_gwas` (wheat integrated compute) |
| `gwas/gwas_plot.nf` | `plot_gwas_association`, `plot_gwas` |

**Include paths:** from `subworkflows/local/*/` use `../../../modules/local/...`; from `modules/local/genotype/*` use `../../infra/<file>.nf`. Partial reruns: launch `subworkflows/local/entry/partial_router.nf` with `--partial_task`, or `include` the smallest `partial/*.nf` / module sub-file.

---

## Important `params` (from `nextflow.config`)

| Parameter | When to set | Notes |
| --- | --- | --- |
| `process_dir` | Reuse prebuilt per-chromosome or merged test pfiles | When merged `A_test.plink2` … `Others_test.plink2` exist under `{process_dir}`, test processors **skip** thin/merge and rebuild basic info + LD only. See frozen test paths in **`doc/TODO.md`** §2. |
| `camp` | **`test_camp` only** | Path to cohort map TSV (e.g. `camp_vmap4_map.tsv`). Required for that mod. |
| `chr` | Optional | Comma-separated chromosome filter after job config builds the VCF channel. |
| `wheat_table_input`, `wheat_gwas_*`, `wheat_kgwas_*` | **`params.mod` starts with `wheat_`** | See **Wheat integrated workflows** below. Requires `user_dir` for the `stats` Conda env (same as genotype stats processes). |

Hard-filter / common-thin tuning for `test_common_thin` uses `params.hf_maf`, `params.hf_geno`, `params.hf_thin_rate` (see `nextflow.config`; retune and log when cohort policy changes).

---

## Output layout (test jobs)

For test modes (e.g. `params.job` = `test_plink` while `params.mod` is `test_thin` / `test_camp` / `test_common_thin`), **`params.job`** is the **top-level** folder under `params.output_dir`. Genotype outputs go under **`{output_dir}/{job}/process/{mod}/`** (sample, variant, info, logs, etc.). Stats outputs go under **`{output_dir}/{job}/stats/{mod}/`** (plots, thresholds, info, logs). Assess debug (`partial_router.nf --partial_task assess_plink`) uses the same trees: compute artefacts under **`process/{mod}/`**, plots and summary tables under **`stats/{mod}/`**. Wheat integrated plots/tables from merged PLINK2 (`partial_router.nf --partial_task wheat_from_plink` / `WHEAT_STUDY_FROM_PLINK`) publish under **`{output_dir}/{job}/integrated/{plink_source_mod}/{wheat_task_mod}/`** (e.g. `integrated/test_thin/wheat_pca_tsne/{info,plots}`). Table-only `wheat_*` modes from **`main.nf`** still use **`{output_dir}/{job}/integrated/{params.mod}/`**. Do not encode `mod` by changing `params.job` alone—the pipeline uses **`params.mod`** under those trees where applicable.

---

## Composed subworkflows (`subworkflows/local/`)

| Folder | Files | Role |
| --- | --- | --- |
| `entry/` | `genetics_helpers.nf`, `partial_router.nf` | `main.nf` helpers; partial rerun router (`--partial_task`) |
| `plink/` | `plink_genotype_modes.nf`, `analysis_extensions.nf` | PLINK router workflows; router-gap module index |
| `wheat/` | `wheat_integrated_study.nf` | `RUN_WHEAT_INTEGRATED`, `WHEAT_STUDY_FROM_PLINK`, `RUN_WHEAT_FROM_PLINK` |
| `upstream/` | `raw_data_upstream.nf` | `RUN_ALIGN_USB_TRANSFER`, `run_FastCall3` |
| `partial/` | `partial_assess.nf`, `partial_stats.nf` | Assess + stats partial workflow bodies |

---

## Auxiliary scripts (`subworkflows/tmp/`)

**Partial reruns** — one entry, no per-task duplicate scripts under `workflow/Genetics/tmp/` (removed 2026-06-09):

```bash
nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/local/entry/partial_router.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --partial_task ld_redraw --mod test_thin --output_dir … --job test_plink …
```

| `--partial_task` | Workflow | Typical `params.mod` / notes |
| --- | --- | --- |
| `assess_plink` | `RUN_ASSESS_PLINK_DEBUG` | `test_thin` / `test_common_thin` |
| `assess_vcf` | `RUN_ASSESS_VCF_DEBUG` | any mod with `export/*.debug.vcf.gz` |
| `ld_redraw` | `RUN_LD_PLOTS_REDRAW` | e.g. `test_thin` |
| `mac_stats` | `RUN_MAC_STATS_FROM_GCOUNT` | test mods or `test_rebulld_lib` |
| `mac_dist_redraw` | `RUN_MAC_DIST_LOG_REDRAW` | multi-mod batch |
| `chr_counts` | `RUN_CHR_VARIANT_COUNTS` | per `params.mod` |
| `chr_compare` | `RUN_CHR_VARIANT_COMPARE` | thin vs common-thin |
| `rebuild_lib_stats` | `RUN_REBUILLD_LIB_STATS` | `test_rebulld_lib` only |
| `wheat_from_plink` | `RUN_WHEAT_FROM_PLINK` | `test_thin` / `test_common_thin` + `wheat_integrated_mod` |

**Ops** (`subworkflows/tmp/ops/`) — self-contained FTP upload helpers (not genetics analysis):

| Script | Purpose |
| --- | --- |
| `ops/gsa_ftp_upload.nf` | GSA FTP upload (`submit.big.ac.cn`); set `GSA_FTP_PASSWORD`. |
| `ops/gwh_ftp_upload.nf` | GWH FTP upload; optional `--verify_md5` / `--verify_only`. Set `GSA_FTP_PASSWORD`. |

**Partial task notes** (workflow bodies in `subworkflows/local/partial/`):

| `--partial_task` | Deliverables / behaviour |
| --- | --- |
| `ld_redraw` | Re-run LD decay + cross-chr plots from existing `process/<mod>/variant/*.vcor` (no full processor). |
| `assess_plink` | `plink2_assess_debug_slice` → `process/<mod>/`; `assess_plink_debug_plots` → `stats/<mod>/` via `assess_slice.py`. |
| `assess_vcf` | Legacy VCF assess on `process/<mod>/export/*.debug.vcf.gz`; prefer PLINK2 when pfiles exist. |
| `chr_counts` | `process/<mod>/info/<mod>.chr_variant_counts.tsv` + `.by_ref.tsv`. |
| `chr_compare` | Genome-wide density panels + thin/common fraction bars under `stats/thin_common_compare/`. |
| `rebuild_lib_stats` | `TEST_PLINK_STATS` on `process/test_rebulld_lib/` chr002 tables; LD plots skipped. |
| `mac_stats` | MAC bar (0–100) + het-fraction heatmap; `*.variant.mac.mac0.info.tsv` audit for MAC=0. |
| `mac_dist_redraw` | Re-draw MAC distribution log plots from existing gcount tables. |
| `wheat_from_plink` | `RUN_WHEAT_FROM_PLINK` on existing merged test pfiles (`wheat_integrated_mod` required). |

Historic runs logged under `doc/NF_CMD.md` may still reference removed `workflow/Genetics/tmp/*.nf` paths; use `partial_router.nf` for new runs.

---

## Wheat integrated workflows (`modules/local/integrated/`)

Wheat integrated modes are run from **`main.nf`** when **`params.mod`** starts with **`wheat_`** (VCF `check_input` is skipped). PLINK2 compute lives in **`processor/processor_plink2.nf`**; plots/reports in **`stats/stats_integrated.nf`**; study wiring in **`subworkflows/local/wheat/wheat_integrated_study.nf`**.

Python matches genotype stats style: inline `#!/usr/bin/env python` imports from **`genetics.genomics.*`** and **`genetics.gwas.association_plot`** (no `python -m` CLI). **`genetics.wheat`** is a deprecated shim only — do not add new imports there.

### Routing summary

When `params.mod` starts with `wheat_`, the workflow runs `modules/local/integrated/integrated_wheat.nf`.

### Active `params.mod` values (wheat)

| `params.mod` | Input params | Python entry |
| --- | --- | --- |
| `wheat_snp_qc` | `wheat_table_input` | `genetics.wheat.snp_qc.run_snp_qc` |
| `wheat_pca_tsne` | `wheat_table_input` | `genetics.wheat.population_structure.run_population_structure` |
| `wheat_tagsnp` | `wheat_table_input` | `genetics.wheat.tagsnp.run_tagsnp_selection` |
| `wheat_hapmap` | `wheat_table_input` | `genetics.wheat.hapmap.run_hapmap_build` |
| `wheat_cnv` | `wheat_table_input` | `genetics.wheat.cnv.run_cnv_calling` |
| `wheat_genetic_map` | `wheat_table_input` | `genetics.wheat.genetic_map.run_genetic_map` |
| `wheat_gwas` | `wheat_gwas_genotype`, `wheat_gwas_phenotype` | `genetics.wheat.wheat_gwas.run_gwas` |
| `wheat_kgwas` | `wheat_kgwas_kmer_matrix`, `wheat_kgwas_phenotype` | `genetics.wheat.kgwas.run_kgwas` |

Tuning knobs (defaults in `conf/base.config`): `wheat_snp_qc_*`, `wheat_pca_*`, `wheat_tsne_*` (t-SNE on PLINK PCs in `wheat_pca_tsne`), `wheat_tagsnp_*`, `wheat_hapmap_window_size`, `wheat_cnv_*`, `wheat_gwas_trait`, `wheat_kgwas_trait`.

### Required launch params (wheat)

- `output_dir`, `job`, `user_dir` (Conda `stats` env, same pattern as `modules/local/genotype/stats/stats.nf`)
- `mod` as one of the table rows above
- Paths per mod (e.g. `--wheat_table_input` for single-table tasks)

`home_dir` / `src_dir` are not used by wheat processes but may still be set for a uniform CLI with other mods.

### Wheat output layout

Artifacts publish under **`{output_dir}/{job}/integrated/{mod}/info`** (TSV) and **`.../plots`** (PNG when figures exist). For per–subgenome **`wheat_pca_tsne`** from PLINK2, expect **`*.pca.tsv`**, **`*.eigen.tsv`**, **`*.tsne.tsv`**, **`*.pca.png`**, **`*.variance.png`**, **`*.tsne.png`** (prefix = subgenome id). Table-only entry modes use **`{mod}`** as the prefix (e.g. `wheat_gwas.gwas.tsv`).

Topic folders under `modules/local/integrated/` (`gwas/`, `snp_qc/`, …) are organizational only; the executable surface is **`integrated_wheat.nf`** plus **`main.nf`** routing.

### Example (`wheat_snp_qc`)

```text
nextflow run /data/home/tusr1/git/script/workflow/Genetics/main.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --user_dir /data/home/tusr1 \
  --output_dir /path/to/out \
  --job demo_wheat \
  --mod wheat_snp_qc \
  --wheat_table_input /path/to/snp_summary.tsv
```

Log production runs in **`doc/NF_CMD.md`** per repo convention.

---

## Where and how to run (workstation policy)

**Do not run from the git repo.** Use the vmap4 project workspace with a **two-level** folder layout:

| Level | Path pattern | Example |
| --- | --- | --- |
| Project root | `/data/home/tusr1/01projects/vmap4/` | All runtime work lives here |
| Task module | `<NN><area>.<topic>/` | `08stats.genome/` — stats, assess, LD/MAC/chr partials, wheat-from-plink |
| Run attempt | `<NN>run` or `<NN>run_<slug>/` | `57run_mac_stats_test_thin`, `23run_assess_plink2_debug_stub` |

**Before each new run:**

1. `cd /data/home/tusr1/01projects/vmap4/<module>/`
2. List existing `*run*` folders; set `NN` to max + 1 (zero-padded: `01`, `02`, … `63`, …).
3. `mkdir -p <NN>run_<descriptive_slug>` and `cd` into it.
4. Launch with absolute repo paths, e.g.:
   - Full router: `…/workflow/Genetics/main.nf`
   - Partial: `…/subworkflows/local/entry/partial_router.nf --partial_task <name>`
   - Always: `-c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config`
5. Append **cwd + full command** to **`doc/NF_CMD.md`**.

**Conventions**

- **Do not reuse** an old run folder for a fresh attempt — create the next `NNrun_…`.
- **Conda `run`** for Nextflow; **`stats`** for Python stats processes inside NF.
- **`screen`** for long production runs; `-preview` / small debug exempt.
- **Read-only inputs:** `/data1/dazheng_tusr1/vmap4.VCF.v1` (see **`doc/TODO.md`** for prepared `test_plink/process` paths).
- **Ephemeral preview dirs** (`00nf_preview*`) under `vmap4/` are OK for syntax smoke tests.

**Partial-task → typical slug examples**

| `--partial_task` | Example run folder under `08stats.genome/` |
| --- | --- |
| `assess_plink` | `24run_assess_plink2_full_test_thin` |
| `ld_redraw` | `16run_ld_plots_redraw_test_thin` |
| `mac_stats` | `57run_mac_stats_test_thin` |
| `chr_compare` | `51run_chr_variant_compare_plots` |
| `wheat_from_plink` | `43run_wheat_pca_grp_test_thin` |

Authoritative guardrails: **`.cursor/rules/workstation-core.mdc`** and **`.cursor/rules/workstation-nextflow.mdc`**.

---

## Downstream tools

Processes invoke **`src/python`**, **`src/r`**, and **`src/java`** as configured in each `process` block (typically conda `stats` under `params.user_dir` for Python stats). Keep processor vs stats boundaries as in **`workstation-nextflow.mdc`** (heavy `plink2` in `modules/local/genotype/processor/`; plotting / summarization under `modules/local/genotype/stats/`).

---

## Further reading

- **`doc/NF_CMD.md`** — chronological full command lines and cwd notes.
- **`doc/TODO.md`** — checklist, router gaps, frozen paths.
- **`.cursor/rules/workstation-*.mdc`** — workstation and pipeline policy.
