# Genetics Nextflow workflow (`workflow/Genetics/`)

**Documentation policy (`workflow/`):** Under `workflow/`, keep **exactly one** README—**this file** (`workflow/Genetics/README.md`). Do **not** add `README.md` elsewhere under `workflow/` (including `tmp/`, `integrated/`, or deeper subtrees). New workflow topics, tmp entry conventions, and wheat integrated modes belong **here** as new sections or bullets.

This directory holds the **Genotype** branch of the pipeline: PLINK/PLINK2 preprocessing, per-mod processors, stats (Python/R via `src/`), small **auxiliary** entry scripts under `tmp/`, and **wheat** table/matrix analytics under `integrated/`. This README focuses on what **`main.nf`** runs and how auxiliary scripts are launched.

---

## Entry point and config

| Item | Role |
| --- | --- |
| `main.nf` | Top-level router: for non-`wheat_*` mods, builds `ch_vcf` in `check_input`, then dispatches on `params.mod`. For `params.mod` starting with `wheat_`, skips VCF input and runs `integrated_wheat`. |
| `nextflow.config` | Default `params`, resource labels, conda hints. Pass with **absolute** `-c` when launching from arbitrary cwd (required for `tmp/*.nf`). |

Typical invocation shape (full examples and chronological runs: **`doc/NF_CMD.md`**):

```text
nextflow run <path-to>/workflow/Genetics/main.nf \
  -c <absolute-path-to>/workflow/Genetics/nextflow.config \
  --home_dir … --user_dir … --src_dir … --output_dir … --mod <mod> --job <job> [optional flags]
```

---

## Active `params.mod` values (`main.nf` router)

These branches are wired in `workflow { }`. PLINK tracks below run end-to-end genotype processors plus stats; `wheat_*` modes run the integrated table/matrix workflow only (no VCF channel). Details: **Wheat integrated workflows** (section below).

| `params.mod` | Processor workflow | Stats |
| --- | --- | --- |
| `v1_plink` | `PLINK_PROCESSOR` (`plink_processor`) | `PLINK_STATS` (`plink_stats`) |
| `test_plink` or `test_thin` | `TEST_PLINK_PROCESSOR` (`test_plink_processor`) | `TEST_PLINK_STATS` (`test_plink_stats`) |
| `test_plink_camp` or `test_camp` | `TEST_PLINK_CAMP_PROCESSOR` + **`params.camp`** TSV | `TEST_PLINK_STATS` |
| `test_common_thin` | `TEST_COMMON_THIN_PROCESSOR` (hard-filter + thin from `hf_*` params) | `TEST_PLINK_STATS` |
| `wheat_*` | *(none — table/matrix integrated track)* | `integrated_wheat` (`integrated/integrated_wheat.nf`) |

Single-table wheat modes use `wheat_table_input`; `wheat_gwas` / `wheat_kgwas` use dedicated genotype or k-mer + phenotype params (see **Wheat integrated workflows** below).

Any other `params.mod` falls through to the “no workflow module was chosen” path (no-op aside from input check).

### Included but not routed (known gap)

`genotype/database.nf`, `dynamic/kinship.nf`, `dynamic/ps.nf`, `static/gwas/gwas.nf`, and `genotype/hail.nf` are **included** at the top of `main.nf` but have **no** `params.mod` branch yet. See **`doc/TODO.md`** §2 (“Router gap”) for the backlog to expose `HAIL`, kinship, population structure, GWAS, and `database` from the main router.

---

## Important `params` (from `nextflow.config`)

| Parameter | When to set | Notes |
| --- | --- | --- |
| `process_dir` | Reuse prebuilt per-chromosome or merged test pfiles | When merged `A_test.plink2` … `Others_test.plink2` exist under `{process_dir}`, test processors **skip** thin/merge and rebuild basic info + LD only. See frozen test paths in **`doc/TODO.md`** §2. |
| `camp` | **`test_plink_camp` / `test_camp` only** | Path to cohort map TSV (e.g. `camp_vmap4_map.tsv`). Required for those mods. |
| `chr` | Optional | Comma-separated chromosome filter after job config builds the VCF channel. |
| `wheat_table_input`, `wheat_gwas_*`, `wheat_kgwas_*` | **`params.mod` starts with `wheat_`** | See **Wheat integrated workflows** below. Requires `user_dir` for the `stats` Conda env (same as genotype stats processes). |

Hard-filter / common-thin tuning for `test_common_thin` uses `params.hf_maf`, `params.hf_geno`, `params.hf_thin_rate` (see `nextflow.config`; retune and log when cohort policy changes).

---

## Output layout (test jobs)

For test modes (e.g. `job` = `test_plink`), **`params.job`** is the **top-level** folder under `params.output_dir`. Genotype outputs go under **`{output_dir}/{job}/process/{mod}/`** (sample, variant, logs, etc.). Stats outputs go under **`{output_dir}/{job}/stats/{mod}/`** (plots, thresholds, info, logs). Assess mini-workflows publish under **`{output_dir}/{job}/assess/{mod}/`**. Wheat integrated modes publish under **`{output_dir}/{job}/integrated/{mod}/`** (`info/` and `plots/`). Do not encode `mod` by changing `params.job` alone—the pipeline uses **`params.mod`** under those trees.

---

## Auxiliary entry scripts (`tmp/`)

Small workflows in **`workflow/Genetics/tmp/`** reuse the same **`workflow/Genetics/nextflow.config`** as `main.nf`. Pass the config with an **absolute** path so it works no matter what the launch directory is:

```bash
nextflow run /data/home/tusr1/git/script/workflow/Genetics/tmp/ld_plots_redraw.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  ...
```

Entry scripts use `include { ... } from '../genotype/...'` paths relative to `tmp/`.

| Script | Purpose |
| --- | --- |
| `tmp/ld_plots_redraw.nf` | Re-run LD decay + cross-chr plot processes from existing `process/<mod>/variant/*.vcor` without rerunning the full processor. |
| `tmp/assess_plink_debug.nf` | Tier-1 assess on `test_thin` / `test_common_thin`: PLINK2 `--freq` / `--missing` on a representative-chr slice of each `*_test.plink2`, MAF-bin TSV from `.afreq`, then Python plots via `assess_slice.py` under `assess/<mod>/plots` and tables under `assess/<mod>/info`. |

---

## Wheat integrated workflows (`integrated/`)

Table- and matrix-based analyses in **`src/python/genetics/wheat/`** are run from **`main.nf`** when **`params.mod`** starts with **`wheat_`** (VCF `check_input` is skipped). Nextflow module: **`integrated/integrated_wheat.nf`** (`workflow integrated_wheat`).

Python matches genotype stats style: each process uses `#!/usr/bin/env python` and **imports** `run_*` from `genetics.wheat.*` (no `python -m` CLI).

### Routing summary

When `params.mod` starts with `wheat_`, the workflow runs `integrated/integrated_wheat.nf`.

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

Tuning knobs (defaults in `nextflow.config`): `wheat_snp_qc_*`, `wheat_pca_*`, `wheat_tagsnp_*`, `wheat_hapmap_window_size`, `wheat_cnv_*`, `wheat_gwas_trait`, `wheat_kgwas_trait`.

### Required launch params (wheat)

- `output_dir`, `job`, `user_dir` (Conda `stats` env, same pattern as `genotype/stats.nf`)
- `mod` as one of the table rows above
- Paths per mod (e.g. `--wheat_table_input` for single-table tasks)

`home_dir` / `src_dir` are not used by wheat processes but may still be set for a uniform CLI with other mods.

### Wheat output layout

Artifacts publish under **`{output_dir}/{job}/integrated/{mod}/info`** (TSV) and **`.../plots`** (PNG when figures exist). Output file prefixes are **`{job}.{mod}`** (e.g. `myrun.wheat_gwas.gwas.tsv`).

Topic folders under `integrated/` (`gwas/`, `snp_qc/`, …) are organizational only; the executable surface is **`integrated_wheat.nf`** plus **`main.nf`** routing.

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

- Launch Nextflow and other long pipeline commands from a **task-specific run folder** under **`/data/home/tusr1/01projects/vmap4/...`** (not from the git repo root). Use monotonically increasing numeric run directories where your team convention applies (e.g. `08stats.genome/NNrun_…`).
- Use conda env **`run`** for Nextflow (`source ~/.bashrc && conda activate run`, or equivalent).
- For **long** production-style runs, use **`screen`** (or equivalent); short debug / `-preview` runs are exempt.
- Treat **`/data1/dazheng_tusr1/vmap4.VCF.v1`** test assets as **read-only**; use prepared paths documented in **`doc/TODO.md`**.

Authoritative guardrails: **`.cursor/rules/workstation-core.mdc`** and **`.cursor/rules/workstation-nextflow.mdc`**.

---

## Downstream tools

Processes invoke **`src/python`**, **`src/r`**, and **`src/java`** as configured in each `process` block (typically conda `stats` under `params.user_dir` for Python stats). Keep processor vs stats boundaries as in **`workstation-nextflow.mdc`** (heavy `plink2` in `processor.nf`; plotting / summarization in `stats.nf`).

---

## Further reading

- **`doc/NF_CMD.md`** — chronological full command lines and cwd notes.
- **`doc/TODO.md`** — checklist, router gaps, frozen paths.
- **`.cursor/rules/workstation-*.mdc`** — workstation and pipeline policy.
