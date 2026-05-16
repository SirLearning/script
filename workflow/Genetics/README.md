# Genetics Nextflow workflow (`workflow/Genetics/`)

This directory holds the **Genotype** branch of the pipeline: PLINK/PLINK2 preprocessing, per-mod processors, stats (Python/R via `src/`), and small **auxiliary** entry scripts under `tmp/`. Phenotype, GWAS product flows, and other themes live elsewhere in the repo; this README focuses on what **`main.nf`** actually runs today.

---

## Entry point and config

| Item | Role |
| --- | --- |
| `main.nf` | Top-level router: builds `ch_vcf` in `check_input`, then dispatches on `params.mod`. |
| `nextflow.config` | Default `params`, resource labels, conda hints. Pass with **absolute** `-c` when launching from arbitrary cwd (required for `tmp/*.nf`). |

Typical invocation shape (full examples and chronological runs: **`doc/NF_CMD.md`**):

```text
nextflow run <path-to>/workflow/Genetics/main.nf \
  -c <absolute-path-to>/workflow/Genetics/nextflow.config \
  --home_dir … --user_dir … --src_dir … --output_dir … --mod <mod> --job <job> [optional flags]
```

---

## Active `params.mod` values (`main.nf` router)

These branches are wired in `workflow { }` and execute end-to-end genotype + stats for PLINK test/production tracks:

| `params.mod` | Processor workflow | Stats |
| --- | --- | --- |
| `v1_plink` | `PLINK_PROCESSOR` (`plink_processor`) | `PLINK_STATS` (`plink_stats`) |
| `test_plink` or `test_thin` | `TEST_PLINK_PROCESSOR` (`test_plink_processor`) | `TEST_PLINK_STATS` (`test_plink_stats`) |
| `test_plink_camp` or `test_camp` | `TEST_PLINK_CAMP_PROCESSOR` + **`params.camp`** TSV | `TEST_PLINK_STATS` |
| `test_common_thin` | `TEST_COMMON_THIN_PROCESSOR` (hard-filter + thin from `hf_*` params) | `TEST_PLINK_STATS` |

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

Hard-filter / common-thin tuning for `test_common_thin` uses `params.hf_maf`, `params.hf_geno`, `params.hf_thin_rate` (see `nextflow.config`; retune and log when cohort policy changes).

---

## Output layout (test jobs)

For test modes (e.g. `job` = `test_plink`), **`params.job`** is the **top-level** folder under `params.output_dir`. Genotype outputs go under **`{output_dir}/{job}/process/{mod}/`** (sample, variant, logs, etc.). Stats outputs go under **`{output_dir}/{job}/stats/{mod}/`** (plots, thresholds, info, logs). Assess mini-workflows publish under **`{output_dir}/{job}/assess/{mod}/`**. Do not encode `mod` by changing `params.job` alone—the pipeline uses **`params.mod`** under those trees.

---

## Auxiliary entry scripts (`tmp/`)

Small partial reruns live next to `main.nf` but **must** load the same config with an absolute `-c` path. See **`tmp/README.md`** for the include convention.

| Script | Purpose |
| --- | --- |
| `tmp/ld_plots_redraw.nf` | Re-run LD decay + cross-chr plot processes from existing `process/<mod>/variant/*.vcor` without rerunning the full processor. |
| `tmp/assess_plink_debug.nf` | Tier-1 assess on `test_thin` / `test_common_thin`: PLINK2 `--freq` / `--missing` on a representative-chr slice of each `*_test.plink2`, MAF-bin TSV from `.afreq`, then Python plots via `assess_slice.py` under `assess/<mod>/plots` and tables under `assess/<mod>/info`. |

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
- **`tmp/README.md`** — `-c` absolute path reminder for `tmp/*.nf`.
