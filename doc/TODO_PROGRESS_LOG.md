# TODO progress log

Append-only audit log for completed TODO / ops work (per `.cursor/rules/workstation-core.mdc`). **Do not delete or rewrite** existing body text; add new material at the **end** only (including `SUPP` sections).

**Section headings:** `## YYYY-MM-DD — <short English title>` (and `## YYYY-MM-DD — SUPP — <topic>` for supplements). No opaque IDs (e.g. `OPS-NF-001`); the date plus title should be enough to find an entry.

Each entry should read like an **engineering report**: goals, what changed and why, validation (**pass / fail / blocked**), root causes and fixes, next steps, and how outputs are laid out (formats, `publishDir`, prefixes). **Full command lines** belong in **`doc/NF_CMD.md`**; here, cite that file by heading and record **working directory(ies)** per run. English; prose or tables are both fine.

---

## 2026-05-13 — Establish master TODO + progress log in repo

| Field | Detail |
| --- | --- |
| **Date** | 2026-05-13 |
| **Summary** | Import vmap4 workflow methodology/process/assess/filter checklist |
| **Files changed** | `doc/TODO.md` (created), `doc/TODO_PROGRESS_LOG.md` (created) |
| **Validation** | Manual: checklist matches user-supplied snapshot; open `[ ]` and done `[x]` preserved. |
| **Outcome** | Master list versioned in repo; future completions should flip boxes in `doc/TODO.md` and append a new section here. |
| **Risks / follow-ups** | LogRef dates inside checklist are narrative references, not git commits. When closing substantive items, cite PR/commit or `main.nf` example tag where applicable. |

---

## 2026-05-13 — Variation library TODO: English + expanded scope

| Field | Detail |
| --- | --- |
| **Date** | 2026-05-13 |
| **Summary** | Rewrite `doc/TODO.md` in English; merge new items from methodology / `3.variation-library` / `2.population-genetics` / `5.wheat-WGS-technology` |
| **Files changed** | `doc/TODO.md` |
| **Validation** | Manual diff vs user checklist: all new bullets present; prior `[x]` / `[ ]` states preserved; wiki links `[[4.process]]`, `[[5.assess]]`, `[[6.filter]]` and paths unchanged. |
| **Outcome** | Single English master list under “Variation library”; nested QC / AF / LD / IBS / filter sections expanded per user input. |
| **Risks / follow-ups** | Internal note IDs (`2.population-genetics`, etc.) are prose references; link to actual notes or commits when those docs exist in-repo. |

---

## 2026-05-13 — TODO.md supplemented from repo survey

| Field | Detail |
| --- | --- |
| **Date** | 2026-05-13 |
| **Summary** | Align `doc/TODO.md` with `workflow/Genetics` implementation (main router, stats/processor/hail, config gaps) |
| **Files changed** | `doc/TODO.md` |
| **Validation** | Manual: reviewed `main.nf`, `genotype/processor.nf`, `genotype/stats.nf`, `genotype/hail.nf`, `genotype/assess.nf`, `nextflow.config`; new bullets match current wiring. |
| **Outcome** | Added “Pipeline & repository snapshot” block; noted `test_common_thin`, test-only LD stats vs `v1_plink`, HAIL/kinship/PS router gap, README/assess parity, provenance + CI follow-ups. |
| **Risks / follow-ups** | Router and mod names may change; keep this section updated when `main.nf` branch logic changes. |

---

## 2026-05-13 — Reorder `doc/TODO.md` by Nextflow execution path

| Field | Detail |
| --- | --- |
| **Date** | 2026-05-13 |
| **Summary** | Structure TODO list along `main.nf` → `processor.nf` (`test_plink*`, `test_common_thin`, `v1_plink`) → `test_plink_stats` / `plink_stats` |
| **Files changed** | `doc/TODO.md` |
| **Validation** | Manual: compared section order to `workflow/Genetics/main.nf` branches and `test_plink_processor` / `test_plink_camp` / `test_common_thin_processor` / `plink_processor` + `test_plink_stats` / `plink_stats`; prior checklist items relocated, not dropped. |
| **Outcome** | Numbered sections (methodology → entry → preprocess → processor by mod → stats → production → filter → parallel modules → extended backlog → future); wiki anchors `[[4.process]]`, `[[5.assess]]`, `[[6.filter]]` kept on relevant sections. |
| **Risks / follow-ups** | If `stats.nf` process order changes, update **§5** ordering note. |

---

## 2026-05-13 — Re-run `test_thin` and `test_common_thin` with merged `process_dir` reuse

| Field | Detail |
| --- | --- |
| **Date** | 2026-05-13 |
| **Summary** | Restats for `test_thin` / `test_common_thin` using `--process_dir` with existing `A_test.plink2` … `Others_test.plink2` (skip thin/merge; rebuild `mk_plink_basic_info` + LD + `test_plink_stats`) |
| **Files changed** | `doc/NF_CMD.md`, `doc/TODO_PROGRESS_LOG.md` (this entry) |
| **Validation** | Ran Nextflow from conda env `run`: (1) `08stats.genome/11run_test_thin_restats` — exit 0, 52 succeeded, 21m 50s, log shows merged reuse; (2) `08stats.genome/12run_test_common_thin_restats` — exit 0, 52 succeeded, 1m 38s. |
| **Outcome** | Stats/plots republished under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/stats/test_thin` and `.../stats/test_common_thin` per existing `publishDir` wiring. |
| **Risks / follow-ups** | Nextflow warned Graphviz missing for DAG render (optional). Reuse same launch dir with `-resume` for retries. |

---

## 2026-05-13 — Mini-workflow `tmp/ld_plots_redraw.nf` for LD-only stats

| Field | Detail |
| --- | --- |
| **Date** | 2026-05-13 |
| **Summary** | Add `workflow/Genetics/tmp/ld_plots_redraw.nf` to re-run `variant_ld_decay_plot` + `variant_ld_crosschr_plot` from existing `process/<mod>/variant/*.vcor`; document commands in `doc/NF_CMD.md`; agent rules: prefer mini-NF over ad-hoc Python for partial reruns |
| **Files changed** | `workflow/Genetics/tmp/ld_plots_redraw.nf`, `workflow/Genetics/tmp/README.md`, `.cursor/rules/workstation-core.mdc`, `.cursor/rules/workstation-nextflow.mdc`, `.cursor/rules/workstation-python.mdc`, `doc/NF_CMD.md`, `doc/TODO_PROGRESS_LOG.md` (this entry); removed top-level `workflow/Genetics/ld_plots_redraw.nf` in favor of `tmp/` |
| **Validation** | `test_thin`: from `08stats.genome/14run_ld_plots_redraw`, `conda activate run`, `nextflow run .../ld_plots_redraw.nf ... --mod test_thin` — **exit 0**, **8** succeeded, **~15m 22s** (original run used repo root `ld_plots_redraw.nf`; equivalent command uses `tmp/ld_plots_redraw.nf`). Repeat from a **new** numbered folder for `--mod test_common_thin` (not run in this session). |
| **Outcome** | Partial reruns reuse the same `publishDir` layout as full `test_plink_stats`; eight LD-related tasks per mod run in parallel across subgenomes. |
| **Risks / follow-ups** | LD decay is CPU-heavy; avoid piping Nextflow stdout through `tail` if live monitoring is needed. |

---

## 2026-05-13 — Relocate LD mini-workflow to `tmp/`; chronological `NF_CMD.md`; rules sync

| Field | Detail |
| --- | --- |
| **Date** | 2026-05-13 |
| **Summary** | Move `ld_plots_redraw.nf` to `workflow/Genetics/tmp/`; document `-c` to `nextflow.config`; restructure `doc/NF_CMD.md` as time-ordered run log only; update `.cursor/rules` and `main.nf` help |
| **Files changed** | `workflow/Genetics/tmp/ld_plots_redraw.nf`, `workflow/Genetics/tmp/README.md`, removed `workflow/Genetics/ld_plots_redraw.nf`, `doc/NF_CMD.md`, `doc/TODO_PROGRESS_LOG.md` (this entry + LD redraw entry heading/files refresh), `.cursor/rules/workstation-core.mdc`, `.cursor/rules/workstation-nextflow.mdc`, `.cursor/rules/workstation-python.mdc`, `workflow/Genetics/main.nf` (help text) |
| **Validation** | Repo grep: no stale `workflow/Genetics/ld_plots_redraw.nf` references in rules; `include` path in `tmp/ld_plots_redraw.nf` points to `../genotype/stats.nf`. Nextflow not re-executed in this session after the move. |
| **Outcome** | Auxiliary entry lives under `tmp/`; config still loaded via absolute `-c .../workflow/Genetics/nextflow.config`; agent docs match. |
| **Risks / follow-ups** | Re-run `nextflow run .../tmp/ld_plots_redraw.nf` once from a vmap4 run folder to smoke-test after pull. |

---

## 2026-05-13 — Tier-1 assess debug mini-workflow (`tmp/assess_plink_debug.nf`)

| Field | Detail |
| --- | --- |
| **Date** | 2026-05-13 |
| **Summary** | §9 backlog tier-1: export narrow VCF from `*_test.plink2` for `test_thin` / `test_common_thin`, scrub PLINK `##chrSet` header for bcftools, run `quick_count` + `bcftools_qc_assess` + MAF-bin TSV; fix `assess.nf` conda + publishDir by `mod` + PLINK-safe `+fill-tags` / GQ handling |
| **Files changed** | `workflow/Genetics/tmp/assess_plink_debug.nf`, `workflow/Genetics/tmp/README.md`, `workflow/Genetics/genotype/assess.nf`, `doc/TODO.md`, `doc/NF_CMD.md`, `doc/TODO_PROGRESS_LOG.md` (this entry), `.cursor/rules/workstation-{core,nextflow,python}.mdc` |
| **Validation** | `conda activate run`: `21run_assess_debug_test_thin` — **exit 0**, **16** succeeded, **~2m 42s**; `22run_assess_debug_test_common_thin` — **exit 0**, **16** succeeded. |
| **Outcome** | Published under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/assess/test_thin/` and `.../assess/test_common_thin/` (`export/*.debug.vcf.gz`, `*.counts.tsv`, `*.maf_missing.tsv`, `*.gq_summary.tsv`, `info/*.mac_site_histogram.tsv`). |
| **Risks / follow-ups** | Slice is one representative chromosome per subgenome (not full genome); `dumpnice` still optional / script path may be absent; extend or wire into `main.nf` when router work is ready. |

---

## 2026-05-14 — Assess debug: PLINK2-native slice + Python infra plots

| Field | Detail |
| --- | --- |
| **Date** | 2026-05-14 |
| **Summary** | Replace bcftools/VCF export path in `tmp/assess_plink_debug.nf` with PLINK2 `--freq`/`--missing` on pfiles; MAF bins from `.afreq`; plots via `assess_slice.py` + `infra.utils.graph`; agent rules for assess + plotting |
| **Files changed** | `workflow/Genetics/genotype/assess.nf` (`plink2_assess_debug_slice`), `workflow/Genetics/tmp/assess_plink_debug.nf`, `workflow/Genetics/tmp/README.md`, `src/python/genetics/genomics/variant/assess_slice.py`, `src/python/infra/utils/graph.py` (`plot_bar_chart` dpi/labels), `.cursor/rules/workstation-core.mdc`, `.cursor/rules/workstation-python.mdc`, `doc/TODO.md`, `doc/NF_CMD.md`, `doc/TODO_PROGRESS_LOG.md` (this entry) |
| **Validation** | `python3` AST parse of `assess_slice.py`; `nextflow run .../tmp/assess_plink_debug.nf -preview` with real `output_dir`/`job`/`mod` from `/data/home/tusr1/01projects/vmap4/08stats.genome/23run_assess_plink2_debug_stub` — **exit 0**. Full `nextflow run` not re-executed (requires existing `process/<mod>/*_test.plink2`). |
| **Outcome** | Debug assess no longer depends on bcftools for the tier-1 slice; tables + PNGs follow stats-style infra helpers; workstation rules document PLINK2-first assess and plot layout under `assess/<mod>/`. |
| **Risks / follow-ups** | Re-run full assess for `test_thin` and `test_common_thin` from a numbered vmap4 run folder once processor outputs are present; remove or archive stale `export/*.vcf.gz` from older runs if disk hygiene matters. |

---

## 2026-05-14 — SUPP — Narrative supplement (2026-05-13 batch: restats, LD redraw, assess, rules)

**Purpose of this supplement:** Earlier table rows (above) captured outcomes briefly. This block adds **goals, rationale, working directories, pass/fail, output layout, and follow-ups** in prose. **Full command lines** for every Nextflow run are in **`doc/NF_CMD.md`** under the matching `### 2026-05-13 …` / `### 2026-05-14 …` headings—do not duplicate them here.

### Restats with merged `process_dir`

**Goal:** Refresh `test_thin` and `test_common_thin` stats (and downstream LD in that track) **without** re-thinning or re-merging, by pointing `--process_dir` at already merged `*_test.plink2` trees so the processor short-circuits to `mk_plink_basic_info` + LD + stats.

**Why:** Saves wall time and keeps test panels stable while validating the reuse branch used in production-like reruns.

**Runs (cwd only; commands → `doc/NF_CMD.md`):**
- `/data/home/tusr1/01projects/vmap4/08stats.genome/11run_test_thin_restats` — **pass** (exit 0, 52 processes, ~22 min); merged reuse visible in log.
- `/data/home/tusr1/01projects/vmap4/08stats.genome/12run_test_common_thin_restats` — **pass** (exit 0, 52 processes, ~2 min).

**Outputs / format:** Republished under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/stats/<mod>/` (plots, thresholds, info TSVs, logs) per existing `test_plink_stats` `publishDir` contract; same tree shape as a full stats run.

**Issues / next steps:** Optional Graphviz warning for DAG image only. Use `-resume` from the same launch dir on transient failures.

### LD-only redraw (`tmp/ld_plots_redraw.nf`)

**Goal:** Regenerate LD decay and cross-chromosome LD plots from existing `.vcor` inputs when plotting code or output prefixes change, without rerunning genotype.

**Why:** LD steps are expensive; isolating them in a mini-workflow avoids invalidating upstream work.

**Runs (cwd; commands → `doc/NF_CMD.md`):**
- `.../14run_ld_plots_redraw` — `test_thin` — **pass** (8 tasks, ~15 m).
- `.../16run_ld_plots_redraw_test_thin` — second `test_thin` pass after decay basename/prefix tweak — **pass** (8 tasks, ~16 m).
- `.../17run_ld_plots_redraw_test_common_thin` — `test_common_thin` — **pass** (8 tasks).

*(Note: the original LD-only redraw table row was conservative about `test_common_thin`; `doc/NF_CMD.md` records the `17run…` run as above.)*

**Outputs / format:** Same publish layout as stats LD processes (`stats/<mod>/plots` etc.). **Interpretation:** `publishDir` copies by pattern; confirm opened PNG basenames match the **current** run prefix to avoid stale files.

### Legacy tier-1 assess (VCF slice + bcftools)

**Goal:** Fast per-subgenome QC on one representative chromosome each for `test_thin` / `test_common_thin`, using a narrow VCF export so bcftools could fill/query `INFO` fields.

**Why (superseded):** bcftools path was needed for `+fill-tags` / per-site tables; later superseded by the 2026-05-14 PLINK2-native assess entry where pfiles suffice.

**Runs (cwd; commands → `doc/NF_CMD.md`):**
- `.../21run_assess_debug_test_thin` — **pass** (16 tasks, ~3 m).
- `.../22run_assess_debug_test_common_thin` — **pass** (16 tasks, ~1 m).

**Outputs / format:** Under `.../test_plink/assess/<mod>/`: `export/*.debug.vcf.gz`, root-level `*.counts.tsv`, `*.maf_missing.tsv`, `*.gq_summary.tsv`, `info/*.mac_site_histogram.tsv`. GQ summary often a placeholder when FORMAT/GQ absent in export.

### PLINK2-native assess + Python plots

**Goal:** Remove bcftools/VCF export from the debug assess mini-workflow where possible; surface MAF vs missing and MAF-bin summaries with **`infra.utils.graph`**-style figures.

**Why:** Aligns with PLINK2-first policy and reuses project plotting conventions (`assess_slice.py`).

**Code touched (high level):** `genotype/assess.nf` (`plink2_assess_debug_slice`), `tmp/assess_plink_debug.nf`, `assess_slice.py`, `plot_bar_chart` dpi/label behavior; rules for assess + Python.

**Validation:** Static parse of `assess_slice.py`; **Nextflow `-preview`** only from cwd `/data/home/tusr1/01projects/vmap4/08stats.genome/23run_assess_plink2_debug_stub` — **pass** (exit 0; preview does not execute tasks). **Full non-preview run:** **not done in that session** (needs existing `process/<mod>/*_test.plink2` under chosen `--output_dir` / `--job`).

**Next step:** Full run from a fresh numbered vmap4 folder, **without** `-preview`; log commands in `doc/NF_CMD.md` and narrative here.

**Expected outputs after full run:** `assess/<mod>/` with `.afreq`, `.vmiss`, TSVs, `plots/*.png`, `info/*.tsv`, `logs/*.log` — no new `export/*.vcf.gz` from this workflow version.

### Relocate LD mini-workflow + rules

**Goal:** Keep auxiliary entry scripts under `workflow/Genetics/tmp/` and keep agent docs aligned.

**Validation:** Static (grep, include paths). **No mandatory NF smoke** in that session.

**Next step:** Optional one-off `nextflow run .../tmp/ld_plots_redraw.nf` after large pulls.

---

## 2026-05-14 — Progress log policy (narrative vs `NF_CMD.md`)

**Goal:** `doc/TODO_PROGRESS_LOG.md` should carry **intent, diffs, rationale, test outcomes, failure analysis, remediation, and how artefacts are structured**; **`doc/NF_CMD.md`** remains the single place for **verbatim** runnable Nextflow/bash command blocks.

**What changed:** Header and `SUPP` in this file; `workstation-core.mdc`, `workstation-nextflow.mdc`, `todo-drive-close/SKILL.md` updated so agents record **cwd** per run in the progress log and **cite** `doc/NF_CMD.md` instead of duplicating multi-line commands.

**Validation:** Doc-only; no pipeline execution.

---

## 2026-05-14 — Refresh `workflow/Genetics/README.md` (todo-drive-close batch)

**Goal:** Close **`doc/TODO.md`** §2 “Docs drift” by replacing the outdated one-line `README` with a concise operator-facing description of the real **`main.nf`** router, optional `process_dir` / `camp`, test output tree, `tmp/` auxiliary entries, and workstation run policy (cwd under vmap4 projects, conda `run`, `screen` for long jobs, read-only vmap4 inputs).

**Changes:** Rewrote **`workflow/Genetics/README.md`** (English tables + sections); flipped the matching checklist item in **`doc/TODO.md`** to `[x]` with `LogRef: 2026-05-14`. No Nextflow execution in this session—**no new block in `doc/NF_CMD.md`**.

**Validation:** Manual review of `main.nf` dispatch (`v1_plink`, `test_plink`/`test_thin`, camp mods, `test_common_thin`) and `nextflow.config` parameter names against the new text; no automated test run.

**Outputs / effect:** README now points readers to **`doc/NF_CMD.md`**, **`doc/TODO.md`**, **`tmp/README.md`**, and **`.cursor/rules/workstation-*.mdc`** for commands and policy instead of duplicating long command blocks.

**Risks / follow-ups:** Router-gap mods (`HAIL`, kinship, PS, GWAS, `database`) remain documented only as backlog in README + **`doc/TODO.md`** §2; refresh again when those branches land in `workflow { }`.

---

## 2026-05-15 — Full `assess_plink_debug.nf` runs (`test_thin`, `test_common_thin`)

**Goal:** Close the validation gap left after the 2026-05-14 assess refactor (preview + AST only): execute the PLINK2-native assess mini-workflow **without** `-preview` for both test mods so `assess_slice.py` plots and TSVs are materialized under the real `publishDir` tree.

**Why:** Confirms merged `*_test.plink2` under `process/<mod>/` are readable end-to-end and that Python plotting processes complete under conda `stats`.

**Runs (cwd; full commands in `doc/NF_CMD.md` → `### 2026-05-15 — Assess debug …`):**
1. `…/08stats.genome/24run_assess_plink2_full_test_thin` — `--mod test_thin` — **pass** (exit 0, 8 processes succeeded; wall ~65 s from launcher; Nextflow `succeedDuration` aggregate ~7m 25s in `.nextflow.log`).
2. `…/08stats.genome/25run_assess_plink2_full_test_common_thin` — `--mod test_common_thin` — **pass** (exit 0, 8 succeeded, succeedDuration ~1m 25s in trace).

**Outputs / format:** Under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/assess/<mod>/`: per-subgenome `*.assess.afreq`, `*.assess.vmiss`, `*.counts.tsv`, `*.mac_site_histogram.tsv`, `*.gq_summary.tsv` (placeholder); `plots/` PNGs (`*.maf_vs_fmiss.scatter.png`, `*.maf.dist.png`, `*.fmiss.dist.png`, `*.maf_bins.bar.png` per `test_*__{A,B,D,Others}` prefix); `info/` joined TSVs (`*.maf_miss.info.tsv`, `*.maf_miss.th.tsv`, `*.counts.echo.tsv`); `logs/` from plot processes.

**Issues:** None blocking. Optional Graphviz missing for DAG render (warn only).

**Risks / follow-ups:** `main.nf` still does not dispatch assess (§8 item remains “router integration optional”); tier-1 slice remains one representative chr per subgenome. Large §9 science backlog unchanged.

---

## 2026-05-16 — Consolidate `workflow/` documentation into `workflow/Genetics/README.md`

**Goal:** Keep **one** README under `workflow/` (`workflow/Genetics/README.md`); merge content from `workflow/Genetics/integrated/README.md` and `workflow/Genetics/tmp/README.md`, delete those files, and encode the rule in **`.cursor/rules/workstation-core.mdc`**, **`.cursor/rules/workstation-nextflow.mdc`**, and cross-references (`main.nf`, `nextflow.config`, `tmp/ld_plots_redraw.nf`, `.cursor/skills/todo-drive-close/SKILL.md`).

**Validation:** Grep confirms no remaining `README.md` under `workflow/` except `Genetics/README.md`; doc-only.

---

## 2026-05-16 — Narrow VCF/PLINK `params.mod` list in docs and schema

**Goal:** Align operator docs, CLI help, JSON parameter schema, and engineering checklists with **`workflow/Genetics/main.nf`**: the VCF/PLINK branch now dispatches only **`v1_plink`**, **`test_thin`**, **`test_camp`**, and **`test_common_thin`** (the `wheat_*` integrated branch is unchanged).

**Changes:** Updated **`workflow/Genetics/subworkflows/local/genetics_helpers.nf`** (`helpMessage`), **`workflow/Genetics/docs/GENETICS_WORKFLOW.md`** (router table + `camp` row + output-layout note), **`workflow/Genetics/nextflow_schema.json`** (`mod` description), **`doc/TODO.md`** §2 router bullets + §5 intro line, and **`doc/NF_CMD.md`** intro (historic command blocks kept; note explains **`test_plink`** / **`test_plink_camp`** → **`test_thin`** / **`test_camp`** for new runs).

**Validation:** Doc-only; manual grep of `main.nf` branches vs edited tables. No Nextflow execution in this session.

**Risks / follow-ups:** Frozen paths and `params.job` examples still use the **`test_plink`** folder name under `output_dir`; that is independent of `--mod` string.

---

## 2026-05-17 — Tier-1 assess mini-workflow revalidation (`todo-drive-close`, assess scope)

**Goal:** Re-run **`workflow/Genetics/tmp/assess_plink_debug.nf`** for **`test_thin`** and **`test_common_thin`** after recent repo activity so the PLINK2 slice + **`assess_slice.py`** path remains green. **`doc/TODO.md`** §8 assess debug and the tier-1 bullet under §9 were already `[x]`; this session did **not** flip additional checklist boxes (§9 singletons / depth / LD backlog remains open).

**Why:** `todo-drive-close` hygiene: smallest proof that merged `*_test.plink2` under `process/<mod>/` still feed representative-chr `--freq` / `--missing` and Python plots without code changes.

**Code touched:** None.

**Validation & runs:** Conda env **`run`**; full commands in **`doc/NF_CMD.md`** (`### 2026-05-17 — Assess debug …` for each mod).

| Mod | cwd | Result |
| --- | --- | --- |
| `test_thin` | `/data/home/tusr1/01projects/vmap4/08stats.genome/28run_assess_revalidate_test_thin` | **pass** — exit 0, 8 succeeded; wall ~1m 14s; aggregate `succeedDuration` ~9m 19s |
| `test_common_thin` | `/data/home/tusr1/01projects/vmap4/08stats.genome/29run_assess_revalidate_test_common_thin` | **pass** — exit 0, 8 succeeded; wall ~23s; aggregate `succeedDuration` ~1m 9s |

**Outputs:** Same publish layout as the **2026-05-15 full assess runs**: `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/assess/<mod>/` with `.afreq` / `.vmiss` / bin TSVs at assess root; `plots/*.png`, `info/*.tsv`, `logs/*.log` under the mod tree.

**Issues:** Graphviz missing for DAG (warn only).

**Follow-ups:** Optional **`main.nf`** router hook for assess remains backlog (§8); science items in §9 (singletons, depth CI, LD interpretation, etc.) unchanged.

---

## 2026-05-17 — §9 singleton / MAC tables from PLINK2 `.acount` (`todo-drive-close`, §9 backlog)

**Goal:** Close **`doc/TODO.md`** §9 “Singletons: counts and distribution” with an **implementable** tier: extend the existing assess mini-workflow to use PLINK2 **`--freq counts`** (`.acount` with `ALT_CTS` / `OBS_CT`), derive **MAC** and **MAF**, and publish **singleton fraction / MAC bucket** tables plus bar + histogram figures next to the existing MAF–missing assess plots—still on the **representative-chromosome slice** per subgenome (not genome-wide).

**Why:** Frequency-only `.afreq` cannot stratify exact MAC (needed for singleton counts); `--freq counts` stays on the pfile-first assess policy.

**Changes (repo):**
- **`workflow/Genetics/modules/local/genotype/assess.nf`:** `plink2_assess_debug_slice` now runs `--freq counts`, publishes `*.acount`, MAF-bin awk reads `ALT_CTS` / `OBS_CT`.
- **`src/python/genetics/genomics/variant/assess_slice.py`:** `ana_assess_plink_debug_slice(..., acount_path=...)`, merged **`REF_CTS`**, **`MAC`**, **`MAF`** with **`F_MISS`**; new artefacts `*.singleton_mac.summary.tsv`, `*.mac_category_counts.tsv`, plots `*.mac_category.bar.png`, `*.mac.dist.png`; MAF-bin plot title notes `.acount`.
- **`workflow/Genetics/tmp/assess_plink_debug.nf`:** Wire `.acount` into the plot process.
- **Docs / policy:** `workflow/Genetics/docs/GENETICS_WORKFLOW.md` (aux table), **`doc/TODO.md`** §8–§9 (singleton parent + tier-1 bullet text, LogRef 2026-05-17), **`.cursor/rules/workstation-core.mdc`** (assess + `--freq counts`).

**Validation & runs:** Conda **`run`**; verbatim commands in **`doc/NF_CMD.md`** (`### 2026-05-17 — Assess debug ... PLINK2 \`--freq counts\` + MAC / singleton summaries` for each mod).

| Mod | cwd | Result |
| --- | --- | --- |
| `test_thin` | `/data/home/tusr1/01projects/vmap4/08stats.genome/30run_assess_singleton_mac_test_thin` | **pass** — 8 succeeded |
| `test_common_thin` | `/data/home/tusr1/01projects/vmap4/08stats.genome/31run_assess_singleton_mac_test_common_thin` | **pass** — 8 succeeded |

**Outputs / layout:** Under `.../test_plink/assess/<mod>/`: root **`*.assess.acount`**, **`*.assess.vmiss`**, **`*.counts.tsv`**, **`*.mac_site_histogram.tsv`**; **`info/`** adds **`*.singleton_mac.summary.tsv`**, **`*.mac_category_counts.tsv`**, updated **`*.maf_miss.info.tsv`** (extra MAC columns), thresholds **`*.maf_miss.th.tsv`**; **`plots/`** adds **`*.mac_category.bar.png`**, **`*.mac.dist.png`**. Older **`*.assess.afreq`** files may still appear in the same folder from **pre-change** runs (not removed by this workflow).

**Issues:** Graphviz missing (warn only).

**Follow-ups (still open in §9):** “Effect on LD”, caller/`lib.gz` singleton **science**, genome-wide singleton stats, **`MAC < 10` het** vs missing, etc.

---

## 2026-05-17 — Subgenome-first artefact names; integrated path encodes PLINK source mod

**Goal:** Remove **`job`** / redundant **`mod`** prefixes from **assess** and **wheat-from-PLINK** file basenames (subgenome id first, e.g. `A.assess.*`, `A.pca.png`). Keep **folder** semantics: **`assess/<plink_mod>/`**, and **`integrated/<plink_source_mod>/<wheat_task_mod>/`** so two PLINK test mods no longer share one integrated directory.

**Changes:** `tmp/assess_plink_debug.nf` (plot channel id only); `wheat_integrated_study.nf` (prefix = subgenome id; `RUN_WHEAT_INTEGRATED` table-mode prefix = `params.mod` only); `utils.nf` `listMergedSubgenomeSnpQcPlotTuples(process_dir)`; `stats.nf` publishDir for **`plot_plink2_population_structure`** / **`report_plink2_tagsnp`**; `processor.nf` **`plink2_pca`** `--pca approx` (fixes invalid legacy `alleles vzs` on PLINK 2.0a; avoids deterministic GRM NaN on high-missing pairs); **`GENETICS_WORKFLOW.md`** output layout text.

**Ops:** Cleared **`…/test_plink/integrated/`** and assess+integrated targets for **`test_thin`** / **`test_common_thin`** before regenerate.

**Validation:** Conda **`run`**; commands in **`doc/NF_CMD.md`** (`### 2026-05-17 — Output naming refactor…`). Assess **pass** (32–33); wheat **pass** (39–40), **8** processes each.

**Outputs:** Assess e.g. `A.assess.maf.dist.png` under `assess/test_thin/`; wheat PCA e.g. `integrated/test_thin/wheat_pca_tsne/plots/A.pca.png` (and parallel `test_common_thin` tree).

---

## 2026-05-20 — t-SNE on PLINK2 PCs in `plot_population_structure`

**Goal:** Restore **t-SNE** alongside PCA for **`wheat_pca_tsne`**: run **sklearn** `TSNE` on the first *K* PLINK PC columns (same idea as before job/mod prefixes were dropped), save **`{subgenome}.tsne.tsv`** + **`{subgenome}.tsne.png`**, without replacing PLINK2 PCA itself.

**Changes:** `src/python/genetics/genomics/sample/pca_structure.py` (perplexity scaled to sample size; optional `n_jobs`); **`plot_plink2_population_structure`** gains **`cpus 8`** and passes **`wheat_tsne_*`** + **`task.cpus`**; **`conf/base.config`** / **`nextflow_schema.json`** for `wheat_tsne_n_input_pcs`, `wheat_tsne_max_iter`, `wheat_tsne_random_state`; **`GENETICS_WORKFLOW.md`** / **`.cursor/rules/workstation-python.mdc`** (t-SNE on PCs allowed for visualisation).

**Validation:** Conda **`stats`**, **`PYTHONPATH=src/python`**, smoke call on real **`A.pca.eigenvec`/`eigenval`** (`test_common_thin`), `tsne_max_iter=300` — **pass** (writes `*.tsne.png` / `*.tsne.tsv`). Full Nextflow **`wheat_integrated_from_plink`** re-run not repeated in this batch (same process entry; operator can `-resume` or rerun when convenient).

**Outputs:** Per subgenome under **`…/integrated/<plink_mod>/wheat_pca_tsne/{info,plots}/`**: adds **`*.tsne.tsv`**, **`*.tsne.png`** next to existing PCA artifacts.

---

## 2026-05-20 — SUPP — Full wheat PCA/t-SNE NF rerun (`test_thin`, `test_common_thin`)

**Goal:** Regenerate **`wheat_pca_tsne`** integrated artefacts after t-SNE code landed.

**Runs:** Conda **`run`**; commands in **`doc/NF_CMD.md`** (`### 2026-05-20 — Wheat PCA/t-SNE …`).

| Mod | cwd | Result |
| --- | --- | --- |
| `test_thin` | `…/08stats.genome/41run_wheat_pca_tsne_test_thin` | **pass** — 8 succeeded, ~12m 36s |
| `test_common_thin` | `…/08stats.genome/42run_wheat_pca_tsne_test_common_thin` | **pass** — 8 succeeded, ~2m 42s |

**Outputs:** `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/integrated/<mod>/wheat_pca_tsne/` — per subgenome **`A|B|D|Others`**: `*.pca.png`, `*.variance.png`, **`*.tsne.png`**, matching **`info/*.tsne.tsv`**.

---

## 2026-05-20 — Sample `Group` on PCA/t-SNE plots

**Goal:** Color PCA and t-SNE scatter plots by germplasm **`Group`** (same `anno_group` path as `ref_ibs` / `test_plink_stats`), and include `Group` in exported TSVs.

**Changes:** `pca_structure.py` (`group_file`, `_attach_sample_groups`); `graph.py` `plot_scatter_with_thresholds(..., group_col=)`; `anno.py` `save_tsv=False` for in-memory annotate; `results_io.load_plink2_eigenvec` maps `#IID` → `Sample`; `stats.nf` passes `${params.output_dir}/sample_group.txt`.

**Validation:** Conda **`stats`** smoke on `A.pca.eigenvec` — **pass** (7 group categories). Nextflow **`run`** env:

| Mod | cwd | Result |
| --- | --- | --- |
| `test_common_thin` | `…/44run_wheat_pca_grp_test_common_thin` | **pass** — 8 succeeded, ~13m 43s (`doc/NF_CMD.md` §2026-05-20 group coloring, `test_common_thin`) |
| `test_thin` | `…/43run_wheat_pca_grp_test_thin` | **in progress** — `screen wheat43`; prior agent interrupt left no run dirs |

**Outputs (when complete):** `…/integrated/<mod>/wheat_pca_tsne/info/*.pca.tsv` with columns `Sample`, `PC*`, `Group`; plots `*.pca.png` / `*.tsne.png` with legend by group.

---

## 2026-06-04 — Progress log headings: date + English title only

**Goal:** Drop opaque section IDs (`OPS-NF-001`, `TODO-ASSESS-009`, etc.) from `doc/TODO_PROGRESS_LOG.md` and agent rules; use `## YYYY-MM-DD — <short English title>` (and `SUPP` supplements) so entries are self-describing.

**Changes:** Retrofitted existing log headings and table **`Summary`** rows; updated **`.cursor/rules/workstation-core.mdc`** and **`.cursor/skills/todo-drive-close/SKILL.md`**. Checklist **`doc/TODO.md`** still uses **`LogRef:`** calendar dates only (unchanged).

**Validation:** Doc-only.

---

## 2026-06-09 — Remove `assess.nf` shim after processor/stats split

**Goal:** Delete the deprecated `workflow/Genetics/modules/local/genotype/assess.nf` re-export shim now that assess compute and plots live in `processor.nf` and `stats.nf`, and drop stale references from agent rules and the checklist.

**Files changed:** removed `modules/local/genotype/assess.nf`; updated `.cursor/rules/workstation-core.mdc`, `doc/TODO.md` (§8 assess bullet + closed progress_overview §2.3.3 parity item).

**Validation:** `nextflow run …/tmp/assess_plink_debug.nf -preview` (`test_thin`) — **pass** (exit 0); includes resolve to `processor.nf` + `stats.nf` only. No `include` of `assess.nf` remains in active workflows.

**Outcome:** Assess debug remains available via `tmp/assess_plink_debug.nf` and `tmp/assess_vcf_debug.nf`; publish layout unchanged under `assess/<mod>/`.

---

## 2026-06-09 — Move dynamic selection plots out of `plots.nf`

**Goal:** Remove the standalone `dynamic/plots.nf` orchestrator; colocate selection-statistic plot processes with the dynamic theme and align naming/publish paths with other `dynamic/` modules.

**Changes:** Deleted `modules/local/dynamic/plots.nf` (including its `workflow plots { }` block). Added `modules/local/dynamic/selection.nf` with `plot_fst_boxplot`, `plot_pi_distribution`, `plot_tajima_d_distribution`, and `plot_ideogram_density` (lowercase process names; `publishDir` under `{output_dir}/{job}/population_genetics/selection/plots/…`). Closed **`doc/TODO.md`** §8 `plots.nf` backlog item.

**Validation:** No active workflow `include`d `plots.nf`; **skipped** full Nextflow run (processes not yet wired from `main.nf` router).

**Outcome:** Call plot processes directly from future selection compute workflows via `include { … } from '…/selection.nf'`; no separate plots workflow file.

---

## 2026-06-09 — Merge assess publish paths into process/stats trees

**Goal:** Drop the separate `{job}/assess/{mod}/` output tree; publish assess compute like other genotype work and assess plots/tables like other stats (test vs production distinguished by dataset/`job`, not by folder name).

**Changes:** Updated `publishDir` in `processor.nf` (`plink2_assess_debug_slice`, `quick_count`, `bcftools_qc_assess`, `mk_vcftools_basic_info` logs) → `process/<mod>/{variant,info,logs}`; `stats.nf` (`assess_plink_debug_plots`, `dumpnice_vcf_qc_assess`) → `stats/<mod>/{plots,info,logs}`; `tmp/assess_vcf_debug.nf` default VCF export dir → `process/<mod>/export`. Docs/rules: `GENETICS_WORKFLOW.md`, `workstation-core.mdc`, `workstation-python.mdc`, `todo-drive-close` skill, `doc/TODO.md` §9 tier-1 bullet.

**Validation:** `nextflow run …/tmp/assess_plink_debug.nf -preview` for `test_thin` — **pass** (exit 0).

**Outcome:** Existing `/data1/.../test_plink/assess/` on disk is legacy; safe to remove after re-run. New assess debug runs write alongside main `process/` and `stats/` artefacts for the same `mod`.

---

## 2026-06-09 — Repo hygiene and pipeline guardrails (no env/packaging change)

**Goal:** Apply the agreed optimization batch without altering `environment_*.yml` or Python/conda dependency wiring (one-step `conda env create -f …` remains the install path; Hail and heavy tools stay in conda).

**Changes:**
- Replaced broken `resources` Windows symlink with `resources/README.md` + `transposon/` placeholder layout for WeaTE local assets.
- Rewrote root `README.md` as a portal (correct `modules/local/` paths, link to `GENETICS_WORKFLOW.md`; clarified `pip install -e .` registers `python_script` only).
- `main.nf`: unknown `params.mod` now **exit 1** with supported-mod list; `GENETICS_WORKFLOW.md` updated to match.
- Python: added `infra.utils.errors`; `io.py` loaders and `mac.py` raise on missing/invalid inputs (Nextflow sees non-zero exit); `graph.combine_plots` default output is relative `combined_plot.png`.
- Tests: `tests/test_io.py`, `tests/test_mac.py`, `pytest.ini`; CI `.github/workflows/ci.yml` (conda `stats` env from existing yml + `pytest`/`ruff`).

**Validation:** `conda run -n stats pytest tests/ -q` — **pass** (11 tests). `ruff check` on touched Python — **pass**.

**Outcome:** No change to conda YAML dependency lists; collaborators still use `conda env create -f environment_stats.yml` then `pip install -e .` as before.

---

## 2026-06-09 — P3: stats.nf split, tmp governance, wheat shim warnings

**Goal:** Continue repo optimization without changing conda/Python dependency wiring: make large Nextflow stats modules navigable, document `tmp/` auxiliary scripts, and surface `genetics.wheat` deprecation at import time.

**Changes:**
- Split `modules/local/genotype/stats.nf` (1013 lines) into themed libraries: `stats_variant.nf`, `stats_sample.nf`, `stats_assess.nf`, `stats_integrated.nf`, `stats_chr_report.nf`; `stats.nf` now holds composed workflows only (~125 lines).
- Updated `tmp/*.nf` and `wheat_integrated_study.nf` includes to reference the smallest sub-library; fixed duplicate `ch_gcount` / redundant `variant_mac_stats` call in `tmp/test_rebulld_lib_stats.nf`.
- Added `workflow/Genetics/tmp/README.md` (script index, lifecycle, promotion backlog).
- `genetics.wheat`: `_deprecation.py` + import-time `DeprecationWarning` on all shim modules; `GENETICS_WORKFLOW.md` wheat section now cites `genetics.genomics.*` as canonical.
- `variant_utils.py`: removed dead `if df is None` branches (loaders now raise).

**Validation:** `nextflow run …/tmp/assess_plink_debug.nf -preview` (`test_thin`) — **pass** (exit 0). `conda run -n stats pytest tests/ -q` — **pass** (11 tests).

**Outcome:** Partial reruns import targeted stats sub-files; full pipelines still use `stats.nf` workflows. `genetics.wheat` imports remain working but warn once per module.

---

## 2026-06-09 — Move genotype stats modules into `genotype/stats/`

**Goal:** Keep `modules/local/genotype/` navigable by grouping all stats workflows and process libraries under one subfolder.

**Changes:** Moved `stats.nf`, `stats_variant.nf`, `stats_sample.nf`, `stats_assess.nf`, `stats_integrated.nf`, `stats_chr_report.nf` → `modules/local/genotype/stats/`; updated `stats/stats.nf` to `include` from `../utils.nf`. Refreshed `include` paths in `plink_genotype_modes.nf`, `wheat_integrated_study.nf`, all `tmp/*.nf`, `GENETICS_WORKFLOW.md`, `tmp/README.md`, `doc/TODO.md`, and `.cursor/rules` (nextflow/python/core). Added `genotype/stats/README.md`.

**Validation:** `nextflow run …/tmp/assess_plink_debug.nf -preview` (`test_thin`) — **pass** (exit 0).

**Outcome:** Canonical stats path is `genotype/stats/`; composed workflows remain in `genotype/stats/stats.nf`.

---

## 2026-06-09 — Consolidate docs: one repo README, Genetics narrative in docs/

**Goal:** Enforce single portal `README.md` at repo root; all `workflow/Genetics/` operator text lives in `docs/GENETICS_WORKFLOW.md` only.

**Changes:** Deleted `workflow/Genetics/tmp/README.md`, `modules/local/genotype/stats/README.md`, and `resources/README.md`; merged their content into `GENETICS_WORKFLOW.md` (Genetics operator text) and root `README.md` (`resources/` layout). Updated `workstation-nextflow.mdc` to forbid nested READMEs under `workflow/Genetics/`.

**Validation:** N/A (doc-only).

**Outcome:** `tmp/` and `genotype/stats/` have no standalone README; use `GENETICS_WORKFLOW.md` for operator narrative.

---

## 2026-06-09 — Split genotype processor into `genotype/processor/`

**Goal:** Mirror the `stats/` refactor: move ~1500-line monolithic `processor.nf` into themed sub-libraries under `modules/local/genotype/processor/`, keep composed DSL2 workflows in a thin barrel file, and preserve the processor vs stats compute/plot boundary.

**Changes:** Removed flat `genotype/processor.nf`; added `processor/processor.nf` (workflows) plus `processor_vcf.nf`, `processor_test.nf`, `processor_plink2.nf`, `processor_legacy.nf`, `processor_depth.nf`, `processor_filter.nf`, `processor_assess.nf`. Updated `include` paths in `plink_genotype_modes.nf`, `wheat_integrated_study.nf` (direct `processor_plink2.nf`), `static/gwas.nf` (`processor_vcf.nf`), `tmp/assess_*.nf` (`processor_assess.nf`). Refreshed `GENETICS_WORKFLOW.md`, `.cursor/rules` (nextflow/core/python), `todo-drive-close` skill, and `doc/TODO.md` §3–§4 / boundary bullet.

**Validation:** `nextflow run …/tmp/assess_plink_debug.nf -preview` (`test_thin`, cwd `/data/home/tusr1/01projects/vmap4/00nf_preview`) — **pass** (exit 0). `nextflow run …/main.nf -preview` (`test_thin`, cwd `/data/home/tusr1/01projects/vmap4/00nf_preview2`) — **pass** (exit 0); processor + stats processes resolve from new paths.

**Outcome:** Canonical processor path is `genotype/processor/`; subworkflows import `processor/processor.nf` for full workflows; partial reruns and cross-theme imports use the smallest `processor_*.nf` file.

---

## 2026-06-09 — Move shared DSL helpers to `modules/local/infra/`

**Goal:** Relocate cross-cutting `genotype/utils.nf` (~729 lines) into a theme-neutral `infra/` folder under `modules/local/`, matching the `processor/` and `stats/` split pattern.

**Changes:** Removed `genotype/utils.nf`; added `infra/utils.nf` (barrel) plus `infra_tools.nf`, `infra_job_config.nf`, `infra_tiger.nf`, `infra_ref_v1.nf`, `infra_plink_reuse.nf`. Updated `include` paths in `processor/processor.nf`, `processor_vcf.nf`, `processor_depth.nf`, `stats/stats.nf`, `genetics_helpers.nf`, `wheat_integrated_study.nf`, `caller.nf`. Documented layout in `GENETICS_WORKFLOW.md` and `workstation-nextflow.mdc`; closed `doc/TODO.md` modules-index bullet.

**Validation:** `nextflow run …/tmp/assess_plink_debug.nf -preview` (`test_thin`, cwd `/data/home/tusr1/01projects/vmap4/00nf_preview3`) — **pass** (exit 0). `nextflow run …/main.nf -preview` (`test_thin`, cwd `/data/home/tusr1/01projects/vmap4/00nf_preview4`) — **pass** (exit 0).

**Outcome:** Shared helpers live under `modules/local/infra/`; genotype modules import the smallest `infra_*.nf` file. Python parity for ref maps remains in `src/python/infra/wheat/ref_v1.py`.

---

## 2026-06-09 — Split upstream calling/align and static GWAS modules

**Goal:** Continue module governance after `processor/` / `stats/` / `infra/` — relocate large standalone genotype entry scripts and align `static/gwas` with the documented `gwas/` layout.

**Changes:** Moved `caller.nf` (852 lines) → `genotype/calling/` (`caller.nf` + `caller_prep.nf` + `caller_fastcall.nf`); `align.nf` (390 lines) → `genotype/align/` (`align.nf` + `align_md5.nf` + `align_bwa.nf` + `align_transfer.nf`); flat `static/gwas.nf` → `static/gwas/` (`gwas.nf` workflow barrel + `gwas_legacy.nf`, `gwas_plink2.nf`, `gwas_plot.nf`). Consolidated duplicate `include` lines in `plink_genotype_modes.nf`. Updated `wheat_integrated_study.nf`, rules, and `GENETICS_WORKFLOW.md`.

**Validation:** `nextflow run …/main.nf -preview` (`test_thin`, cwd `/data/home/tusr1/01projects/vmap4/00nf_preview6`) — **pass** (exit 0). `nextflow run …/calling/caller.nf -preview --help` — **pass** (exit 0).

**Outcome:** Upstream FastCall3 and align tracks are under themed subfolders; wheat GWAS imports target `gwas_plink2.nf` / `gwas_plot.nf` directly. Composing `align` + `calling` into a `subworkflows/local/raw_data.nf` router remains open (`doc/TODO.md` §2).

---

## 2026-06-09 — Hail/database split and raw-data subworkflow

**Goal:** Finish genotype root cleanup (`hail.nf`, `database.nf` still flat) and compose upstream align + calling under one subworkflow file for future `main.nf` routing.

**Changes:** Split `hail.nf` → `genotype/hail/` (`hail.nf`, `hail_io.nf`, `hail_stats.nf`, `hail_gwas.nf`); `database.nf` → `genotype/database/` (`database.nf`, `database_build.nf`, `database_annotate.nf`). Extracted `RUN_ALIGN_USB_TRANSFER` workflow in `align/align.nf`. Added `subworkflows/local/raw_data_upstream.nf` (includes `RUN_ALIGN_USB_TRANSFER`, `run_FastCall3`). Updated `plink_genotype_modes.nf`, `GENETICS_WORKFLOW.md`, `doc/TODO.md` §2.

**Validation:** `nextflow run …/main.nf -preview` (`test_thin`, cwd `/data/home/tusr1/01projects/vmap4/00nf_preview9`) — **pass** (exit 0).

**Outcome:** `genotype/` root has no monolithic NF modules left except none—all themed subfolders. `static/database.nf` (DBone / phenotype ops) merge with genotype database remains open (`doc/TODO.md` §2).

---

## 2026-06-09 — Promote tmp partial reruns to subworkflows/local

**Goal:** Move channel-wiring logic out of `tmp/*.nf` into reusable named workflows under `subworkflows/local/`; slim `plink_genotype_modes.nf` to PLINK-only includes.

**Changes:** Added `partial_assess.nf` (`RUN_ASSESS_PLINK_DEBUG`, `RUN_ASSESS_VCF_DEBUG`), `partial_stats.nf` (LD/MAC/chr/rebuild partial workflows), `analysis_extensions.nf` (router-gap `DATABASE`, `KINSHIP`, `PS`, `GWAS`, `HAIL`). Rewrote eight `tmp/*.nf` files as thin `workflow { RUN_*() }` wrappers. Removed unused includes from `plink_genotype_modes.nf`; consolidated `main.nf` and `genetics_helpers.nf` includes. Updated `GENETICS_WORKFLOW.md`, `workstation-nextflow.mdc`, `doc/TODO.md` §2.

**Validation:** `nextflow run …/main.nf -preview` (`test_thin`) — **pass**. `tmp/assess_plink_debug.nf -preview` and `tmp/ld_plots_redraw.nf -preview` — **pass** (exit 0).

**Outcome:** `tmp/` is launch-only; partial-rerun logic is testable/reusable from `subworkflows/local/partial_*.nf` without duplicating module `include` paths.

---

## 2026-06-09 — Remove `workflow/Genetics/tmp/`; layer `subworkflows/local/`

**Goal:** Keep `modules/local/`; relocate scratch/ops scripts beside `local/`; delete redundant genetics tmp entry scripts; organize composed subworkflows into themed subfolders.

**Changes:**

- Deleted **`workflow/Genetics/tmp/`** genetics scripts (assess, LD/MAC/chr/rebuild, wheat-from-plink wrappers) whose logic already lives in **`subworkflows/local/partial/`** and **`subworkflows/local/wheat/`**.
- Moved FTP ops to **`subworkflows/tmp/ops/`** (`gsa_ftp_upload.nf`, `gwh_ftp_upload.nf`).
- Layered **`subworkflows/local/`** into `entry/`, `plink/`, `wheat/`, `upstream/`, `partial/`.
- Added **`entry/partial_router.nf`** as the single partial-rerun launcher (`--partial_task`).
- Updated **`main.nf`** includes, **`GENETICS_WORKFLOW.md`**, **`nextflow_schema.json`** (`partial_task`), **`.cursor/rules`** (core/nextflow/python), **`todo-drive-close`** skill, **`doc/TODO.md`** §2/§8/§9, and **`doc/NF_CMD.md`** migration note (historic `tmp/` blocks retained).

**Validation:** cwd `/data/home/tusr1/01projects/vmap4/00nf_preview11` — `nextflow run …/main.nf -preview` (`test_thin`) **pass**; `partial_router.nf -preview --partial_task assess_plink` **pass** (exit 0).

**Outcome:** Partial reruns use one router entry; `subworkflows/local/` layout matches functional areas; only FTP helpers remain under `subworkflows/tmp/`.
