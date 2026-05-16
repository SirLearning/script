# TODO progress log

Append-only audit log for completed TODO / ops work (per `.cursor/rules/workstation-core.mdc`). **Do not delete or rewrite** existing body text; add new material at the **end** only (including `SUPP` sections). Each entry should read like an **engineering report**: goals, what changed and why, validation (**pass / fail / blocked**), root causes and fixes, next steps, and how outputs are laid out (formats, `publishDir`, prefixes). **Full command lines** belong in **`doc/NF_CMD.md`**; here, cite that file by heading and record **working directory(ies)** per run. English; prose or tables are both fine.

---

## 2026-05-13 — TODO-INIT-001 — Establish master TODO + progress log in repo

| Field | Detail |
| --- | --- |
| **Date** | 2026-05-13 |
| **TODO ID / title** | TODO-INIT-001 — Import vmap4 workflow methodology/process/assess/filter checklist |
| **Files changed** | `doc/TODO.md` (created), `doc/TODO_PROGRESS_LOG.md` (created) |
| **Validation** | Manual: checklist matches user-supplied snapshot; open `[ ]` and done `[x]` preserved. |
| **Outcome** | Master list versioned in repo; future completions should flip boxes in `doc/TODO.md` and append a new section here. |
| **Risks / follow-ups** | LogRef dates inside checklist are narrative references, not git commits. When closing substantive items, cite PR/commit or `main.nf` example tag where applicable. |

---

## 2026-05-13 — TODO-DOC-002 — Variation library TODO: English + expanded scope

| Field | Detail |
| --- | --- |
| **Date** | 2026-05-13 |
| **TODO ID / title** | TODO-DOC-002 — Rewrite `doc/TODO.md` in English; merge new items from methodology / `3.variation-library` / `2.population-genetics` / `5.wheat-WGS-technology` |
| **Files changed** | `doc/TODO.md` |
| **Validation** | Manual diff vs user checklist: all new bullets present; prior `[x]` / `[ ]` states preserved; wiki links `[[4.process]]`, `[[5.assess]]`, `[[6.filter]]` and paths unchanged. |
| **Outcome** | Single English master list under “Variation library”; nested QC / AF / LD / IBS / filter sections expanded per user input. |
| **Risks / follow-ups** | Internal note IDs (`2.population-genetics`, etc.) are prose references; link to actual notes or commits when those docs exist in-repo. |

---

## 2026-05-13 — TODO-DOC-003 — TODO.md supplemented from repo survey

| Field | Detail |
| --- | --- |
| **Date** | 2026-05-13 |
| **TODO ID / title** | TODO-DOC-003 — Align `doc/TODO.md` with `workflow/Genetics` implementation (main router, stats/processor/hail, config gaps) |
| **Files changed** | `doc/TODO.md` |
| **Validation** | Manual: reviewed `main.nf`, `genotype/processor.nf`, `genotype/stats.nf`, `genotype/hail.nf`, `genotype/assess.nf`, `nextflow.config`; new bullets match current wiring. |
| **Outcome** | Added “Pipeline & repository snapshot” block; noted `test_common_thin`, test-only LD stats vs `v1_plink`, HAIL/kinship/PS router gap, README/assess parity, provenance + CI follow-ups. |
| **Risks / follow-ups** | Router and mod names may change; keep this section updated when `main.nf` branch logic changes. |

---

## 2026-05-13 — TODO-DOC-004 — Reorder `doc/TODO.md` by Nextflow execution path

| Field | Detail |
| --- | --- |
| **Date** | 2026-05-13 |
| **TODO ID / title** | TODO-DOC-004 — Structure TODO list along `main.nf` → `processor.nf` (`test_plink*`, `test_common_thin`, `v1_plink`) → `test_plink_stats` / `plink_stats` |
| **Files changed** | `doc/TODO.md` |
| **Validation** | Manual: compared section order to `workflow/Genetics/main.nf` branches and `test_plink_processor` / `test_plink_camp` / `test_common_thin_processor` / `plink_processor` + `test_plink_stats` / `plink_stats`; prior checklist items relocated, not dropped. |
| **Outcome** | Numbered sections (methodology → entry → preprocess → processor by mod → stats → production → filter → parallel modules → extended backlog → future); wiki anchors `[[4.process]]`, `[[5.assess]]`, `[[6.filter]]` kept on relevant sections. |
| **Risks / follow-ups** | If `stats.nf` process order changes, update **§5** ordering note. |

---

## 2026-05-13 — OPS-NF-001 — Re-run `test_thin` and `test_common_thin` with merged `process_dir` reuse

| Field | Detail |
| --- | --- |
| **Date** | 2026-05-13 |
| **TODO ID / title** | OPS-NF-001 — Restats for `test_thin` / `test_common_thin` using `--process_dir` with existing `A_test.plink2` … `Others_test.plink2` (skip thin/merge; rebuild `mk_plink_basic_info` + LD + `test_plink_stats`) |
| **Files changed** | `doc/NF_CMD.md`, `doc/TODO_PROGRESS_LOG.md` (this entry) |
| **Validation** | Ran Nextflow from conda env `run`: (1) `08stats.genome/11run_test_thin_restats` — exit 0, 52 succeeded, 21m 50s, log shows merged reuse; (2) `08stats.genome/12run_test_common_thin_restats` — exit 0, 52 succeeded, 1m 38s. |
| **Outcome** | Stats/plots republished under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/stats/test_thin` and `.../stats/test_common_thin` per existing `publishDir` wiring. |
| **Risks / follow-ups** | Nextflow warned Graphviz missing for DAG render (optional). Reuse same launch dir with `-resume` for retries. |

---

## 2026-05-13 — OPS-NF-002 — Mini-workflow `tmp/ld_plots_redraw.nf` for LD-only stats

| Field | Detail |
| --- | --- |
| **Date** | 2026-05-13 |
| **TODO ID / title** | OPS-NF-002 — Add `workflow/Genetics/tmp/ld_plots_redraw.nf` to re-run `variant_ld_decay_plot` + `variant_ld_crosschr_plot` from existing `process/<mod>/variant/*.vcor`; document commands in `doc/NF_CMD.md`; agent rules: prefer mini-NF over ad-hoc Python for partial reruns |
| **Files changed** | `workflow/Genetics/tmp/ld_plots_redraw.nf`, `workflow/Genetics/tmp/README.md`, `.cursor/rules/workstation-core.mdc`, `.cursor/rules/workstation-nextflow.mdc`, `.cursor/rules/workstation-python.mdc`, `doc/NF_CMD.md`, `doc/TODO_PROGRESS_LOG.md` (this entry); removed top-level `workflow/Genetics/ld_plots_redraw.nf` in favor of `tmp/` |
| **Validation** | `test_thin`: from `08stats.genome/14run_ld_plots_redraw`, `conda activate run`, `nextflow run .../ld_plots_redraw.nf ... --mod test_thin` — **exit 0**, **8** succeeded, **~15m 22s** (original run used repo root `ld_plots_redraw.nf`; equivalent command uses `tmp/ld_plots_redraw.nf`). Repeat from a **new** numbered folder for `--mod test_common_thin` (not run in this session). |
| **Outcome** | Partial reruns reuse the same `publishDir` layout as full `test_plink_stats`; eight LD-related tasks per mod run in parallel across subgenomes. |
| **Risks / follow-ups** | LD decay is CPU-heavy; avoid piping Nextflow stdout through `tail` if live monitoring is needed. |

---

## 2026-05-13 — OPS-NF-003 — Relocate LD mini-workflow to `tmp/`; chronological `NF_CMD.md`; rules sync

| Field | Detail |
| --- | --- |
| **Date** | 2026-05-13 |
| **TODO ID / title** | OPS-NF-003 — Move `ld_plots_redraw.nf` to `workflow/Genetics/tmp/`; document `-c` to `nextflow.config`; restructure `doc/NF_CMD.md` as time-ordered run log only; update `.cursor/rules` and `main.nf` help |
| **Files changed** | `workflow/Genetics/tmp/ld_plots_redraw.nf`, `workflow/Genetics/tmp/README.md`, removed `workflow/Genetics/ld_plots_redraw.nf`, `doc/NF_CMD.md`, `doc/TODO_PROGRESS_LOG.md` (this entry + OPS-NF-002 heading/files refresh), `.cursor/rules/workstation-core.mdc`, `.cursor/rules/workstation-nextflow.mdc`, `.cursor/rules/workstation-python.mdc`, `workflow/Genetics/main.nf` (help text) |
| **Validation** | Repo grep: no stale `workflow/Genetics/ld_plots_redraw.nf` references in rules; `include` path in `tmp/ld_plots_redraw.nf` points to `../genotype/stats.nf`. Nextflow not re-executed in this session after the move. |
| **Outcome** | Auxiliary entry lives under `tmp/`; config still loaded via absolute `-c .../workflow/Genetics/nextflow.config`; agent docs match. |
| **Risks / follow-ups** | Re-run `nextflow run .../tmp/ld_plots_redraw.nf` once from a vmap4 run folder to smoke-test after pull. |

---

## 2026-05-13 — OPS-ASSESS-001 — Tier-1 assess debug mini-workflow (`tmp/assess_plink_debug.nf`)

| Field | Detail |
| --- | --- |
| **Date** | 2026-05-13 |
| **TODO ID / title** | OPS-ASSESS-001 — §9 backlog tier-1: export narrow VCF from `*_test.plink2` for `test_thin` / `test_common_thin`, scrub PLINK `##chrSet` header for bcftools, run `quick_count` + `bcftools_qc_assess` + MAF-bin TSV; fix `assess.nf` conda + publishDir by `mod` + PLINK-safe `+fill-tags` / GQ handling |
| **Files changed** | `workflow/Genetics/tmp/assess_plink_debug.nf`, `workflow/Genetics/tmp/README.md`, `workflow/Genetics/genotype/assess.nf`, `doc/TODO.md`, `doc/NF_CMD.md`, `doc/TODO_PROGRESS_LOG.md` (this entry), `.cursor/rules/workstation-{core,nextflow,python}.mdc` |
| **Validation** | `conda activate run`: `21run_assess_debug_test_thin` — **exit 0**, **16** succeeded, **~2m 42s**; `22run_assess_debug_test_common_thin` — **exit 0**, **16** succeeded. |
| **Outcome** | Published under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/assess/test_thin/` and `.../assess/test_common_thin/` (`export/*.debug.vcf.gz`, `*.counts.tsv`, `*.maf_missing.tsv`, `*.gq_summary.tsv`, `info/*.mac_site_histogram.tsv`). |
| **Risks / follow-ups** | Slice is one representative chromosome per subgenome (not full genome); `dumpnice` still optional / script path may be absent; extend or wire into `main.nf` when router work is ready. |

---

## 2026-05-14 — OPS-ASSESS-002 — Assess debug: PLINK2-native slice + Python infra plots

| Field | Detail |
| --- | --- |
| **Date** | 2026-05-14 |
| **TODO ID / title** | OPS-ASSESS-002 — Replace bcftools/VCF export path in `tmp/assess_plink_debug.nf` with PLINK2 `--freq`/`--missing` on pfiles; MAF bins from `.afreq`; plots via `assess_slice.py` + `infra.utils.graph`; agent rules for assess + plotting |
| **Files changed** | `workflow/Genetics/genotype/assess.nf` (`plink2_assess_debug_slice`), `workflow/Genetics/tmp/assess_plink_debug.nf`, `workflow/Genetics/tmp/README.md`, `src/python/genetics/genomics/variant/assess_slice.py`, `src/python/infra/utils/graph.py` (`plot_bar_chart` dpi/labels), `.cursor/rules/workstation-core.mdc`, `.cursor/rules/workstation-python.mdc`, `doc/TODO.md`, `doc/NF_CMD.md`, `doc/TODO_PROGRESS_LOG.md` (this entry) |
| **Validation** | `python3` AST parse of `assess_slice.py`; `nextflow run .../tmp/assess_plink_debug.nf -preview` with real `output_dir`/`job`/`mod` from `/data/home/tusr1/01projects/vmap4/08stats.genome/23run_assess_plink2_debug_stub` — **exit 0**. Full `nextflow run` not re-executed (requires existing `process/<mod>/*_test.plink2`). |
| **Outcome** | Debug assess no longer depends on bcftools for the tier-1 slice; tables + PNGs follow stats-style infra helpers; workstation rules document PLINK2-first assess and plot layout under `assess/<mod>/`. |
| **Risks / follow-ups** | Re-run full assess for `test_thin` and `test_common_thin` from a numbered vmap4 run folder once processor outputs are present; remove or archive stale `export/*.vcf.gz` from older runs if disk hygiene matters. |

---

## 2026-05-14 — SUPP — Narrative supplement (OPS-NF-001, OPS-NF-002, OPS-ASSESS-001, OPS-ASSESS-002, OPS-NF-003)

**Purpose of this supplement:** Earlier table rows (above) captured outcomes briefly. This block adds **goals, rationale, working directories, pass/fail, output layout, and follow-ups** in prose. **Full command lines** for every Nextflow run are in **`doc/NF_CMD.md`** under the matching `### 2026-05-13 …` / `### 2026-05-14 …` headings—do not duplicate them here.

### OPS-NF-001 — Restats with merged `process_dir`

**Goal:** Refresh `test_thin` and `test_common_thin` stats (and downstream LD in that track) **without** re-thinning or re-merging, by pointing `--process_dir` at already merged `*_test.plink2` trees so the processor short-circuits to `mk_plink_basic_info` + LD + stats.

**Why:** Saves wall time and keeps test panels stable while validating the reuse branch used in production-like reruns.

**Runs (cwd only; commands → `doc/NF_CMD.md`):**
- `/data/home/tusr1/01projects/vmap4/08stats.genome/11run_test_thin_restats` — **pass** (exit 0, 52 processes, ~22 min); merged reuse visible in log.
- `/data/home/tusr1/01projects/vmap4/08stats.genome/12run_test_common_thin_restats` — **pass** (exit 0, 52 processes, ~2 min).

**Outputs / format:** Republished under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/stats/<mod>/` (plots, thresholds, info TSVs, logs) per existing `test_plink_stats` `publishDir` contract; same tree shape as a full stats run.

**Issues / next steps:** Optional Graphviz warning for DAG image only. Use `-resume` from the same launch dir on transient failures.

### OPS-NF-002 — LD-only redraw (`tmp/ld_plots_redraw.nf`)

**Goal:** Regenerate LD decay and cross-chromosome LD plots from existing `.vcor` inputs when plotting code or output prefixes change, without rerunning genotype.

**Why:** LD steps are expensive; isolating them in a mini-workflow avoids invalidating upstream work.

**Runs (cwd; commands → `doc/NF_CMD.md`):**
- `.../14run_ld_plots_redraw` — `test_thin` — **pass** (8 tasks, ~15 m).
- `.../16run_ld_plots_redraw_test_thin` — second `test_thin` pass after decay basename/prefix tweak — **pass** (8 tasks, ~16 m).
- `.../17run_ld_plots_redraw_test_common_thin` — `test_common_thin` — **pass** (8 tasks).

*(Note: the original OPS-NF-002 table row was conservative about `test_common_thin`; `doc/NF_CMD.md` records the `17run…` run as above.)*

**Outputs / format:** Same publish layout as stats LD processes (`stats/<mod>/plots` etc.). **Interpretation:** `publishDir` copies by pattern; confirm opened PNG basenames match the **current** run prefix to avoid stale files.

### OPS-ASSESS-001 — Legacy tier-1 assess (VCF slice + bcftools)

**Goal:** Fast per-subgenome QC on one representative chromosome each for `test_thin` / `test_common_thin`, using a narrow VCF export so bcftools could fill/query `INFO` fields.

**Why (superseded):** bcftools path was needed for `+fill-tags` / per-site tables; later superseded by PLINK2-native slice (OPS-ASSESS-002) where pfiles suffice.

**Runs (cwd; commands → `doc/NF_CMD.md`):**
- `.../21run_assess_debug_test_thin` — **pass** (16 tasks, ~3 m).
- `.../22run_assess_debug_test_common_thin` — **pass** (16 tasks, ~1 m).

**Outputs / format:** Under `.../test_plink/assess/<mod>/`: `export/*.debug.vcf.gz`, root-level `*.counts.tsv`, `*.maf_missing.tsv`, `*.gq_summary.tsv`, `info/*.mac_site_histogram.tsv`. GQ summary often a placeholder when FORMAT/GQ absent in export.

### OPS-ASSESS-002 — PLINK2-native assess + Python plots

**Goal:** Remove bcftools/VCF export from the debug assess mini-workflow where possible; surface MAF vs missing and MAF-bin summaries with **`infra.utils.graph`**-style figures.

**Why:** Aligns with PLINK2-first policy and reuses project plotting conventions (`assess_slice.py`).

**Code touched (high level):** `genotype/assess.nf` (`plink2_assess_debug_slice`), `tmp/assess_plink_debug.nf`, `assess_slice.py`, `plot_bar_chart` dpi/label behavior; rules for assess + Python.

**Validation:** Static parse of `assess_slice.py`; **Nextflow `-preview`** only from cwd `/data/home/tusr1/01projects/vmap4/08stats.genome/23run_assess_plink2_debug_stub` — **pass** (exit 0; preview does not execute tasks). **Full non-preview run:** **not done in that session** (needs existing `process/<mod>/*_test.plink2` under chosen `--output_dir` / `--job`).

**Next step:** Full run from a fresh numbered vmap4 folder, **without** `-preview`; log commands in `doc/NF_CMD.md` and narrative here.

**Expected outputs after full run:** `assess/<mod>/` with `.afreq`, `.vmiss`, TSVs, `plots/*.png`, `info/*.tsv`, `logs/*.log` — no new `export/*.vcf.gz` from this workflow version.

### OPS-NF-003 — Relocate LD mini-workflow + rules

**Goal:** Keep auxiliary entry scripts under `workflow/Genetics/tmp/` and keep agent docs aligned.

**Validation:** Static (grep, include paths). **No mandatory NF smoke** in that session.

**Next step:** Optional one-off `nextflow run .../tmp/ld_plots_redraw.nf` after large pulls.

---

## 2026-05-14 — DOC-PROGLOG-001 — Progress log policy (narrative vs `NF_CMD.md`)

**Goal:** `doc/TODO_PROGRESS_LOG.md` should carry **intent, diffs, rationale, test outcomes, failure analysis, remediation, and how artefacts are structured**; **`doc/NF_CMD.md`** remains the single place for **verbatim** runnable Nextflow/bash command blocks.

**What changed:** Header and `SUPP` in this file; `workstation-core.mdc`, `workstation-nextflow.mdc`, `todo-drive-close/SKILL.md` updated so agents record **cwd** per run in the progress log and **cite** `doc/NF_CMD.md` instead of duplicating multi-line commands.

**Validation:** Doc-only; no pipeline execution.

---

## 2026-05-14 — TODO-DOC-README — Refresh `workflow/Genetics/README.md` (todo-drive-close batch)

**Goal:** Close **`doc/TODO.md`** §2 “Docs drift” by replacing the outdated one-line `README` with a concise operator-facing description of the real **`main.nf`** router, optional `process_dir` / `camp`, test output tree, `tmp/` auxiliary entries, and workstation run policy (cwd under vmap4 projects, conda `run`, `screen` for long jobs, read-only vmap4 inputs).

**Changes:** Rewrote **`workflow/Genetics/README.md`** (English tables + sections); flipped the matching checklist item in **`doc/TODO.md`** to `[x]` with `LogRef: 2026-05-14`. No Nextflow execution in this session—**no new block in `doc/NF_CMD.md`**.

**Validation:** Manual review of `main.nf` dispatch (`v1_plink`, `test_plink`/`test_thin`, camp mods, `test_common_thin`) and `nextflow.config` parameter names against the new text; no automated test run.

**Outputs / effect:** README now points readers to **`doc/NF_CMD.md`**, **`doc/TODO.md`**, **`tmp/README.md`**, and **`.cursor/rules/workstation-*.mdc`** for commands and policy instead of duplicating long command blocks.

**Risks / follow-ups:** Router-gap mods (`HAIL`, kinship, PS, GWAS, `database`) remain documented only as backlog in README + **`doc/TODO.md`** §2; refresh again when those branches land in `workflow { }`.

---

## 2026-05-15 — OPS-ASSESS-003 — Full `assess_plink_debug.nf` runs (`test_thin`, `test_common_thin`)

**Goal:** Close the validation gap left after OPS-ASSESS-002 (preview + AST only): execute the PLINK2-native assess mini-workflow **without** `-preview` for both test mods so `assess_slice.py` plots and TSVs are materialized under the real `publishDir` tree.

**Why:** Confirms merged `*_test.plink2` under `process/<mod>/` are readable end-to-end and that Python plotting processes complete under conda `stats`.

**Runs (cwd; full commands in `doc/NF_CMD.md` → `### 2026-05-15 — Assess debug …`):**
1. `…/08stats.genome/24run_assess_plink2_full_test_thin` — `--mod test_thin` — **pass** (exit 0, 8 processes succeeded; wall ~65 s from launcher; Nextflow `succeedDuration` aggregate ~7m 25s in `.nextflow.log`).
2. `…/08stats.genome/25run_assess_plink2_full_test_common_thin` — `--mod test_common_thin` — **pass** (exit 0, 8 succeeded, succeedDuration ~1m 25s in trace).

**Outputs / format:** Under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/assess/<mod>/`: per-subgenome `*.assess.afreq`, `*.assess.vmiss`, `*.counts.tsv`, `*.mac_site_histogram.tsv`, `*.gq_summary.tsv` (placeholder); `plots/` PNGs (`*.maf_vs_fmiss.scatter.png`, `*.maf.dist.png`, `*.fmiss.dist.png`, `*.maf_bins.bar.png` per `test_*__{A,B,D,Others}` prefix); `info/` joined TSVs (`*.maf_miss.info.tsv`, `*.maf_miss.th.tsv`, `*.counts.echo.tsv`); `logs/` from plot processes.

**Issues:** None blocking. Optional Graphviz missing for DAG render (warn only).

**Risks / follow-ups:** `main.nf` still does not dispatch assess (§8 item remains “router integration optional”); tier-1 slice remains one representative chr per subgenome. Large §9 science backlog unchanged.

