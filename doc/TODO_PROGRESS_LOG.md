# TODO progress log

Append-only log for completed TODO items (per `.cursor/rules/workstation-core.mdc`). One entry per completion batch or item.

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
