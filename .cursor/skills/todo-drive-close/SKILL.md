---
name: todo-drive-close
description: >-
  Drives completion of open items in `doc/TODO.md`: scope unchecked work, implement
  in-repo changes, run validation from the vmap4 project run roots with conda `run`,
  debug, then append `doc/TODO_PROGRESS_LOG.md` (narrative: goals, rationale, outcomes,
  cwd per run; no duplicate of full commands in `doc/NF_CMD.md`) and `doc/NF_CMD.md`
  when applicable. Use when the user asks to work through the TODO checklist, close
  backlog items, batch-complete `doc/TODO.md`, or mentions TODO drive / progress log /
  NF command log.
---

# TODO drive — close checklist items

Use this skill when the user **explicitly** asks to work through the TODO checklist, close backlog items, or mentions TODO drive / progress log / NF command log.

**Do not** open `doc/TODO.md`, `doc/TODO_PROGRESS_LOG.md`, or `doc/NF_CMD.md` for general task planning — use **`doc/KNOWLEDGE_README.md`** and **`doc/project_knowledge/`** instead (workstation-core guardrail 14).

## Non-negotiables (this repo)

Follow **`.cursor/rules/workstation-core.mdc`** and **`.cursor/rules/workstation-nextflow.mdc`** in full. In particular:

- **Repo vs run:** Edit source only under `/data/home/tusr1/git/script`. Run Nextflow and vmap4-tied analysis from **`/data/home/tusr1/01projects/vmap4/<module>/<NNrun_<slug>>/`** (two-level layout; see workstation-core guardrail 9). Example module: `08stats.genome/` with runs like `57run_mac_stats_test_thin`. Create a **new** numbered folder per attempt; never reuse an old run dir.
- **Conda:** Use conda env **`run`** for Nextflow and related runtime (`source ~/.bashrc && conda activate run`).
- **Long Nextflow in `screen`:** Start real pipeline runs inside **`screen` from conda `run`** (not system `/usr/bin/screen`). Pattern: activate `run` → `screen -dmS … bash -lc 'source ~/.bashrc && conda activate run && …'`. Re-activate `run` after interactive attach. Small tests / debug runs are exempt from `screen`.
- **vmap4 inputs:** Treat `/data1/dazheng_tusr1/vmap4.VCF.v1` as **read-only**. For prepared test data, use **`/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process`** and existing prefixes (e.g. `A_test.`).
- **Assess track:** Compute in **`processor/processor_assess.nf`**, plots in **`stats/stats_assess.nf`**; launch via **`partial_router.nf --partial_task assess_plink`** (PLINK2) or **`assess_vcf`** (VCF). Prefer **PLINK2 on pfiles** over bcftools; assess **plots/tables** must use **`src/python/genetics/...`** + **`infra.utils.io` / `infra.utils.graph`** and publish under `stats/<mod>/{plots,info,logs}`; compute under `process/<mod>/` (see `.cursor/rules/workstation-core.mdc` “Genotype assess”).
- **Partial NF reruns:** Use `subworkflows/local/entry/partial_router.nf` with `--partial_task`; pass **absolute** `-c` to `nextflow.config` (see `GENETICS_WORKFLOW.md`).
- **Docs language:** New or edited checklist text in `doc/TODO.md` stays **English**. **`doc/TODO_PROGRESS_LOG.md` is append-only** — never delete or rewrite past sections.
- **File edits:** Use editor/patch tools for source and docs; **do not** create or modify tracked files with shell redirection (`cat`/`echo`/`sed`/`tee` + `>`/`>>`).

## Workflow

### 1) Inventory open work

1. Read **`doc/TODO.md`** and list bullets that are still **`- [ ]`** (nested items count independently).
2. If the user named a **section, progress-log date/title, or bullet**, scope to that. Otherwise propose **one coherent batch** (single feature, single NF area, or single doc pass) and confirm implicit scope is reasonable; avoid mixing unrelated science + infra in one batch without user intent.

### 2) Plan the batch

- Map each targeted `[ ]` to **concrete deliverables** (files to touch, processes to add, tests to run).
- Identify whether validation is **unit/script**, **Nextflow**, or **manual review** only.

### 3) Implement

- **Read before edit** for shared infra called out in workstation-core (e.g. `src/python/infra/utils/io.py`, `src/python/infra/utils/graph.py`).
- Keep diffs **minimal** and aligned with existing naming and patterns.
- For plotting helpers, **reuse** `plot_*` in `src/python/infra/utils/graph.py` per workstation-core; add only **generic** helpers when nothing fits.

### 4) Validate (iterate until pass or blocked)

- Run the **smallest** commands that prove the change (tests, `-resume` NF slice, lint).
- Execute from the correct **`/data/home/tusr1/01projects/vmap4/<module>/<NNrun_<slug>>`** folder. If none exists yet for this attempt, create it first (`mkdir -p`, pick next `NN` from sibling listing).
- On failure: inspect logs/work dirs, fix, re-run. Document **partial validation** honestly in the progress log if blocked.

### 5) Close the checklist (same session as the shipped work)

When one or more targeted items are **actually done**:

1. **`doc/TODO.md`:** Flip matching **`[ ]` → `[x]`**. Preserve `LogRef`, dates, paths, `[[wiki]]` links, and checkbox nesting. Do **not** bulk-uncheck unrelated items.
2. **`doc/TODO_PROGRESS_LOG.md`:** Append a **new** section at the **end** (English) with heading `## YYYY-MM-DD — <short English title>` (no opaque IDs). Write a **report-style narrative** (or a table with a **`Summary`** row if needed) covering: **goal** of the work; **what changed** in the repo and **why**; each **test/run** (working directory only—full commands go to **`doc/NF_CMD.md`**); **whether it passed**; if not, **root cause** and **next fix**; **final outcome** and **how outputs are shaped** (paths, formats, prefixes). Cite `doc/NF_CMD.md` by date heading for command detail. Do **not** paste duplicate multi-line bash blocks. For thin older entries, append `## YYYY-MM-DD — SUPP — <topic>` at the end—**do not rewrite** earlier text.
3. **`doc/NF_CMD.md`:** If **any** `nextflow` (or `screen ... nextflow`) command was executed for this repo during the session, append a **new block at the end** in chronological order: `### YYYY-MM-DD` heading, optional one-line cwd/outcome, then a fenced **`bash`** block with the **full** command line(s). English prose only.

### 6) User did not ask for doc-only work

Do **not** edit `doc/TODO.md`, `doc/TODO_PROGRESS_LOG.md`, or `doc/NF_CMD.md` unless this skill session (or explicit user request) includes **closing** items or **logging** runs. Avoid drive-by doc rewrites.

## Quick reference — progress log (narrative + cwd; commands in `NF_CMD.md`)

Append to **`doc/TODO_PROGRESS_LOG.md`** (English). Example shape:

```markdown
---

## YYYY-MM-DD — <Short English title>

**Goal:** …

**Changes (what / why):** …

**Validation & runs:** For each attempt: **cwd** `…`; link **commands** to `doc/NF_CMD.md` (`### …`); **result** (pass/fail, task counts if useful); if fail, **cause** + **remediation** or **next step**.

**Outputs:** Where files land; PNG/TSV/etc.; how to read them.

**Risks / follow-ups:** …
```

If an earlier entry was table-only, add **`## YYYY-MM-DD — SUPP — …`** at the file end with the missing narrative (still **no** full command duplication if NF_CMD already has them).

## Invocation hint for the user

In chat, attach **`@doc/TODO.md`** and optionally **`@.cursor/skills/todo-drive-close/SKILL.md`**, then state which section or bullet to drive (or "pick the next unfinished infra item").
