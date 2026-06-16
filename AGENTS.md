# Agent instructions (`script`)

This repository is the **Genetics Analysis Pipeline & Library** (`script`): Nextflow workflows under `workflow/Genetics/`, shared Python/R/Java libraries under `src/`, and engineering docs under `doc/`.

**`AGENTS.md` is a portal only.** Canonical agent rules, conventions, and skills live under [`.cursor/`](.cursor/). Do not treat this file as a second copy of those policies—read the linked files below.

## How to use this portal

| Agent host | What to load |
|------------|--------------|
| **Cursor** | Rules and skills under `.cursor/` are applied automatically when configured in the IDE. |
| **Codex, Claude Code, or other agents** | At session start, read every file listed in **Rules** and any **Skills** relevant to the task. For `.mdc` rules, the markdown body below the YAML frontmatter is binding. |

When instructions conflict, prefer the **most specific** rule for the path or task (e.g. Nextflow edits → `workstation-nextflow.mdc`; Python under `src/python/` → `workstation-python.mdc`). Global guardrails always apply via `workstation-core.mdc`.

## Rules (`.cursor/rules/`)

| File | Scope |
|------|--------|
| [`.cursor/rules/workstation-core.mdc`](.cursor/rules/workstation-core.mdc) | **Always** — role, edit guardrails, repo map, vmap4 run policy, plotting/TODO conventions |
| [`.cursor/rules/workstation-nextflow.mdc`](.cursor/rules/workstation-nextflow.mdc) | `workflow/**/*.nf`, `workflow/Genetics/docs/**` |
| [`.cursor/rules/workstation-python.mdc`](.cursor/rules/workstation-python.mdc) | `src/python/**/*.py` |

## Skills (`.cursor/skills/`)

| File | When |
|------|------|
| [`.cursor/skills/todo-drive-close/SKILL.md`](.cursor/skills/todo-drive-close/SKILL.md) | User asks to work through `doc/TODO.md`, close backlog items, or update progress / NF run logs |

Add new skills under `.cursor/skills/` and link them here—keep `AGENTS.md` as an index, not the skill body.

## Operator docs (human + agent context)

These are not agent rules, but agents should cite them when running or extending pipelines:

- [`workflow/Genetics/docs/GENETICS_WORKFLOW.md`](workflow/Genetics/docs/GENETICS_WORKFLOW.md) — canonical Genetics Nextflow operator narrative
- [`doc/TODO.md`](doc/TODO.md) — checklist; [`doc/TODO_PROGRESS_LOG.md`](doc/TODO_PROGRESS_LOG.md) — append-only engineering log; [`doc/NF_CMD.md`](doc/NF_CMD.md) — Nextflow command log

## Minimal agent behavior

- Follow **workstation-core** for language (English in versioned text), minimal diffs, and operational guardrails.
- **Repo vs run vs publish:** edit code under `/data/home/tusr1/git/script` only (no `work/`, `.nextflow*`, `pipeline_info/`). Run Nextflow from `/data/home/tusr1/01projects/vmap4/<module>/<NNrun_*>/` (e.g. `10stats.genome/01run_main_raw_popdepth`). Published artefacts and frozen refs live under `/data1/dazheng_tusr1/vmap4.VCF.v1/…` via `params.output_dir` / `publishDir`.
- **Pipeline layout:** `modules/local/` (process libs); `subworkflows/local/{entry,plink,wheat,upstream,partial}/`; partial reruns via `partial_router.nf --partial_task`; ops FTP under `subworkflows/tmp/ops/`. No `workflow/Genetics/tmp/` or repo `resources/`.
- **Conda:** `run` for Nextflow; `stats` for Python stats / pytest; immutable inputs under `/data1/dazheng_tusr1/vmap4.VCF.v1`.
- Do not `git commit` unless the user explicitly asks.
- Before editing [`.cursor/rules/`](.cursor/rules/) or [`.cursor/skills/`](.cursor/skills/), confirm the user wants agent policy changed—unless the user explicitly requests a policy update (as in “update AI config”).
