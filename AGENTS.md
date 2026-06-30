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
| [`.cursor/rules/workstation-core.mdc`](.cursor/rules/workstation-core.mdc) | **Always** — role, edit guardrails, repo map, vmap4 run policy, plotting conventions |
| [`.cursor/rules/progress-logging.mdc`](.cursor/rules/progress-logging.mdc) | **Always** — daily progress log, NF replay blocks, write triggers |
| [`.cursor/rules/workstation-nextflow.mdc`](.cursor/rules/workstation-nextflow.mdc) | `workflow/**/*.nf`, `workflow/Genetics/docs/**` |
| [`.cursor/rules/workstation-python.mdc`](.cursor/rules/workstation-python.mdc) | `src/python/**/*.py` |

Project domain knowledge (taxaBamMap, vmap4 run registry) lives in [`doc/KNOWLEDGE_README.md`](doc/KNOWLEDGE_README.md) and [`doc/project_knowledge/`](doc/project_knowledge/) — not in rules.

## Agent skills

### Issue tracker

GitHub Issues on `SirLearning/script` (via `gh` CLI). See [`doc/agents/issue-tracker.md`](doc/agents/issue-tracker.md).

### Triage labels

Default mattpocock/skills vocabulary (`needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`). See [`doc/agents/triage-labels.md`](doc/agents/triage-labels.md).

### Domain docs

Single-context layout: `CONTEXT.md` and `doc/adr/` at repo root (created lazily by `/grill-with-docs`). See [`doc/agents/domain.md`](doc/agents/domain.md).

## Skills (`.cursor/skills/`)

### Project-specific

| Skill | When |
|-------|------|
| [todo-drive-close](.cursor/skills/todo-drive-close/SKILL.md) | Drive GitHub Issues / scoped work batches; append daily progress and NF run logs |
| [session-to-agent-docs](.cursor/skills/session-to-agent-docs/SKILL.md) | Summarize a conversation into project knowledge and today's `doc/progress/` |

### [mattpocock/skills](https://github.com/mattpocock/skills) — setup (run once)

| Skill | When |
|-------|------|
| [setup-matt-pocock-skills](.cursor/skills/setup-matt-pocock-skills/SKILL.md) | Reconfigure issue tracker, triage labels, or domain doc layout |

### Engineering (daily)

| Skill | When |
|-------|------|
| [grill-with-docs](.cursor/skills/grill-with-docs/SKILL.md) | **Before new features** — align on design, update `CONTEXT.md` and ADRs |
| [tdd](.cursor/skills/tdd/SKILL.md) | Red-green-refactor test-driven development |
| [diagnose](.cursor/skills/diagnose/SKILL.md) | Structured debugging for hard bugs |
| [improve-codebase-architecture](.cursor/skills/improve-codebase-architecture/SKILL.md) | Find deepening opportunities; run periodically on the codebase |
| [to-prd](.cursor/skills/to-prd/SKILL.md) | Turn conversation into a PRD (GitHub issue) |
| [to-issues](.cursor/skills/to-issues/SKILL.md) | Break a plan into vertical-slice GitHub issues |
| [triage](.cursor/skills/triage/SKILL.md) | Triage issues through the label state machine |
| [zoom-out](.cursor/skills/zoom-out/SKILL.md) | Higher-level context on unfamiliar code |
| [prototype](.cursor/skills/prototype/SKILL.md) | Throwaway prototype to flesh out a design |

### Productivity

| Skill | When |
|-------|------|
| [grill-me](.cursor/skills/grill-me/SKILL.md) | Grilling session for non-code plans |
| [caveman](.cursor/skills/caveman/SKILL.md) | Ultra-compressed communication (~75% fewer tokens) |
| [handoff](.cursor/skills/handoff/SKILL.md) | Compact conversation into a handoff doc |
| [teach](.cursor/skills/teach/SKILL.md) | Multi-session teaching workspace |
| [write-a-skill](.cursor/skills/write-a-skill/SKILL.md) | Author new skills with proper structure |

Add new skills under `.cursor/skills/` and link them here—keep `AGENTS.md` as an index, not the skill body.

## Operator docs (human archive + agent primary source)

| Doc | Agents |
|-----|--------|
| [`doc/KNOWLEDGE_README.md`](doc/KNOWLEDGE_README.md) + [`doc/project_knowledge/`](doc/project_knowledge/) | **Read first** for paths, runs, resource semantics |
| [`workflow/Genetics/docs/GENETICS_WORKFLOW.md`](workflow/Genetics/docs/GENETICS_WORKFLOW.md) | **Read** for Nextflow operator narrative |
| [`doc/progress/`](doc/progress/), [`doc/PROGRESS_README.md`](doc/PROGRESS_README.md) | **Do not read** for task planning (archive only; often stale). **Append** today's progress file when rules/skills require logging. Load only if the user explicitly asks for replay or log archaeology. |

## Minimal agent behavior

- Follow **workstation-core** for language (English in versioned text), minimal diffs, and operational guardrails.
- **Repo vs run vs publish:** edit code under `/data/home/tusr1/git/script` only (no `work/`, `.nextflow*`, `pipeline_info/`). Run Nextflow from `/data/home/tusr1/01projects/vmap4/<module>/<NNrun_*>/` (e.g. `10stats.genome/01run_main_raw_popdepth`). Published artefacts and frozen refs live under `/data1/dazheng_tusr1/vmap4.VCF.v1/…` via `params.output_dir` / `publishDir`.
- **Pipeline layout:** `modules/local/` (process libs); `subworkflows/local/{entry,plink,wheat,upstream,partial}/`; partial reruns via `partial_router.nf --partial_task`; ops FTP under `subworkflows/tmp/ops/`. No `workflow/Genetics/tmp/` or repo `resources/`.
- **Conda:** `run` for Nextflow; `stats` for Python stats / pytest; immutable inputs under `/data1/dazheng_tusr1/vmap4.VCF.v1`.
- Do not `git commit` unless the user explicitly asks.
- **Never run destructive shell** (`rm`, `find -delete`, publish-tree cleanup, etc.) without the user's explicit same-message approval; see **workstation-core** guardrail 5.
- Before editing [`.cursor/rules/`](.cursor/rules/) or [`.cursor/skills/`](.cursor/skills/), confirm the user wants agent policy changed—unless the user explicitly requests a policy update (as in “update AI config”).
