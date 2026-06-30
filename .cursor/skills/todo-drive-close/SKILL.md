---
name: todo-drive-close
description: >-
  Drives scoped engineering work from GitHub Issues or user-named batches: implement
  in-repo changes, run validation from vmap4 project run roots with conda `run`,
  debug, then append today's doc/progress/YYYY-MM-DD.md (narrative + NF replay blocks).
  Use when the user asks to drive an issue, close a backlog item, batch-complete scoped
  work, or mentions issue drive / progress log. Task backlog is not in repo.
---

# Issue drive — scoped work batches

Use when the user **explicitly** asks to drive a GitHub Issue, named work batch, or mentions issue drive / progress log.

**Do not** open `doc/progress/` for general task planning — use **`doc/KNOWLEDGE_README.md`** and **`doc/project_knowledge/`** instead (workstation-core guardrail 14). Task backlog: **GitHub Issues** (`SirLearning/script`) or vault.

## Non-negotiables (this repo)

Follow **`.cursor/rules/workstation-core.mdc`**, **`progress-logging.mdc`**, and **`.cursor/rules/workstation-nextflow.mdc`** in full. See those files for repo vs run, conda `run`, `screen`, vmap4 inputs, assess track, and partial NF reruns.

**Progress:** Append to **`doc/progress/YYYY-MM-DD.md`** for **any progress**. **Append-only.** After Nextflow runs, include an **`#### NF replay`** block in the **same** daily file (see **`progress-logging.mdc`**).

## Workflow

### 1) Inventory scope

Read the **GitHub issue** (`gh issue view`) or confirm the **batch scope** with the user.

### 2–4) Plan, implement, validate

Same as before: minimal diffs, correct vmap4 run cwd, iterate on failures.

### 5) Log progress (same session)

1. **`doc/progress/YYYY-MM-DD.md`:** `## YYYY-MM-DD — title` with bullets; add **`<!-- nf-replay -->` … `#### NF replay` … `<!-- /nf-replay -->`** when commands were run.
2. **`PROGRESS_README.md`:** prepend wikilink if new daily file.
3. **GitHub Issues:** close/comment when requested and verified.

### 6) No drive-by doc edits unless logging this session.

## Example progress entry

```markdown
## YYYY-MM-DD — Assess debug revalidation

- Goal: smoke-test assess_plink on test_thin …
- Outcome: exit 0, 32 succeeded

<!-- nf-replay -->
#### NF replay — assess_plink test_thin

Working directory: `/data/home/tusr1/01projects/vmap4/08stats.genome/23run_…`. Outcome: exit 0.

```bash
cd … && source ~/.bashrc && conda activate run && nextflow run …
```
<!-- /nf-replay -->
```

## Invocation

Cite a **GitHub issue** or scoped batch; optionally attach **`@.cursor/skills/todo-drive-close/SKILL.md`**.
