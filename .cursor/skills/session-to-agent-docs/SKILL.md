---
name: session-to-agent-docs
description: >-
  Summarize a completed agent conversation into durable repo docs: append
  doc/progress/YYYY-MM-DD.md, update doc/project_knowledge/*.yaml
  for run/resource registry changes, and link doc/KNOWLEDGE_README.md. Use when
  the user asks to preserve session learnings, capture runs/paths into project
  knowledge, summarize dialogue into docs, or "write this conversation into
  agent config". Registry = resources and layout only — not routine stats/plot
  outputs or analysis parameters. Do not put project facts in .cursor/rules —
  see boundary table in this skill.
---

# Session → agent docs

Use when a **conversation produced durable operational facts** that should outlive the chat — new run folders, canonical paths, domain policies, replay commands, or engineering narrative.

**Complements** [todo-drive-close](.cursor/skills/todo-drive-close/SKILL.md): that skill drives scoped issue/work batches; this skill **captures registry + logs** for any session with durable output.

## Where things live — boundary table

| Location | Purpose | Put here | Do **not** put here |
|----------|---------|----------|---------------------|
| **`doc/project_knowledge/`** | Structured **project registry** (YAML) | Run-folder inventory, module index, canonical **resource** paths (BAM, maps, refs), publish **layout** conventions, domain policies (e.g. taxaBamMap), `nf_cmd_ref` pointers | Pipeline **output** catalogues (plots, info TSVs under `{job}/stats/…`), analysis parameters, plot recipes, scientific conclusions, full bash transcripts |
| **`doc/KNOWLEDGE_README.md`** | Human + agent **index** to project_knowledge | Overview, YAML kinds, table of files, scope / out-of-scope | Long prose, investigation detail |
| **`.cursor/rules/`** | **AI behaviour** only | Edit guardrails, conda/screen policy, coding conventions, when to read/update knowledge (workstation-core §14) | taxaBamMap facts, run lists, BAM paths, task-specific playbooks |
| **`doc/progress/YYYY-MM-DD.md`** | **Engineering narrative + NF replay** (append-only) | Intent, repo diffs, pass/fail, blockers; **`#### NF replay`** with full commands | Registry rows (use YAML) |
| **`doc/PROGRESS_README.md`** | **Index** to daily progress files | Recent-day wikilinks only | Session body text |
| **User notes / published artefacts** | Science & QC interpretation | Personal notes, PNG/TSV under data1 or run dirs | `project_knowledge` or rules |

**Rule of thumb:** Registry facts → **`project_knowledge`**. Replay commands → **`#### NF replay`** in today's progress file. Behaviour → **`.cursor/rules/`**.

## Agent read policy

| Source | Read for task work? |
|--------|---------------------|
| `KNOWLEDGE_README.md` + `project_knowledge/*.yaml` | **Yes — primary** |
| `.cursor/rules/`, `GENETICS_WORKFLOW.md`, repo source | **Yes** |
| `doc/progress/`, `PROGRESS_README.md` | **No** unless user explicitly requests replay / log archaeology |

Still **write** to archive logs when this skill or **`progress-logging.mdc`** require post-session logging.

## What is **not** project knowledge

Do **not** create or extend YAML to catalogue:

- **Routine pipeline outputs** — PNG/TSV under `{job}/stats/{mod}/plots|info|…` when layout is already in `domain/data_publish_tree.yaml`
- **Per-analysis deliverable lists** — heatmap names, hybrid merge policies, colormap conventions, GAM/LOESS parameters
- **Scientific or QC conclusions** — sample QC flags, plot interpretations, correlation summaries

Those belong in **published artefacts** or **daily progress** — not the registry.

## `doc/project_knowledge` layout

Start at **`doc/KNOWLEDGE_README.md`**, then **`doc/project_knowledge/manifest.yaml`**.

```
doc/project_knowledge/
  manifest.yaml          # registry of all YAML files + operational_logs pointers
  domain/                # stable semantics (paths, policies, layout)
    taxa_bam_map.yaml
    directory_layout.yaml
    data_publish_tree.yaml
  workspace/             # vmap4 modules and numbered run folders
    vmap4_root.yaml
    vmap4_00data.yaml
    vmap4_10stats_genome.yaml
    …
```

| `kind` (in YAML) | Layer | Example content |
|------------------|-------|-----------------|
| `domain_concept` | `domain/` | taxaBamMap membership, three-root layout, `00data` / data1 trees |
| `workspace_registry` | `workspace/` | Top-level vmap4 module list |
| `workspace_module` | `workspace/` | One module path + `runs[]` (folder, one-line role, `artefacts`, `nf_cmd_ref`) |

**Registry entry shape** (minimal):

```yaml
- folder: 17run_example_slug
  summary: One-line role only (no analysis numbers)
  artefacts: [key/output/paths/relative/to/run/or/publish]
  nf_cmd_ref: doc/progress/YYYY-MM-DD.md — NF replay: slug
```

Cross-link logs; do not paste commands or findings into YAML.

## Non-negotiables

Follow **`.cursor/rules/workstation-core.mdc`** and **`progress-logging.mdc`** (especially **guardrail 14**):

- **English** in all versioned text.
- **`doc/progress/YYYY-MM-DD.md`:** append-only; `##` narrative + optional **`#### NF replay`** (`<!-- nf-replay -->` wrapper).
- **Repo vs run:** code under `/data/home/tusr1/git/script`; runs under `/data/home/tusr1/01projects/vmap4/<module>/<NNrun_*>/`.
- **`.cursor/rules/`:** change only for **AI policy** updates — user request required (this skill counts). **Never** migrate project facts from YAML into rules.
- **New run / new canonical path:** update **`project_knowledge` in the same session** as progress logging when possible — not rules.
- **Preserve existing project knowledge:** do **not** lightly rewrite or restructure existing `doc/project_knowledge/*.yaml`. Read the full target file first; prefer **append** (new `runs[]`, new keys) over renaming fields or editing unrelated entries. Use a **minimal diff** when fixing a prior error. If terminology is unclear, read canonical sources (`workflow/Genetics/docs/GENETICS_WORKFLOW.md`, existing YAML) — e.g. **`mod`** = Nextflow `params.mod` (`test_thin`, `test_common_thin`), **not** subgenome `A`/`B`/`D`/`Others`. Do not cascade bulk renames across knowledge files without explicit user approval.

## Workflow

### 1) Extract from the conversation

Split facts by sink:

| Registry (`project_knowledge`) | Logs / elsewhere |
|-------------------------------|------------------|
| New `NNrun_*` folder name and module | Full commands → `#### NF replay` in today's progress file |
| Stable **resource** path or publish **layout** convention | Engineering story → today's progress file |
| Domain policy change (e.g. new pop map) | Plot/stat outputs under data1 → artefacts only |
| Key **non-routine** resource paths (checkpoints, frozen refs) | Analysis params, plot recipes, QC conclusions → skip registry |

Drop: dead ends, typos, unverified speculation, verbatim chat.

### 2) Choose sinks

| Sink | When |
|------|------|
| **`doc/project_knowledge/*.yaml`** | New run folder, new module, new stable **resource** path, or publish **layout** / domain policy — **not** routine stats/plot outputs |
| **`doc/progress/YYYY-MM-DD.md`** | Any narrative, investigation, validation, or replay-worthy shell |
| **`doc/PROGRESS_README.md`** | Prepend wikilink when creating a **new** daily file |
| **`manifest.yaml` + `KNOWLEDGE_README.md`** | New YAML file or new top-level section |
| **`.cursor/rules/*.mdc`** | **Rare:** AI behaviour change only — **not** project facts |
| **`AGENTS.md`** | Portal index row if new skill or doc entry |

### 3) Update project knowledge

1. Read **`manifest.yaml`**.
2. Edit the matching file:

| Change | File | Action |
|--------|------|--------|
| taxaBamMap, publish layout, `00data` | `domain/*.yaml` | Extend `policies`, `entities`, `examples` |
| New run in existing module | `workspace/vmap4_{module}.yaml` | Append `runs[]` |
| New vmap4 top-level module | `workspace/vmap4_root.yaml` + new workspace YAML | Add module + manifest entry |
| Analysis / QC conclusions | — | **Skip** project_knowledge |

3. Register in **`manifest.yaml`** if new file.
4. Add **`KNOWLEDGE_README.md`** table row if new file.

**When editing existing YAML:** change only what the session requires; leave unrelated keys, paths, and naming conventions intact unless the user explicitly asks for a registry-wide correction.

### 4) Append logs

**Progress (today's file):** narrative bullets; add **`#### NF replay`** when commands were run (cwd, outcome, fenced `bash`).

Reference new YAML in bullets, e.g. `doc/project_knowledge/workspace/vmap4_10stats_genome.yaml`.

### 5) Update indexes

- `manifest.yaml`, `KNOWLEDGE_README.md`, `PROGRESS_README.md` (if new day file), `AGENTS.md` (if needed).
- **Do not** duplicate registry tables into `.cursor/rules/`.

### 6) Confirm with user

List: files updated; what went to **registry vs logs vs neither** (and why).

## Checklist

- [ ] Boundary respected: **facts → project_knowledge**, **behaviour → rules**
- [ ] Existing project_knowledge edited with **minimal diff** only (no unrelated renames/restructure)
- [ ] New run / path → workspace or domain YAML + `nf_cmd_ref`
- [ ] `#### NF replay` present when runs occurred
- [ ] No routine pipeline outputs, plot catalogues, or analysis params in YAML
- [ ] `manifest.yaml` / `KNOWLEDGE_README` updated if new YAML
- [ ] English throughout
- [ ] Progress appended to today's daily file when session had any durable narrative

## Example trigger phrases

- “Summarize this conversation into agent config”
- “Write runs and paths into project knowledge”
- “Preserve this session — update progress log and registry”
- “把对话总结到 agent 配置” (translate versioned content to English)
