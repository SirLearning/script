# Project knowledge (`doc/project_knowledge`)

Structured, machine-readable **registry** for the **Genetics/script** pipeline and **vmap4** workspace — what files mean, where they live, and which run folders exist.

**Scope (what belongs here):** paths, naming conventions, module/run inventory, domain policies (e.g. taxaBamMap membership), and **resource / input file semantics** (BAM, refs, maps, frozen publish layout).

**Out of scope (do not store here):**

- **Pipeline outputs** already implied by publish layout (`data_publish_tree.yaml`, `directory_layout.yaml`) — e.g. per-run PNG plots, info TSVs, thresholds under `{job}/stats/{mod}/plots|info|…`
- **Analysis task catalogues** — which plots were produced, hybrid merge recipes, plotting parameters, colormap choices, GAM/LOESS settings
- Detailed analysis, root-cause narratives, or scientific conclusions → user notes or published artefacts (not the registry)
- Plot interpretations or figure-derived findings → published PNG/TSV under data1 or run dirs
- Agent operational policy → [`.cursor/rules/`](../.cursor/rules/) (`workstation-core`, `workstation-nextflow`, `workstation-python`, `codegraph`)

YAML entries may optionally record `nf_cmd_ref` as a human pointer; agents must **not** load archive logs to resolve facts (see **Agent read policy** below).

## How to use

1. Read [`project_knowledge/manifest.yaml`](project_knowledge/manifest.yaml) — registry of all knowledge files.
2. Load YAML under `domain` or `workspace` for your task.
3. If a fact is missing from YAML, **verify on disk or in workflow source** — do not infer from archive logs.

### Agent read policy (canonical vs archive)

| Source | Agents **read** for task work? | Purpose |
|--------|----------------------------------|---------|
| **`KNOWLEDGE_README.md` + `project_knowledge/*.yaml`** | **Yes — primary** | Run folders, resource paths, layout, domain policies |
| **`.cursor/rules/`**, **`workflow/Genetics/docs/GENETICS_WORKFLOW.md`**, **repo source** | **Yes** | Behaviour, pipeline wiring, parameters |
| **`doc/progress/`**, **`doc/PROGRESS_README.md`** | **No** (unless user asks for log archaeology / command replay) | Human append-only daily log (narrative + `#### NF replay` commands) |

**Write:** Agents still **append** to today's **`doc/progress/YYYY-MM-DD.md`** when rules or skills require logging — that is **output**, not a read source for planning. Task backlog: GitHub Issues — not in repo.

Agents: update YAML via [session-to-agent-docs](../.cursor/skills/session-to-agent-docs/SKILL.md) when adding runs or domain paths. That skill defines the **boundary** between `project_knowledge`, `.cursor/rules`, and append-only logs.

## YAML kinds

| `kind` | Purpose | Example |
|--------|---------|---------|
| `domain_concept` | Stable file/path/population semantics | `taxa_bam_map`, `directory_layout` |
| `workspace_registry` | Index of vmap4 top-level modules | `vmap4_root` |
| `workspace_module` | One module + numbered `runs[]` (folder, role, artefact paths) | `vmap4_10stats_genome` |

Common fields: `id`, `title`, `path`, `summary`, `nf_cmd_ref`, `artefacts`.

## Index (from manifest)

### Domain

| File | Summary |
|------|---------|
| [`domain/taxa_bam_map.yaml`](project_knowledge/domain/taxa_bam_map.yaml) | taxaBamMap population membership, paths, job concat scope |
| [`domain/directory_layout.yaml`](project_knowledge/domain/directory_layout.yaml) | Code repo vs vmap4 run cwd vs data1 publish |
| [`domain/data_publish_tree.yaml`](project_knowledge/domain/data_publish_tree.yaml) | `00data` folders and data1 `{job}/process|stats` layout |
| [`domain/main_raw_variant_assets.yaml`](project_knowledge/domain/main_raw_variant_assets.yaml) | main_raw VCF/plink per chr, chr032 paths, PopDep vs VCF gate |

### Workspace (vmap4)

| File | Summary |
|------|---------|
| [`workspace/vmap4_root.yaml`](project_knowledge/workspace/vmap4_root.yaml) | Top-level module index under `01projects/vmap4` |
| [`workspace/vmap4_00data.yaml`](project_knowledge/workspace/vmap4_00data.yaml) | BAM, refs, taxaBamMap under `00data` |
| [`workspace/vmap4_10stats_genome.yaml`](project_knowledge/workspace/vmap4_10stats_genome.yaml) | PopDep / `main_raw_popdepth` runs (`01`–`16`) |
| [`workspace/vmap4_08stats_genome.yaml`](project_knowledge/workspace/vmap4_08stats_genome.yaml) | Legacy assess / LD / wheat / MAC runs |
| [`workspace/vmap4_04runScreens.yaml`](project_knowledge/workspace/vmap4_04runScreens.yaml) | FastCall scan screens, CS_mp run folders |
| [`workspace/vmap4_11aneuploidy.yaml`](project_knowledge/workspace/vmap4_11aneuploidy.yaml) | Mosdepth aneuploidy benchmark module + publish layout |

### Operational logs (Markdown, append-only — agents do not read for facts)

| File | Role |
|------|------|
| [`PROGRESS_README.md`](PROGRESS_README.md) | Index to daily progress (`progress/YYYY-MM-DD.md`; includes NF replay blocks) |

## Adding knowledge

- **New vmap4 run:** append a `runs[]` entry in `project_knowledge/workspace/*.yaml` (folder name, one-line role, key **resource** paths). Optionally add **`#### NF replay`** in today's progress file; registry facts belong in YAML first.
- **New domain fact:** extend `domain/*.yaml` for **resources and layout only**; register in `manifest.yaml`.
- **Stats/plots/info TSV outputs:** covered by `data_publish_tree` — **do not** add per-task YAML registries.
- **Analysis / QC conclusions:** personal notes or artefacts — **not** `project_knowledge`.
- **Do not** add task playbooks to `.cursor/rules/` — rules are for AI behaviour only.
