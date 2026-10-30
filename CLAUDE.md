# Claude Code project configuration

This file only wires Claude Code to the repo's canonical agent policy; it holds no
policy of its own. Canonical sources are [`AGENTS.md`](AGENTS.md) and
[`.cursor/rules/`](.cursor/rules/).

Claude Code defaults to `claude-md-or-agents-md`: when this file exists, `AGENTS.md`
is **not** read automatically. The import below keeps the portal active.

## Auto-loaded

@AGENTS.md
@.cursor/rules/workstation-core.mdc
@.cursor/rules/progress-logging.mdc

## Load on demand (path-scoped)

Read these before the first edit in the matching area:

- `src/python/**/*.py` → [`.cursor/rules/workstation-python.mdc`](.cursor/rules/workstation-python.mdc)
- `workflow/**/*.nf`, `workflow/Genetics/docs/**` → follow
  [`workflow/Genetics/docs/GENETICS_WORKFLOW.md`](workflow/Genetics/docs/GENETICS_WORKFLOW.md);
  `.cursor/rules/workstation-nextflow.mdc` is referenced by `AGENTS.md` but is not
  present in this repo.

## CodeGraph

[`.cursor/rules/codegraph.mdc`](.cursor/rules/codegraph.mdc) applies only when the
CodeGraph MCP server (`codegraph_*` tools) is configured for the current host. It is
not configured for Claude Code today, so use native search/read instead. If the MCP
server is added later, load that rule then.

## Project skills

`.cursor/skills/*` are exposed to Claude Code through symlinks in
[`.claude/skills/`](.claude/skills/). Invoke them by name (e.g. `/todo-drive-close`);
their authored source stays under `.cursor/skills/`.
