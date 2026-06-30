# Progress Log

Engineering records — narrative **and** Nextflow/shell command replay — are archived by day under [[progress/2026-06-30]] and sibling files in `doc/progress/`.

**Log any progress** — agents append to today's `doc/progress/YYYY-MM-DD.md` (create the file if missing). No need to finish a task or close an issue before writing.

**NF replay:** verbatim commands go in the same daily file under **`#### NF replay`** (wrapped in `<!-- nf-replay -->` … `<!-- /nf-replay -->` for nightly sync). See **`.cursor/rules/progress-logging.mdc`**.

Nightly sync: `doc/progress/YYYY-MM-DD.md` → local cache → progress_overview review inbox.

## Historical notes

**Router (VCF/PLINK branch, `main.nf`):** As of **2026-05-16**, only `--mod v1_plink`, `test_thin`, `test_camp`, and `test_common_thin` are dispatched (plus any `wheat_*` integrated mod). Older replay blocks below May 2026 may still use `--mod test_plink` or `test_plink_camp`; for **new** full-pipeline runs, substitute **`test_thin`** and **`test_camp`** respectively.

**PublishDir:** Stats processes often use `publishDir` with `mode: 'copy'`. Changing `output_prefix` produces **new** filenames under `stats/.../plots` but **does not delete** older PNG/TSV names. If plots “look unchanged,” confirm you are opening files whose basenames match the **current** prefix.

Migrated monolithic command blocks from the former `doc/NF_CMD.md` live under matching dates in **`doc/progress/`** (section **`## YYYY-MM-DD — NF replays`** or inline **`#### NF replay`** subsections).

## Recent entries

- [[progress/2026-06-30]]
- [[progress/2026-06-28]]
- [[progress/2026-06-25]]
- [[progress/2026-06-23]]
- [[progress/2026-06-22]]
- [[progress/2026-06-10]]
- [[progress/2026-06-09]]
- [[progress/2026-06-04]]

Older daily files remain under `doc/progress/` but are omitted here when outside the ~14-day window.

## Related

| File | Role |
|------|------|
| [`KNOWLEDGE_README.md`](KNOWLEDGE_README.md) | Structured project registry index |

Task backlog: GitHub Issues (`SirLearning/script`) — not stored in this repo.
