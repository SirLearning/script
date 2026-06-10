# Genetics Analysis Pipeline & Script Library

Nextflow workflows, Python/R/Java libraries, and engineering docs for large-scale crop genetics (VMap4 / wheat variation library).

**Canonical operator guide:** [`workflow/Genetics/docs/GENETICS_WORKFLOW.md`](workflow/Genetics/docs/GENETICS_WORKFLOW.md)  
**Backlog:** [`doc/TODO.md`](doc/TODO.md) · **Run log:** [`doc/NF_CMD.md`](doc/NF_CMD.md) · **Agent rules:** [`AGENTS.md`](AGENTS.md)

---

## Project layout

```
.
├── environment_*.yml          # Conda envs (one-step create per task)
├── setup.py                   # Editable install of src/python (python_script)
├── src/
│   ├── python/                # genetics, infra, WeaTE
│   ├── r/
│   └── java/
├── workflow/Genetics/
│   ├── main.nf                # Entry router (--mod)
│   ├── modules/local/         # Active process libraries (genotype/stats/ for stats workflows)
│   ├── subworkflows/
│   │   ├── local/             # Composed workflows (entry/, plink/, wheat/, partial/, …)
│   │   └── tmp/ops/           # Ops-only scripts (FTP upload)
│   ├── conf/                  # nextflow.config fragments
│   ├── docs/                  # Operator narrative (GENETICS_WORKFLOW.md)
│   └── tests/                 # Nextflow test harness config
├── doc/                       # TODO, progress log, NF command log
└── note/                      # Analysis notes
```

Legacy theme trees may remain under `workflow/Genetics/old/` for reference only.

---

## Environment setup

Environments are created **once from the repo-root YAML files** (no extra dependency manifests required).

| Env | File | Role |
| --- | --- | --- |
| `run` | `environment_run.yml` | Nextflow, screen |
| `stats` | `environment_stats.yml` | Python 3.12, PLINK/PLINK2, bcftools, Hail (conda + pip), stats stack |
| `stats_r` | `environment_stats_r.yml` | R plotting / BioConductor |
| `tiger` | `environment_tiger.yml` | Alignment / variant calling |
| `dbone` | `environment_dbone.yml` | DBone / basic ops |

```bash
conda env create -f environment_stats.yml
conda activate stats
pip install -e .    # registers python_script; conda yml already provides runtime deps
```

`pip install -e .` only installs the **in-repo package** so Nextflow and scripts can `import genetics` / `import infra`. Heavy tools (Hail, PLINK2, bcftools, compilers) stay in the conda YAML by design.

---

## Running the Genetics pipeline

Entry point: `workflow/Genetics/main.nf`. Always pass an **absolute** path to `nextflow.config` when launching from a project run folder.

```bash
conda activate run
nextflow run /data/home/tusr1/git/script/workflow/Genetics/main.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --home_dir … --user_dir … --src_dir … \
  --output_dir … --mod test_thin --job myjob
```

| `params.mod` | Track |
| --- | --- |
| `v1_plink` | Production PLINK processor + stats |
| `test_thin` / `test_camp` / `test_common_thin` | Test PLINK processors + stats |
| `wheat_*` | Integrated wheat table/matrix analytics (no VCF channel) |

Partial reruns use `workflow/Genetics/subworkflows/local/entry/partial_router.nf` with `--partial_task`; ops FTP scripts live under `subworkflows/tmp/ops/` — see **GENETICS_WORKFLOW.md**.

### vmap4 run workspace (not the git repo)

All production Nextflow runs and project-scale analysis use a **two-level** layout under `/data/home/tusr1/01projects/vmap4/`:

```text
01projects/vmap4/
  08stats.genome/                      # task module (stats, assess, partial genetics)
    01run/                             # full-pipeline attempts
    23run_assess_plink2_debug_stub/    # partial_router --partial_task assess_plink
    57run_mac_stats_test_thin/
    63run_mac_dist_log_redraw/
    …                                  # NN = monotonic; slug describes the attempt
```

1. Pick the **task module** folder (`08stats.genome` for stats/assess/LD/MAC/chr/wheat partial work).
2. List siblings, choose the next **`NN`**, create **`NNrun_<slug>`**, `cd` there.
3. Run with **absolute** paths to `main.nf` or `partial_router.nf` and `nextflow.config`.
4. Log cwd + command in `doc/NF_CMD.md`.

Authoritative policy: `.cursor/rules/workstation-core.mdc` and `workflow/Genetics/docs/GENETICS_WORKFLOW.md` § “Where and how to run”.

---

## Python library (`src/python`)

| Package | Role |
| --- | --- |
| `genetics.genomics` | Sample/variant QC, LD, PLINK result I/O, assess plots |
| `genetics.gwas` | Association result plotting |
| `infra.utils.io` / `graph` | Shared parsing and plotting (reuse before adding task wrappers) |
| `infra.server` | File copy, integrity checks |
| `WeaTE` | Wheat transposon analysis |

`genetics.wheat` is a deprecated re-export shim; new code should import from `genetics.genomics.*`.

**WeaTE transposon inputs:** Scripts under `src/python/WeaTE/` expect a project-local `transposon/` tree (e.g. `TEcode`, `np.stats`, `vu_reads_depth/`) relative to the run working directory—not a repo-level `resources/` folder.

---

## Tests

Lightweight unit tests live under `tests/` (no extra env files — use the existing `stats` conda env):

```bash
conda activate stats
pip install pytest ruff
pytest
ruff check src/python tests
```

CI runs the same checks via `.github/workflows/ci.yml`.
