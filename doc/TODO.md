# Variation library — project TODO

Master checklist for the variation-library / vmap4 genetics workstream. Completed items stay checked for audit. Log completions in `doc/TODO_PROGRESS_LOG.md`.

**Structure:** Sections follow the **runtime order** in `workflow/Genetics/main.nf` → `genotype/processor.nf` → `genotype/stats.nf`: entry router → shared preprocess → per-`mod` processor workflows → stats. Upstream science (methodology) sits in **§1**; backlog analyses not yet wired into NF sit in **§9**.

---

## 1. Upstream methodology & design

Cross-cutting science and QC policy before or beside the PLINK/Hail machinery.

- [x] Workflow study / references
	- [x] `funFastCall3.nf`: SNP calling flow
	- [x] Cohort-scale processing practices (peer notes)
- [x] Baseline setup
	- [x] Python package layout (`environment.yml`, `setup.py`)
	- [x] Global `nextflow` configuration
- [ ] Sequencing QC — 📅 2026-03-27
	- [ ] Noise from short-read alignment
	- [ ] **Depth distribution and statistical modeling** (from `5.wheat-WGS-technology`)
		- [ ] Check Mahalanobis assumptions: normality required? Prefer or add a direct non-parametric alternative
		- [ ] K-mer uniform null model: quantify deviation of observed reads from ideal uniformity
		- [ ] CV distribution: why still poor after FastCall filters? High-depth sites with variance ≫ mean
		- [ ] Negative-binomial QQ-plot: assess model fit
		- [ ] Wheat WGS depth threshold: filter unreliable variants
	- [ ] **Missing rate and mapping quality** (from `5.wheat-WGS-technology`)
		- [ ] Test whether high-depth sites map to repeats; both low- and high-depth extremes are problematic
		- [ ] Tune FastCall3 parameters
		- [ ] High-depth, low-MQ sites (polyploidy): drop sites with MQ < 20
		- [ ] MQ vs missing-rate scatter and correlation
		- [ ] Variant distribution among low-MQ sites
- [ ] Large-cohort variation library build-out
- [x] Depth-only thresholds rejected; use joint QC on Depth / MQ / MAF / Missing (LogRef: 2026-04-22)
- [ ] Encode multi-metric thresholds as reusable parameter templates (per A / B / D / Others)

---

## 2. `main.nf` — entry router, inputs, frozen paths — [[4.process]]

`workflow { }` selects `params.mod` after `check_input` builds `ch_vcf` (optional `--chr` filter).

- [x] `main.nf` dispatches: `v1_plink` \| `test_plink` / `test_thin` \| `test_plink_camp` / `test_camp` \| `test_common_thin`
- [x] Freeze `test_plink` / `test_plink_camp` entry points and command baseline (`main.nf` example 20260401.1)
- [x] Fixed test input dir: `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process` (read-only) (LogRef: 2026-04-22)
- [x] Fixed VCF root: `/data1/dazheng_tusr1/vmap4.VCF.v1` (read-only) (LogRef: 2026-04-22)
- [x] A / B / D / Others + `A_test.` prefix conventions frozen (LogRef: 2026-04-22)
- [x] `plink` version pinned
- [ ] **Router gap:** `HAIL`, `KINSHIP`, `POPULATION_STRUCTURE`, `GWAS`, genotype `database` are `include`d but **not** branched in `workflow { }` — add `params.mod` (or meta-workflow) + examples
- [ ] **Docs drift:** refresh `workflow/Genetics/README.md` (active mods, `process_dir` / `--camp` TSV, output layout, run-folder / `run` conda / `screen` policy)
- [ ] **Provenance:** emit `run_manifest.tsv` (or JSON): `params` subset, plink2 version, conda envs, Tiger JAR tag, optional git commit
- [ ] `hail` version + gnomAD-style reference; expose `--mod hail` (or equivalent) in `main.nf` + resource profile

---

## 3. Shared preprocess — `plink_preprocess` (`processor.nf`)

Used by test and `v1` tracks unless `process_dir` short-circuits pfile/VCF paths.

- [x] `params.process_dir`: reuse prebuilt per-chr `.plink2` at `{process_dir}/{id}.plink2` and matching `sample/` + `variant/` artifacts for `v1_plink` when set

---

## 4. Processor — test mods (`genotype/processor.nf`)

Sub-workflows run **in this order** inside each test path: **input pfiles** → **(optional thin / hard-filter)** → **group A/B/D/Others** → **`merge_subgenome_test_pfile`** → **`mk_plink_basic_info*`** → **`calc_plink_ld_unphased`** → **`calc_plink_ld_crosschr_random`**.

### 4.1 `test_plink` / `test_thin` → `test_plink_processor`

- [x] `process_dir` branch OR `plink_preprocess` → `subsampling_pfile_for_test` (`params.thin_rate`)
- [x] Map chr → subgenome (A / B / D / Others); `groupTuple`; merge per subgenome
- [x] `mk_plink_basic_info` → emits `smiss`, `vmiss`, `scount`, `gcount`, `afreq`, `hardy`
- [x] `calc_plink_ld_unphased` + `calc_plink_ld_crosschr_random` (outputs feed **§5**)
- [ ] Runtime / summary comparison table by A / B / D / Others (trace or `work/` timing)

### 4.2 `test_plink_camp` / `test_camp` → `test_plink_camp`

- [x] Same pipeline shape as **§4.1**, but `mk_plink_basic_info_camp_pop_with_filter` + `params.camp` TSV (`camp_vmap4_map.tsv` in examples)
- [ ] (Open) Document / validate CAMP map columns vs `psam` IDs whenever cohort map changes

### 4.3 `test_common_thin` → `test_common_thin_processor`

- [x] `subsampling_common_variant_pfile_for_test` after preprocess: `--geno` / `--maf` / `--thin` from `params.hf_geno`, `params.hf_maf`, `params.hf_thin_rate`
- [x] Same merge → `mk_plink_basic_info` → LD pair as **§4.1**
- [ ] **Config debt:** `nextflow.config` `hf_*` “TODO: rerun with new thresholds” — retune, log values in `TODO_PROGRESS_LOG.md`

---

## 5. Stats — `test_plink_stats` (`genotype/stats.nf`)

Runs **after** any of **§4.1–4.3** (`main.nf` workflows `test_plink`, `test_plink_camp`, `test_common_thin`). Order in code: **sample** pipelines → **variant** stats → **LD** plots.

### 5.1 Sample-level (Python `stats` env)

- [x] `sample_missing_stats` → `sample_coverage_stats` (with taxaBamMap) → `sample_heterozygosity_stats` → `sample_mapping_rate_stats` → `sample_ref_ibs_stats` → `sample_germplasm_dedup`
- [ ] Threshold / report polish driven by cohort QC needs (tie to **§9** depth / het anomalies)

### 5.2 Variant-level

- [x] `variant_missing_stats`, `variant_maf_stats`
- [ ] Extend plots / tables as needed for MAF–missing joint interpretation (see **§9**)

### 5.3 LD summaries (from processor outputs)

- [x] `variant_ld_decay_plot`, `variant_ld_crosschr_plot` consume **§4** `.vcor` outputs
- [ ] Promote same LD stats path to **§6** `v1_plink` / `plink_stats`, or document test-only rationale

---

## 6. Production — `v1_plink` → `plink_processor` + `plink_stats`

`main.nf` → `vmap4_v1_plink`: **PLINK_PROCESSOR** then **PLINK_STATS**.

### 6.1 `plink_processor`

- [x] `plink_preprocess` → per-chr pfile/VCF
- [x] Either reuse `process_dir` artifacts **or** `mk_plink_basic_info` + Tiger `calc_population_depth` → `popdep`
- [ ] Keep processor/stats boundary: all `plink2` heavy steps stay in `processor.nf` (**workstation-nextflow** rule)

### 6.2 `plink_stats`

- [x] `variant_missing_stats`, `variant_maf_stats`, `variant_popdep_mahalanobis` (NB / Mahalanobis in `src/python/genetics/genomics/variant/popdep.py`)
- [ ] **Parity:** port missing `test_plink_stats` sample-level diagnostics where appropriate for production runs
- [ ] No LD decay/cross-chr plots in this branch yet — align with **§5.3** open item

---

## 7. Filtering & reliable-set delivery — [[6.filter]]

Downstream of caller + library merge; parameters touch `nextflow.config` / future Hail.

- [ ] `MAF`, missing rate, `HWE` (v1: MAF + missing first; watch `HWE`)
- [x] Versioned filter provenance (params / software / sample set)
- [ ] **Heterozygosity and filtering** (from `3.variation-library`)
	- [ ] Recompute and validate median inbreeding coefficient 0.333
	- [ ] 10× heterozygosity bound too loose? (filters good sites when allele count < 3)
	- [ ] Heatmap: high-heterozygosity sites at low MAC
	- [ ] Wheat CNV / polyploidy impact on heterozygosity calls
- [ ] **Missing-rate filter tuning** (from `2.population-genetics`)
	- [ ] Is missing < 0.30 too strict? (removed entire AA group)
	- [ ] Validate new missing-rate threshold 0.75
- [ ] `hail` filtering: try annotating filter criteria as `INFO` only (no hard drop) and benchmark speed
	- [ ] Reliable set definition
	- [ ] `gVCF` merge across samples
		- [ ] Embed `vlib` for speed?
- [ ] Ship `VCF filter params v1` (A / B / D / Others retention and behavior)
- [ ] **Integrated filter / QC** (cross-cutting)
	- [ ] Hard-filter tuning: inbreeding, heterozygosity, chromosome context
	- [ ] Per-chromosome QC (SV differences)
	- [ ] Check independence assumptions when merging samples by chromosome

---

## 8. Parallel / legacy modules (not on current router)

- [x] **Assess (debug):** `genotype/assess.nf` wired for `test_thin` / `test_common_thin` via `tmp/assess_plink_debug.nf` (export slice + `quick_count` + `bcftools_qc_assess` + MAF bins); full router integration still optional (LogRef: 2026-05-13)
- [x] Hail scaffold: `genotype/hail.nf` + `src/python/genetics/hail/*` — await **§2** router + mod docs

---

## 9. Extended assessment backlog — [[5.assess]]

Science and plots **beyond** the current `test_plink_stats` / `plink_stats` wiring; keep for library interpretation.

- [x] Assess inputs: fixed test dir only; do not rebuild test VCFs (LogRef: 2026-04-22)
- [x] Core plot: variance vs depth (aberrant sites / repeats)
- [x] Tier-1 assess debug (both mods): `workflow/Genetics/tmp/assess_plink_debug.nf` — per-subgenome PLINK2 export (A=chr1, B=3, D=5, Others=0), strip `##chrSet` for bcftools, `quick_count` + `bcftools_qc_assess` + MAF-bin table under `assess/<mod>/` (LogRef: 2026-05-13)
- [ ] Singletons: counts and distribution
	- [x] Debug slice: MAF-bin counts + per-site `maf_missing.tsv` / `gq_summary.tsv` (no FORMAT/GQ in export → placeholder summary)
	- [ ] Effect on LD
	- [ ] Sample size vs variant discovery
	- [ ] Why singletons are scarce
		- [ ] Many variants with MAC < 10?
		- [ ] `blib` step had no filter params — why still insufficient?
		- [ ] **(from `3.variation-library`)** Investigate singleton count anomaly (still 6); confirm FastCall `blib` filtering behavior
- [ ] **Allele frequency** (from `3.variation-library`)
	- [x] Representative-chromosome MAF + `F_MISSING` table from merged test pfiles (`assess_plink_debug.nf`)
	- [ ] MAC < 10 heterozygous sites: link to per-site missing rate?
	- [ ] MAF curve rises then falls: likely unmerged samples — rebuild `lib.gz` and verify
- [ ] Rebuild `lib.gz` and compare MAC ≥ 1 vs MAC ≥ 2
- [ ] MAF
	- [ ] Distribution: mechanism behind rise-then-fall
	- [ ] **(from `3.variation-library`)** MAF vs LD association: biological vs artifact?
- [ ] Depth
	- [x] Empirical depth distribution
	- [ ] Statistics: confidence intervals
	- [ ] **Per-sample depth anomalies** (from `2.population-genetics`)
		- [ ] `CS_mp_2018_8X`: metadata 8× but mean depth ~5.44 — targeted sequencing?
		- [ ] `SAMEA10944181` depth zero — handle or exclude
- [ ] Dimensionality reduction
	- [ ] PCA for population layout (flag outliers first; do not drop samples preemptively)
		- [ ] **(from `2.population-genetics`)** Heterozygosity vs missing-rate regression line looks wrong — biological interpretation?
		- [ ] **(from `2.population-genetics`)** PCA on hexaploid-only panel vs Schulthess et al. 2022 reference
	- [ ] UMAP
- [ ] LD (biology / maps beyond **§5.3** plots)
	- [ ] Decay curves: cohort-wide LD distribution; genetic drivers
	- [ ] **LD deep dive** (from `3.variation-library`)
		- [ ] Chromosome-wide maps; flag centromeric high-LD regions
		- [ ] D-subgenome LD anomaly vs population structure
	- [ ] LD vs marker density
- [ ] Relatedness
	- [ ] K matrix
		- [x] Poor fit for selfing wheat (large negative values)
		- [ ] Explain bundled / streaky patterns
		- [ ] **(from `2.population-genetics`)** Kinship vs IBS0 scatter: slope clusters (−1500, −600, −420, …); flag low-heterozygosity samples
	- [ ] IBS
		- [ ] Define taxa; merge BAMs if needed
		- [ ] **IBS anomalies** (from `2.population-genetics`)
			- [ ] No IBS = 1 pairs despite duplicate samples — why?
			- [ ] Verify IBS implementation
			- [ ] Chinese Spring v1.0 reference: gaps and mis-mapping
			- [ ] ABD cohort IBS “smoothing”; locate duplicate pairs in IBS space
			- [ ] Validate AB–AB bimodal pattern
	- [ ] GRM
		- [ ] Z-score computation
- [ ] Population structure: Q matrix

---

## 10. Future

- [ ] Scale GBS sample size
- [ ] Step-3 direct scan: runtime / logic differences?
- [ ] Optional CI smoke: `test_plink` dry run or `-stub-run` when CI paths exist
