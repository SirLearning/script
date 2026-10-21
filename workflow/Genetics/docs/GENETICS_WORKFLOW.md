# Genetics Nextflow workflow (`workflow/Genetics/`)

**Documentation:** This file (`docs/GENETICS_WORKFLOW.md`) is the **only** narrative operator guide for this tree. Previous theme trees may live under `old/` for reference; **active** process libraries are **`modules/local/{genotype,dynamic,static,integrated}/`**.

This directory holds the **Genotype** branch of the pipeline: PLINK/PLINK2 preprocessing, per-mod processors, stats (Python/R via `src/`), **composed subworkflows** under `subworkflows/local/`, ops scripts under `subworkflows/tmp/ops/`, and **wheat** table/matrix analytics.

---

## Entry point and config

| Item | Role |
| --- | --- |
| `main.nf` | Anonymous top-level `workflow { }`: for non-`wheat_*` mods, builds `ch_vcf` via `check_input` from `subworkflows/local/entry/genetics_helpers.nf`, then dispatches on `params.mod`. For `params.mod` starting with `wheat_`, skips VCF input and runs `RUN_WHEAT_INTEGRATED` (`subworkflows/local/wheat/wheat_integrated_study.nf`). |
| `subworkflows/local/entry/genetics_helpers.nf` | `header()`, `helpMessage()`, and `workflow check_input` (VCF channel + summaries). |
| `subworkflows/local/plink/plink_genotype_modes.nf` | Composed workflows `RUN_V1_PLINK`, `RUN_TEST_PLINK`, `RUN_TEST_PLINK_CAMP`, `RUN_TEST_COMMON_THIN` (processor + stats only). Router-gap modules: `plink/analysis_extensions.nf`. |
| `subworkflows/local/entry/partial_router.nf` | Partial reruns via `--partial_task` (assess, stats redraw, wheat-from-plink). |
| `nextflow.config` | Default `params`, resource labels, conda hints. Pass with **absolute** `-c` when launching auxiliary scripts from arbitrary cwd. |

Typical invocation shape (full examples in **`doc/progress/`** — search `#### NF replay`):

```text
nextflow run <path-to>/workflow/Genetics/main.nf \
  -c <absolute-path-to>/workflow/Genetics/nextflow.config \
  --home_dir … --user_dir … --src_dir … --output_dir … --mod <mod> --job <job> [optional flags]
```

---

## Active `params.mod` values (`main.nf` router)

These branches are wired in the **anonymous** top-level workflow in `main.nf`. PLINK tracks below run end-to-end genotype processors plus stats; `wheat_*` modes run the integrated table/matrix workflow only (no VCF channel). Details: **Wheat integrated workflows** (section below).

| `params.mod` | Processor workflow | Stats |
| --- | --- | --- |
| `v1_plink` | `PLINK_PROCESSOR` (`plink_processor`) | `PLINK_STATS` (`plink_stats`) |
| `test_thin` | `TEST_PLINK_PROCESSOR` (`test_plink_processor`) | `TEST_PLINK_STATS` (`test_plink_stats`) |
| `test_camp` | `TEST_PLINK_CAMP_PROCESSOR` + **`params.camp`** TSV | `TEST_PLINK_STATS` |
| `test_common_thin` | `TEST_COMMON_THIN_PROCESSOR` (hard-filter + thin from `hf_*` params) | `TEST_PLINK_STATS` |
| `test_common_inter` | `TEST_COMMON_INTER_PROCESSOR` (exclude all-group-rare then thin; reuse merged pfiles when present) | `TEST_PLINK_STATS` |
| `test_common_inter_refchr` | Recode `test_common_inter` subgenome pfiles onto linear vmap4 chromosomes (chr3A = PLINK 13+14, …), republish `{A,B,D,Others}_test.plink2`, then full stats + along-chr r2=0.2 for `ld_along_chr_chromosomes` (default chr1–chr7 × A/B/D, one stacked PNG per subgenome). Source `--process_dir` defaults to `process/test_common_inter` (read-only). | `TEST_PLINK_STATS` + `variant_ld_r2_decay_along_subgenome_plot` |
| `test_wap_refchr` | Recode mixed-MAF `test_wap` subgenome pfiles onto linear vmap4 chromosomes (same maps as `test_common_inter_refchr`). Republish `{A,B,D,Others}_test.plink2` under `process/test_wap_refchr` — does not overwrite `test_wap` or `test_wap_common_refchr`. Rebuilds basic info + windowed unphased LD, then overall LD decay plus rare-rare / rare-common / common-common class plots. Source `--process_dir` defaults to `process/test_wap` (read-only). | `variant_ld_decay_plot` + `variant_ld_decay_rare_common_plot` |
| `test_wap_thin_refchr` | Random `--thin` of existing `test_wap_refchr` pfiles to **2–3M ABD sites** (`test_set_target_n` / `test_set_max_n`). One keep rate from the A+B+D total, applied to every chromosome (A/B stay larger, ~1M each; D smaller). Does not overwrite `test_wap` or `test_wap_refchr`. Then basic info + windowed LD + overall / class LD plots. Source `--process_dir` defaults to `process/test_wap_refchr`. | `variant_ld_decay_plot` + `variant_ld_decay_rare_common_plot` |
| `test_hex_thin` | Same 2–3M ABD target and one global keep rate on existing `test_hex` pfiles. Does not overwrite `test_hex`. Source `--process_dir` defaults to `process/test_hex`. | `variant_ld_decay_plot` + `variant_ld_decay_rare_common_plot` |
| `test_wap_common_refchr` | From `main_raw`: `--keep` WAP, `--maf test_wap_common_maf` (default 0.05), `--thin test_wap_common_thin_rate` (default 0.01895, ~2M sites), merge subgenomes, recode onto linear vmap4 chromosomes. Does not overwrite `test_wap`. Source `--process_dir` defaults to `process/main_raw`. Along-chr r2=0.2 for `ld_along_chr_chromosomes` (default chr1–chr7 × A/B/D; one stacked PNG per subgenome). | `variant_ld_r2_decay_along_subgenome_plot` (no full `TEST_PLINK_STATS`) |
| `test_wap` | `TEST_WAP_PROCESSOR` (`--keep` WAP sample_group, `--mac 2`, `--thin test_wap_thin_rate`; input `--process_dir` = `main_raw`) | `TEST_PLINK_STATS` |
| `test_wap_exon` | From `main_raw`: `--keep` WAP, `--extract bed0` of merged Lulab GFF **exon** intervals (`test_wap_exon_gff_type`), `--mac 2`, **no thin**. Merge A/B/D only (no organelle/unplaced exons). Basic info, no LD/MQ/popdep. Does not overwrite `test_wap`. Source `--process_dir` defaults to being required (`process/main_raw`). | none (VEP / SnpEff / burden via partials) |
| `wheat_*` | *(none — table/matrix integrated track)* | `integrated_wheat` (`modules/local/integrated/integrated_wheat.nf`) |

Single-table wheat modes use `wheat_table_input`; `wheat_gwas` / `wheat_kgwas` use dedicated genotype or k-mer + phenotype params (see **Wheat integrated workflows** below).

Any other `params.mod` causes `main.nf` to **exit 1** with a list of supported mods (no silent no-op).

### Included but not routed (known gap)

`analysis_extensions.nf` collects router-gap modules (`database`, `kinship`, `ps`, `GWAS`, `HAIL`) for future `main.nf` branches — they are **not** loaded by the active PLINK router. Track router-gap work in GitHub Issues.

### Process libraries (`modules/local/`)

**Active** Nextflow process definitions live under **`workflow/Genetics/modules/local/genotype/`**, **`.../dynamic/`**, **`.../static/`**, and **`.../integrated/`**. Shared DSL helpers live under **`.../infra/`**. Earlier copies or experiments may remain under **`workflow/Genetics/old/`**; new work should `include` from **`modules/local/...`** only. Cross-theme imports (e.g. GWAS pulling `format_vcf_plink` from genotype) use paths relative to the importing file under `modules/local/`.

**Infra library** (`modules/local/infra/`) — cross-cutting Groovy helpers (no `process` blocks):

| File | Contents |
| --- | --- |
| `infra/utils.nf` | Index / compose-only barrel (do not `include { def }` from here; use `infra_*.nf` directly) |
| `infra/infra_tools.nf` | Java setup, software/TIGER/samtools paths, `validatePaths` |
| `infra/infra_job_config.nf` | `getJobConfig`, calling/pop/chr config, `getVcfIdFromPath` |
| `infra/infra_tiger.nf` | TIGER JAR config, FastCall server population lists |
| `infra/infra_ref_v1.nf` | vmap4 IWGSC v1 chr maps, offsets, lengths, taxaBam paths |
| `infra/infra_plink_reuse.nf` | Merged test pfile reuse tuples (`hasMergedSubgenomeTestPfiles`, …) |

**Genotype processor library** (`modules/local/genotype/processor/`) — heavy `plink2` / TIGER / bcftools compute; workflows only in the barrel file:

| File | Contents |
| --- | --- |
| `processor/processor.nf` | Composed workflows `test_plink_processor`, `test_plink_camp`, `test_common_thin_processor`, `plink_processor`, `plink_preprocess`, `vcf_arrange_merge`, `vcf_bcftools_filter` |
| `processor/processor_vcf.nf` | VCF format, arrange, merge, `format_vcf_plink` |
| `processor/processor_test.nf` | Test subsampling and subgenome merge |
| `processor/processor_plink2.nf` | PLINK2 basic info, PCA, tag-SNP prune, CNV awk, LD matrices |
| `processor/processor_legacy.nf` | Legacy vcftools basic info |
| `processor/processor_depth.nf` | Population depth (TIGER): `calc_population_depth` (PopDepFull, per chr), `calc_population_depth_crosschr` (PopDepCrossChr), shared `popdep_tiger_gz_to_bgzip_tabix` (gzip → BGZF + tabix, label `popdep_bgz_tabix`) |
| `processor/processor_filter.nf` | Sample/variant PLINK filters; VCF v0 filters |
| `processor/processor_assess.nf` | Assess compute (`plink2_assess_debug_slice`, `quick_count`, `bcftools_qc_assess`) |

**Genotype stats library** (`modules/local/genotype/stats/`) — plotting and summarization only; no heavy `plink2` compute:

| File | Contents |
| --- | --- |
| `stats/stats.nf` | Composed workflows `test_plink_stats`, `plink_stats` (includes all sub-libraries below) |
| `stats/stats_variant.nf` | Variant QC, MAC, popdep, LD plot processes |
| `stats/stats_sample.nf` | Sample QC, kinship/IBS, PCA helper processes |
| `stats/stats_assess.nf` | Assess debug plot processes (`assess_plink_debug_plots`, `dumpnice_vcf_qc_assess`) |
| `stats/stats_integrated.nf` | Wheat integrated PLINK2 plot/report processes |
| `stats/stats_chr_report.nf` | Chr variant count tables and thin vs common-thin compare plots |

**Upstream genotype libraries** (not routed from `main.nf` yet — standalone entry scripts; composed in `subworkflows/local/upstream/raw_data_upstream.nf`):

| Path | Entry | Contents |
| --- | --- | --- |
| `genotype/calling/` | `calling/caller.nf` | FastCall3 workflows (`run_FastCall3`, `load_lib_files`) + `caller_prep.nf` / `caller_fastcall.nf`. Default jar **`TIGER_F3_20261012.jar`**. Logs on stderr; scan/scan2/scan3/scan4 VCF on stdout. Mods: `disc` → `blib` → `scan`; optional `vlib` / `clib`. `--mod scan_iac` starts from an existing `.lib.gz` (same entry as `scan_only`) and runs Java `-mod scan_iac`: same pileup parameters as `scan`, writes `{gen}/indiCounts/<taxon>/<chrom>_<binStart>_<binEnd>.iac.gz`, then stops (no genotype call, no VCF). Those `.iac.gz` files feed `scan2`, `scan3`, and `scan4`. `scan2` rewrites VCF from existing `.iac` (whole bin in RAM); `scan4` same job in 256-site windows; `scan3` sums source-taxa `.iac` AD. VCF filter: `--scan_midr` (default 0.2). This jar has no MaDR (`-n` is unrecognized) and no absolute MDC. Scan/scan2/scan3/scan4 pipe Java stdout through `bgzip -c` to `{gen}/VCF/chrNNN.vcf.gz`. `scan_iac` uses the same resources as `fastcall3_scan`. scan2/scan3: 32 CPUs, 640 GB + 128 GB/retry. scan4: 8 CPUs, 32 GB + 16 GB/retry (local executor 128 CPUs → up to 16 concurrent scan4). |
| `genotype/align/` | `align/align.nf` | `RUN_ALIGN_USB_TRANSFER` workflow; BAM USB transfer + MD5 (`align_transfer.nf`, `align_md5.nf`); BWA-MEM2 (`align_bwa.nf`) |

**Hail library** (`genotype/hail/`) — included as `HAIL` workflow (router gap; not wired in `main.nf`):

| File | Contents |
| --- | --- |
| `hail/hail.nf` | Composed workflow `HAIL` |
| `hail/hail_io.nf` | `vcf_to_mt_hail`, `filter_hail` |
| `hail/hail_stats.nf` | QC, kinship, PCA processes |
| `hail/hail_gwas.nf` | `hail_normal_gwas` |

**Genotype database library** (`genotype/database/`) — VCF → Hail MT + bcftools annotation helpers:

| File | Contents |
| --- | --- |
| `database/database.nf` | Composed workflow `database` |
| `database/database_build.nf` | `build_hail_from_vcf` |
| `database/database_annotate.nf` | bcftools annotate / tab / BED helpers |

**Static GWAS library** (`modules/local/static/gwas/`):

| File | Contents |
| --- | --- |
| `gwas/gwas.nf` | Composed workflow `GWAS` (legacy PLINK1 track) |
| `gwas/gwas_legacy.nf` | PCA, pheno/covar prep, PLINK1 / GAPIT / rMVP processes |
| `gwas/gwas_plink2.nf` | `plink2_gwas_glm`, `gcta_gwas` (wheat integrated compute) |
| `gwas/gwas_plot.nf` | `plot_gwas_association`, `plot_gwas` |

**Include paths:** from `subworkflows/local/*/` use `../../../modules/local/...`; from `modules/local/genotype/*` use `../../infra/<file>.nf`. Partial reruns: launch `subworkflows/local/entry/partial_router.nf` with `--partial_task`, or `include` the smallest `partial/*.nf` / module sub-file.

---

## Important `params` (from `nextflow.config`)

| Parameter | When to set | Notes |
| --- | --- | --- |
| `process_dir` | Reuse prebuilt per-chromosome or merged test pfiles | When merged `A_test.plink2` … `Others_test.plink2` exist under `{process_dir}`, most test processors **skip** thin/merge and rebuild basic info + LD only. **`test_common_inter_refchr`** reads merged pfiles as source (default `process/test_common_inter`). **`test_wap_refchr`** recodes mixed-MAF `test_wap` (default `process/test_wap`) onto linear chroms and publishes under `process/test_wap_refchr`. **`test_wap_thin_refchr`** thins `test_wap_refchr` to 2–3M ABD sites with one keep rate (default source `process/test_wap_refchr`). **`test_hex_thin`** does the same from `process/test_hex`. **`test_wap_common_refchr`** extracts from per-chr `main_raw` pfiles (default `process/main_raw`; `--keep` WAP, `--maf`, `--thin`) and publishes recoded pfiles under `process/test_wap_common_refchr` — it does not overwrite `test_wap` or `main_raw`. See **`doc/project_knowledge/domain/data_publish_tree.yaml`**. |
| `mq_dir` | Site MQ reference (partial only) | Frozen **`{mq_dir}/reference/chrNNN.site_mq.{ref,calls}.*`** under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/abstract_mq_50_bams`. Built once via **`--partial_task abstract_mq_50_bams --job test_plink`**; skipped when refs exist (`mq_force_rerun=false`). Downstream Python: `genetics.genomics.variant.mq.site_mq_ref_path()`. |
| `popdep_dir` | Population depth reference (partial + test annotate) | Frozen **`{popdep_dir}/variant/chrNNN.popdep.txt.bgz`** (+ `.tbi`) under `/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/main_raw`. Legacy plain `*.popdep.txt` is still accepted for reuse/skip. Built once via **`--partial_task main_raw_popdepth --mod main_raw_popdepth --job test_plink`** (default **`PopDepCrossChr`** in **`TIGER_PD_20260619.jar`**: one pass over **`all.ALL.taxaBamMap.txt`** (FastCall-style 3-column map; per-chr **`all.{A,B,D,ALL}.taxaBamMap.txt`** for PopDepFull); length file **`Chr\\tLength\\tnTaxa`**; **`popdep_crosschr_threads`** default 16, **`popdep_crosschr_memory_gb`** default **640** (relativeDepth in 20260619); **`popdep_tiger_memory_gb`** default **640** for PopDepFull; **`popdep_crosschr_taxa_order`** (`interleave` \| `shuffle` \| `sorted`, TIGER `-o`); optional **`popdep_min_mapq`** (TIGER `-g` / samtools `depth -q`, default **0**; **`TIGER_PD_20260818.jar`**); optional **`popdep_crosschr_checkpoint_dir`** + **`popdep_crosschr_checkpoint_interval`** (TIGER `-k` / `-ci`, resume with same paths); override **`--popdep_tiger_app PopDepFull`** for per-chr scans); skipped when refs exist (`popdep_force_rerun=false`). Test modes extract **`process/{mod}/variant/{A,B,D,Others}.popdep.info.tsv`** via `annotate_subgenome_variant_popdep` (tabix lookup). Python: `genetics.genomics.variant.popdep.popdep_chr_ref_path()`. Bench: **`popdep_taxa_bam_file`**, **`popdep_chr_list`**, **`popdep_tiger_max_forks`** — see **`10stats.genome/10run_popdep_head2head/`**. |
| `test_site_target_bp` | **`test_site` mod** (tmp build) | Random genome site panel size (default **10_000_000** bp). Build via **`subworkflows/tmp/test_site_build.nf`**; popdep-only stats via **`test_site_popdep_stats.nf`**. Panel + popdep sidecar under **`process/test_site/sites/`**; variant pfiles under **`process/test_site/variant/chrNNN.site.plink2`**. Stats plots highlight **`Is_Variant`** sites (grey vs red). |
| `test_popdep_target_bp` | **`test_popdep` mod** (tmp build) | Random **100_000_000** popdep sites (not variants). Build + class-coloured popdep plots via **`subworkflows/tmp/test_popdep_build.nf`** only — no pfile extract, no variant stats. Classes: **gap** (AGP BED), **repeat** (IWGSC v1.0 CLARITE TE GFF), **gc** (extreme-GC 1 kb windows from the v1.0 FASTA; IWGSC does not ship a GC track), **cnv** (1 Mb mapq_bin_popdep bins with `cnv_frac ≥ 0.2`; IWGSC v1.0 has no CS CNV track), else **normal**. Publish **`process/test_popdep/{annotation,sites}/`** and **`stats/test_popdep/`**. Plots: mean-depth + relative-mean histograms; Depth and RelativeDepth SD/CV vs mean JointGrids. To invert June CrossChr **zero-padded nTaxa** on those site tables (44 chr → n=500; chr032 unchanged) and redraw plots, use **`subworkflows/tmp/test_popdep_ntaxa_rescale.nf`**. Plot-only redraw from existing sites: **`subworkflows/tmp/test_popdep_class_replot.nf`**. |
| `camp` | **`test_camp` only** | Path to cohort map TSV (e.g. `camp_vmap4_map.tsv`). Required for that mod. |
| `chr` | Optional | Comma-separated chromosome filter after job config builds the VCF channel. |
| `wheat_table_input`, `wheat_gwas_*`, `wheat_kgwas_*` | **`params.mod` starts with `wheat_`** | See **Wheat integrated workflows** below. Requires `user_dir` for the `stats` Conda env (same as genotype stats processes). |

Hard-filter / common-thin tuning for `test_common_thin` uses `params.hf_maf`, `params.hf_geno`, `params.hf_thin_rate` (see `nextflow.config`; retune and log when cohort policy changes).

---

## Output layout (test jobs)

For test modes (e.g. `params.job` = `test_plink` while `params.mod` is `test_thin` / `test_camp` / `test_common_thin`), **`params.job`** is the **top-level** folder under `params.output_dir`. Genotype outputs go under **`{output_dir}/{job}/process/{mod}/`** (sample, variant, info, logs, etc.). Stats outputs go under **`{output_dir}/{job}/stats/{mod}/`** (plots, thresholds, info, logs). Assess debug (`partial_router.nf --partial_task assess_plink`) uses the same trees: compute artefacts under **`process/{mod}/`**, plots and summary tables under **`stats/{mod}/`**. Wheat integrated plots/tables from merged PLINK2 (`partial_router.nf --partial_task wheat_from_plink` / `WHEAT_STUDY_FROM_PLINK`) publish under **`{output_dir}/{job}/integrated/{plink_source_mod}/{wheat_task_mod}/`** (e.g. `integrated/test_thin/wheat_pca_tsne/{info,plots}`). Table-only `wheat_*` modes from **`main.nf`** still use **`{output_dir}/{job}/integrated/{params.mod}/`**. Do not encode `mod` by changing `params.job` alone—the pipeline uses **`params.mod`** under those trees where applicable.

---

## Composed subworkflows (`subworkflows/local/`)

| Folder | Files | Role |
| --- | --- | --- |
| `entry/` | `genetics_helpers.nf`, `partial_router.nf` | `main.nf` helpers; partial rerun router (`--partial_task`) |
| `plink/` | `plink_genotype_modes.nf`, `analysis_extensions.nf` | PLINK router workflows; router-gap module index |
| `wheat/` | `wheat_integrated_study.nf` | `RUN_WHEAT_INTEGRATED`, `WHEAT_STUDY_FROM_PLINK`, `RUN_WHEAT_FROM_PLINK` |
| `upstream/` | `raw_data_upstream.nf` | `RUN_ALIGN_USB_TRANSFER`, `run_FastCall3` |
| `partial/` | `partial_assess.nf`, `partial_stats.nf` | Assess + stats partial workflow bodies |

---

## Auxiliary scripts (`subworkflows/tmp/`)

**Partial reruns** — one entry, no per-task duplicate scripts under `workflow/Genetics/tmp/` (removed 2026-06-09):

```bash
nextflow run /data/home/tusr1/git/script/workflow/Genetics/subworkflows/local/entry/partial_router.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --partial_task ld_redraw --mod test_thin --output_dir … --job test_plink …
```

| `--partial_task` | Workflow | Typical `params.mod` / notes |
| --- | --- | --- |
| `assess_plink` | `RUN_ASSESS_PLINK_DEBUG` | `test_thin` / `test_common_thin` |
| `assess_vcf` | `RUN_ASSESS_VCF_DEBUG` | any mod with `export/*.debug.vcf.gz` |
| `ld_redraw` | `RUN_LD_PLOTS_REDRAW` | e.g. `test_thin` |
| `ld_rare_common` | `RUN_LD_RARE_COMMON` | mixed-MAF mods (`test_thin`, `test_wap_refchr`, `test_wap_thin_refchr`, `test_hex_thin`) — rare-rare / rare-common / common-common LD r2 |
| `mac_stats` | `RUN_MAC_STATS_FROM_GCOUNT` | test mods or `test_rebulld_lib` |
| `mac_het_excess` | `RUN_MAC_HET_EXCESS` | `test_rare_only` (gcount) + `het_source_mod` common pfiles for `--het` |
| `mac_dist_redraw` | `RUN_MAC_DIST_LOG_REDRAW` | multi-mod batch |
| `chr_counts` | `RUN_CHR_VARIANT_COUNTS` | per `params.mod` |
| `chr_compare` | `RUN_CHR_VARIANT_COMPARE` | thin vs common-thin |
| `rebuild_lib_stats` | `RUN_REBUILLD_LIB_STATS` | `test_rebulld_lib` only |
| `mq_missing_reg` | `RUN_MQ_MISSING_REG` | `test_thin` / `test_common_thin`; needs `*.mq.info.tsv` |
| `popdep_missing_reg` | `RUN_POPDEP_MISSING_REG` | `test_thin` / `test_common_thin`; needs `*.popdep.info.tsv` |
| `abstract_mq_50_bams` | `RUN_ABSTRACT_MQ_50_BAMS` | `--mod abstract_mq_50_bams`; one-off site MQ build |
| `main_raw_popdepth` | `RUN_MAIN_RAW_POPDEPTH` | `--mod main_raw_popdepth`; one-off TIGER popdepth per chr → `main_raw/variant/` |
| `aneuploidy_mapq` | `RUN_ANEUPLOIDY_MAPQ` | `--mod aneuploidy_check`; full-cohort BAM MAPQ + 1 Mb bin depth TSVs |
| `aneuploidy_breachr` | `RUN_BREACHR_FLAGGED` | `--mod aneuploidy_check`; breachr segment calls on flagged samples with existing MAPQ TSVs + mosdepth reconciliation |
| `wheat_from_plink` | `RUN_WHEAT_FROM_PLINK` | `test_thin` / `test_common_thin` + `wheat_integrated_mod` |
| `sample_manifest` | `RUN_SAMPLE_MANIFEST` | `--mod sample_manifest`; WAP assoc keep + WGS-core keep + RefSeq v1.0 inventory |
| `cohort_freq` | `RUN_COHORT_FREQ` | merged test pfile mods; PLINK2 `--freq` per cohort keep list |
| `annotate_vep` | `RUN_ANNOTATE_VEP` | `--vep_from_mod <mod>` (sites VCFs from pfile IDs) or `--vep_vcf <one VCF>`; conda env **`anno`**, local GFF, no Ensembl cache |
| `annotate_snpeff` | `RUN_ANNOTATE_SNPEFF` | Same sites VCFs as VEP; reuses `{home_dir}/00data/03ref/07snpeff` (`wheat_lulab_v1.1`, SnpEff **4.3t** jar) when present; otherwise builds from GFF+FASTA. Conda **`anno`** supplies Java/bgzip |
| `site_anno` | `RUN_SITE_ANNO` | `--site_anno_from_mod <mod>`; streams `00data/11epi_anno/chr*_featHighRes.bed.gz` for conservation / epigenome features per variant |
| `gene_expr` | `RUN_GENE_EXPR` | `--gene_expr_from_dir <dir>` or `--gene_expr_tpm <files>`; TPM atlases normalized onto v1.1 gene IDs |
| `variant_sets` | `RUN_VARIANT_SETS` | `--variant_sets_from_mod <mod>` (per-subgenome fan-out) or a published `*.coding_mask.tsv` + `*.dual_af.tsv` pair; optional `site_anno` / `gene_expr` tables add more sets |
| `assoc_wap` | `RUN_ASSOC_WAP` | merged test pfile mods + WAP phenotype; `--assoc_wap_burden` switches to the gene burden track |
| `phenotype_wap` | `RUN_PHENOTYPE_WAP` | `--phenotype_wap_input <wide table>`; aligns traits to `assoc.keep` |

### Cohort-to-association chain

The five tasks above form one ordered chain; each consumes the previous step's
published artefacts, so a rerun starts wherever inputs changed.

```text
sample_manifest ──> cohort_freq ──> (dual_af)  ─┐
                                                ├──> variant_sets ──> assoc_wap
annotate_vep ──> (coding_mask) ─────────────────┤                        ▲
annotate_snpeff ─ (comparison table; not burden)┤                        │
site_anno ──────────────────────────────────────┤                        │
gene_expr ──────────────────────────────────────┘                        │
phenotype_wap ───────────────────────────────────────────────────────────┘
```

`annotate_vep` is the only task in Conda env **`anno`** (`environment_anno.yml`).
Its GFF run yields no SIFT/PolyPhen columns — wheat has no Ensembl cache — so
`coding_mask` resolves missense rows to `protein_altering` until a deleteriousness
predictor is added.

`site_anno` and `gene_expr` are optional: `variant_sets` runs without them and
simply skips the sets that need them, logging which ones it dropped.

### Ref-chr chromosome naming

The ref-chr recode merges each chromosome's two published arms into one linear
chromosome. Until 2026-09-03 it labelled the result with the **first arm's numeric
id** — chr2A appeared as `7` — which is unambiguous but reads as "chromosome 7" on
every plot axis. It now writes **chromosome names**: `chr1A`…`chr7D`, plus `chrUn`,
`Mit`, `Chl` for the `Others` subgenome.

The variant `ID` column is untouched by the recode and keeps its
`{original arm}-{local pos on that arm}` form, so VEP, `coding_mask`, `dual_af`,
`site_anno` and `variant_sets` tables join across the change.

Both codings still exist on disk, and readers accept either:
`load_pvar_ref_variants` detects names versus arm numbers, and `site_anno_extract`
matches a pvar row on either token. `test_wap_refchr`, `test_wap_common_refchr` and
`test_wap_thin_refchr` are renamed; `test_common_inter_refchr` is still numeric.

To rename an already-merged mod without rebuilding it — important for the common and
thin mods, whose random subsamples a rebuild would redraw — use
**`subworkflows/tmp/refchr_rename.nf`**, which republishes over the same mod with POS
and IDs unchanged.

June 2026 `PopDepCrossChr` grids for the N500 GAM-credible panel wrote FastCall map
sizes into `chromosomeLength` nTaxa (A=7675, B=7594, D=7738, Others=8285) while
scanning 500 BAMs, so 44 chromosomes were diluted by zero-padded taxa. chr032 used
`PopDepFull` with nTaxa=500. **`subworkflows/tmp/test_popdep_ntaxa_rescale.nf`**
inverts mean/SD (population second-moment identity) on published
`process/test_popdep/sites` tables to n=500 and redraws `stats/test_popdep` plots.
Class histograms are pre-binned (default mean-depth 0–40 / relative 0–5, 100 bins)
and scatter JointGrids (Depth and RelativeDepth SD/CV vs mean) use a 2M stratified
subsample; the 100M-row panel is not concatenated. Plot-only redraw:
**`subworkflows/tmp/test_popdep_class_replot.nf`**.

### Rare-variant gene burden

`--partial_task assoc_wap --assoc_wap_burden --assoc_wap_burden_set <set>` collapses
every variant a gene carries in one `variant_sets` set into a single per-sample score
and regresses the trait on it. Single-variant GLM has no power on a variant three
accessions carry; aggregating the gene's variants is the standard answer.

The four design choices the track had left open are settled as follows.

* **Grouping** is the VEP `Gene` column of `stats/{mod}/info/{id}.variant_sets.tsv`.
  The `*.extract` lists carry variant IDs only, so the gene map has to come from the
  set table — both are resolved per subgenome, and a run fans out over A/B/D/Others.
* **Homoeologs are not collapsed.** A/B/D live in separate pfiles and no triad map is
  wired in, so a triad's three copies are three independent genes.
* **Collapsing** is `--assoc_wap_burden_collapse count` (sum of minor alleles across
  the gene, the default) or `carrier` (0/1 indicator).
* **No frequency weighting.** Madsen-Browning or beta(1,25) weights only pay off when
  a gene holds many variants across a range of frequencies.

plink2 applies the keep list and the extract list and writes a variant-major
additive `.traw` (`--export Av`); the collapsing and the regression run in
`genetics.gwas.burden`, because plink2 has no gene-collapsing test. Each `.traw`
row is one site, so frequency, the major-ALT flip, and the per-gene sum happen
in a single scan — a set like `protein_altering_rare_both` (~2.6M sites × the
WAP panel) never loads the full dosage matrix. Covariates are the same shared structure PCs the
single-variant track uses. Output is
`process/{mod}/assoc/{id}.{trait}.{set}.burden.tsv`. Gene-level Manhattan and QQ
go to `stats/{mod}/plots/{id}.{trait}.{set}.{manhattan,qq}.png`; the dashed line
is Bonferroni `0.05 / n_genes`, not single-variant `5e-8`. Coordinates come from
`params.ref_v1_gff_hc` (or `--assoc_wap_burden_gff`). To redraw existing tables
without re-exporting genotypes, use `subworkflows/tmp/burden_replot.nf`.

Two properties of the score are worth knowing when reading that table. plink2
`--export A` counts the ALT allele, which for a rare pLoF site is usually the
*major* one, so those columns are flipped to count minor alleles — and the flip
happens before missing calls are zero-filled, since `2 - NaN` must stay missing
rather than become a homozygote. Missing calls then score as non-carrier, which
biases a poorly called gene toward the null; `mean_missing_rate` is reported per
gene so that bias is visible.

**Check the aggregation factor before trusting a burden run.** The test only buys
power when genes hold several variants. On `test_wap_refchr` / `plof_rare_wap` the
sets average **1.01 variants per gene**, so the run is a single-variant scan wearing
a burden coat. Widening the mask (`protein_altering_rare_both`) only reaches 1.34.
The `[Info] collapsed …` log line reports this figure per run.

**Ops** (`subworkflows/tmp/ops/`) — self-contained FTP upload helpers (not genetics analysis):

| Script | Purpose |
| --- | --- |
| `ops/gsa_ftp_upload.nf` | GSA FTP upload (`submit.big.ac.cn`); set `GSA_FTP_PASSWORD`. |
| `ops/gwh_ftp_upload.nf` | GWH FTP upload; optional `--verify_md5` / `--verify_only`. Set `GSA_FTP_PASSWORD`. |

**Plot-only replots** (read published TSVs; no processor rerun):

| Script | Purpose |
| --- | --- |
| `breachr_reconcile_replot.nf` | Replot breachr mosdepth reconciliation summary charts from `info/breachr_check/breachr_reconcile_all.tsv`. |
| `mosdepth_aneuploidy_replot.nf` | Replot mosdepth rel_depth distributions from `info/chr_depth/*.depth.tsv`. |
| `aneuploidy_chromosome_plot.nf` | Stacked bar of confirmed whole-chromosome aneuploidy (nullisomy / monosomy / trisomy / tetrasomy) per chromosome from `info/breachr_check/breachr_reconcile_all.tsv`. Params: `--aneuploidy_breachr_input_dir`, `--aneuploidy_plot_id`, `--aneuploidy_plot_title`. |
| `ld_rare_common_decay_replot.nf` | Replot `*.variant.ld_decay.rare_common.png` from published `*.variant.ld_decay.rare_common.info.tsv` (no vcor). Rare-common line uses alpha so overlap with opaque rare-rare is visible. |
| `ld_pplus_vs_q_overall_replot.nf` | Overall \(P(D'\approx +1 \mid |D'|=1)\) vs background MAF \(q\) from `*.maf_scatter.info.tsv`, restricted to `--pplus_vs_q_max_kb` (default 1000 = 1 Mb). Writes `*.dprime.signed_pplus_vs_q.{suffix}.png` (default suffix `1mb`). Does not recompute dosages or overwrite the default wider-window P+ vs q figure. |
| `ld_maf_scatter_r2max_replot.nf` | Replot `*.maf_scatter.png` / `*.maf_scatter_0.2.png` plus rare (both MAF < 0.05) and ultra-rare (both MAF < 0.01) zooms from `*.maf_scatter.info.tsv` with \(r_{\max}^2=0.8,0.6,0.2\) contour pairs. Recomputes sign-aware `|D'|=\|D\|/D_{\mathrm{den}}\) from pairwise-complete pfile dosages (`D` / `D_PRIME` / `DISTANCE_KB`) when `D_PRIME` is absent, and writes `*.dprime.png` / `*.dprime.scatter.png` / `*.dprime.dist.png` / `*.dprime.dist.log.png` / `*.dprime.info.tsv`. Density lines are fine-bin frequency polygons (400 bins on [0, 1]), one curve per pair class. The log panel plots \(\log_{10}(\mathrm{density}+1)\) (empty bins stay at 0). Signed \(D'=D/D_{\mathrm{den}}\) (minor-allele coding, not recovered from \(|D'|\)) is classified with \(\varepsilon=10^{-6}\) into \(\pm 1\), interior \(+/-\), and 0. Figures: `*.dprime.signed_prop.png` (\(P_{+1},P_{-1},P_{\mathrm{other}}\) vs distance, rare-rare / rare-common facets), `*.dprime.signed_boundary_stack.png` (100% stack of \(+1\) vs \(-1\) among \(|D'|=1\)), `*.dprime.signed_interior_dist.png` (signed \(D'\) density excluding the \(\pm 1\) boundary). Do **not** recover |D'| from r2. |
| `mapq_bin_popdep.nf` | Cohort 1 Mb popdep grid from published `info/mapq_check/*.mapq_depth_bin.tsv`. Default **unfiltered** depth (`bin_mean_depth`, `-c conf/mapq_bin_popdep_all.config`, prefix `allmapq_1mb`); legacy MQ30+ approx via `-c conf/mapq_bin_popdep.config`. Publishes `info/mapq_popdep/` and `plots/mapq_popdep/`. `--mapq_popdep_plot_only true` replots from existing grid TSV. |
| `lof_vs_gene_length_replot.nf` | Per-gene LoF count vs Lulab CDS / gene span from a LoF site table (`--lof_table`) plus `--lof_ancestry_png`. Writes `stats/<mod>/plots/all.gene.lof_vs_{cds,cds.nonzero,span.nonzero,cds.joint}.png` and `all.gene.lof_ancestry.class.png`. |
| `rare_vs_gene_length_replot.nf` | Per-gene `rare_wap` (WAP MAF<1%) count vs Lulab CDS / gene span from published `{A,B,D}.variant_sets.tsv`. Writes `stats/<mod>/plots/all.gene.rare_{vs_cds,per_kb_vs_cds,vs_cds.nonzero,vs_span.nonzero,vs_cds.joint}.png`. Optional `--variant_sets_dir`. |

**Compute tmp** (processor-style annotation from a vmap4 run folder):

| Script | Purpose |
| --- | --- |
| `popdep_annotate.nf` | Tabix-annotate merged `*_test.plink2` variants from frozen `chrNNN.popdep.txt.bgz` grids (`--popdep_dir`, typically `_popdep_n500_gam_credible`) into `process/<mod>/variant/*.popdep.info.tsv`. |
| `mosdepth_single_bam.nf` | Per-BAM `mosdepth -n` (optional `--mosdepth_mapq`, `--mosdepth_cpus`, `--mosdepth_max_forks`). Publish under `--depth_publish_dir` in `00data/04depth/` (new tree for MAPQ-filtered runs; do not overwrite unfiltered `{cohort}/`). Analysis plots stay on data1. |
| `ld_r2_decay_along_chr.nf` | In each `ld_along_chr_spatial_bin_mb` (default 5 Mb) window, sample `ld_along_chr_n_focal_per_bin` (default 100) focals. Compute unphased r2 against partners in the full `ld_along_chr_max_kb` window (not limited to the 5 Mb bin). Plot the *mean* distance where mean r2 falls to `ld_along_chr_r2` (default 0.2). Non-default r2 writes `*.ld_r2_decay_along.r2_*` (plus `*.chr_summary.info.tsv`) so default 0.2 artefacts stay. Default chromosomes are chr1–chr7 × A/B/D; one stacked PNG per subgenome (chr1 at the top through chr7). Unreached bins are drawn at the window cap in crimson with dashed connectors. `--ld_along_chr_plot_only true` replots from existing `variant/*.focal.ld.vcor`. Split-arm pfiles are merged in the work dir; `test_common_inter_refchr` / `test_wap_common_refchr` pfiles are already linear. |
| `ld_decay_from_pfile.nf` | From merged `*_test.plink2`, run `calc_plink_ld_unphased` + `calc_plink_ld_crosschr_random`, then `variant_ld_decay_plot` / `variant_ld_crosschr_plot`. For mods that skipped all-pairs LD (e.g. `test_wap_common_refchr`). Decay PNG: one coloured mean-r2 curve per reference chromosome (chr1A…chr7A) plus the overall mean as a thick dashed line. `--ld_decay_plot_only true` replots from existing `variant/*.info.ld.vcor`. |
| `ld_decay_r2_threshold.nf` | Fine-bin (default 1 kb, 0–100 kb) all-pairs mean r2 from existing `variant/{A,B,D}.info.ld.vcor` and the per-chromosome distance where that mean first falls to `ld_decay_r2_threshold` (default 0.5). Writes `*.ld_decay.r2_0p5.1kb.{png,by_chr.info.tsv,chr_summary.info.tsv}` — does not overwrite the default 50 kb `*.ld_decay.*` artefacts. |
| `ld_homoeolog_compare.nf` | From published `*.variant.ld_decay.by_chr.info.tsv`, plot A/B/D parallel-coordinate compare (one colour per homoeolog number) of pair-weighted mean r2 and of within-subgenome rank. |
| `ld_maf_bin_decay.nf` | Within-bin MAF LD decay from `process/<mod>/variant/{A,B,D}.info.ld.vcor` + `*.info.afreq`. Each curve uses only pairs whose two sites share one interval: (0, 0.01], (0.01, 0.05], (0.05, 0.1], (0.1, 0.2], (0.2, 0.5]. Writes `stats/<mod>/plots/{A,B,D}.variant.ld_decay.maf_bin.png` plus pooled A+B+D `all.variant.ld_decay.maf_bin.png` (does not overwrite `*.rare_common.png`). `--ld_maf_bin_plot_only true` replots from published `*.maf_bin.info.tsv`. |
| `circos_ld_map.nf` | Three-ring Circos: inner rare and middle common density from `--mod` (default `test_wap_thin_refchr`); outer LD blocks from `--circos_ld_mod` (e.g. `test_wap_common_refchr`; defaults to `--mod`). Blocks are adjacent-r2 ≥ 0.8 runs, then merged when junction Q1(\|D′\|) > 0.98. Uses PLINK 1.9 `--r2 dprime` on the LD-mod A/B/D bfiles (not the 40+ GB vcor). Plot step scans each `*.info.ld_block.ld` in chunks and processes A/B/D in parallel (`task.cpus`, label `cpus_16`). The outer ring paints `circos_bin_bp` (1 Mb) windows whose unioned block coverage is ≥ `circos_ld_block_coverage_frac` (default 0.5) and merges adjacent passing windows. Writes `stats/<mod>/plots/all.variant.circos_ld_map.png` plus `info/all.variant.circos_ld_map.{density,blocks,blocks_1mb}.info.tsv`. `--circos_ld_map_plot_only true` replots from the published density/blocks tables (no PLINK, no `.ld` scan). |
| `vep_atlas_plots.nf` | Unique-site VEP atlas from published `{A,B,D}.vep.tsv`. Genome-wide `--mod` (default `test_wap_refchr`) supplies common (MAF≥5%) and rare (MAC≥2 and MAF<1%) window counts. Rare predicted LoF uses `--vep_atlas_exon_mod` (`test_wap_exon`; GFF exon including UTR) plus that mod's afreq/acount/pvar. Collapses multi-transcript rows to the most severe consequence. Writes `stats/<mod>/plots/all.variant.vep_atlas.{composition,freq_function,circos_windows,mirror,enrichment,windows.violin}.png`. Counts are per reference window, not per callable Mb. |
| `vep_atlas_circos_replot.nf` | Replot `all.variant.vep_atlas.circos_windows.png` from published `windows.tsv`. Thin outer ideogram (no cytoband legend). Next: LD blocks, gene counts (yellow), predicted LoF (purple), common (red), rare (blue). Centre legend; chromosome names sit just outside the Mb ticks. |
| `filter_cascade_popdep_density.nf` | Stream all n500 (`_popdep_n500_gam_credible_ntaxa500`) RelativeDepth mean-vs-SD cells (250×150 on [0,5]×[0,5]). Trust the densest 60% and 80% of non-empty cells. Keep `test_thin` variants in the 80% cells (re-annotate in the work dir; does not overwrite `process/test_thin` popdep TSVs), then wheat excess-het (`heatmap_mac_h.filter` / `P_excess>=1e-6`), MAC≥5, WAP MAC≥2. Writes step counts and a genome estimate scaled by `thin_rate`. |
| `popdep_cv_reliable_lib.nf` | Stream all n500 (`_popdep_n500_gam_credible_ntaxa500`) 1 bp RelMean vs RelCV cells. Linear **[0, 6] × [0, 25]**, **500×500 = 250_000** cells. Keep densest non-empty cells until **60% and 80% of in-range sites** remain. Writes two merged BED libraries (`process/_popdep_n500_cv_reliable/keep{60,80}/` plus genome `all.site.popdep.cv_reliable.keep{60,80}.bed.gz`; Chrom is the PopDep chr_id). Heatmap and site-coverage curve go to `stats/_popdep_n500_cv_reliable/`. |
| `filter_cascade_cv_reliable.nf` | Intersect `test_thin` variants with those keep-80 and keep-60 BEDs, each as its own cascade: PopDep → excess het (`P>=1e-6`) → MAC≥5 → WAP MAC≥2. Writes per-step counts and a full-library estimate (`n / thin_rate`) to `stats/test_thin/info/all.variant.cv_reliable.keep{60,80}.filter_steps_*.tsv`. |
| `wap_cv_reliable_vep_atlas.nf` | PopDep keep 80% and keep 60% on `--mod` (default `test_wap_refchr`) and `--vep_atlas_exon_mod` (`test_wap_exon`), then excess het from that mod's own `plink2 --het`, then MAC≥`--filter_cascade_min_mac` (default 5). Arm ID in the pvar, not linear `#CHROM`/`POS`. Does not rewrite pfiles. Step counts, MAF histograms, and the kept-site table go to `stats/<mod>/`. The VEP atlas is redrawn from those IDs as `stats/<mod>/plots/all.variant.vep_atlas.cv_reliable.keep{60,80}.*` (existing `all.variant.vep_atlas.*` stay). |
| `wap_mac1_testsets.nf` | Rebuild WAP pfiles at `--test_wap_min_mac 1` into new mods. `--wap_mac1_stage arm` writes `process/test_wap_mac1` (`--keep` WAP, `--mac 1`, `--thin` `test_wap_thin_rate`, seed `ld_decay_random_seed`). `refchr` recodes that merge onto linear chroms (`test_wap_mac1_refchr`) plus basic info and `*.assoc.acount`, no LD. `exon` extracts GFF exons from `main_raw` with `--mac 1` and no thin (`test_wap_mac1_exon`). Does not overwrite `test_wap`, `test_wap_refchr`, or `test_wap_exon`. |
| `vep_atlas_rare_common_share.nf` | From published `all.variant.vep_atlas.cv_reliable.keep{60,80}.counts.tsv`, one 100% stacked bar per subgenome (common under rare). Denominator is `n_common + n_rare` on that subgenome. Writes `*.rare_common_share.png`. |
| `test_popdep_density_keep_replot.nf` | Colour published `test_popdep` scatter (~2M) by those genome 60%/80% cells. Kept blue, dropped grey. Depth and RelativeDepth mean-vs-SD JointGrids, linear and log, stacked x/y marginals. Writes `stats/test_popdep/plots/test_popdep.popdep.density.keep{60,80}.scatter_*.png`. |
| `test_popdep_log_density.nf` | log10 RelativeDepth mean-vs-SD histogram on all `test_popdep` sites (~100M). 90_000 equal-log cells (300×300); axis range is the observed min/max (no `[0,5]` clip). Keep densest non-empty cells until **60%/80%/90% of sites** remain; write the implied cell-top fraction. Also writes retained-site % vs nonempty-cell top % (linear and log-y). Heatmap counts every site; keep/drop JointGrids use the published scatter (point alpha 0.2). Linear RelDepth keep 60/80/90 scatters use axes [0, 10]; log JointGrids keep the observed log10 range. `--popdep_log_counts_npz` reuses a previous histogram. `--popdep_log_plot_keep_drop false` skips JointGrids. Does not scan n500 grids. |
| `test_popdep_linear_density.nf` | Linear RelativeDepth mean-vs-SD histogram on all `test_popdep` sites (~100M). Range **[0.02, 6]×[0.02, 6]** (no log). 90_000 equal-linear cells (300×300). Keep densest non-empty cells until **60/80/90/95/97/98/99% of in-range sites** remain. Heatmap and keep/drop JointGrids (axes 0.02–6, point alpha 0.2). `--popdep_lin_counts_npz` reuses a previous histogram. `--popdep_lin_plot_keep_drop false` skips JointGrids. Does not overwrite `log_density` artefacts. |
| `test_popdep_cv_density.nf` | Linear RelativeDepth mean-vs-CV histogram on all `test_popdep` sites (~100M). RelMean **[0, 6]**; CV is the observed min/max among in-range positive-mean sites (no CV clip). 90_000 equal-linear cells (300×300). Keep densest non-empty cells until **60/70/80/85/90/95/97/98/99% of in-range sites** remain. Heatmap and keep/drop JointGrids (mean 0–6, CV free, point alpha 0.2) plus RelDepth `cv_vs_mean_by_class` (mean 0–6). `--popdep_cv_counts_npz` reuses a previous histogram. `--popdep_cv_plot_keep_drop false` skips keep/drop JointGrids. Does not overwrite `log_density` / `lin_density` artefacts. |
| `species_ld_plot.nf` | Comparative LD decay distance vs genome size from `00data/10meta_data/species_LD.csv`. Canvas 2:1; log genome-size y-axis. Writes log-x and linear-x pairs (`species_ld_vs_genome{,.labeled}.png` and `species_ld_vs_genome.linear_x{,.labeled}.png`), including *Homo sapiens*. |
| `subgenome_maf_tajima.nf` | Stacked A/B/D MAF histograms (linear y and log y) from `process/<mod>/variant/{A,B,D}.info.afreq`, plus stacked Tajima's D histograms. Tajima's D: PLINK2 `--export vcf-4.2 bgz` (SNPs only, work dir only) then `vcftools --TajimaD` (`tajima_d_window_bp`, default 100 kb). Publishes `stats/<mod>/plots/all.variant.{maf,tajima_d}.subgenome.dist*.png` and per-subgenome `*.Tajima.D`. `--tajima_d_plot_only true` replots from published Tajima tables. |
| `coding_mask_chr_pilot.nf` | Classify a VEP-like TSV into project masks (`plof_high`, `mis_severe`, `mis_consensus`, `protein_altering`). Does not run VEP. `--coding_vep_table` required. |
| `dual_af_from_afreq.nf` | Join existing WGS-core and WAP `*.afreq` tables on variant ID; rare/discordant flags. `--wgs_afreq` and `--wap_afreq`. |
| `singer_arg_common_inter.nf` | WAP China subset (~1300; `r1_id` plus `PxxxxA` aliases) of `test_common_inter_refchr`: 64 spread 2 Mb windows (A/B/D, skip 0–2 Mb and centromeres), rewrite GT `/` to `\|`, Beagle-impute, SINGER ARG per window. Publishes under `process/test_common_inter_2mb_wap_china_refchr/` (`WAP_China.keep`, `singer_windows.tsv`, `export/`, `singer/`). Does not overwrite `test_common_inter` or `test_common_inter_refchr`. `--singer_from_beagle true` reruns only windows missing `singer_*_nodes_0.txt`, from published Beagle VCFs, with our random `-seed` values (not `singer_master -debug`, not 42). |
| `singer_arg_preview.nf` | Tables and plots from published SINGER texts (no new ARG sampling). Writes `stats/<mod>/info/all.variant.singer_arg.window.info.tsv` plus status / Ne / TMRCA / recombination figures. Default `--mod test_common_inter_2mb_wap_china_refchr`. |
| `singer_arg_relatedness.nf` | Local relatedness ranking from published SINGER ARGs: mean of four haplotype TMCRAs per sample pair, averaged over successful windows and a grid of positions (default 11). Writes per-sample top-K neighbors and the closest pairs (`*.neighbor.info.tsv`, `*.pair_rank.info.tsv`). Not a pedigree. |
| `sample_identity_hom_snp_ibs.nf` | Extra IBS on merged A/B/D bfiles: biallelic SNPs only, score a site only when both samples are homozygous. Publishes `process/{mod}/ibs_hom_snp/{A,B,D}.ibs_hom_snp.mibs` beside the existing PLINK `process/{mod}/ibs/*.ibs.mibs` (does not overwrite). Pair/neighbor compare tables and PLINK-vs-hom scatter go to `stats/{mod}/`. `--mods` default `test_thin,test_common_inter`. |
| `sample_identity_hom_eq1.nf` | From existing hom-SNP and PLINK square IBS matrices, extract pairs with hom-SNP IBS = 1 and plot the PLINK IBS of those pairs (`*.eq1.pairs.tsv`, `*.eq1.plink.dist.png`, `*.eq1.slope.png`). |
| `sample_identity_ibs_pair_dist.nf` | Full-range PLINK and hom-SNP IBS histograms from the published 200k proportional pair sample (`*.hom_snp.vs_plink.pairs.tsv`). Writes `*.ibs_pair.{plink,hom}.dist.png` plus log-count versions. Does not reuse the hom-SNP `*.dist.png` that only shows pairs above 0.99. |
| `sample_identity_ploidy_ibs_dist.nf` | Within-ploidy PLINK and hom-SNP IBS pair histograms. Previous genome formula from `sample_group.txt`; breachr full-cohort `breachr_segments_all.tsv` drops a claimed subgenome when all 7 homoeologues are `subgenome_absent` or `subgenome_reference_mismatch` (hexaploid-labeled D mismatch → AABB; AB B mismatch → AA). Expected nulls (AB D-null, diploid off-subgenomes) are not corrections. Subgenome matrices: A → AA/AABB/AABBDD; B → SS/AABB/AABBDD; D → DD/AABBDD. Faceted histograms (linear + log count) and overlay density. `--mods` default `test_thin,test_common_inter`. |
| `sample_identity_ploidy_group_ibs_dist.nf` | Within-class IBS using each ploidy's own subgenomes together: AA on A, SS on B, DD on D; AABB = mean of A+B; AABBDD = mean of A+B+D. Two IBS types (PLINK and hom-SNP). Breachr-corrected labels as in `sample_identity_ploidy_ibs_dist.nf`. Publishes `all.sample.identity.ploidy_group_ibs.*` (does not overwrite per-subgenome `*.ploidy_ibs.*`). `--mods` default `test_thin,test_common_inter`. |

**Partial task notes** (workflow bodies in `subworkflows/local/partial/`):

| `--partial_task` | Deliverables / behaviour |
| --- | --- |
| `ld_redraw` | Re-run LD decay + cross-chr plots from existing `process/<mod>/variant/*.vcor` (no full processor). |
| `ld_rare_common` | Re-run three-colour LD decay / stacked r2 dist / scatter from existing `*.info.ld.vcor` + `*.info.afreq` on mixed-MAF mods (e.g. `test_thin`, `test_wap_refchr`, `test_wap_thin_refchr`). Rare = MAF < `hf_maf` (0.05). Class plots include rare-rare, rare-common, and common-common. Sampled MAF-vs-MAF scatter uses **pairwise-complete MAF** from `{mod}/{id}_test.plink2`. Pair table: `*.maf_scatter.info.tsv` (adds signed `D`, dosage-based `D_PRIME`, `DISTANCE_KB`). `|D'|` uses the sign-aware Lewontin denominator from dosages, not \(r^2/r_{\max}^2\). Decay / scatter / density: `*.dprime.png` + `*.dprime.scatter.png` + `*.dprime.dist.png` + `*.dprime.dist.log.png` (`log_{10}(density+1)`) + signed-D' `*.dprime.signed_prop.png` / `*.dprime.signed_boundary_stack.png` / `*.dprime.signed_interior_dist.png` + `*.dprime.info.tsv`. r2 decay bins remain `*.info.tsv`. Prefix `*.variant.ld_decay.rare_common.*`. |
| `assess_plink` | `plink2_assess_debug_slice` → `process/<mod>/`; `assess_plink_debug_plots` → `stats/<mod>/` via `assess_slice.py`. |
| `assess_vcf` | Legacy VCF assess on `process/<mod>/export/*.debug.vcf.gz`; prefer PLINK2 when pfiles exist. |
| `chr_counts` | `process/<mod>/info/<mod>.chr_variant_counts.tsv` + `.by_ref.tsv`. |
| `chr_compare` | Genome-wide density panels + thin/common fraction bars under `stats/thin_common_compare/`. |
| `rebuild_lib_stats` | `TEST_PLINK_STATS` on `process/test_rebulld_lib/` chr002 tables; LD plots skipped. |
| `mac_stats` | MAC bar (0–100) + het-fraction heatmap; `*.variant.mac.mac0.info.tsv` audit for MAC=0. |
| `mac_het_excess` | Wheat selfing-model excess-het sweep: PLINK2 `--het` median F on `het_source_mod` (default `test_common_only`), then `P(H>=H_obs\|MAC,F)` retain/drop counts on gcount under `mac_het_variant_mod` or `params.mod` (e.g. `test_rare_only`). Does **not** write a filtered pfile. Heatmaps: MAC vs H and MAC vs H/MAC for unfiltered / retained / dropped sites at MAC 0-100 / 0-200 / 0-500 / 0-1000, plus `heatmap_mac_h.filter.0_{100,200,500}.png` (white-to-red PowerNorm scale, excess-het H cutoff, MAC=5 line). Popdep JointGrid scatters (`*.variant.popdep.het_filter.scatter_{Depth,RelativeDepth}_{SD,CV}_vs_*_Mean.png` plus `*.log.png` log-log axes): mean vs SD / CV and relative mean vs relative SD / CV, coloured by pass/fail at `mac_het_retain_alpha`, with 100% fill x/y marginals (within-bin pass/fail count fraction; fail on top; 250×150 bins) (needs filled `*.popdep.info.tsv` under `process/{mod}/variant`; MAPQ30 grids: annotate with `popdep_annotate.nf --mod test_rare_only_mapq30 --popdep_dir …/_popdep_n500_gam_credible_mapq30` and `--mac_het_variant_mod test_rare_only` so plots do not overwrite MAPQ0 artefacts). `--mac_het_rh_diag` classifies `R_H=H/MAC` into bottom / middle / diagonal bands, re-estimates site-moment F after dropping `H=MAC` (`MAC>=5`) sites, and joins unfiltered plus optional `--mac_het_popdep_mapq30_mod test_rare_only_mapq30` depth, site MQ, and missingness. `--mac_het_rh_sample` (default on with rh_diag) runs PLINK2 `--het` on the all-het extract to test whether the diagonal is concentrated in few samples. Flags: `mac_het_sweep`, `mac_het_heatmap`, `mac_het_popdep_scatter`, `mac_het_rh_diag`, `mac_het_rh_sample`. |
| `mac_dist_redraw` | Re-draw MAC distribution log plots from existing gcount tables. |
| `mq_missing_reg` | MQ vs missing regression plots from `*.mq.info.tsv` + vmiss. |
| `popdep_missing_reg` | Depth vs missing regression plots from `*.popdep.info.tsv` + vmiss. |
| `abstract_mq_50_bams` | Full-chr bcftools mpileup site MQ → `{mq_dir}/reference/`. |
| `main_raw_popdepth` | Full-chr TIGER popdepth → `{popdep_dir}/variant/chrNNN.popdep.txt.bgz` (+ tabix `.tbi`). |
| `aneuploidy_mapq` | Per-sample `samtools view -F 1796` scan → `{output_dir}/{job}/stats/aneuploidy_check/info/mapq_check/{sample}.mapq_depth_bin.tsv` (1 Mb bins: MAPQ histogram + `bin_mean_depth`). Skips samples whose TSV already exists when `mapq_skip_existing=true`. Optional `-c conf/mosdepth_aneuploidy_mapq.config` sets `mapq_max_forks=16`, `mapq_plot=false`. Sample list default: `{home_dir}/00data/05taxaBamMap/final.ALL.taxaBamMap.txt`. |
| `aneuploidy_breachr` | Reads `info/flagged_samples.tsv` (`aneuploid_flag=true`), runs `breachr` on existing `info/mapq_check/*.mapq_depth_bin.tsv`, publishes `info/breachr_check/` (`breachr_segments_all.tsv`, `cnv_misclassified_as_aneuploidy.tsv`, `breachr_reconcile_summary.tsv`). Optional `-c conf/breachr_flagged.config`. |
| `wheat_from_plink` | `RUN_WHEAT_FROM_PLINK` on existing merged test pfiles (`wheat_integrated_mod` required). |
| `sample_manifest` | Join `sample_group.txt` to taxaBamMap coverage; write `assoc.keep` (WAP) and `wgs_core.keep` (A,AB,ABD,WAP,HZNU,Nature); QC tables/plots; RefSeq v1.0 FASTA vs missing GFF/triad inventory. |
| `cohort_freq` | PLINK2 `--freq` / `--freq counts` / `--missing` per (subgenome, cohort) → `process/{mod}/variant/{sg}.{cohort}.{afreq,acount,vmiss}`. Cohorts from `cohort_freq_cohorts`; keep files from `cohort_keep_dir` (default the `sample_manifest` mod). |
| `annotate_vep` | Sort/bgzip/tabix the IWGSC v1.1 GFF, then `vep --gff --fasta --tab` → `process/{mod}/variant/{id}.vep.tsv`. With `--vep_from_mod` it first builds one sites VCF per subgenome from the pvar **ID** column, which carries `<segment>-<segment_pos>` coordinates matching the GFF and FASTA even when `#CHROM`/`POS` are linear (`*_refchr` mods). Segments missing from the GFF are dropped, and segment 0 is renamed (`vep_unplaced_name`) because VEP reads the chromosome name `"0"` as Perl false and aborts. Narrow a rerun with `--vep_subgenomes`; skip the GFF sort with `--vep_prepared_gff`. |
| `annotate_snpeff` | Annotate the same sites VCFs as VEP. Default genome is the prebuilt 4.3 Lulab DB at `{home_dir}/00data/03ref/07snpeff` (numeric chroms 0–42, longest transcript; conda `anno` is 5.4c and cannot load that bin). Builds from GFF+FASTA only when that tree is missing. Reuses `process/{mod}/annotate/*.sites.vcf.gz` when present. FastCall ALT `I`/`D` are recoded for SnpEff 4.3 (`REF+N` / `<DEL>`) and restored in the published VCF. Writes `{id}.snpeff.vcf.gz` and a VEP-shaped `{id}.snpeff.tsv`. Burden still uses the VEP coding_mask. |
| `site_anno` | Stream one `00data/11epi_anno/chr*_featHighRes.bed.gz` per pfile chromosome and keep the rows a variant sits on → `stats/{mod}/info/{sg}.site_anno.tsv` with degeneracy, sequence element, conserved and regulatory classes. Those BEDs use RefSeq v1.0 **linear** coordinates, the same as `*_refchr` `#CHROM`/`POS`, so this join needs no lift-over — unlike `annotate_vep`, which must parse segment coordinates out of the variant ID. They cover chr1A–chr7D only, so `Others` is skipped. The `degeneration` code is **2 = 0-fold, 1 = 4-fold, 0 = 2/3-fold** (verified against reference codons, see `site_anno.py`); 4-fold sites are the neutral reference. |
| `gene_expr` | Normalize a TPM atlas onto v1.1 gene IDs → `stats/{mod}/info/*.gene_expr.tsv` with peak/mean TPM, per-group means, expression breadth and tau. v1.0 IDs map through the `previous_id` attribute on the GFF gene rows; LC genes are dropped rather than suffix-stripped, because HC and LC are numbered independently in v1.0. |
| `variant_sets` | Join coding masks to WGS-core/WAP frequencies → `stats/{mod}/info/*.variant_sets.tsv` plus one `process/{mod}/variant/*.{set}.extract` list per set for `plink2 --extract`. `--variant_sets_from_mod` derives coding masks and dual AF per subgenome from that mod's `{id}.vep.tsv` and `{id}.{cohort}.afreq` files, and picks up `{id}.site_anno.tsv` from the same mod's stats tree when present. With site annotation it adds `neutral_4fold`, `conserved_noncoding_rare_wap`, `regulatory_rare_wap` and `plof_conserved_rare_wap`; with `--variant_sets_gene_expr` it adds `plof_expressed_rare_wap`. On a WAP-only mod the two cohorts resolve to the same samples, so the dual-AF axis is degenerate — use an all-sample mod for a real WGS-core frequency. |
| `assoc_wap` | `plink2 --glm` on WAP only (`assoc.keep`), optionally restricted by a `variant_sets` extract list → `process/{mod}/assoc/`, Manhattan/QQ under `stats/{mod}/plots`. Structure PCs are computed once (`--assoc_wap_pca_source`, default `A`) and shared by every subgenome, since structure is a property of the individual; `--assoc_wap_covar` overrides them. `--assoc_wap_burden --assoc_wap_burden_set <set>` switches to the gene burden track (see below). |
| `phenotype_wap` | Align a wide phenotype table to `assoc.keep` → `process/{mod}/phenotype/*.pheno.tsv` (plink2 `--pheno`) plus per-trait coverage QC. No BLUE/BLUP; adjust replicates upstream. |

Historic runs in **`doc/progress/`** may still reference removed `workflow/Genetics/tmp/*.nf` paths; use `partial_router.nf` for new runs.

---

## Wheat integrated workflows (`modules/local/integrated/`)

Wheat integrated modes are run from **`main.nf`** when **`params.mod`** starts with **`wheat_`** (VCF `check_input` is skipped). PLINK2 compute lives in **`processor/processor_plink2.nf`**; plots/reports in **`stats/stats_integrated.nf`**; study wiring in **`subworkflows/local/wheat/wheat_integrated_study.nf`**.

Python matches genotype stats style: inline `#!/usr/bin/env python` imports from **`genetics.genomics.*`** and **`genetics.gwas.association_plot`** (no `python -m` CLI). **`genetics.wheat`** is a deprecated shim only — do not add new imports there.

### Routing summary

When `params.mod` starts with `wheat_`, the workflow runs `modules/local/integrated/integrated_wheat.nf`.

### Active `params.mod` values (wheat)

| `params.mod` | Input params | Python entry |
| --- | --- | --- |
| `wheat_snp_qc` | `wheat_table_input` | `genetics.wheat.snp_qc.run_snp_qc` |
| `wheat_pca_tsne` | `wheat_table_input` | `genetics.wheat.population_structure.run_population_structure` |
| `wheat_tagsnp` | `wheat_table_input` | `genetics.wheat.tagsnp.run_tagsnp_selection` |
| `wheat_hapmap` | `wheat_table_input` | `genetics.wheat.hapmap.run_hapmap_build` |
| `wheat_cnv` | `wheat_table_input` | `genetics.wheat.cnv.run_cnv_calling` |
| `wheat_genetic_map` | `wheat_table_input` | `genetics.wheat.genetic_map.run_genetic_map` |
| `wheat_gwas` | `wheat_gwas_genotype`, `wheat_gwas_phenotype` | `genetics.wheat.wheat_gwas.run_gwas` |
| `wheat_kgwas` | `wheat_kgwas_kmer_matrix`, `wheat_kgwas_phenotype` | `genetics.wheat.kgwas.run_kgwas` |

Tuning knobs (defaults in `conf/base.config`): `wheat_snp_qc_*`, `wheat_pca_*`, `wheat_tsne_*` (t-SNE on PLINK PCs in `wheat_pca_tsne`), `wheat_tagsnp_*`, `wheat_hapmap_window_size`, `wheat_cnv_*`, `wheat_gwas_trait`, `wheat_kgwas_trait`.

### Required launch params (wheat)

- `output_dir`, `job`, `user_dir` (Conda `stats` env, same pattern as `modules/local/genotype/stats/stats.nf`)
- `mod` as one of the table rows above
- Paths per mod (e.g. `--wheat_table_input` for single-table tasks)

`home_dir` / `src_dir` are not used by wheat processes but may still be set for a uniform CLI with other mods.

### Wheat output layout

Artifacts publish under **`{output_dir}/{job}/integrated/{mod}/info`** (TSV) and **`.../plots`** (PNG when figures exist). For per–subgenome **`wheat_pca_tsne`** from PLINK2, expect **`*.pca.tsv`**, **`*.eigen.tsv`**, **`*.tsne.tsv`**, **`*.pca.png`**, **`*.variance.png`**, **`*.tsne.png`** (prefix = subgenome id). Table-only entry modes use **`{mod}`** as the prefix (e.g. `wheat_gwas.gwas.tsv`).

Topic folders under `modules/local/integrated/` (`gwas/`, `snp_qc/`, …) are organizational only; the executable surface is **`integrated_wheat.nf`** plus **`main.nf`** routing.

### Example (`wheat_snp_qc`)

```text
nextflow run /data/home/tusr1/git/script/workflow/Genetics/main.nf \
  -c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config \
  --user_dir /data/home/tusr1 \
  --output_dir /path/to/out \
  --job demo_wheat \
  --mod wheat_snp_qc \
  --wheat_table_input /path/to/snp_summary.tsv
```

Log production runs in today's **`doc/progress/YYYY-MM-DD.md`** per **`progress-logging.mdc`**.

---

## Where and how to run (workstation policy)

**Do not run from the git repo.** Use the vmap4 project workspace with a **two-level** folder layout:

| Level | Path pattern | Example |
| --- | --- | --- |
| Project root | `/data/home/tusr1/01projects/vmap4/` | All runtime work lives here |
| Task module | `<NN><area>.<topic>/` | `12stats.genome/` — current stats (10stats.genome ended at 99) |
| Run attempt | `<NN>run` or `<NN>run_<slug>/` | `57run_mac_stats_test_thin`, `23run_assess_plink2_debug_stub` |

**Before each new run:**

1. `cd /data/home/tusr1/01projects/vmap4/<module>/`
2. List existing `*run*` folders; set `NN` to max + 1 (zero-padded: `01`, `02`, … `63`, …).
3. `mkdir -p <NN>run_<descriptive_slug>` and `cd` into it.
4. Launch with absolute repo paths, e.g.:
   - Full router: `…/workflow/Genetics/main.nf`
   - Partial: `…/subworkflows/local/entry/partial_router.nf --partial_task <name>`
   - Always: `-c /data/home/tusr1/git/script/workflow/Genetics/nextflow.config`
5. Append **cwd + full command** to today's **`doc/progress/YYYY-MM-DD.md`** (`#### NF replay`).

**Conventions**

- **Do not reuse** an old run folder for a fresh attempt — create the next `NNrun_…`. When a module's highest `NNrun_*` is **99**, open successor **`NN+2<name>/`** (e.g. `10stats.genome` → `12stats.genome`); do not create `100run_*`.
- **Conda `run`** for Nextflow and **`screen`** (use `screen` from env `run` via `conda activate run`; not system `/usr/bin/screen`); **`stats`** for Python stats processes inside NF.
- **`screen`** for long production runs (`conda activate run` before `screen`; re-activate inside the session); `-preview` / small debug exempt.
- **Read-only inputs:** `/data1/dazheng_tusr1/vmap4.VCF.v1` (prepared `test_plink/process` paths in **`doc/KNOWLEDGE_README.md`** / **`doc/project_knowledge/`**).
- **Ephemeral preview dirs** (`00nf_preview*`) under `vmap4/` are OK for syntax smoke tests.

**Partial-task → typical slug examples**

| `--partial_task` | Example run folder under `08stats.genome/` |
| --- | --- |
| `assess_plink` | `24run_assess_plink2_full_test_thin` |
| `ld_redraw` | `16run_ld_plots_redraw_test_thin` |
| `ld_rare_common` | `80run_ld_rare_common_test_thin` under `10stats.genome/` |
| `mac_stats` | `57run_mac_stats_test_thin` |
| `mac_het_excess` | `01run_popdep_het_filter_fill_test_rare_only` under `12stats.genome/` (earlier `79run_*` under `10stats.genome/`) |
| `chr_compare` | `51run_chr_variant_compare_plots` |
| `wheat_from_plink` | `43run_wheat_pca_grp_test_thin` |
| `sample_manifest` | `57run_sample_manifest` under `12stats.genome/` |
| `cohort_freq` / `annotate_vep` / `variant_sets` / `assoc_wap` / `phenotype_wap` | `59run_wap_ph_gwas_chain`, `60run_wap_ph_assoc_final`, `61run_wap_annotate_vep`, `62run_wap_variant_sets` under `12stats.genome/` (first full plant-height run on `test_wap_refchr`) |
| `site_anno` / `gene_expr` | `63run_wap_site_anno`, `64run_wap_gene_expr` under `12stats.genome/` |
| `aneuploidy_mapq` | `15run_mapq_full_cohort` under `11aneuploidy/` |
| `aneuploidy_breachr` | `16run_breachr_flagged_test` under `11aneuploidy/` |

Authoritative guardrails: **`.cursor/rules/workstation-core.mdc`** and **`.cursor/rules/workstation-nextflow.mdc`**.

---

## Downstream tools

Processes invoke **`src/python`**, **`src/r`**, and **`src/java`** as configured in each `process` block (typically conda `stats` under `params.user_dir` for Python stats). Keep processor vs stats boundaries as in **`workstation-nextflow.mdc`** (heavy `plink2` in `modules/local/genotype/processor/`; plotting / summarization under `modules/local/genotype/stats/`).

---

## Further reading

- **`doc/progress/`** — daily engineering narrative and **`#### NF replay`** command blocks.
- **`doc/PROGRESS_README.md`** — index to recent daily files.
- **`doc/KNOWLEDGE_README.md`** — structured project registry (paths, runs, layout).
- **`.cursor/rules/workstation-*.mdc`** — workstation and pipeline policy.
