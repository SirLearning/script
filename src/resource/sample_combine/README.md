# sample_combine

Canonical snapshots for germplasm/sample merging. From 2026-10-27 the
`curation_germplasm_combine_name_groups*.tsv` tables are the merge basis.
`sample.dedup.checked.info.tsv` is kept as a reference only.

## Files

| File | Source | md5 | mtime | Logical rows |
|---|---|---|---|---|
| `curation_germplasm_combine_name_groups.tsv` | `/data/home/tusr1/git/DBone/src/main/resources/curation_dict/` | `84a02584842e17c88a25b46456f5cded` | 2026-10-27 20:44 | 3665 |
| `curation_germplasm_combine_name_groups_ipk_hzau_watkins.tsv` | `/data/home/tusr1/git/DBone/src/main/resources/curation_dict/` | `a25d55391b7ed3bb3b0ffce7bd0650b2` | 2026-10-28 | 2108 (adds `watkins_samc`) |
| `sample.dedup.checked.info.tsv` | `/data/home/tusr1/git/script/src/resource/` | `b3c7264f6f345746490e52274b8254d5` | 2026-08-26 | 7590 samples |
| `combine_groups_not_merged_hom_ibs.tsv` | derived: `test_plink/stats/test_common_inter/info/all.sample.combine_ibs.pairs.tsv`, `hom_ibs > 0.995` | `b9c33ac5a9075e1c8254ebee5bcd2909` | 2026-10-28 11:06 | 228 groups |
| `combine_pairs_decisions_hom_ibs.tsv` | every WAP_PQ/VMap3 same-ploidy pair with a pair-level status (pass / fail / unjudged_low_coverage) | `4f2ba271d302cfce457b0e4114377b80` | 2026-10-28 11:06 | 878 pairs |
| `combine_groups_support_hom_ibs.tsv` | group summary of those statuses (all_pairs_supported / supported_except_low_coverage / contains_fail) | `52d4f9ea59d5f9fb774074133a3d766b` | 2026-10-28 11:06 | 632 groups |
| `combine_pairs_decisions_ipk_hzau_watkins_hom_ibs.tsv` | same, IPK/HZAU/Watkins cohort | `bb2da84fb338cd31a546cea8e8aa4d33` | 2026-10-28 11:06 | 4783 pairs |
| `combine_groups_support_ipk_hzau_watkins_hom_ibs.tsv` | same, IPK/HZAU/Watkins cohort | `166a255fdf2d37034f5eb53b5b0d681f` | 2026-10-28 11:06 | 1053 groups |
| `combine_pairs_not_merged_hom_ibs.tsv` | same pairs table, `hom_ibs > 0.99` | `2a358d9d192d5f1e3d63277f0e9db269` | 2026-10-28 11:06 | 277 pairs |
| `combine_groups_not_merged_ipk_hzau_watkins_hom_ibs.tsv` | derived: `…/all.sample.combine_ibs_ipk_hzau_watkins.pairs.tsv`, `hom_ibs > 0.99` | `a050613dfe1c79c28cf211069ba89388` | 2026-10-28 11:06 | 5 groups |
| `combine_pairs_not_merged_ipk_hzau_watkins_hom_ibs.tsv` | same IPK/HZAU/Watkins pairs table, `hom_ibs > 0.99` (carries one `keep_merge` decision) | `ae64f0aa396dda510327804887bcdef6` | 2026-10-28 11:06 | 5 pairs |
| `combine_pairs_not_merged_ipk_hzau_watkins.tsv` | same, plink_ibs judged — diagnostic only | `c39bb7f63917e6c4b027397d4ca9c53d` | 2026-10-27 16:44 | 214 pairs |
| `combine_map.tsv` | derived: pairs tables + curation tables (see below) | `70c991ff1a5538ef666474446b6ba129` | 2026-10-28 | 6196 germplasm rows |
| `final.ALL.taxaMergeMap.txt` | derived: `combine_map.tsv` + `final.ALL.taxaBamMap.txt` | `906772bece4f45f6f6c552a9a9cfa99d` | 2026-10-28 | 6196 lines = 8288 taxa |
| `final.A.taxaMergeMap.txt` | same, `final.A.taxaBamMap.txt` | `9e0365e8efdc5910d69906a62ed21576` | 2026-10-28 | 5588 lines = 7677 taxa |
| `final.B.taxaMergeMap.txt` | same, `final.B.taxaBamMap.txt` | `59636234d7c6f5e228fa94be6f530e8c` | 2026-10-28 | 5507 lines = 7596 taxa |
| `final.D.taxaMergeMap.txt` | same, `final.D.taxaBamMap.txt` | `b70b3a41f10f0c4a4dfcf133c775f601` | 2026-10-28 | 5656 lines = 7741 taxa |
| `fastcall3_scan3_merge_map.ALL.tsv` | multi-member subset of `final.ALL.taxaMergeMap.txt` | `ece2314166e9626928d7623cac8e97dc` | 2026-10-28 | 1361 merge lines |
| `fastcall3_scan3_merge_map.A.tsv` | same, A | `a60d25537ffc1e4c5ac2aa42afcfb190` | 2026-10-28 | 1358 merge lines |
| `fastcall3_scan3_merge_map.B.tsv` | same, B | `a60d25537ffc1e4c5ac2aa42afcfb190` | 2026-10-28 | 1358 merge lines |
| `fastcall3_scan3_merge_map.D.tsv` | same, D | `2d3426c3ea53d0660144fa93f683c612` | 2026-10-28 | 1354 merge lines |

The two combine tables are frozen snapshots. Do not edit them in place; if the
DBone sources change, copy new snapshots here and record their md5/mtime.

## Columns

`curation_germplasm_combine_name_groups.tsv`:

- `germplasm_id`, `germplasm_name`, `germplasm_name_cn`,
  `germplasm_accession_number`, `accession_numbers`, `synonyms`
- `wap_pq_combine_names` / `wap_pq_combine_count`
- `vmap3_combine_names` / `vmap3_combine_count`
- `total_combine_count`, `multiple_combine_names`

`curation_germplasm_combine_name_groups_ipk_hzau_watkins.tsv`:

- `germplasm_id`, `germplasm_name`, `germplasm_name_cn`,
  `germplasm_accession_number`, `accession_numbers`, `synonyms`
- `ipk_combine_names` / `ipk_combine_count`
- `hzau_combine_names` / `hzau_combine_count`
- `watkins_combine_names` / `watkins_combine_count`
- `total_combine_count`, `multiple_combine_names`

Within a name column, names are pipe-separated (`P0001|P0002`). A row is one
germplasm/merge group; the names in its cohort columns are the candidate merged
units.

## Merge policy

- Use the combine-group rows as the merge basis.
- Exclude names whose `status` is `excluded` in the DBone control tables:
  - `curation_germplasm_combine_sample_untrusted_excluded.tsv`
  - `curation_germplasm_combine_sample_conflict_ck_excluded.tsv`
- The dedup table is reference only, not the merge basis.
- No merged IDs are assigned here. Groups are reviewed one by one; the merged
  ID scheme is designed separately.

## Coverage annotation

The `*.with_coverage.tsv` files keep every original column and append:

- `coverage_total` — sum of `Coverage-Of-All-Bams` over the unique genotype
  names in the row (P#### falls back to P####A when needed)
- `coverage_n_names`, `coverage_n_mapped`, `coverage_n_missing`
- `coverage_missing_names`, `coverage_mapped_names`
- `coverage_lt_1`, `coverage_incomplete`, `coverage_flag`
  (`ok`, `lt_1`, `incomplete`, `lt_1_incomplete`, `empty`)

Coverage is read from `00data/05taxaBamMap/final.ALL.taxaBamMap.txt` and the
cohort taxaBamMap files. `final.ALL` has 8288 taxa and adds `ABD_0915`,
`ABD_0938`, and `D_0132`; coverage for the 8285 common taxa is identical to
`all.ALL`. Names absent from the taxaBamMaps are not counted in
`coverage_total`; those rows are marked incomplete and the total is a lower
bound.

The IBS pair tables take sample depth from the dedup snapshot, which omits some
genotyped taxa (25 Watkins samples are reported `NotInWatkinsExcel` in the
2026-08-26 dedup run). `add_dedup_fields()` backfills those depths by sample ID
from `final.ALL.taxaBamMap.txt`, so `coverage1` / `coverage2` / `coverage_min`
are defined for every pair whose two samples are in the map.

- `low_coverage_rows.tsv` — rows with at least one name and `coverage_total < 1`.
  It also carries `coverage_mapped_names`, `bam_paths`, `bam_basenames`,
  `n_bams`, and `bam_missing_names`. A blank `coverage_missing_names` means all
  names were found in the taxaBamMaps; it does not mean the BAM is missing.
- `coverage_incomplete_rows.tsv` — rows with names missing from the taxaBamMaps
- `coverage_summary.tsv` — per-table counts

On the 2026-10-27 20:44 snapshots all 34 low-coverage rows belong to
IPK/HZAU/Watkins and every one resolves to at least one BAM (`rows_with_bam=34`,
`rows_without_bam=0`); WAP_PQ/VMap3 has none. `coverage_incomplete_rows.tsv` is
now header-only — the 83 names `济麦22` previously had missing from the
taxaBamMaps all resolve against the current maps, so no row is incomplete.

## Final combine map

`combine_map.tsv` is the merge decision in the form the release needs: **one row
per germplasm fragment**, with every genotype ID that merges into it in a single
`members` column (`|`-separated). Fill in `final_name` and that value becomes the
sample name in the released VCF.

Rules applied (decided 2026-10-28):

- a pair whose samples are all **Watkins** merges unconditionally;
- every other pair (IPK / HZAU / WAP_PQ / VMap3) merges only when its hom-SNP IBS
  is **> 0.999**;
- components are per curated germplasm, so a blocked pair splits that germplasm
  into separate rows instead of vetoing it wholesale.

Columns: `germplasm_id`, `germplasm_name`, `germplasm_name_cn`, `final_name`
(blank, curator's choice), `members`, `n_members`, `sources` (e.g.
`WAP:2|VMap3:1`), `n_pairs`, `n_pairs_merged`, `n_pairs_blocked`,
`blocked_pairs` (each with its IBS), `n_low_coverage_members`,
`low_coverage_members` (coverage < 2x, flagged not removed) and `merge_basis`
(`watkins_all` or `hom_ibs>0.999`).

`final_name` is filled automatically:

- Watkins rows take the curated `watkins_samc` value (1042 rows, all filled);
- other multi-member rows group the IDs by class prefix and list the numbers once:
  `ABD_0451|ABD_1063|ABD_1180|P0361` becomes `ABD_0451_1063_1180_P0361`, and
  `CS_sg_2014_3X|CS_sg_2017_60X` collapses to `CS_2014_2017`; a class with a
  single member keeps its original ID (including the trailing `A`);
- single-member rows keep their own ID.

Only names above 30 characters are rewritten in base 36 (`0749 -> KT`); today
that is a single row, the five-member `ABD_0749_1090_CS_2014_2017_P1359A`, which
becomes `ABD_KT_UA_CS_1JY_1K1_P11RA` (26). The 24-character
`ABD_0451_1063_1180_P0361` is left in its readable form — 26 characters is an
accepted ceiling. No name is blank or duplicated.

Merged block first, then every sample that formed no pair as a singleton row
(`merge_basis = singleton`), so the table lists all 8288 genotype IDs of
`final.ALL.taxaBamMap.txt`.

Current output: 6196 rows — 2104 merged/fragment rows from 1687 germplasms
(417 split by blocked pairs) plus 4092 singleton rows; 1361 rows carry more than
one member. `ABD_0915`, `ABD_0938` and `D_0132` were missing until 2026-10-28:
the name universe used to be `all.sample.identity.ploidy_group_ibs.samples.tsv`,
which predates their addition to `final.*`. `build_combine_map.py --taxa-map` now
takes the taxaBamMap as the universe, so new taxa cannot be dropped again. Watkins rows (1042) have 0 blocked pairs; 783 rows are affected by
the 0.999 rule; 379 rows contain at least one low-coverage member.

`combine_map.pair_decisions.tsv` is the supporting edge list: one row per pair
with its IBS, the sources involved and whether it merged.

## Taxa merge map

`final.<A|B|D|ALL>.taxaMergeMap.txt` is the authoritative merge definition: one
line per **final taxon**, in TIGER merge-map syntax
(`mergedTaxon <TAB> coverage <TAB> source1 [<TAB> source2 ...]`). A merged
germplasm lists its members; an unmerged sample maps to itself as a single
source. `coverage = .` makes TIGER sum the source coverages from the `-b`
taxaBamMap, and lines starting with `#` are skipped as headers.

The file mirrors the taxaBamMap row for row — the multiset of members is exactly
the taxon set of the matching `final.<sub>.taxaBamMap.txt`:

| file | lines | of which merged | taxa covered |
|---|---|---|---|
| `final.ALL.taxaMergeMap.txt` | 6196 | 1361 | 8288 |
| `final.A.taxaMergeMap.txt` | 5588 | 1358 | 7677 |
| `final.B.taxaMergeMap.txt` | 5507 | 1358 | 7596 |
| `final.D.taxaMergeMap.txt` | 5656 | 1354 | 7741 |

Generated by `build_scan3_merge_map.py` in
`/data/home/tusr1/01projects/vmap4/16stats.genome/40run_combine_map/`.

One file per sub-genome, because TIGER checks every source taxon against the
`-b` map of the chromosome being scanned. That is also why A/B hold 3 fewer
merged groups than ALL (the D-only groups `D_0043_0254`, `D_0036_0064`,
`D_0031_0135`) and D holds 7 fewer (the AB-only groups `AB_0089_0164`,
`AB_0091_0163`, `AB_0087_0141`, `AB_0084_0134`, `AB_0023_0170`, `AB_0038_0174`,
`AB_0032_0176`); their samples map to themselves in that sub-genome's table
because they do not exist there.

### scan3 subset (old jar only)

`fastcall3_scan3_merge_map.<A|B|D|ALL>.tsv` is the same table minus the
single-source lines (1361 / 1358 / 1358 / 1354). It exists for
`TIGER_F3_20261012.jar`, whose `ScanGenotypeF3.addMergeMapLine` throws
`Merged taxon P0332 needs at least two source taxa: P0332 . P0332`
(`ScanGenotypeF3.java:380`; runtime-confirmed on chr44, 2026-10-28).
`TIGER_F3_20261028.jar` dropped that rule — a source list may hold a single
taxon — so the default jar takes `final.<sub>.taxaMergeMap.txt` directly and the
subset is only needed to rerun on the older jar.

Usage — chromosomes `0`, `43`, `44` map to the `ALL` sub-genome in
`getChrConfig`, A/B/D chromosomes use the matching table. `TIGER_F3_20261028.jar`
is the pipeline default (`conf/base.config`), registered in
`infra/infra_tiger.nf`:

```
java -Xmx<g>g -jar TIGER_F3_20261028.jar -app FastCall3 -mod scan3 \
  -a <ref> -b 00data/05taxaBamMap/final.<A|B|D|ALL>.taxaBamMap.txt \
  -c <vLib/NN_1_len.lib.gz> -d <chr> -h 0.05 -j 32 -k <gen3> -m 0.2 \
  -o 04runScreens/final -p final.<sub>.taxaMergeMap.txt
```

Verified on chr44 through the pipeline (2026-10-28,
`04runScreens/22run_scan3_chr44_taxamerge`, new default jar):

- 6196 merge groups parsed, 0 `missing source iac` warnings, 8288 taxa read from
  `-b`, 6124 sites, **6196 sample columns** — every final taxon in one VCF;
- the 1361 merged columns are GT:AD:GL-identical to the run that used the
  multi-source subset on `TIGER_F3_20261012.jar` (8,334,764 cells, 0 differences);
- singleton columns are literal passthrough: the merged
  `gen3/indiCounts/<taxon>/44_1_134546.iac.gz` is byte-identical to
  `final/indiCounts/<taxon>/44_1_134546.iac.gz`, so no AD is summed twice;
- VCF column order is the merged name's dictionary order, not the file's line
  order (per the jar notes).

Earlier, `fastcall3_scan3_merge_map.ALL.tsv` on the old jar produced the same
6124 sites with only the 1361 merged columns
(`16stats.genome/40run_combine_map/scan3_smoke_chr44`). scan3 sets its sample
list from the merge-map keys, so the complete map is what makes the one-VCF
release possible; with the subset the release would need the `scan` VCF minus
the 3453 merged-away sources plus the 1361 scan3 columns.

## Pair-level decisions

`combine_pairs_decisions_*_hom_ibs.tsv` replaces the whole-group verdict with a
status per pair:

- `pass` — `hom_ibs` clears the cohort cutoff
- `fail` — below the cutoff with both samples at or above the coverage floor
- `unjudged_low_coverage` — at least one sample is under 2x coverage or has no
  depth at all, so the genotype call cannot be trusted either way

Each row also carries `above_0_999` so a stricter reading needs no rerun, and a
blank `decision` column for annotation. The companion
`combine_groups_support_*_hom_ibs.tsv` describes each group by the mix of its
pairs instead of rejecting it wholesale. Coverage dominates the split: 3,505 of
4,783 IPK pairs (73%) are unjudged for coverage (1,278 pass, none fail) even
though only 22% of IPK samples sit under 2x, because those samples hang off the
large multi-sample groups. Those are 679 fully supported groups and 374 blocked
only by coverage. WAP_PQ/VMap3 splits 552 pass / 323 fail / 3 unjudged, with 405
groups fully supported, 225 containing a real fail, and 2 blocked only by
coverage.

## IBS merge audit (hom-SNP IBS)

A combine group is treated as merging only when every genotype pair inside it
exceeds the cohort cutoff — the hom-SNP IBS is the merge basis; plain PLINK
IBS is **not** a merge criterion and survives only as a diagnostic (the one
`*_plink_ibs`-judged file still here, `combine_pairs_not_merged_ipk_hzau_watkins.tsv`)
and should not drive decisions. Cutoffs differ per cohort:

| Cohort | Group cutoff | Pairs cutoff |
|---|---|---|
| WAP_PQ/VMap3 | `hom_ibs > 0.995` | `hom_ibs > 0.99` |
| IPK/HZAU/Watkins | `hom_ibs > 0.99` | `hom_ibs > 0.99` |

The WAP_PQ/VMap3 pairs table is deliberately the looser 0.99 list — it is the
annotation worklist, a superset of what the 0.995 group cutoff rejects.
Derived from the published pair tables, not hand-curated — regenerate after the
pair table changes.

- `combine_groups_not_merged.tsv` — one row per failing group, sorted by the
  weakest pair: `germplasm_name`, `germplasm_id`, `group_key`,
  `check_status` (`no_pair_above_threshold` when no pair clears the cutoff,
  `partial_below_threshold` when only some do), `n_pairs`,
  `n_pairs_above_threshold`, `frac_pairs_above_threshold`, min/max
  `plink_ibs` and `hom_ibs`, the member `samples` and `ploidies`, and the
  weakest pair inline (`worst_*`, including `worst_coverage_min`).
- `combine_pairs_*_hom_ibs.tsv` — every pair at or below the cohort cutoff:
  samples, raw combine names, `pair_type`, ploidies, `same_ploidy`,
  `plink_ibs`, `hom_ibs`, coverage, dedup class, tier, the PLINK random
  class, and a blank `decision` column.

The `decision` column in the two `combine_pairs_*_hom_ibs.tsv` files is for
hand annotation of which pairs are dropped from the merge (suggested values
`keep_merge` / `discard` / `review`; blank means not yet judged). Pass
`--carry-decisions` when regenerating: it re-attaches the previous file's
values by `germplasm_id` + `sample1` + `sample2` and reports how many it
carried. Without the flag the column comes back blank.

`*_hom_ibs.tsv` covers the WAP_PQ/VMap3 cohort; the
`_ipk_hzau_watkins_hom_ibs` files cover the IPK/HZAU/Watkins cohort, generated
with `--name-columns ipk_combine_names,hzau_combine_names,watkins_combine_names`.

Results on the 2026-10-28 rebuilt pair tables, judged on `hom_ibs`:

| Cohort | Cutoff | Groups failing | Pairs ≤ cutoff |
|---|---|---|---|
| WAP_PQ/VMap3 | 0.995 | 228 / 631 | 329 / 882 |
| IPK/HZAU/Watkins | 0.99 | 5 / 1053 | 5 / 4783 |

At a uniform 0.99 the WAP_PQ/VMap3 cohort fails 190 of 631 groups; the 0.995
cutoff used for its groups adds the 38 groups whose weakest pair lands between
the two lines.

The same tables judged on `plink_ibs` fail far more (266/631 and 77/1053),
which is why the metric matters: PLINK IBS mixes heterozygous sites in, so a
true duplicate can sit below 0.99 there while its hom-SNP IBS is clean.
Example: `WATDE0278` (1128 pairs) has 11 pairs below 0.99 on plink_ibs but
none on hom_ibs, so it passes the merge criterion.

The 4 cross-ploidy pairs in the WAP_PQ/VMap3 table are flagged
`same_ploidy=no` — a merge decision should weigh the same-ploidy pairs.

`group_key` uses the curated `germplasm_name`, falling back to
`germplasm_id`: 40 failing groups have no `germplasm_name` in the source
table and are identified by UUID only. Every pair is mapped back to its
curated row through its two raw combine names, which also recovers the 122
pairs whose pair-table `germplasm` field is empty (631 groups instead of 560
if grouped by that field alone).

Regenerate with `extract_not_merged.py` from
`/data/home/tusr1/01projects/vmap4/16stats.genome/04run_combine_ibs_merge_threshold/`:

```bash
cd /data/home/tusr1/01projects/vmap4/16stats.genome/04run_combine_ibs_merge_threshold && \
/data/home/tusr1/miniconda3/envs/stats/bin/python -I extract_not_merged.py \
  --pairs /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/stats/test_common_inter/info/all.sample.combine_ibs.pairs.tsv \
  --combine /data/home/tusr1/git/script/src/resource/sample_combine/curation_germplasm_combine_name_groups.tsv \
  --metric hom_ibs --threshold 0.995 --pairs-threshold 0.99 \
  --suffix _hom_ibs --outdir /data/home/tusr1/git/script/src/resource/sample_combine \
  --carry-decisions
```

The IPK/HZAU/Watkins cohort uses its own pair table, combine table and name
columns; write it to a separate `--outdir` so the two cohorts never overwrite
each other:

```bash
cd /data/home/tusr1/01projects/vmap4/16stats.genome/04run_combine_ibs_merge_threshold && \
/data/home/tusr1/miniconda3/envs/stats/bin/python -I extract_not_merged.py \
  --pairs /data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/stats/test_common_inter/info/all.sample.combine_ibs_ipk_hzau_watkins.pairs.tsv \
  --combine /data/home/tusr1/git/script/src/resource/sample_combine/curation_germplasm_combine_name_groups_ipk_hzau_watkins.tsv \
  --name-columns ipk_combine_names,hzau_combine_names,watkins_combine_names \
  --metric hom_ibs --threshold 0.99 --pairs-threshold 0.99 \
  --suffix _ipk_hzau_watkins_hom_ibs \
  --outdir /data/home/tusr1/git/script/src/resource/sample_combine \
  --carry-decisions
```

`--suffix` keeps the two cohorts from overwriting each other; `--threshold`
judges groups, `--pairs-threshold` builds the annotation worklist (WAP keeps it
looser at 0.99). `--metric plink_ibs` produces diagnostic-only files.
