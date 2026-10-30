"""Validate combine-name duplicate pairs against test-set IBS matrices.

The combine-name table links WAP_PQ and VMap3 names for the same germplasm.
This module maps those names to genotype sample IDs, removes curation-excluded
control slots, scores each pair from the published square IBS matrices, and
groups the AABBDD pairs into low/mid/high IBS tiers.

It is a library called from Nextflow: paths are passed in, and no publish-tree
paths or command-line entry points live here.
"""

from __future__ import annotations

import hashlib
import logging
import mmap
import os
import re
from collections.abc import Sequence
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from genetics.genomics.sample.sample_utils import load_ids_file
from infra.utils.errors import fail
from infra.utils.graph import (
    plot_connected_scatter,
    plot_grouped_histogram,
    plot_kde_density_lines,
    plot_scatter_with_thresholds,
)
from infra.utils.io import save_df_to_tsv

logger = logging.getLogger(__name__)

PLOIDY_SUBGENOMES = {
    "AABBDD": "ABD",
    "AABB": "AB",
    "DD": "D",
    "AA": "A",
    "SS": "B",
}
TIER_ORDER = ("low", "mid", "high")
TIER_LOW_MAX = 0.975
TIER_HIGH_MIN = 0.990
TIER_PALETTE = {
    "low": "#2a78d6",
    "mid": "#eb6834",
    "high": "#1baf7a",
}
RANDOM_PALETTE = {
    "random": "#2a78d6",
    "duplicate": "#eb6834",
}
DEDUP_PALETTE = {
    "random": "#2a78d6",
    "same_dedup_group": "#eb6834",
    "different_dedup_group": "#1baf7a",
}
RANDOM_HUE_ORDER = ("random", "duplicate")
DEDUP_HUE_ORDER = ("random", "same_dedup_group", "different_dedup_group")
P_A_NAME_RE = re.compile(r"^P\d+A$")
# Depth bins widen with coverage: coverage clusters on a few discrete
# sequencing depths, so equal-width or narrow bins follow those clusters
# instead of the trend. These edges separate the two low-depth clusters
# (<0.15x, 0.15-0.3x), the 0.3-0.75x block, and the callability transition.
LOW_COVERAGE_BINS = (0.05, 0.15, 0.30, 0.50, 0.75, 1.5, 3.0, 20.0)
LOW_COVERAGE_POINT_COLOR = "#2a78d6"
LOW_COVERAGE_HIGHLIGHT_COLOR = "#eb6834"
HIGHLIGHT_PAIR_COLOR = "#1a1a1a"


def parse_name_list(value: str) -> list[str]:
    """Split a pipe-separated combine-name cell into stripped names."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    return [part.strip() for part in str(value).split("|") if part.strip()]


def map_to_genotype_name(name: str, genotype_names: set[str]) -> str | None:
    """Map a combine name to its genotype ID, trying the trailing-A form.

    WAP combine names omit the trailing ``A`` that the genotype and dedup
    tables use for the P0030A-P2096A series. Exact matches win; otherwise the
    ``name + 'A'`` form is used when present.
    """
    if name in genotype_names:
        return name
    candidate = f"{name}A"
    if candidate in genotype_names:
        return candidate
    return None


def load_excluded_controls(paths: list[str | Path] | None) -> set[str]:
    """Collect combine names whose curation status is ``excluded``."""
    excluded: set[str] = set()
    for path in paths or []:
        if not path:
            continue
        path = Path(path)
        if not path.exists():
            fail(f"Excluded-control table not found: {path}")
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        if "germplasm_combine_name_sample" not in df.columns:
            fail(
                f"{path} is missing germplasm_combine_name_sample; "
                f"columns={list(df.columns)}"
            )
        if "status" in df.columns:
            mask = df["status"].str.lower() == "excluded"
        else:
            mask = pd.Series(True, index=df.index)
        names = df.loc[mask, "germplasm_combine_name_sample"]
        excluded.update(name.strip() for name in names if name.strip())
    return excluded


def load_genotype_meta(samples_path: str | Path) -> pd.DataFrame:
    """Load sample ploidy and subgenome assignments."""
    path = Path(samples_path)
    if not path.exists():
        fail(f"Genotype sample table not found: {path}")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"Sample", "ploidy", "subgenomes"}
    missing = required.difference(df.columns)
    if missing:
        fail(f"{path} is missing columns: {sorted(missing)}")
    return df[["Sample", "ploidy", "subgenomes"]].drop_duplicates("Sample")


def _pair_type_from_cohorts(first: set[str], second: set[str]) -> str:
    shared = first & second
    if shared:
        if len(first) == 1 and len(second) == 1:
            return f"within_{next(iter(first))}"
        return "other"
    return "cross"


def build_combine_pairs(
    combine_path: str | Path,
    genotype_names: set[str],
    excluded_names: set[str] | None = None,
    *,
    name_columns: Sequence[str] | None = None,
    cohort_labels: Sequence[str] | None = None,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Build genotype-available name pairs from a combine table.

    ``name_columns`` lists the cohort name columns (default: WAP_PQ and VMap3).
    ``cohort_labels`` gives one label per column; labels are used in
    ``pair_type`` (``cross`` / ``within_<cohort>``). Excluded control names are
    removed before mapping, so the Jimai22 control slots contribute no pairs.
    """
    path = Path(combine_path)
    if not path.exists():
        fail(f"Combine-name table not found: {path}")
    excluded_names = excluded_names or set()
    columns = list(name_columns or ("wap_pq_combine_names", "vmap3_combine_names"))
    labels = list(cohort_labels or ("wap", "vmap3"))
    if len(columns) != len(labels):
        fail("name_columns and cohort_labels must have the same length")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"germplasm_name", *columns}
    missing = required.difference(df.columns)
    if missing:
        fail(f"{path} is missing columns: {sorted(missing)}")

    stats = {
        "germplasm_rows": int(len(df)),
        "single_name_rows_raw": 0,
        "multi_name_rows_raw": 0,
        "rows_with_zero_names_after_exclusion": 0,
        "rows_with_one_name_after_mapping": 0,
        "rows_with_ge2_genotype_names": 0,
        "excluded_name_occurrences": 0,
    }
    records: dict[tuple[str, str], dict] = {}

    for _, row in df.iterrows():
        column_names = [
            parse_name_list(row.get(column, "")) for column in columns
        ]
        raw_names = [name for names in column_names for name in names]
        n_raw = len(raw_names)
        if n_raw == 1:
            stats["single_name_rows_raw"] += 1
        elif n_raw > 1:
            stats["multi_name_rows_raw"] += 1
        stats["excluded_name_occurrences"] += sum(
            1 for name in raw_names if name in excluded_names
        )

        kept_by_column = [
            [name for name in names if name not in excluded_names]
            for names in column_names
        ]
        if not any(kept_by_column):
            stats["rows_with_zero_names_after_exclusion"] += 1
            continue

        cohorts_by_sample: dict[str, set[str]] = {}
        raw_by_sample: dict[str, str] = {}
        for names, label in zip(kept_by_column, labels):
            for name in names:
                genotype = map_to_genotype_name(name, genotype_names)
                if genotype:
                    cohorts_by_sample.setdefault(genotype, set()).add(label)
                    raw_by_sample.setdefault(genotype, name)

        union = sorted(cohorts_by_sample)
        if len(union) < 2:
            stats["rows_with_one_name_after_mapping"] += 1
            continue
        stats["rows_with_ge2_genotype_names"] += 1
        union_n = len(union)
        germplasm = str(row.get("germplasm_name", "")).strip()

        for sample1, sample2 in combinations(union, 2):
            key = (sample1, sample2)
            pair_type = _pair_type_from_cohorts(
                cohorts_by_sample[sample1], cohorts_by_sample[sample2]
            )
            name1_raw = raw_by_sample.get(sample1, "")
            name2_raw = raw_by_sample.get(sample2, "")
            p_has_a = int(
                bool(P_A_NAME_RE.match(sample1))
                or bool(P_A_NAME_RE.match(sample2))
            )
            current = records.get(key)
            if current is None:
                records[key] = {
                    "sample1": sample1,
                    "sample2": sample2,
                    "name1_raw": name1_raw,
                    "name2_raw": name2_raw,
                    "pair_type": pair_type,
                    "p_has_A": p_has_a,
                    "union_n": union_n,
                    "germplasm": {germplasm} if germplasm else set(),
                }
            else:
                if pair_type == "cross":
                    current["pair_type"] = "cross"
                current["union_n"] = max(current["union_n"], union_n)
                if germplasm:
                    current["germplasm"].add(germplasm)

    rows = []
    for record in records.values():
        record = dict(record)
        record["germplasm"] = "|".join(sorted(record["germplasm"]))
        rows.append(record)
    pairs = pd.DataFrame(rows)
    if not pairs.empty:
        pairs = pairs.sort_values(["sample1", "sample2"]).reset_index(drop=True)
    return pairs, stats


def add_dedup_fields(
    pairs: pd.DataFrame,
    dedup_path: str | Path,
    coverage_map: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Join sample-level dedup group, status, coverage, and name agreement.

    ``coverage_map`` (a taxaBamMap loaded with :func:`load_taxa_coverage`)
    fills sample depths that the dedup snapshot lacks: the dedup input list can
    omit genotyped taxa, and a missing depth would otherwise leave the pair
    without ``coverage_min``. Filled values are formatted like the dedup ones.
    """
    path = Path(dedup_path)
    if not path.exists():
        fail(f"Dedup table not found: {path}")
    dedup = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"Sample", "Dup_Group_ID", "Dup_Status"}
    missing = required.difference(dedup.columns)
    if missing:
        fail(f"{path} is missing columns: {sorted(missing)}")

    records = {row["Sample"]: row for _, row in dedup.iterrows()}
    out = pairs.copy()
    n_filled = 0
    for index, row in out.iterrows():
        first = records.get(row["sample1"], {})
        second = records.get(row["sample2"], {})
        group1 = first.get("Dup_Group_ID", "")
        group2 = second.get("Dup_Group_ID", "")
        if group1 and group2:
            dedup_class = (
                "same_dedup_group"
                if group1 == group2
                else "different_dedup_group"
            )
        else:
            dedup_class = "not_in_dedup"
        out.at[index, "dedup_class"] = dedup_class
        out.at[index, "status1"] = first.get("Dup_Status", "")
        out.at[index, "status2"] = second.get("Dup_Status", "")
        out.at[index, "dup_group1"] = group1
        out.at[index, "dup_group2"] = group2
        cov1 = pd.to_numeric(first.get("Coverage", ""), errors="coerce")
        cov2 = pd.to_numeric(second.get("Coverage", ""), errors="coerce")
        coverage1 = first.get("Coverage", "")
        coverage2 = second.get("Coverage", "")
        if coverage_map:
            if pd.isna(cov1) and row["sample1"] in coverage_map:
                cov1 = coverage_map[row["sample1"]]
                coverage1 = f"{cov1:.4f}"
                n_filled += 1
            if pd.isna(cov2) and row["sample2"] in coverage_map:
                cov2 = coverage_map[row["sample2"]]
                coverage2 = f"{cov2:.4f}"
                n_filled += 1
        out.at[index, "coverage1"] = coverage1
        out.at[index, "coverage2"] = coverage2
        out.at[index, "coverage_min"] = (
            float(min(cov1, cov2)) if pd.notna(cov1) and pd.notna(cov2) else np.nan
        )
        cn1 = first.get("Clean_CN", "")
        cn2 = second.get("Clean_CN", "")
        acc1 = first.get("Clean_Acc", "")
        acc2 = second.get("Clean_Acc", "")
        out.at[index, "clean_cn1"] = cn1
        out.at[index, "clean_cn2"] = cn2
        out.at[index, "clean_acc1"] = acc1
        out.at[index, "clean_acc2"] = acc2
        out.at[index, "cn_status"] = _agreement_status(cn1, cn2)
        out.at[index, "acc_status"] = _agreement_status(acc1, acc2)
    if coverage_map:
        logger.info(
            "Filled %s pair sample depths from the taxaBamMap", n_filled
        )
    return out


def _agreement_status(first: str, second: str) -> str:
    if not first or not second:
        return "NA"
    return "match" if first == second else "mismatch"


def add_ploidy_fields(pairs: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """Attach ploidy and common subgenomes for each pair."""
    meta_by_sample = meta.set_index("Sample").to_dict("index")
    out = pairs.copy()
    for index, row in out.iterrows():
        first = meta_by_sample.get(row["sample1"], {})
        second = meta_by_sample.get(row["sample2"], {})
        ploidy1 = first.get("ploidy", "")
        ploidy2 = second.get("ploidy", "")
        sub1 = first.get("subgenomes", "") or PLOIDY_SUBGENOMES.get(ploidy1, "")
        sub2 = second.get("subgenomes", "") or PLOIDY_SUBGENOMES.get(ploidy2, "")
        out.at[index, "ploidy1"] = ploidy1
        out.at[index, "ploidy2"] = ploidy2
        out.at[index, "ploidy"] = f"{ploidy1}/{ploidy2}"
        out.at[index, "subgenomes"] = sub1 + "/" + sub2
        out.at[index, "common_subgenomes"] = "".join(
            sg for sg in "ABD" if sg in sub1 and sg in sub2
        )
    return out


class _TextMatrix:
    """Lazy reader for a PLINK square IBS text matrix.

    It memory-maps the file and parses only requested sample rows. This avoids
    loading a 500 MB matrix per subgenome and metric into memory.
    """

    def __init__(self, matrix_path: str | Path, id_path: str | Path):
        self.matrix_path = str(matrix_path)
        self.ids = {
            name: index
            for index, name in enumerate(load_ids_file(str(id_path)))
        }
        self._handle = open(self.matrix_path, "rb")
        self._mmap = mmap.mmap(
            self._handle.fileno(), 0, access=mmap.ACCESS_READ
        )
        offsets = [0]
        position = 0
        while True:
            newline = self._mmap.find(b"\n", position)
            if newline < 0:
                break
            offsets.append(newline + 1)
            position = newline + 1
        if offsets[-1] != len(self._mmap):
            offsets.append(len(self._mmap))
        self._offsets = offsets
        self._cache: dict[str, np.ndarray] = {}

    def value(self, first: str, second: str) -> float | None:
        left = self.ids.get(first)
        right = self.ids.get(second)
        if left is None or right is None:
            return None
        row = self._cache.get(first)
        if row is None:
            line = self._mmap[
                self._offsets[left] : self._offsets[left + 1]
            ].rstrip(b"\n")
            row = np.fromstring(
                line.replace(b"\t", b" "), dtype=np.float32, sep=" "
            )
            self._cache[first] = row
        if right >= len(row):
            return None
        return float(row[right])

    def close(self) -> None:
        self._mmap.close()
        self._handle.close()


def add_ibs_fields(
    pairs: pd.DataFrame,
    ibs_dir: str | Path,
    hom_dir: str | Path,
) -> pd.DataFrame:
    """Score each pair from the PLINK and hom-SNP square IBS matrices."""
    out = pairs.copy()
    n_rows = len(out)
    plink_values: dict[tuple[int, str], float] = {}
    hom_values: dict[tuple[int, str], float] = {}
    subgenomes = ("A", "B", "D")

    for sg in subgenomes:
        wanted = [
            index
            for index, row in out.iterrows()
            if sg in str(row.get("common_subgenomes", ""))
        ]
        if not wanted:
            continue
        plink_matrix = _TextMatrix(
            Path(ibs_dir) / f"{sg}.ibs.mibs",
            Path(ibs_dir) / f"{sg}.ibs.mibs.id",
        )
        hom_matrix = _TextMatrix(
            Path(hom_dir) / f"{sg}.ibs_hom_snp.mibs",
            Path(hom_dir) / f"{sg}.ibs_hom_snp.mibs.id",
        )
        try:
            for index in wanted:
                row = out.loc[index]
                first = row["sample1"]
                second = row["sample2"]
                plink = plink_matrix.value(first, second)
                hom = hom_matrix.value(first, second)
                if plink is not None:
                    plink_values[(index, sg)] = plink
                if hom is not None:
                    hom_values[(index, sg)] = hom
        finally:
            plink_matrix.close()
            hom_matrix.close()
        logger.info(
            "IBS %s: %s pairs requested, %s PLINK, %s hom-SNP values",
            sg,
            len(wanted),
            sum(1 for key in plink_values if key[1] == sg),
            sum(1 for key in hom_values if key[1] == sg),
        )

    plink_mean: list[float] = []
    hom_mean: list[float] = []
    n_sg: list[int] = []
    plink_by_sg: list[str] = []
    hom_by_sg: list[str] = []
    for index in range(n_rows):
        plink_parts = []
        hom_parts = []
        for sg in subgenomes:
            plink = plink_values.get((index, sg))
            hom = hom_values.get((index, sg))
            if plink is not None and hom is not None:
                plink_parts.append(f"{sg}:{plink:.4f}")
                hom_parts.append(f"{sg}:{hom:.4f}")
        plink_mean.append(
            float(np.mean([plink_values[(index, sg)] for sg in subgenomes
                           if (index, sg) in plink_values]))
            if any((index, sg) in plink_values for sg in subgenomes)
            else np.nan
        )
        hom_mean.append(
            float(np.mean([hom_values[(index, sg)] for sg in subgenomes
                           if (index, sg) in hom_values]))
            if any((index, sg) in hom_values for sg in subgenomes)
            else np.nan
        )
        n_sg.append(len(plink_parts))
        plink_by_sg.append(";".join(plink_parts))
        hom_by_sg.append(";".join(hom_parts))

    out["plink_ibs"] = plink_mean
    out["hom_ibs"] = hom_mean
    out["n_subgenomes"] = n_sg
    out["plink_by_sg"] = plink_by_sg
    out["hom_by_sg"] = hom_by_sg
    logger.info("Scored %s pairs with at least one common subgenome", int(
        (out["n_subgenomes"] > 0).sum()
    ))
    return out


def assign_ibs_tiers(
    pairs: pd.DataFrame,
    *,
    metric: str = "plink_ibs",
    low_max: float = TIER_LOW_MAX,
    high_min: float = TIER_HIGH_MIN,
) -> pd.DataFrame:
    """Assign low/mid/high tiers to AABBDD pairs from one IBS metric.

    Defaults are the PLINK cutoffs (< 0.975 low, 0.975-0.990 mid, >= 0.990
    high). For the hom-SNP matrix pass ``metric='hom_ibs'`` and its own
    cutoffs; the hom-SNP distribution sits higher than the PLINK one.
    """
    if metric not in pairs.columns:
        fail(f"IBS metric column not found: {metric}")
    out = pairs.copy()
    out["tier"] = "NA"
    mask = (
        (out["ploidy1"] == "AABBDD")
        & (out["ploidy2"] == "AABBDD")
        & out[metric].notna()
    )
    values = pd.to_numeric(out.loc[mask, metric], errors="coerce")
    tier = pd.Series("high", index=values.index)
    tier.loc[values < float(high_min)] = "mid"
    tier.loc[values < float(low_max)] = "low"
    out.loc[mask, "tier"] = tier
    logger.info(
        "IBS tiers (%s, cutoffs %.4g/%.4g): %s",
        metric,
        low_max,
        high_min,
        out.loc[mask, "tier"].value_counts().to_dict(),
    )
    return out


def _numeric_median(df: pd.DataFrame, column: str) -> float:
    if column not in df.columns:
        return float("nan")
    values = pd.to_numeric(df[column], errors="coerce").dropna()
    return float(values.median()) if not values.empty else float("nan")


def _random_group_stats(
    random_pairs_path: str | Path,
) -> tuple[pd.DataFrame, dict[str, dict]]:
    """Return the random pair table plus per-ploidy sorted values and stats."""
    path = Path(random_pairs_path)
    if not path.exists():
        fail(f"Random-pair table not found: {path}")
    random_pairs = pd.read_csv(path, sep="\t")
    required = {"ploidy", "plink_ibs", "hom_ibs"}
    missing = required.difference(random_pairs.columns)
    if missing:
        fail(f"{path} is missing columns: {sorted(missing)}")
    by_group: dict[str, dict] = {}
    for group, sub in random_pairs.groupby("ploidy"):
        plink = (
            pd.to_numeric(sub["plink_ibs"], errors="coerce")
            .dropna()
            .to_numpy(dtype=float)
        )
        hom = (
            pd.to_numeric(sub["hom_ibs"], errors="coerce")
            .dropna()
            .to_numpy(dtype=float)
        )
        if plink.size == 0:
            continue
        by_group[str(group)] = {
            "n": int(plink.size),
            "plink_sorted": np.sort(plink),
            "hom_sorted": np.sort(hom),
            "plink_median": float(np.median(plink)),
            "plink_p95": float(np.quantile(plink, 0.95)),
            "plink_p99": float(np.quantile(plink, 0.99)),
            "hom_median": float(np.median(hom)) if hom.size else float("nan"),
            "hom_p95": float(np.quantile(hom, 0.95)) if hom.size else float("nan"),
            "hom_p99": float(np.quantile(hom, 0.99)) if hom.size else float("nan"),
            "plink_A_median": _numeric_median(sub, "plink_ibs_A"),
            "plink_B_median": _numeric_median(sub, "plink_ibs_B"),
            "plink_D_median": _numeric_median(sub, "plink_ibs_D"),
            "hom_A_median": _numeric_median(sub, "hom_ibs_A"),
            "hom_B_median": _numeric_median(sub, "hom_ibs_B"),
            "hom_D_median": _numeric_median(sub, "hom_ibs_D"),
        }
    if not by_group:
        fail("Random-pair table has no usable ploidy groups")
    return random_pairs, by_group


def _random_reference(
    random_pairs_path: str | Path,
) -> tuple[pd.DataFrame, dict[str, float], pd.DataFrame, dict[str, dict]]:
    random_pairs, by_group = _random_group_stats(random_pairs_path)
    if "AABBDD" not in by_group:
        fail("Random-pair table has no AABBDD rows")
    stats_abd = by_group["AABBDD"]
    stats = {
        "random_n": stats_abd["n"],
        "random_plink_median": stats_abd["plink_median"],
        "random_plink_p95": stats_abd["plink_p95"],
        "random_plink_p99": stats_abd["plink_p99"],
        "random_hom_median": stats_abd["hom_median"],
        "random_hom_p95": stats_abd["hom_p95"],
        "random_hom_p99": stats_abd["hom_p99"],
    }
    abd = random_pairs.loc[random_pairs["ploidy"] == "AABBDD"].copy()
    return abd, stats, random_pairs, by_group


def assign_random_classes(
    pairs: pd.DataFrame,
    by_group: dict[str, dict],
) -> pd.DataFrame:
    """Add per-ploidy random percentiles and p95/p99 classes to every pair."""
    out = pairs.copy()
    out["random_group"] = "NA"
    out["random_n"] = 0
    for prefix in ("plink", "hom"):
        out[f"{prefix}_percentile"] = np.nan
        out[f"{prefix}_above_p95"] = False
        out[f"{prefix}_above_p99"] = False
        out[f"{prefix}_random_class"] = "NA"
    mask = (
        (out["ploidy1"] == out["ploidy2"])
        & out["ploidy1"].astype(str).isin(by_group)
    )
    for group, index in out.loc[mask].groupby("ploidy1").groups.items():
        stats = by_group[str(group)]
        sub = out.loc[index]
        for prefix, column in (
            ("plink", "plink_ibs"),
            ("hom", "hom_ibs"),
        ):
            values = pd.to_numeric(
                sub[column], errors="coerce"
            ).to_numpy(dtype=float)
            valid = np.isfinite(values)
            if not valid.any():
                continue
            sorted_values = stats[f"{prefix}_sorted"]
            p95 = float(stats[f"{prefix}_p95"])
            p99 = float(stats[f"{prefix}_p99"])
            percentiles = (
                np.searchsorted(sorted_values, values[valid], side="left")
                / float(stats["n"])
            )
            selected = values[valid]
            classes = np.where(
                selected >= p99,
                "high",
                np.where(selected >= p95, "mid", "low"),
            )
            target = sub.index[valid]
            out.loc[target, f"{prefix}_percentile"] = percentiles
            out.loc[target, f"{prefix}_above_p95"] = selected >= p95
            out.loc[target, f"{prefix}_above_p99"] = selected >= p99
            out.loc[target, f"{prefix}_random_class"] = classes
        out.loc[sub.index, "random_group"] = str(group)
        out.loc[sub.index, "random_n"] = int(stats["n"])
    return out


def _subgenome_median(series: pd.Series, subgenome: str) -> float:
    values: list[float] = []
    for text in series.dropna():
        parts = {
            part.split(":", 1)[0]: part.split(":", 1)[1]
            for part in str(text).split(";")
            if ":" in part
        }
        if subgenome in parts:
            values.append(float(parts[subgenome]))
    return float(np.median(values)) if values else float("nan")


def summarize_by_ploidy(
    pairs: pd.DataFrame,
    by_group: dict[str, dict],
) -> pd.DataFrame:
    """One-row-per-ploidy-group summary, using the group's own subgenomes."""
    order = ["AABBDD", "AABB", "DD", "AA", "SS"]
    present = sorted(
        set(pairs.loc[pairs["ploidy1"] == pairs["ploidy2"], "ploidy1"])
    )
    groups = [group for group in order if group in present]
    groups += [group for group in present if group not in groups]
    rows: list[dict] = []
    for group in groups:
        sub = pairs[
            (pairs["ploidy1"] == group) & (pairs["ploidy2"] == group)
        ].copy()
        if sub.empty:
            continue
        plink = pd.to_numeric(sub["plink_ibs"], errors="coerce").dropna()
        hom = pd.to_numeric(sub["hom_ibs"], errors="coerce").dropna()
        stats = by_group.get(group)
        subgenomes = (
            str(sub["common_subgenomes"].mode().iloc[0])
            if not sub["common_subgenomes"].empty
            else ""
        )
        row = {
            "ploidy_group": group,
            "n_pairs": int(len(sub)),
            "subgenomes": subgenomes,
            "dup_plink_median": float(plink.median()) if not plink.empty else np.nan,
            "dup_plink_p95": float(plink.quantile(0.95)) if not plink.empty else np.nan,
            "dup_plink_p99": float(plink.quantile(0.99)) if not plink.empty else np.nan,
            "dup_hom_median": float(hom.median()) if not hom.empty else np.nan,
            "dup_hom_p95": float(hom.quantile(0.95)) if not hom.empty else np.nan,
            "dup_hom_p99": float(hom.quantile(0.99)) if not hom.empty else np.nan,
            "random_n": int(stats["n"]) if stats else 0,
            "random_plink_median": stats["plink_median"] if stats else np.nan,
            "random_plink_p95": stats["plink_p95"] if stats else np.nan,
            "random_plink_p99": stats["plink_p99"] if stats else np.nan,
            "random_hom_median": stats["hom_median"] if stats else np.nan,
            "random_hom_p95": stats["hom_p95"] if stats else np.nan,
            "random_hom_p99": stats["hom_p99"] if stats else np.nan,
        }
        if stats:
            row["frac_plink_above_p95"] = float(
                (plink >= stats["plink_p95"]).mean()
            )
            row["frac_plink_above_p99"] = float(
                (plink >= stats["plink_p99"]).mean()
            )
            row["frac_hom_above_p95"] = float(
                (hom >= stats["hom_p95"]).mean()
            )
            row["frac_hom_above_p99"] = float(
                (hom >= stats["hom_p99"]).mean()
            )
            row["class_plink_high"] = int((plink >= stats["plink_p99"]).sum())
            row["class_plink_mid"] = int(
                ((plink >= stats["plink_p95"]) & (plink < stats["plink_p99"])).sum()
            )
            row["class_plink_low"] = int((plink < stats["plink_p95"]).sum())
            row["class_hom_high"] = int((hom >= stats["hom_p99"]).sum())
            row["class_hom_mid"] = int(
                ((hom >= stats["hom_p95"]) & (hom < stats["hom_p99"])).sum()
            )
            row["class_hom_low"] = int((hom < stats["hom_p95"]).sum())
        for sg in ("A", "B", "D"):
            row[f"dup_hom_{sg}_median"] = _subgenome_median(
                sub["hom_by_sg"], sg
            )
            row[f"random_hom_{sg}_median"] = (
                stats.get(f"hom_{sg}_median", np.nan) if stats else np.nan
            )
        rows.append(row)
    return pd.DataFrame(rows)


def summarize_pairs(
    pairs: pd.DataFrame,
    stats: dict[str, int],
    random_stats: dict[str, float],
    *,
    combine_path: str | Path,
    excluded_names: set[str],
    ibs_dir: str | Path,
    hom_dir: str | Path,
    name_columns: Sequence[str] | None = None,
    cohort_labels: Sequence[str] | None = None,
    metric: str = "plink_ibs",
    low_max: float = TIER_LOW_MAX,
    high_min: float = TIER_HIGH_MIN,
) -> pd.DataFrame:
    """Build the long-form metric/value summary table."""
    combine_path = Path(combine_path)
    metrics: list[tuple[str, object]] = [
        ("combine_file", str(combine_path)),
        (
            "combine_mtime",
            pd.Timestamp(os.path.getmtime(combine_path), unit="s").isoformat(),
        ),
        ("combine_md5", _md5(combine_path)),
        ("combine_name_columns", ",".join(name_columns or ())),
        ("combine_cohort_labels", ",".join(cohort_labels or ())),
        ("excluded_names", len(excluded_names)),
        ("tier_metric", metric),
        ("tier_low_max", float(low_max)),
        ("tier_high_min", float(high_min)),
    ]
    metrics.extend(sorted(stats.items()))
    metrics.extend(
        [
            ("ibs_dir", str(ibs_dir)),
            ("hom_ibs_dir", str(hom_dir)),
            ("genotype_pairs", int(len(pairs))),
            (
                "paired_samples",
                int(len(set(pairs["sample1"]) | set(pairs["sample2"])))
                if not pairs.empty
                else 0,
            ),
        ]
    )
    if not pairs.empty:
        for ploidy, count in pairs["ploidy"].value_counts().items():
            metrics.append((f"pairs_ploidy_{ploidy}", int(count)))
        for pair_type, count in pairs["pair_type"].value_counts().items():
            metrics.append((f"pairs_type_{pair_type}", int(count)))
        for dedup_class, count in pairs["dedup_class"].value_counts().items():
            metrics.append((f"pairs_dedup_{dedup_class}", int(count)))
        for tier, count in pairs["tier"].value_counts().items():
            metrics.append((f"pairs_tier_{tier}", int(count)))
        unique_unique = (
            (pairs["status1"] == "Unique") & (pairs["status2"] == "Unique")
        ).sum()
        metrics.append(("pairs_unique_unique", int(unique_unique)))
        for tier in TIER_ORDER:
            tier_pairs = pairs.loc[pairs["tier"] == tier]
            if tier_pairs.empty:
                continue
            plink = pd.to_numeric(
                tier_pairs["plink_ibs"], errors="coerce"
            ).dropna()
            hom = pd.to_numeric(
                tier_pairs["hom_ibs"], errors="coerce"
            ).dropna()
            metrics.extend(
                [
                    (f"tier_{tier}_plink_median", float(plink.median())),
                    (f"tier_{tier}_hom_median", float(hom.median())),
                    (
                        f"tier_{tier}_frac_plink_gt_p95",
                        float((plink > random_stats["random_plink_p95"]).mean()),
                    ),
                    (
                        f"tier_{tier}_frac_plink_gt_p99",
                        float((plink > random_stats["random_plink_p99"]).mean()),
                    ),
                    (
                        f"tier_{tier}_frac_hom_gt_p95",
                        float((hom > random_stats["random_hom_p95"]).mean()),
                    ),
                    (
                        f"tier_{tier}_frac_hom_gt_p99",
                        float((hom > random_stats["random_hom_p99"]).mean()),
                    ),
                ]
            )
    metrics.extend(sorted(random_stats.items()))
    return pd.DataFrame(metrics, columns=["metric", "value"])


def _ibs_threshold_lines(
    value: float | Sequence[float] | None,
) -> list[dict] | None:
    """Vertical reference-line specs for one or more IBS cutoffs, or ``None``."""
    if value is None:
        return None
    return [
        {
            "value": float(cutoff),
            "color": "black",
            "linestyle": "--",
            "label": f"IBS {float(cutoff):g}",
        }
        for cutoff in (value if isinstance(value, (list, tuple)) else [value])
    ]


def _plot_distribution(
    data: pd.DataFrame,
    hue_order: Sequence[str],
    palette_map: dict[str, str],
    hue_labels: dict[str, str],
    title: str,
    stem: str,
    xlabel: str,
    xlim: tuple[float, float],
    *,
    make_kde: bool,
    make_hist: bool,
    bins: int = 80,
    thresholds_v: list[dict] | None = None,
) -> None:
    """Emit the KDE line and/or histogram-bar version of one distribution."""
    if make_kde:
        plot_kde_density_lines(
            data=data,
            x_col="value",
            hue_col="group",
            hue_order=list(hue_order),
            palette_map=palette_map,
            hue_labels=hue_labels,
            title=title,
            filename=f"{stem}.png",
            xlabel=xlabel,
            xlim=xlim,
            bw_method=0.15,
        )
    if make_hist:
        plot_grouped_histogram(
            data=data,
            x_col="value",
            hue_col="group",
            hue_order=list(hue_order),
            palette_map=palette_map,
            hue_labels=hue_labels,
            title=title,
            filename=f"{stem}.hist.png",
            xlabel=xlabel,
            xlim=xlim,
            bins=bins,
            thresholds_v=thresholds_v,
        )


def plot_combine_ibs(
    pairs: pd.DataFrame,
    random_abd: pd.DataFrame,
    output_prefix: str | Path,
    *,
    metric: str = "plink_ibs",
    make_kde: bool = True,
    make_hist: bool = True,
    bins: int = 80,
    ibs_threshold: float | Sequence[float] | None = None,
) -> None:
    """Draw duplicate-vs-random, dedup-class, and tier figures."""
    if metric not in pairs.columns or metric not in random_abd.columns:
        fail(f"IBS metric column not found: {metric}")
    prefix = str(output_prefix)
    thresholds_v = _ibs_threshold_lines(ibs_threshold)
    metric_label = {
        "plink_ibs": "PLINK IBS (AABBDD)",
        "hom_ibs": "hom-SNP IBS (AABBDD)",
    }.get(metric, f"{metric} (AABBDD)")
    abd = pairs.loc[
        (pairs["ploidy1"] == "AABBDD")
        & (pairs["ploidy2"] == "AABBDD")
        & pairs[metric].notna()
    ].copy()
    random_values = pd.to_numeric(
        random_abd[metric], errors="coerce"
    ).dropna()
    duplicate_values = pd.to_numeric(
        abd[metric], errors="coerce"
    ).dropna()

    vs_random = pd.concat(
        [
            pd.DataFrame({"value": random_values, "group": "random"}),
            pd.DataFrame({"value": duplicate_values, "group": "duplicate"}),
        ],
        ignore_index=True,
    )
    _plot_distribution(
        data=vs_random,
        hue_order=RANDOM_HUE_ORDER,
        palette_map=RANDOM_PALETTE,
        hue_labels={
            "random": "random AABBDD pairs",
            "duplicate": "combine-name pairs",
        },
        title="Combine-name pairs vs random AABBDD pairs",
        stem=f"{prefix}.vs_random",
        xlabel=metric_label,
        xlim=(0.94, 1.0),
        make_kde=make_kde,
        make_hist=make_hist,
        bins=bins,
        thresholds_v=thresholds_v,
    )

    dedup_frames = [
        pd.DataFrame({"value": random_values, "group": "random"})
    ]
    for dedup_class in ("same_dedup_group", "different_dedup_group"):
        values = pd.to_numeric(
            abd.loc[abd["dedup_class"] == dedup_class, metric],
            errors="coerce",
        ).dropna()
        dedup_frames.append(
            pd.DataFrame({"value": values, "group": dedup_class})
        )
    _plot_distribution(
        data=pd.concat(dedup_frames, ignore_index=True),
        hue_order=DEDUP_HUE_ORDER,
        palette_map=DEDUP_PALETTE,
        hue_labels={
            "random": "random AABBDD pairs",
            "same_dedup_group": "same dedup group",
            "different_dedup_group": "different dedup group",
        },
        title="Combine-name pairs by dedup agreement",
        stem=f"{prefix}.by_dedup",
        xlabel=metric_label,
        xlim=(0.94, 1.0),
        make_kde=make_kde,
        make_hist=make_hist,
        bins=bins,
        thresholds_v=thresholds_v,
    )

    tier_frames = []
    for tier in TIER_ORDER:
        values = pd.to_numeric(
            abd.loc[abd["tier"] == tier, metric], errors="coerce"
        ).dropna()
        tier_frames.append(pd.DataFrame({"value": values, "group": tier}))
    _plot_distribution(
        data=pd.concat(tier_frames, ignore_index=True),
        hue_order=TIER_ORDER,
        palette_map=TIER_PALETTE,
        hue_labels={
            "low": "low-IBS tier",
            "mid": "mid-IBS tier",
            "high": "high-IBS tier",
        },
        title="Combine-name pair IBS tiers",
        stem=f"{prefix}.tiers",
        xlabel=metric_label,
        xlim=(0.94, 1.0),
        make_kde=make_kde,
        make_hist=make_hist,
        bins=bins,
        thresholds_v=thresholds_v,
    )


def plot_combine_ibs_by_ploidy(
    pairs: pd.DataFrame,
    random_pairs: pd.DataFrame,
    output_prefix: str | Path,
    *,
    metric: str = "plink_ibs",
    min_duplicates: int = 5,
    make_kde: bool = True,
    make_hist: bool = True,
    bins: int = 60,
    ibs_threshold: float | Sequence[float] | None = None,
) -> None:
    """One duplicate-vs-random distribution per ploidy group."""
    if metric not in pairs.columns or metric not in random_pairs.columns:
        fail(f"IBS metric column not found: {metric}")
    prefix = str(output_prefix)
    thresholds_v = _ibs_threshold_lines(ibs_threshold)
    metric_label = {
        "plink_ibs": "PLINK IBS",
        "hom_ibs": "hom-SNP IBS",
    }.get(metric, metric)
    order = ["AABBDD", "AABB", "DD", "AA", "SS"]
    # Groups come from both tables: a ploidy can have random pairs without any
    # combine-name duplicates (e.g. SS), and its histogram still belongs here.
    present = sorted(
        set(pairs.loc[pairs["ploidy1"] == pairs["ploidy2"], "ploidy1"])
        | set(random_pairs["ploidy"].dropna().astype(str))
    )
    groups = [group for group in order if group in present]
    groups += [group for group in present if group not in groups]
    for group in groups:
        duplicate = pd.to_numeric(
            pairs.loc[
                (pairs["ploidy1"] == group) & (pairs["ploidy2"] == group),
                metric,
            ],
            errors="coerce",
        ).dropna()
        random_values = pd.to_numeric(
            random_pairs.loc[random_pairs["ploidy"] == group, metric],
            errors="coerce",
        ).dropna()
        # Bar panels need no smoothing, so they also draw groups with fewer
        # than ``min_duplicates`` pairs (even none); a KDE needs the full count.
        if len(random_values) < 2 or (make_kde and len(duplicate) < min_duplicates):
            logger.info(
                "Skip %s by-ploidy plot: duplicates=%s random=%s",
                group,
                len(duplicate),
                len(random_values),
            )
            continue
        has_duplicates = not duplicate.empty
        if not has_duplicates:
            logger.info(
                "%s by-ploidy plot has no combine-name pairs; "
                "drawing the random reference only",
                group,
            )
        combined = pd.concat([duplicate, random_values], ignore_index=True)
        low = max(0.0, float(combined.quantile(0.001)) - 0.02)
        frame = pd.concat(
            [
                pd.DataFrame(
                    {"value": random_values, "group": "random"}
                ),
                pd.DataFrame(
                    {"value": duplicate, "group": "duplicate"}
                ),
            ],
            ignore_index=True,
        )
        # Never request more bars than the values support: a small group (SS has
        # only a few dozen random pairs) would otherwise break into single-value
        # spikes instead of a distribution.
        n_values = len(random_values) + len(duplicate)
        group_bins = min(bins, max(10, n_values // 5))
        if group_bins < bins:
            logger.info(
                "%s by-ploidy plot: %s bins instead of %s for %s values",
                group,
                group_bins,
                bins,
                n_values,
            )
        _plot_distribution(
            data=frame,
            hue_order=RANDOM_HUE_ORDER if has_duplicates else ["random"],
            palette_map=RANDOM_PALETTE,
            hue_labels={
                "random": "random pairs",
                "duplicate": "combine-name pairs",
            },
            title=(
                f"{group} combine pairs vs random ({metric_label})"
                if has_duplicates
                else f"{group} random pairs ({metric_label}; "
                "no combine-name pairs)"
            ),
            stem=f"{prefix}.by_ploidy.{group}.vs_random",
            xlabel=f"{metric_label} [{group}]",
            xlim=(low, 1.0),
            make_kde=make_kde,
            make_hist=make_hist,
            bins=group_bins,
            thresholds_v=thresholds_v,
        )


def plot_combine_ibs_pca(
    pairs_path: str | Path,
    eigenvec_path: str | Path,
    output_prefix: str | Path,
    *,
    subgenome: str,
    metric: str = "hom_ibs",
    threshold: float = 0.99,
    suffix: str = "pca",
    allow_empty_edges: bool = False,
    edge_categories: Sequence[str] | None = None,
) -> None:
    """Draw one PCA projection with every combine group wired together.

    Points are the samples of one subgenome PCA; each combine-name pair becomes
    an edge, red when it clears the merge cutoff and blue when it does not, with
    both endpoints of every blue edge labelled.

    ``edge_categories`` restricts the drawing to a subset of those categories
    (``"pass"`` / ``"fail"``), which produces a figure that only shows the pairs
    on one side of the cutoff while the whole population stays in the
    background.
    """
    pairs_file = Path(pairs_path)
    eigenvec_file = Path(eigenvec_path)
    if not pairs_file.exists():
        fail(f"Pair table not found: {pairs_file}")
    if not eigenvec_file.exists():
        fail(f"Eigenvector table not found: {eigenvec_file}")
    pairs = pd.read_csv(pairs_file, sep="\t")
    if metric not in pairs.columns:
        fail(f"IBS metric column not found: {metric}")

    coords = pd.read_csv(eigenvec_file, sep=r"\s+")
    coords = coords.rename(columns={"#IID": "id"})
    coords = coords[["id", "PC1", "PC2"]]

    edges = pairs.loc[
        pairs[metric].notna(), ["sample1", "sample2", metric]
    ].rename(columns={"sample1": "source", "sample2": "target"})
    categories = {
        "pass": f"hom-SNP IBS > {threshold:g}",
        "fail": f"hom-SNP IBS ≤ {threshold:g}",
    }
    edges["category"] = np.where(
        edges[metric] > threshold, categories["pass"], categories["fail"]
    )
    above, below = categories["pass"], categories["fail"]
    missing = (
        ~edges["source"].isin(coords["id"]) | ~edges["target"].isin(coords["id"])
    )
    if missing.any():
        logger.info(
            "%s %s: skip %s pairs whose samples are absent from the %s PCA",
            Path(output_prefix).name,
            subgenome,
            int(missing.sum()),
            subgenome,
        )
        edges = edges.loc[~missing]
    if edges.empty and not allow_empty_edges:
        # A cohort can be single-ploidy (IPK/HZAU/Watkins is AABBDD only), so
        # the other population panels would be background dots with no wiring.
        logger.warning(
            "%s: no combine pairs in the %s view; no figure written",
            Path(output_prefix).name,
            subgenome,
        )
        return
    if edges.empty:
        logger.warning(
            "%s: no combine pairs in the %s view; drawing the population only",
            Path(output_prefix).name,
            subgenome,
        )

    if edge_categories is None:
        selected = edges
        drawn_suffix = suffix
        title_suffix = ""
        # Only the combined view labels the failing endpoints; once pass and
        # fail are drawn apart the split is carried by title and colour, and
        # labels would only add clutter.
        label_categories = [below]
    else:
        keep = {categories[name] for name in edge_categories}
        selected = edges.loc[edges["category"].isin(keep)]
        drawn_suffix = f"{suffix}.{'_'.join(edge_categories)}"
        title_suffix = f" — {'/'.join(edge_categories)} only"
        label_categories = []
    if selected.empty:
        logger.warning(
            "%s: no %s edges in the %s view; no figure written",
            Path(output_prefix).name,
            "/".join(edge_categories or ["selected"]),
            subgenome,
        )
        return

    plot_connected_scatter(
        points=coords,
        edges=selected,
        title=(
            f"Combine groups on {subgenome} PCA — hom-SNP IBS cutoff "
            f"{threshold:g}{title_suffix}"
        ),
        filename=f"{output_prefix}.{drawn_suffix}.{subgenome}.png",
        xlabel=f"{subgenome} PC1",
        ylabel=f"{subgenome} PC2",
        edge_colors={above: "#c44e52", below: "#4c72b0"},
        label_categories=label_categories,
        focus=True,
    )
    logger.info(
        "Wrote %s.%s.%s.png (%s pairs, %s below cutoff)",
        output_prefix,
        drawn_suffix,
        subgenome,
        len(selected),
        int((selected["category"] == below).sum()),
    )


def load_final_merge_edges(
    decisions_path: str | Path,
    *,
    allowed_sources: Sequence[str] = ("VMap3+WAP", "VMap3", "WAP"),
) -> pd.DataFrame:
    """Return merged pairs from the final combine map edge list.

    Only edges whose ``sources`` fall in ``allowed_sources`` are kept, so the
    Watkins and IPK/HZAU merges stay out of the release-facing WAP_PQ/VMap3
    view.
    """
    path = Path(decisions_path)
    if not path.exists():
        fail(f"Pair-decision table not found: {path}")
    decisions = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"sample1", "sample2", "sources", "merges"}
    missing = required.difference(decisions.columns)
    if missing:
        fail(f"{path} is missing columns: {sorted(missing)}")
    merged = decisions["merges"].str.lower().isin(("true", "1", "yes"))
    keep = merged & decisions["sources"].isin(set(allowed_sources))
    edges = decisions.loc[keep, ["sample1", "sample2", "hom_ibs", "sources"]].copy()
    edges = edges.drop_duplicates(["sample1", "sample2"])
    logger.info(
        "Final merge edges from %s: %s merged pairs over %s sources",
        path.name,
        len(edges),
        sorted(edges["sources"].unique()) if not edges.empty else [],
    )
    return edges


def load_taxa_merge_map(path: str | Path) -> list[tuple[str, list[str]]]:
    """Parse ``final.<sub>.taxaMergeMap.txt`` into (taxon, source taxa)."""
    merge_file = Path(path)
    if not merge_file.exists():
        fail(f"Taxa merge map not found: {merge_file}")
    records: list[tuple[str, list[str]]] = []
    with merge_file.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            sources = [part for part in parts[2:] if part]
            if sources:
                records.append((parts[0], sources))
    if not records:
        fail(f"No merge records in {merge_file}")
    return records


def plot_merge_coverage_hist(
    merge_map_path: str | Path,
    taxa_map_path: str | Path,
    output_prefix: str | Path,
    *,
    combine_map_path: str | Path | None = None,
    sources_pattern: str | None = None,
    cohort_label: str | None = None,
    include_singletons: bool = False,
    x_max: float = 30.0,
    bins: int = 100,
) -> dict[str, float]:
    """Grouped histogram of member coverage before vs merged coverage after.

    The merges come from ``final.<sub>.taxaMergeMap.txt``: ``after`` is the sum
    of each line's source coverages from the taxaBamMap, exactly the value
    TIGER uses when the merge map says ``coverage = .``. By default only
    multi-member lines are used, because singleton final taxa are unchanged and
    would otherwise dominate both series; set ``include_singletons`` to keep
    every line. ``sources_pattern`` (e.g. ``WAP|VMap3``) restricts the lines to
    one cohort, using ``combine_map_path`` only to label which taxa belong to
    that cohort. Writes a KDE curve figure (``{prefix}.png``) and a grouped bar
    histogram (``{prefix}.hist.png``).
    """
    merge_file = Path(merge_map_path)
    taxa_file = Path(taxa_map_path)
    for path in (merge_file, taxa_file):
        if not path.exists():
            fail(f"Input not found: {path}")
    records = load_taxa_merge_map(merge_file)
    taxa = pd.read_csv(taxa_file, sep="\t", dtype=str, keep_default_na=False)
    coverage = pd.Series(
        pd.to_numeric(taxa["Coverage-Of-All-Bams"], errors="coerce").to_numpy(),
        index=taxa["Taxa"],
    )

    if sources_pattern:
        if combine_map_path is None:
            fail("sources_pattern requires combine_map_path for cohort labels")
        combine_file = Path(combine_map_path)
        if not combine_file.exists():
            fail(f"Combine map not found: {combine_file}")
        combine = pd.read_csv(
            combine_file, sep="\t", dtype=str, keep_default_na=False
        )
        selected = combine.loc[
            combine["sources"].str.contains(sources_pattern, regex=True)
        ]
        allowed = {
            member
            for cell in selected["members"]
            for member in str(cell).split("|")
            if member
        }
        records = [
            (taxon, sources)
            for taxon, sources in records
            if all(source in allowed for source in sources)
        ]
    if not include_singletons:
        records = [
            (taxon, sources) for taxon, sources in records if len(sources) > 1
        ]

    before: list[float] = []
    after: list[float] = []
    for _, sources in records:
        values = [
            float(coverage[source])
            for source in sources
            if source in coverage.index and pd.notna(coverage[source])
        ]
        if not values:
            continue
        before.extend(values)
        after.append(float(np.sum(values)))
    if not before or not after:
        fail("No merge coverage values were found for the selected rows")

    before_series = pd.Series(before, dtype=float)
    after_series = pd.Series(after, dtype=float)
    stats = {
        "rows": float(len(after_series)),
        "member_samples": float(len(before_series)),
        "before_median": float(before_series.median()),
        "after_median": float(after_series.median()),
        "before_q25": float(before_series.quantile(0.25)),
        "before_q75": float(before_series.quantile(0.75)),
        "after_q25": float(after_series.quantile(0.25)),
        "after_q75": float(after_series.quantile(0.75)),
    }
    plot_df = pd.concat(
        [
            pd.DataFrame({"coverage": before_series, "stage": "before merge"}),
            pd.DataFrame({"coverage": after_series, "stage": "after merge"}),
        ],
        ignore_index=True,
    )
    plot_grouped_histogram(
        data=plot_df,
        x_col="coverage",
        hue_col="stage",
        hue_order=["before merge", "after merge"],
        palette_map={
            "before merge": "#2a78d6",
            "after merge": "#eb6834",
        },
        hue_labels={
            "before merge": "before merge (member samples)",
            "after merge": "after merge (final taxa)",
        },
        title=(
            "Coverage before vs after taxa merge"
            f" — {cohort_label or 'all sources'}\n"
            f"{'all rows (singletons unchanged)' if include_singletons else 'multi-member groups only'}: "
            f"{int(stats['rows']):,} final taxa, "
            f"{int(stats['member_samples']):,} member samples"
        ),
        filename=f"{output_prefix}.hist.png",
        xlabel="Coverage-Of-All-Bams (×)",
        xlim=(0.0, float(x_max)),
        bins=int(bins),
    )
    plot_kde_density_lines(
        data=plot_df,
        x_col="coverage",
        hue_col="stage",
        hue_order=["before merge", "after merge"],
        palette_map={
            "before merge": "#2a78d6",
            "after merge": "#eb6834",
        },
        hue_labels={
            "before merge": "before merge (member samples)",
            "after merge": "after merge (final taxa)",
        },
        title=(
            "Coverage before vs after taxa merge"
            f" — {cohort_label or 'all sources'}\n"
            f"{'all rows (singletons unchanged)' if include_singletons else 'multi-member groups only'}: "
            f"{int(stats['rows']):,} final taxa, "
            f"{int(stats['member_samples']):,} member samples"
        ),
        filename=f"{output_prefix}.png",
        xlabel="Coverage-Of-All-Bams (×)",
        xlim=(0.0, float(x_max)),
        bw_method=0.15,
    )
    logger.info(
        "Wrote %s.png and %s.hist.png: before median %.2fx (n=%s), "
        "after median %.2fx (n=%s)",
        output_prefix,
        output_prefix,
        stats["before_median"],
        int(stats["member_samples"]),
        stats["after_median"],
        int(stats["rows"]),
    )
    return stats


def _robust_limits(values, pad_fraction: float = 0.08) -> tuple[float, float]:
    """1st-99th percentile range plus padding, for a zoomed wiring view."""
    low = float(np.nanquantile(values, 0.01))
    high = float(np.nanquantile(values, 0.99))
    span = high - low
    if span <= 0:
        span = max(abs(high), 1.0)
    pad = span * pad_fraction
    return low - pad, high + pad


def plot_final_merge_wiring(
    decisions_path: str | Path,
    coords_path: str | Path,
    output_prefix: str | Path,
    *,
    axis_label: str,
    view: str,
    edge_color: str = "#c0392b",
) -> int:
    """Draw final WAP_PQ/VMap3 merge edges on one PCA/PCoA coordinate table.

    Returns the number of edges drawn; zero means no figure was written.
    """
    edges = load_final_merge_edges(decisions_path)
    coords_file = Path(coords_path)
    if not coords_file.exists():
        fail(f"Coordinate table not found: {coords_file}")
    coords = pd.read_csv(coords_file, sep=r"\s+")
    coords = coords.rename(columns={"#IID": "id"})
    if "id" not in coords.columns or not {"PC1", "PC2"}.issubset(coords.columns):
        fail(f"{coords_file} needs #IID, PC1, PC2")
    coords = coords[["id", "PC1", "PC2"]].dropna()

    edges = edges.rename(columns={"sample1": "source", "sample2": "target"})
    drawable = (
        edges["source"].isin(coords["id"])
        & edges["target"].isin(coords["id"])
    )
    n_drawn = int(drawable.sum())
    n_missing = int(len(edges) - n_drawn)
    if n_drawn == 0:
        logger.warning(
            "%s: no final merge edges land in %s; no figure written",
            Path(output_prefix).name,
            view,
        )
        return 0
    edges = edges.loc[drawable].copy()
    edges["category"] = "final merge"
    involved = coords.loc[
        coords["id"].isin(set(edges["source"]) | set(edges["target"]))
    ]
    plot_connected_scatter(
        points=coords,
        edges=edges,
        title=(
            f"WAP_PQ/VMap3 final merges on {view}\n"
            f"{n_drawn} merged pairs (hom-SNP IBS > 0.999)"
        ),
        filename=f"{output_prefix}.png",
        xlabel=f"{axis_label} PC1",
        ylabel=f"{axis_label} PC2",
        edge_colors={"final merge": edge_color},
        label_categories=[],
        xlim=_robust_limits(involved["PC1"]),
        ylim=_robust_limits(involved["PC2"]),
    )
    logger.info(
        "Wrote %s.png on %s: %s edges drawn, %s edges dropped "
        "(endpoints absent from the coordinates)",
        output_prefix,
        view,
        n_drawn,
        n_missing,
    )
    return n_drawn


def plot_combine_ibs_edge_distance(
    pairs_path: str | Path,
    eigenvec_path: str | Path,
    output_prefix: str | Path,
    *,
    subgenome: str,
    metric: str = "hom_ibs",
    threshold: float = 0.99,
    bins: int = 60,
    suffix: str = "pcoa_edge_distance",
) -> None:
    """Grouped histogram of PCoA edge lengths, split by the merge cutoff.

    Every combine pair becomes one Euclidean distance over PC1..PC10 in the
    given coordinate table; bars are drawn side by side for pairs above and at
    or below the cutoff, so the two distributions can be compared directly.
    """
    pairs_file = Path(pairs_path)
    eigenvec_file = Path(eigenvec_path)
    if not pairs_file.exists():
        fail(f"Pair table not found: {pairs_file}")
    if not eigenvec_file.exists():
        fail(f"Coordinate table not found: {eigenvec_file}")
    pairs = pd.read_csv(pairs_file, sep="\t")
    if metric not in pairs.columns:
        fail(f"IBS metric column not found: {metric}")
    coords = pd.read_csv(eigenvec_file, sep=r"\s+")
    coords = coords.rename(columns={"#IID": "id"}).set_index("id")
    pcs = [f"PC{i}" for i in range(1, 11) if f"PC{i}" in coords.columns]
    coords = coords[pcs]

    keep = pairs.loc[
        pairs["sample1"].isin(coords.index) & pairs["sample2"].isin(coords.index)
    ]
    first = coords.loc[keep["sample1"], pcs].to_numpy(dtype=float)
    second = coords.loc[keep["sample2"], pcs].to_numpy(dtype=float)
    distance = np.sqrt(((first - second) ** 2).sum(axis=1))

    above = f"hom-SNP IBS > {threshold:g}"
    below = f"hom-SNP IBS ≤ {threshold:g}"
    data = pd.DataFrame(
        {
            "value": distance,
            "group": np.where(keep[metric].to_numpy(float) > threshold, above, below),
        }
    )
    if data.empty:
        logger.warning(
            "%s: no combine pairs land in the %s coordinates; no figure written",
            Path(output_prefix).name,
            subgenome,
        )
        return
    upper = float(np.quantile(data["value"], 0.995))
    plot_grouped_histogram(
        data=data,
        x_col="value",
        hue_col="group",
        hue_order=[above, below],
        palette_map={above: "#c44e52", below: "#4c72b0"},
        hue_labels={above: "pass", below: "fail"},
        title=(
            f"{subgenome} combine pairs — PCoA edge distance "
            f"(cutoff hom-SNP IBS {threshold:g})"
        ),
        filename=f"{output_prefix}.{suffix}.{subgenome}.hist.png",
        xlabel="PCoA distance over PC1–PC10",
        xlim=(0.0, upper),
        bins=bins,
    )
    logger.info(
        "Wrote %s.%s.%s.hist.png (%s pairs, %s below cutoff)",
        output_prefix,
        suffix,
        subgenome,
        len(data),
        int((data["group"] == below).sum()),
    )


def load_taxa_coverage(path: str | Path) -> dict[str, float]:
    """Map taxa IDs to mean coverage from a taxaBamMap table."""
    coverage_path = Path(path)
    if not coverage_path.exists():
        fail(f"taxaBamMap not found: {coverage_path}")
    frame = pd.read_csv(coverage_path, sep="\t", dtype=str, keep_default_na=False)
    required = {"Taxa", "Coverage-Of-All-Bams"}
    missing = required.difference(frame.columns)
    if missing:
        fail(f"{coverage_path} is missing columns: {sorted(missing)}")
    values = pd.to_numeric(frame["Coverage-Of-All-Bams"], errors="coerce")
    return {
        str(taxa): float(value)
        for taxa, value in zip(frame["Taxa"], values)
        if pd.notna(value)
    }


def load_sample_missing_rates(
    paths: Sequence[str | Path],
) -> dict[str, float]:
    """Mean per-sample F_MISS across the given subgenome tables."""
    frames = []
    for path in paths:
        table = Path(path)
        if not table.exists():
            fail(f"Missing-rate table not found: {table}")
        frame = pd.read_csv(table, sep="\t")
        required = {"Sample", "F_MISS"}
        missing = required.difference(frame.columns)
        if missing:
            fail(f"{table} is missing columns: {sorted(missing)}")
        values = pd.to_numeric(frame["F_MISS"], errors="coerce")
        series = pd.Series(values.to_numpy(), index=frame["Sample"].astype(str))
        frames.append(series.groupby(level=0).mean())
    if not frames:
        fail("No missing-rate tables were provided")
    merged = pd.concat(frames, axis=1)
    return merged.mean(axis=1).to_dict()


def _load_highlight_pairs(
    path: str | Path,
    metric: str,
    min_ibs: float,
    pair_type: str | None = None,
) -> set[tuple[str, str]]:
    """Return unordered sample pairs whose ``metric`` exceeds ``min_ibs``.

    ``pair_type`` filters the source table the same way the plotted table is
    filtered, so a highlight selection (for example test_common_inter pairs
    above 0.999) only marks pairs that belong to the plotted cohort.
    """
    frame = pd.read_csv(path, sep="\t")
    required = {"sample1", "sample2", metric}
    missing = required.difference(frame.columns)
    if missing:
        fail(f"{path} is missing columns: {sorted(missing)}")
    if pair_type and pair_type != "all" and "pair_type" in frame.columns:
        frame = frame.loc[frame["pair_type"] == pair_type]
    values = pd.to_numeric(frame[metric], errors="coerce")
    frame = frame.loc[values > float(min_ibs)]
    return {
        tuple(sorted((str(first), str(second))))
        for first, second in zip(frame["sample1"], frame["sample2"])
    }


def plot_combine_ibs_low_coverage(
    pairs_path: str | Path,
    output_prefix: str | Path,
    *,
    pair_type: str | None = "within_watkins",
    metric: str = "hom_ibs",
    cutoff: float = 0.99,
    coverage_floor: float = 3.0,
    bins: Sequence[float] | None = None,
    coverage_map_path: str | Path | None = None,
    additional_cutoffs: Sequence[float] = (0.999,),
    highlight_union_n: int | None = None,
    highlight_pairs_path: str | Path | None = None,
    highlight_min_ibs: float | None = None,
    highlight_label: str = "selected pairs",
    missing_rate_paths: Sequence[str | Path] | None = None,
    label: str | None = None,
    y_min: float | None = None,
) -> pd.DataFrame:
    """Scatter lower-coverage depth against hom-SNP IBS for trusted pairs.

    ``pair_type`` selects the trusted combine pairs (default: both names from
    the Watkins cohort); pass ``None`` or ``"all"`` to keep every pair in the
    table. ``label`` overrides the title/legend cohort name. For each pair
    ``low_coverage`` is the smaller of
    ``coverage1`` and ``coverage2``; pairs missing either depth are dropped
    because they cannot be placed on the x axis.

    The pair table takes sample depth from the dedup snapshot, which does not
    cover every genotyped taxon. Pass ``coverage_map_path`` (a taxaBamMap) to
    fill depths missing from that snapshot by sample ID before taking the
    minimum.

    When ``highlight_union_n`` is set, pairs whose combine group has that many
    members (for example the two 48-member Watkins groups) are drawn in a
    separate colour and legend entry.

    ``highlight_pairs_path`` marks pairs whose ``metric`` exceeds
    ``highlight_min_ibs`` in another pair table (for example the
    test_common_inter table at 0.999) in a third colour, so the same pairs can
    be followed across test sets. Selection is applied to the highlight table,
    not to the plotted one, and it takes precedence over ``highlight_union_n``.

    Depth is the primary x axis: a log figure (``{prefix}.png``) and a linear
    one (``{prefix}.linear.png``). With ``missing_rate_paths`` (the per-subgenome
    ``*.sample.miss.rate.info.tsv`` tables, averaged per sample) two more x
    axes are added: the missing rate of the higher-coverage sample and the sum
    of both samples' missing rates, each in log and linear form. Every figure
    carries the merge cutoff, the ``additional_cutoffs`` review lines (default
    0.999), and - on the depth axis only - the coverage floor. The
    missing-rate-sum axis also carries reference lines at x = 1 and y = 0.997.
    ``bins``
    overrides the default widening depth edges. ``y_min`` fixes the lower y
    limit (default: the data minimum minus 0.001); the upper limit is always
    1.0005. The binned summary is still returned for logging or a companion
    table.
    """
    pairs_file = Path(pairs_path)
    if not pairs_file.exists():
        fail(f"Pair table not found: {pairs_file}")
    pairs = pd.read_csv(pairs_file, sep="\t")
    required = {"sample1", "sample2", "pair_type", "coverage1", "coverage2", metric}
    missing = required.difference(pairs.columns)
    if missing:
        fail(f"{pairs_file} is missing columns: {sorted(missing)}")

    if pair_type and pair_type != "all":
        selected = pairs.loc[pairs["pair_type"] == pair_type].copy()
    else:
        selected = pairs.copy()
    if selected.empty:
        fail(f"No {pair_type or 'usable'} pairs in {pairs_file}")
    for column in ("coverage1", "coverage2", metric):
        selected[column] = pd.to_numeric(selected[column], errors="coerce")
    n_filled = 0
    if coverage_map_path is not None:
        coverage = load_taxa_coverage(coverage_map_path)
        for coverage_col, sample_col in (
            ("coverage1", "sample1"),
            ("coverage2", "sample2"),
        ):
            absent = selected[coverage_col].isna()
            filled = selected.loc[absent, sample_col].map(coverage)
            selected.loc[absent, coverage_col] = filled
            n_filled += int(filled.notna().sum())
        logger.info(
            "Filled %s sample depths from %s", n_filled, coverage_map_path
        )
    selected["low_coverage"] = np.minimum(
        selected["coverage1"], selected["coverage2"]
    )
    if missing_rate_paths:
        missing_rates = load_sample_missing_rates(missing_rate_paths)
        selected["missing_rate1"] = selected["sample1"].map(missing_rates)
        selected["missing_rate2"] = selected["sample2"].map(missing_rates)
        selected["high_missing_rate"] = selected["missing_rate1"].where(
            selected["coverage1"] >= selected["coverage2"],
            selected["missing_rate2"],
        )
        selected["missing_rate_sum"] = (
            selected["missing_rate1"] + selected["missing_rate2"]
        )
        logger.info(
            "Missing rates loaded for %s/%s pair samples",
            int(selected["missing_rate1"].notna().sum()),
            int(selected["missing_rate2"].notna().sum()),
        )
    valid = selected.dropna(subset=["low_coverage", metric])
    n_dropped = int(len(selected) - len(valid))
    if valid.empty:
        fail(f"No {pair_type} pairs have both coverage and {metric}")

    group_col = None
    hue_palette = None
    n_highlight = 0
    pair_set = pd.Series("other pairs", index=valid.index)
    union_label = None
    if highlight_union_n is not None and "union_n" in valid.columns:
        union_n = pd.to_numeric(valid["union_n"], errors="coerce")
        highlighted = union_n == int(highlight_union_n)
        if highlighted.any():
            union_label = f"{int(highlight_union_n)}-member group"
            pair_set[highlighted] = union_label

    threshold_matched = 0
    if highlight_pairs_path is not None and highlight_min_ibs is not None:
        selected = _load_highlight_pairs(
            highlight_pairs_path, metric, float(highlight_min_ibs), pair_type
        )
        if selected:
            keys = pd.Series(
                [
                    tuple(sorted((str(first), str(second))))
                    for first, second in zip(valid["sample1"], valid["sample2"])
                ],
                index=valid.index,
            )
            threshold_mask = keys.isin(selected)
            threshold_matched = int(threshold_mask.sum())
            pair_set[threshold_mask] = highlight_label
        logger.info(
            "Highlighted %s pairs with %s > %s from %s",
            threshold_matched,
            metric,
            highlight_min_ibs,
            highlight_pairs_path,
        )

    if (pair_set != "other pairs").any():
        valid = valid.copy()
        valid["Pair set"] = pair_set.to_numpy()
        group_col = "Pair set"
        hue_palette = {"other pairs": LOW_COVERAGE_POINT_COLOR}
        if union_label is not None:
            hue_palette[union_label] = LOW_COVERAGE_HIGHLIGHT_COLOR
        if (pair_set == highlight_label).any():
            hue_palette[highlight_label] = HIGHLIGHT_PAIR_COLOR
        n_highlight = int((pair_set != "other pairs").sum())

    n_pass = int((valid[metric] > cutoff).sum())
    edges = tuple(bins) if bins is not None else LOW_COVERAGE_BINS
    summary = _low_coverage_bins(valid, metric, cutoff, edges)

    if not label:
        label = (
            pair_type.removeprefix("within_").capitalize()
            if pair_type and pair_type != "all"
            else "Combine"
        )
    metric_label = "hom-SNP IBS" if metric == "hom_ibs" else metric
    y_low = (
        float(y_min)
        if y_min is not None
        else max(0.0, float(valid[metric].min()) - 0.001)
    )
    threshold_values = [float(cutoff)]
    for value in additional_cutoffs:
        if float(value) not in threshold_values:
            threshold_values.append(float(value))
    threshold_lines = [
        {
            "value": value,
            "color": "#c0392b" if value == float(cutoff) else "#8e44ad",
            "linestyle": "--" if value == float(cutoff) else "-.",
            "label": (
                f"merge cutoff {value:g}"
                if value == float(cutoff)
                else f"review line {value:g}"
            ),
        }
        for value in threshold_values
    ]

    x_specs = [
        {
            "col": "low_coverage",
            "label": "Lower-coverage sample depth (×)",
            "stem": "",
            "log_xlim": (float(edges[0]), float(edges[-1])),
            "linear_xlim": (0.0, float(valid["low_coverage"].max()) * 1.05),
            "coverage_floor": True,
            "below_note": True,
        }
    ]
    if missing_rate_paths:
        for column, xlabel, stem in (
            (
                "high_missing_rate",
                "Missing rate of the higher-coverage sample",
                ".high_missing_rate",
            ),
            (
                "missing_rate_sum",
                "Sum of both samples' missing rates",
                ".missing_sum",
            ),
        ):
            values = pd.to_numeric(valid[column], errors="coerce").dropna()
            if values.empty:
                logger.warning(
                    "No %s values; %s figures skipped", column, stem
                )
                continue
            spec = {
                "col": column,
                "label": xlabel,
                "stem": stem,
                "log_xlim": (
                    max(float(values.min()) * 0.9, 1e-6),
                    float(values.max()) * 1.1,
                ),
                "linear_xlim": (0.0, float(values.max()) * 1.05),
                "coverage_floor": False,
                "below_note": False,
            }
            if column == "missing_rate_sum":
                spec["extra_h"] = [
                    {
                        "value": 0.997,
                        "color": "#2e7d32",
                        "linestyle": "--",
                        "label": "IBS 0.997",
                    }
                ]
                spec["extra_v"] = [
                    {
                        "value": 1.0,
                        "color": "#1f4e79",
                        "linestyle": ":",
                        "label": "missing-rate sum 1",
                    }
                ]
            x_specs.append(spec)

    written: list[str] = []
    for spec in x_specs:
        x_col = spec["col"]
        pair_values = (
            valid[[x_col, metric]]
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
        )
        if len(pair_values) < 2:
            logger.warning("Not enough %s values; figure skipped", x_col)
            continue
        rho = float(pair_values.corr(method="spearman").iloc[0, 1])
        title = (
            f"{label} combine pairs — {metric_label} vs {spec['label']}\n"
            f"n = {len(pair_values):,}; Spearman ρ = {rho:.2f}"
        )
        if spec["below_note"]:
            n_below = int((pair_values[x_col] < coverage_floor).sum())
            title += f"; {n_below:,} pairs below {coverage_floor:g}×"
        thresholds_v = []
        if spec["coverage_floor"]:
            thresholds_v.append(
                {
                    "value": float(coverage_floor),
                    "color": "#7f7f7f",
                    "linestyle": ":",
                    "label": f"{coverage_floor:g}× callability floor",
                }
            )
        thresholds_v.extend(spec.get("extra_v", []))
        thresholds_h = threshold_lines + spec.get("extra_h", [])
        for axis_suffix, log_x, xlim in (
            ("", True, spec["log_xlim"]),
            (".linear", False, spec["linear_xlim"]),
        ):
            filename = f"{output_prefix}{spec['stem']}{axis_suffix}.png"
            plot_scatter_with_thresholds(
                data=valid,
                x_col=x_col,
                y_col=metric,
                title=title + ("" if log_x else " (linear x axis)"),
                filename=filename,
                thresholds_h=thresholds_h,
                thresholds_v=thresholds_v or None,
                xlabel=spec["label"],
                ylabel=metric_label,
                color=LOW_COVERAGE_POINT_COLOR,
                group_col=group_col,
                hue_palette=hue_palette,
                alpha=0.35,
                s=14,
                figure_size=(10, 6.5),
                xlim=xlim,
                ylim=(y_low, 1.0005),
                log_x=log_x,
            )
            written.append(Path(filename).name)
    logger.info(
        "Wrote %s: %s %s pairs, %s without coverage dropped, "
        "%s highlighted, %s pass %s",
        ", ".join(written),
        len(valid),
        pair_type,
        n_dropped,
        n_highlight,
        n_pass,
        cutoff,
    )
    logger.info(
        "Low-coverage bins:\n%s",
        summary.to_string(index=False),
    )
    return summary


def _low_coverage_bins(
    frame: pd.DataFrame,
    metric: str,
    cutoff: float,
    edges: Sequence[float],
) -> pd.DataFrame:
    """Median IBS and pass share per fixed depth bin."""
    cut = pd.cut(frame["low_coverage"], bins=list(edges), right=False)
    rows: list[dict] = []
    for interval, sub in frame.groupby(cut, observed=True):
        if sub.empty:
            continue
        rows.append(
            {
                "coverage_low": float(interval.left),
                "coverage_high": float(interval.right),
                "n_pairs": int(len(sub)),
                "median_low_coverage": float(sub["low_coverage"].median()),
                f"median_{metric}": float(sub[metric].median()),
                f"q25_{metric}": float(sub[metric].quantile(0.25)),
                f"q75_{metric}": float(sub[metric].quantile(0.75)),
                "frac_above_cutoff": float((sub[metric] > cutoff).mean()),
            }
        )
    return pd.DataFrame(rows)


def replot_combine_ibs_hist(
    pairs_path: str | Path,
    random_pairs_path: str | Path,
    output_prefix: str | Path,
    *,
    metric: str = "plink_ibs",
    bins: int | None = None,
    ibs_threshold: float | Sequence[float] | None = None,
) -> None:
    """Add histogram-bar distribution plots from published pair tables.

    ``bins`` overrides the per-figure bar counts (80 for the AABBDD panels,
    60 for the by-ploidy panels) when a finer or coarser histogram is wanted;
    ``ibs_threshold`` draws one or more IBS values as dashed reference lines.
    """
    pairs_file = Path(pairs_path)
    random_file = Path(random_pairs_path)
    if not pairs_file.exists():
        fail(f"Pair table not found: {pairs_file}")
    if not random_file.exists():
        fail(f"Random-pair table not found: {random_file}")
    pairs = pd.read_csv(pairs_file, sep="\t")
    random_pairs = pd.read_csv(random_file, sep="\t")
    random_abd = random_pairs.loc[
        random_pairs["ploidy"] == "AABBDD"
    ].copy()
    bin_kwargs = {} if bins is None else {"bins": int(bins)}
    threshold_kwargs = (
        {} if ibs_threshold is None else {"ibs_threshold": ibs_threshold}
    )
    plot_combine_ibs(
        pairs,
        random_abd,
        output_prefix,
        metric=metric,
        make_kde=False,
        make_hist=True,
        **bin_kwargs,
        **threshold_kwargs,
    )
    plot_combine_ibs_by_ploidy(
        pairs,
        random_pairs,
        output_prefix,
        metric=metric,
        make_kde=False,
        make_hist=True,
        **bin_kwargs,
        **threshold_kwargs,
    )


def run_combine_ibs_validate(
    combine_path: str | Path,
    excluded_paths: list[str | Path] | None,
    dedup_path: str | Path,
    samples_path: str | Path,
    ibs_dir: str | Path,
    hom_dir: str | Path,
    random_pairs_path: str | Path,
    output_prefix: str | Path = "all.sample.combine_ibs",
    name_columns: Sequence[str] | None = None,
    cohort_labels: Sequence[str] | None = None,
    ibs_metric: str = "plink_ibs",
    tier_low_max: float = TIER_LOW_MAX,
    tier_high_min: float = TIER_HIGH_MIN,
    coverage_map_path: str | Path | None = None,
) -> pd.DataFrame:
    """Run the full combine-name IBS validation and write its artefacts."""
    prefix = Path(output_prefix)
    meta = load_genotype_meta(samples_path)
    genotype_names = set(meta["Sample"])
    excluded_names = load_excluded_controls(excluded_paths)
    coverage_map = (
        load_taxa_coverage(coverage_map_path)
        if coverage_map_path is not None
        else None
    )
    logger.info(
        "Loaded %s genotype samples and %s excluded control names",
        len(genotype_names),
        len(excluded_names),
    )

    pairs, stats = build_combine_pairs(
        combine_path,
        genotype_names,
        excluded_names,
        name_columns=name_columns,
        cohort_labels=cohort_labels,
    )
    logger.info(
        "Combine rows=%s, genotype pairs=%s, excluded occurrences=%s",
        stats["germplasm_rows"],
        len(pairs),
        stats["excluded_name_occurrences"],
    )
    if pairs.empty:
        fail("No genotype-available combine-name pairs were built")

    pairs = add_ploidy_fields(pairs, meta)
    pairs = add_dedup_fields(pairs, dedup_path, coverage_map)
    pairs = add_ibs_fields(pairs, ibs_dir, hom_dir)
    pairs = assign_ibs_tiers(
        pairs,
        metric=ibs_metric,
        low_max=tier_low_max,
        high_min=tier_high_min,
    )
    (
        random_abd,
        random_stats,
        random_all,
        random_by_group,
    ) = _random_reference(random_pairs_path)
    pairs = assign_random_classes(pairs, random_by_group)
    summary = summarize_pairs(
        pairs,
        stats,
        random_stats,
        combine_path=combine_path,
        excluded_names=excluded_names,
        ibs_dir=ibs_dir,
        hom_dir=hom_dir,
        name_columns=name_columns,
        cohort_labels=cohort_labels,
        metric=ibs_metric,
        low_max=tier_low_max,
        high_min=tier_high_min,
    )

    by_ploidy = summarize_by_ploidy(pairs, random_by_group)

    save_df_to_tsv(pairs, f"{prefix}.pairs.tsv", float_format="%.6g")
    save_df_to_tsv(summary, f"{prefix}.summary.tsv", float_format="%.6g")
    save_df_to_tsv(
        by_ploidy, f"{prefix}.by_ploidy.tsv", float_format="%.6g"
    )
    plot_combine_ibs(pairs, random_abd, prefix, metric=ibs_metric)
    plot_combine_ibs_by_ploidy(
        pairs, random_all, prefix, metric=ibs_metric
    )
    logger.info(
        "Wrote %s.pairs.tsv, %s.summary.tsv, %s.by_ploidy.tsv",
        prefix,
        prefix,
        prefix,
    )
    return summary


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
