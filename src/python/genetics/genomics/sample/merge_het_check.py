"""Heterozygosity check for merged germplasms (FastCall3 scan3 output).

scan3 sums the raw allele depths of every source taxon of a germplasm, so a
merged sample can look heterozygous where its sources disagree. This module
compares the per-sample heterozygosity of merged and unmerged germplasms from a
PLINK2 ``--het`` table and controls for call rate, the obvious confounder:
merged samples sum several sources and therefore genotype more sites, while
low-coverage samples carry error-driven heterozygotes.

Heterozygosity rate is ``(OBS_CT - O(HOM)) / OBS_CT``, i.e. the fraction of
non-missing genotypes that are heterozygous; it is already normalised for call
rate.

Figures are separate files (one figure per question) and go through
:mod:`infra.utils.graph`; tables go through :mod:`infra.utils.io`.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from infra.utils import graph as infra_graph
from infra.utils import io as infra_io

LOG = logging.getLogger(__name__)

HET_COLUMNS = ("IID", "O(HOM)", "E(HOM)", "OBS_CT", "F")
MERGED_LABEL = "merged"
UNMERGED_LABEL = "unmerged"
GROUP_ORDER = (UNMERGED_LABEL, MERGED_LABEL)
GROUP_COLORS = {
    UNMERGED_LABEL: "steelblue",
    MERGED_LABEL: "darkorange",
    "unmerged_D_": "seagreen",
}
HIGHLIGHT_COLOR = "seagreen"


def load_plink_het(het_path: str | Path) -> pd.DataFrame:
    """Read a PLINK2 ``.het`` table (``#IID O(HOM) E(HOM) OBS_CT F``)."""
    frame = infra_io.load_df_from_space_sep(str(het_path))
    frame = frame.rename(columns={"#IID": "IID"})
    missing = [c for c in HET_COLUMNS if c not in frame.columns]
    if missing:
        raise ValueError(f"{het_path}: missing columns {missing}")
    frame["IID"] = frame["IID"].astype(str)
    for column in HET_COLUMNS[1:]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    if frame["OBS_CT"].isna().any():
        raise ValueError(f"{het_path}: non-numeric OBS_CT")
    if (frame["OBS_CT"] <= 0).any():
        raise ValueError(f"{het_path}: samples with zero non-missing genotypes")
    LOG.info("loaded %d samples from %s", len(frame), het_path)
    return frame


def load_merge_map(merge_map_path: str | Path) -> pd.DataFrame:
    """Read a taxaMergeMap (``merged<TAB>coverage<TAB>source...``).

    The number of source columns varies per line, so this is parsed line by line
    rather than with a fixed-column reader; ``#`` lines are headers and ``.`` in
    the coverage column is not a source.
    """
    samples: list[str] = []
    n_sources: list[int] = []
    members: list[str] = []
    with Path(merge_map_path).open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            sources = [f.strip() for f in fields[1:] if f.strip() and f.strip() != "."]
            samples.append(fields[0].strip())
            n_sources.append(len(sources))
            members.append("|".join(sources))
    if not samples:
        raise ValueError(f"{merge_map_path}: no merge rows")
    frame = pd.DataFrame({"IID": samples, "n_sources": n_sources, "members": members})
    if frame["IID"].duplicated().any():
        raise ValueError(f"{merge_map_path}: duplicated merged taxon names")
    LOG.info(
        "loaded %d merge rows from %s (%d merged, %d single-source)",
        len(frame),
        merge_map_path,
        int((frame["n_sources"] > 1).sum()),
        int((frame["n_sources"] == 1).sum()),
    )
    return frame


def add_het_rate(frame: pd.DataFrame) -> pd.DataFrame:
    """Add ``n_het`` and ``het_rate`` (heterozygotes / non-missing genotypes)."""
    frame = frame.copy()
    frame["n_het"] = frame["OBS_CT"] - frame["O(HOM)"]
    frame["het_rate"] = frame["n_het"] / frame["OBS_CT"]
    return frame


def assign_groups(frame: pd.DataFrame, highlight_prefix: str | None = None) -> pd.DataFrame:
    """Label samples merged / unmerged, optionally splitting one ID prefix out.

    ``highlight_prefix`` (e.g. ``D_``) moves those unmerged samples into their own
    series, which is how a cohort-specific cluster inside the unmerged
    distribution is shown without hiding the rest of the group.
    """
    frame = frame.copy()
    group = np.where(frame["n_sources"] > 1, MERGED_LABEL, UNMERGED_LABEL)
    if highlight_prefix:
        label = f"unmerged_{highlight_prefix}"
        group = np.where(
            (frame["n_sources"] <= 1)
            & frame["IID"].str.startswith(highlight_prefix),
            label,
            group,
        )
        GROUP_COLORS.setdefault(label, HIGHLIGHT_COLOR)
    frame["group"] = group
    return frame


def group_order(frame: pd.DataFrame) -> list[str]:
    """Series order for plots: unmerged (and its split-off cohorts) then merged."""
    preferred = [UNMERGED_LABEL, *[g for g in GROUP_COLORS if g.startswith("unmerged_")], MERGED_LABEL]
    present = set(frame["group"].unique())
    ordered = [g for g in preferred if g in present]
    ordered += sorted(present - set(ordered))
    return ordered


def build_het_table(
    het_path: str | Path,
    merge_map_path: str | Path,
    highlight_prefix: str | None = None,
) -> pd.DataFrame:
    """Join the ``--het`` table with the merge map and derive the het rate."""
    het = add_het_rate(load_plink_het(het_path))
    merge = load_merge_map(merge_map_path)
    frame = het.merge(merge, on="IID", how="left")
    unmatched = int(frame["n_sources"].isna().sum())
    if unmatched:
        raise ValueError(f"{unmatched} VCF samples are absent from {merge_map_path}")
    frame["n_sources"] = frame["n_sources"].astype(int)
    frame = assign_groups(frame, highlight_prefix)
    return frame.sort_values("het_rate", ascending=False).reset_index(drop=True)


def summarize_by_group(frame: pd.DataFrame) -> pd.DataFrame:
    """Per-group distribution summary of the heterozygosity rate."""
    rows = []
    for label, subset in (
        (UNMERGED_LABEL, frame[frame["group"] == UNMERGED_LABEL]),
        (MERGED_LABEL, frame[frame["group"] == MERGED_LABEL]),
    ):
        rate = subset["het_rate"]
        rows.append(
            {
                "group": label,
                "n_samples": int(len(subset)),
                "mean_het_rate": float(rate.mean()),
                "median_het_rate": float(rate.median()),
                "sd_het_rate": float(rate.std()),
                "q05_het_rate": float(rate.quantile(0.05)),
                "q95_het_rate": float(rate.quantile(0.95)),
                "mean_non_missing": float(subset["OBS_CT"].mean()),
            }
        )
    for extra in [g for g in group_order(frame) if g.startswith("unmerged_")]:
        subset = frame[frame["group"] == extra]
        rate = subset["het_rate"]
        rows.append(
            {
                "group": extra,
                "n_samples": int(len(subset)),
                "mean_het_rate": float(rate.mean()),
                "median_het_rate": float(rate.median()),
                "sd_het_rate": float(rate.std()),
                "q05_het_rate": float(rate.quantile(0.05)),
                "q95_het_rate": float(rate.quantile(0.95)),
                "mean_non_missing": float(subset["OBS_CT"].mean()),
            }
        )
    for n_sources in sorted(
        frame.loc[frame["n_sources"] > 1, "n_sources"].unique()
    ):
        subset = frame[frame["n_sources"] == n_sources]
        rows.append(
            {
                "group": f"merged_{n_sources}_sources",
                "n_samples": int(len(subset)),
                "mean_het_rate": float(subset["het_rate"].mean()),
                "median_het_rate": float(subset["het_rate"].median()),
                "sd_het_rate": float(subset["het_rate"].std()),
                "q05_het_rate": float(subset["het_rate"].quantile(0.05)),
                "q95_het_rate": float(subset["het_rate"].quantile(0.95)),
                "mean_non_missing": float(subset["OBS_CT"].mean()),
            }
        )
    return pd.DataFrame(rows)


def call_rate_bins(frame: pd.DataFrame, n_bins: int = 10) -> pd.DataFrame:
    """Median heterozygosity per group inside equal-count call-rate bins.

    Removing the call-rate confound is the point of the comparison: within a bin
    both groups genotype a comparable number of sites.
    """
    edges = np.unique(np.quantile(frame["OBS_CT"], np.linspace(0.0, 1.0, n_bins + 1)))
    rows = []
    for left, right in zip(edges[:-1], edges[1:]):
        cell = frame[(frame["OBS_CT"] >= left) & (frame["OBS_CT"] < right)]
        if len(cell) < 20 or cell["group"].nunique() < 2:
            continue
        rows.append(
            {
                "call_rate_million": float(cell["OBS_CT"].median() / 1e6),
                "n_samples": int(len(cell)),
                "n_unmerged": int((cell["group"] != MERGED_LABEL).sum()),
                "n_merged": int((cell["group"] == MERGED_LABEL).sum()),
                "median_het_unmerged": float(
                    cell.loc[cell["group"] != MERGED_LABEL, "het_rate"].median()
                ),
                "median_het_merged": float(
                    cell.loc[cell["group"] == MERGED_LABEL, "het_rate"].median()
                ),
            }
        )
    binned = pd.DataFrame(rows)
    binned["median_diff"] = binned["median_het_merged"] - binned["median_het_unmerged"]
    return binned


def plot_het_distribution(
    frame: pd.DataFrame, filename: str | Path, title: str, quantiles=(0.005, 0.995)
) -> None:
    """Grouped histogram of per-sample heterozygosity, merged vs unmerged."""
    low, high = frame["het_rate"].quantile(list(quantiles))
    order = group_order(frame)
    thresholds = [
        {
            "value": float(frame.loc[frame["group"] == label, "het_rate"].median()),
            "label": f"{label} median",
            "color": GROUP_COLORS.get(label, "gray"),
            "linestyle": "--",
        }
        for label in order
    ]
    hue_labels = {
        label: f"{label} (n={int((frame['group'] == label).sum())})" for label in order
    }
    infra_graph.plot_grouped_histogram(
        frame,
        x_col="het_rate",
        hue_col="group",
        hue_order=order,
        palette_map=GROUP_COLORS,
        hue_labels=hue_labels,
        title=title,
        filename=str(filename),
        xlabel="Heterozygosity rate (het / non-missing)",
        ylabel="Density",
        xlim=(float(low), float(high)),
        bins=60,
        density=True,
        thresholds_v=thresholds,
    )
    LOG.info("wrote %s", filename)


def plot_het_density(
    frame: pd.DataFrame,
    filename: str | Path,
    title: str,
    quantiles=(0.005, 0.995),
    bw_method=0.2,
) -> None:
    """KDE version of the heterozygosity distribution, one curve per group."""
    low, high = frame["het_rate"].quantile(list(quantiles))
    order = group_order(frame)
    infra_graph.plot_kde_density_lines(
        frame,
        x_col="het_rate",
        hue_col="group",
        hue_order=order,
        palette_map=GROUP_COLORS,
        hue_labels={
            label: f"{label} (n={int((frame['group'] == label).sum())})"
            for label in order
        },
        title=title,
        filename=str(filename),
        xlabel="Heterozygosity rate (het / non-missing)",
        ylabel="Density",
        xlim=(float(low), float(high)),
        bw_method=bw_method,
    )
    LOG.info("wrote %s", filename)


def plot_het_vs_call_rate(
    binned: pd.DataFrame, filename: str | Path, title: str
) -> None:
    """Median heterozygosity per group across call-rate bins."""
    infra_graph.plot_multi_line_series(
        binned,
        x_col="call_rate_million",
        y_specs=[
            {
                "y_col": "median_het_unmerged",
                "label": UNMERGED_LABEL,
                "color": GROUP_COLORS[UNMERGED_LABEL],
                "marker": "o",
            },
            {
                "y_col": "median_het_merged",
                "label": MERGED_LABEL,
                "color": GROUP_COLORS[MERGED_LABEL],
                "marker": "s",
            },
        ],
        title=title,
        filename=str(filename),
        x_label="Non-missing genotypes per sample (millions)",
        y_label="Median heterozygosity rate",
    )
    LOG.info("wrote %s", filename)


def run_merge_het_check(
    het_path: str | Path,
    merge_map_path: str | Path,
    output_prefix: str,
    info_dir: str | Path,
    plot_dir: str | Path,
    title_prefix: str = "",
    highlight_prefix: str | None = None,
    exclude_prefix: str | None = None,
) -> dict[str, Path]:
    """Summarise and plot merged vs unmerged heterozygosity.

    Writes ``<prefix>.sample.het.info.tsv``, ``<prefix>.het.summary.info.tsv``,
    ``<prefix>.het.call_rate_bins.info.tsv`` and two figures
    (``<prefix>.sample.het.distribution.png``, ``<prefix>.sample.het.vs_call_rate.png``) into the given
    directories (the process work dir; ``publishDir`` moves them into
    ``stats/<mod>/{info,plots}``). Returns the written paths.
    """
    info = Path(info_dir)
    plots = Path(plot_dir)
    info.mkdir(parents=True, exist_ok=True)
    plots.mkdir(parents=True, exist_ok=True)

    frame = build_het_table(het_path, merge_map_path, highlight_prefix)
    summary = summarize_by_group(frame)
    binned = call_rate_bins(frame)

    samples_path = info / f"{output_prefix}.sample.het.info.tsv"
    summary_path = info / f"{output_prefix}.het.summary.info.tsv"
    bins_path = info / f"{output_prefix}.het.call_rate_bins.info.tsv"
    infra_io.save_df_to_tsv(frame, str(samples_path), float_format="%.6g")
    infra_io.save_df_to_tsv(summary, str(summary_path), float_format="%.6g")
    infra_io.save_df_to_tsv(binned, str(bins_path), float_format="%.6g")

    distribution_path = plots / f"{output_prefix}.sample.het.distribution.png"
    density_path = plots / f"{output_prefix}.sample.het.distribution.kde.png"
    call_rate_path = plots / f"{output_prefix}.sample.het.vs_call_rate.png"
    label = title_prefix or "Scan3"
    plot_het_distribution(
        frame,
        distribution_path,
        f"{label}: merged vs unmerged heterozygosity",
    )
    plot_het_density(
        frame,
        density_path,
        f"{label}: merged vs unmerged heterozygosity (KDE)",
    )
    plot_het_vs_call_rate(
        binned,
        call_rate_path,
        f"{label}: heterozygosity matched by call rate",
    )

    written: dict[str, Path] = {}
    if exclude_prefix:
        # Second comparison with that cohort removed entirely (merged groups of the
        # same prefix go too), so the two remaining series are not stretched by it.
        kept = frame[~frame["IID"].str.startswith(exclude_prefix)].copy()
        kept = assign_groups(kept)
        kept_summary = summarize_by_group(kept)
        kept_summary_path = info / f"{output_prefix}.het.summary.no{exclude_prefix.strip('_')}.info.tsv"
        infra_io.save_df_to_tsv(kept_summary, str(kept_summary_path), float_format="%.6g")
        kept_plot = plots / f"{output_prefix}.sample.het.distribution.no{exclude_prefix.strip('_')}.png"
        plot_het_distribution(
            kept,
            kept_plot,
            f"{label}: merged vs unmerged, {exclude_prefix} samples excluded",
        )
        written["no_excluded_summary"] = kept_summary_path
        written["no_excluded_distribution"] = Path(kept_plot)
        LOG.info(
            "%s excluded: %d samples left (%d merged, %d unmerged)",
            exclude_prefix,
            len(kept),
            int((kept["group"] == MERGED_LABEL).sum()),
            int((kept["group"] == UNMERGED_LABEL).sum()),
        )

    unmerged = summary.loc[summary["group"] == UNMERGED_LABEL].iloc[0]
    merged = summary.loc[summary["group"] == MERGED_LABEL].iloc[0]
    LOG.info(
        "median het rate: unmerged %.5f (n=%d) vs merged %.5f (n=%d)",
        unmerged["median_het_rate"],
        unmerged["n_samples"],
        merged["median_het_rate"],
        merged["n_samples"],
    )
    if not binned.empty:
        LOG.info(
            "call-rate matched median difference: %.5f to %.5f",
            binned["median_diff"].min(),
            binned["median_diff"].max(),
        )
    return {
        **written,
        "samples": samples_path,
        "summary": summary_path,
        "call_rate_bins": bins_path,
        "distribution": Path(distribution_path),
        "density": Path(density_path),
        "vs_call_rate": Path(call_rate_path),
    }


def compare_merged_before_after(
    het_before_path: str | Path,
    het_after_path: str | Path,
    merge_map_path: str | Path,
) -> pd.DataFrame:
    """Per merged germplasm: source-sample heterozygosity vs the merged sample.

    ``het_before`` is the per-sample ``--het`` table of the pre-merge genotypes
    (every source taxon as its own sample), ``het_after`` the table of the scan3
    output where the germplasm is one sample. Requested sources missing from the
    before set are counted, not silently dropped.
    """
    before = add_het_rate(load_plink_het(het_before_path))[["IID", "het_rate", "OBS_CT"]].rename(
        columns={"het_rate": "het_before", "OBS_CT": "obs_ct_before"}
    )
    after = add_het_rate(load_plink_het(het_after_path))[["IID", "het_rate", "OBS_CT"]].rename(
        columns={"het_rate": "het_after", "OBS_CT": "obs_ct_after"}
    )
    merge = load_merge_map(merge_map_path)
    merged = merge[merge["n_sources"] > 1]

    before_by_id = before.set_index("IID")
    after_by_id = after.set_index("IID")
    rows = []
    missing_sources = 0
    for row in merged.itertuples():
        sources = [s for s in row.members.split("|") if s]
        found = [s for s in sources if s in before_by_id.index]
        missing_sources += len(sources) - len(found)
        if row.IID not in after_by_id.index or not found:
            continue
        rates = before_by_id.loc[found, "het_before"].to_numpy(dtype=np.float64)
        after_rate = float(after_by_id.loc[row.IID, "het_after"])
        rows.append(
            {
                "merged_taxon": row.IID,
                "n_sources": int(row.n_sources),
                "n_sources_found": len(found),
                "sources": "|".join(sources),
                "het_before_mean": float(np.mean(rates)),
                "het_before_median": float(np.median(rates)),
                "het_before_min": float(np.min(rates)),
                "het_before_max": float(np.max(rates)),
                "het_after": after_rate,
                "delta_median": after_rate - float(np.median(rates)),
                "delta_mean": after_rate - float(np.mean(rates)),
                "obs_ct_after": float(after_by_id.loc[row.IID, "obs_ct_after"]),
            }
        )
    table = pd.DataFrame(rows)
    if table.empty:
        raise ValueError("no merged germplasm could be matched across both het tables")
    if missing_sources:
        LOG.warning(
            "%d source taxa are absent from %s and were skipped",
            missing_sources,
            het_before_path,
        )
    LOG.info(
        "compared %d merged germplasms (%d sources missing before)",
        len(table),
        missing_sources,
    )
    return table.sort_values("delta_median", ascending=False).reset_index(drop=True)


def source_level_before_after(
    het_before_path: str | Path,
    het_after_path: str | Path,
    merge_map_path: str | Path,
) -> pd.DataFrame:
    """One row per source sample: its own heterozygosity vs its merged sample's.

    ``compare_merged_before_after`` collapses each germplasm to the median of its
    sources, which weights a 2-source group like a 48-source one and hides how
    far apart the sources are. This keeps every source as its own point against
    the same merged value, so a group whose sources disagree shows a spread.
    """
    before = add_het_rate(load_plink_het(het_before_path))
    after = add_het_rate(load_plink_het(het_after_path))
    merge = load_merge_map(merge_map_path)
    before_by_id = before.set_index("IID")
    after_by_id = after.set_index("IID")
    rows = []
    for row in merge[merge["n_sources"] > 1].itertuples():
        if row.IID not in after_by_id.index:
            continue
        merged_rate = float(after_by_id.loc[row.IID, "het_rate"])
        for source in row.members.split("|"):
            if not source or source not in before_by_id.index:
                continue
            rows.append(
                {
                    "merged_taxon": row.IID,
                    "source": source,
                    "n_sources": int(row.n_sources),
                    "het_source": float(before_by_id.loc[source, "het_rate"]),
                    "het_merged": merged_rate,
                }
            )
    table = pd.DataFrame(rows)
    if table.empty:
        raise ValueError("no source sample could be matched across both het tables")
    table["delta"] = table["het_merged"] - table["het_source"]
    LOG.info(
        "source-level comparison: %d sources over %d germplasms; %.1f%% of sources "
        "sit below their merged value",
        len(table),
        table["merged_taxon"].nunique(),
        100 * float((table["delta"] > 0).mean()),
    )
    return table


def plot_source_before_after(
    table: pd.DataFrame, filename: str | Path, title: str
) -> None:
    """Every source sample against its merged germplasm, one point per source."""
    infra_graph.plot_scatter_with_thresholds(
        table,
        x_col="het_source",
        y_col="het_merged",
        title=title,
        filename=str(filename),
        xlabel="Source sample: heterozygosity (before merge)",
        ylabel="Merged germplasm: heterozygosity (after merge)",
        alpha=0.35,
        s=12,
        identity_line=True,
        identity_label=r"$y=x$",
    )
    LOG.info("wrote %s", filename)


def summarize_before_after(table: pd.DataFrame) -> pd.DataFrame:
    """Delta summary overall and per number of sources."""
    rows = []

    def _row(label: str, subset: pd.DataFrame) -> dict:
        return {
            "group": label,
            "n_groups": int(len(subset)),
            "median_het_before": float(subset["het_before_median"].median()),
            "median_het_after": float(subset["het_after"].median()),
            "median_delta": float(subset["delta_median"].median()),
            "q05_delta": float(subset["delta_median"].quantile(0.05)),
            "q95_delta": float(subset["delta_median"].quantile(0.95)),
            "frac_after_higher": float((subset["delta_median"] > 0).mean()),
        }

    rows.append(_row("all_merged", table))
    for n_sources in sorted(table["n_sources"].unique()):
        rows.append(_row(f"merged_{n_sources}_sources", table[table["n_sources"] == n_sources]))
    return pd.DataFrame(rows)


def plot_before_after_scatter(
    table: pd.DataFrame, filename: str | Path, title: str
) -> None:
    """Source-median heterozygosity against the merged sample's heterozygosity."""
    infra_graph.plot_scatter_with_thresholds(
        table,
        x_col="het_before_median",
        y_col="het_after",
        title=title,
        filename=str(filename),
        xlabel="Source samples: median heterozygosity (before merge)",
        ylabel="Merged sample: heterozygosity (after merge)",
        identity_line=True,
        identity_label=r"$y=x$",
    )
    LOG.info("wrote %s", filename)


def plot_delta_distribution(
    table: pd.DataFrame, filename: str | Path, title: str
) -> None:
    """Histogram of per-germplasm heterozygosity change (after - before)."""
    infra_graph.plot_distribution_with_stats(
        data=table,
        col="delta_median",
        title=title,
        filename=str(filename),
        mean_val=float(table["delta_median"].mean()),
        median_val=float(table["delta_median"].median()),
        thresholds=[{"value": 0.0, "label": "no change", "color": "black", "linestyle": "--"}],
        x_label="Heterozygosity change after merge (merged - source median)",
        y_label="Germplasms",
        color="darkorange",
        bins=50,
    )
    LOG.info("wrote %s", filename)


def run_merged_before_after_check(
    het_before_path: str | Path,
    het_after_path: str | Path,
    merge_map_path: str | Path,
    output_prefix: str,
    info_dir: str | Path,
    plot_dir: str | Path,
    title_prefix: str = "",
) -> dict[str, Path]:
    """Compare each merged germplasm's heterozygosity before and after merging."""
    info = Path(info_dir)
    plots = Path(plot_dir)
    info.mkdir(parents=True, exist_ok=True)
    plots.mkdir(parents=True, exist_ok=True)

    table = compare_merged_before_after(het_before_path, het_after_path, merge_map_path)
    summary = summarize_before_after(table)

    table_path = info / f"{output_prefix}.merged.before_after.info.tsv"
    summary_path = info / f"{output_prefix}.merged.before_after.summary.info.tsv"
    infra_io.save_df_to_tsv(table, str(table_path), float_format="%.6g")
    infra_io.save_df_to_tsv(summary, str(summary_path), float_format="%.6g")

    label = title_prefix or "Scan3"
    sources = source_level_before_after(het_before_path, het_after_path, merge_map_path)
    sources_path = info / f"{output_prefix}.merged.before_after.sources.info.tsv"
    infra_io.save_df_to_tsv(sources, str(sources_path), float_format="%.6g")
    sources_plot = plots / f"{output_prefix}.merged.before_after.sources.png"
    plot_source_before_after(
        sources,
        sources_plot,
        f"{label}: each source sample vs its merged germplasm",
    )

    scatter_path = plots / f"{output_prefix}.merged.before_after.scatter.png"
    delta_path = plots / f"{output_prefix}.merged.before_after.delta.png"
    plot_before_after_scatter(
        table, scatter_path, f"{label}: heterozygosity before vs after merging"
    )
    plot_delta_distribution(
        table, delta_path, f"{label}: heterozygosity change after merging"
    )

    overall = summary.iloc[0]
    LOG.info(
        "median delta %.6f (before %.6f -> after %.6f), %.1f%% of germplasms higher after",
        overall["median_delta"],
        overall["median_het_before"],
        overall["median_het_after"],
        100 * overall["frac_after_higher"],
    )
    return {
        "per_germplasm": table_path,
        "summary": summary_path,
        "sources": sources_path,
        "sources_plot": Path(sources_plot),
        "scatter": Path(scatter_path),
        "delta": Path(delta_path),
    }
