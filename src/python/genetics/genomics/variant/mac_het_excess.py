"""Wheat selfing-model excess-heterozygosity filter (probability sweep, no pfile rewrite).

Given per-site MAC and heterozygous sample count H from PLINK2 ``--geno-counts``,
approximate

    H | MAC, F  ~  Binomial(MAC, 1-F)

and compute the upper-tail probability ``P(H >= H_obs | MAC, F)``.  Sites are
never rewritten into a new pfile; callers sweep retain/drop counts across
hypothesis-testing thresholds.

When diploid N is supplied, the success probability uses the rare-variant
correction ``(1 - MAC/(2N))(1-F)``.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence

import numpy as np
import pandas as pd
from scipy.stats import binom

from genetics.genomics.variant.mac import (
    _compute_mac_table,
    _plot_mac_h_heatmap,
    _plot_mac_hetfrac_heatmap,
)
from genetics.genomics.variant.variant_utils import load_df_from_plink_gcount
from infra.utils.errors import fail
from genetics.genomics.variant.popdep import load_df_from_tsv_popdep_info
from infra.utils.graph import (
    SEQUENTIAL_WHITE_RED,
    TWO_CATEGORY_RED,
    plot_bar_chart,
    plot_joint_regression,
    plot_multi_line_series,
)
from infra.utils.io import load_df_from_space_sep, save_df_to_tsv, save_thresholds

logger = logging.getLogger(__name__)

DEFAULT_P_THRESHOLDS: tuple[float, ...] = (
    5e-2,
    1e-2,
    1e-3,
    1e-4,
    1e-5,
    1e-6,
    1e-7,
    1e-8,
    1e-9,
    1e-10,
    1e-12,
)
DEFAULT_F_GRID: tuple[float, ...] = (0.90, 0.95, 0.98, 0.99)
ALL_HET_MAC_MAX = 10
MAC_MAX_ALLOWED_TABLE = 100
RETAIN_ALPHA = 1e-6
HEATMAP_MAC_MAXES = (100, 200, 500, 1000)
HEATMAP_MAC_FILTER_OVERLAYS = (100, 200, 500)
HEATMAP_MAC_MIN_KEEP = 5
HEATMAP_FRACTION_GAMMA = 0.35


def excess_het_cutoff_heatmap_coords(
    mac_max: int,
    f_panel: float,
    alpha: float = RETAIN_ALPHA,
    n_diploid: float | None = None,
) -> tuple[list[float], list[float]]:
    """Heatmap (x, y) polyline for the largest retained H at each MAC.

    ``sns.heatmap`` cell *i* spans ``[i, i+1]``; H is flipped so H=0 is the
    bottom row. The line sits on the top edge of the last retained H cell
    (sites with larger H are dropped).
    """
    xs: list[float] = []
    ys: list[float] = []
    for mac in range(1, int(mac_max) + 1):
        h_cut = max_allowed_het(mac, f_panel, alpha, n_diploid=n_diploid)
        y_edge = float(mac_max) - float(h_cut)
        xs.extend([float(mac), float(mac) + 1.0])
        ys.extend([y_edge, y_edge])
    return xs, ys


_FLOAT_TSV = "%.8g"


def _as_float_array(values) -> np.ndarray:
    return np.asarray(values, dtype=float)


def _clip_f(f: float) -> float:
    f = float(f)
    if not np.isfinite(f):
        fail(f"Inbreeding coefficient F must be finite, got {f}")
    return float(min(max(f, 0.0), 1.0 - 1e-12))


def _het_success_probability(mac: np.ndarray, f: float, n_diploid: np.ndarray | None) -> np.ndarray:
    p = np.full(np.shape(mac), 1.0 - _clip_f(f), dtype=float)
    if n_diploid is None:
        return np.clip(p, 0.0, 1.0)
    two_n = 2.0 * np.asarray(n_diploid, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        freq_term = 1.0 - mac / two_n
    freq_term = np.where((two_n > 0) & np.isfinite(freq_term), freq_term, 1.0)
    return np.clip(p * np.clip(freq_term, 0.0, 1.0), 0.0, 1.0)


def p_excess_het_binomial(
    h_obs,
    mac,
    f: float,
    n_diploid=None,
) -> np.ndarray | float:
    """Upper-tail ``P(H >= h_obs | MAC=m, F)`` under ``H ~ Binomial(m, 1-F)``.

    ``n_diploid`` (optional) applies ``p = (1 - MAC/(2N))(1-F)``.  MAC=0 or
    ``h_obs <= 0`` returns 1 (the site is not an excess-het candidate).
    """
    h = _as_float_array(h_obs)
    m = _as_float_array(mac)
    scalar = np.ndim(h_obs) == 0 and np.ndim(mac) == 0 and (n_diploid is None or np.ndim(n_diploid) == 0)
    n_arr = None if n_diploid is None else _as_float_array(n_diploid)
    p = _het_success_probability(m, f, n_arr)
    n_int = np.clip(np.rint(m), 0, None).astype(int)
    k = h - 1.0
    with np.errstate(invalid="ignore"):
        tail = binom.sf(k, n=n_int, p=p)
    tail = np.where((m <= 0) | (h <= 0), 1.0, tail)
    tail = np.clip(np.nan_to_num(tail, nan=1.0), 0.0, 1.0)
    if scalar:
        return float(np.asarray(tail).reshape(-1)[0])
    return tail


def max_allowed_het(mac: int, f: float, alpha: float, n_diploid: float | None = None) -> int:
    """Largest heterozygous count H still retained at ``P_excess >= alpha``."""
    mac_i = int(mac)
    if mac_i <= 0:
        return 0
    h = np.arange(0, mac_i + 1, dtype=int)
    n = None if n_diploid is None else np.full(h.shape, n_diploid)
    p = p_excess_het_binomial(h, np.full(h.shape, mac_i), f, n_diploid=n)
    ok = np.asarray(p) >= float(alpha)
    if not np.any(ok):
        return 0
    return int(h[ok][-1])


def median_inbreeding_from_het_table(df: pd.DataFrame) -> float:
    """Panel inbreeding coefficient as median of sample-level PLINK2 ``--het`` F."""
    if df is None or df.empty:
        fail("PLINK2 --het table is empty; cannot estimate F_panel")
    work = df.copy()
    work.columns = [str(c).replace("#", "") for c in work.columns]
    if "F" not in work.columns:
        fail(f"PLINK2 --het table missing F column; got {list(work.columns)}")
    f_vals = pd.to_numeric(work["F"], errors="coerce").dropna()
    if f_vals.empty:
        fail("PLINK2 --het F column has no finite values")
    return float(np.nanmedian(f_vals.to_numpy()))


def load_median_inbreeding_from_het(het_path: str) -> float:
    df = load_df_from_space_sep(het_path)
    return median_inbreeding_from_het_table(df)


def _mac_stratum(mac: int) -> str:
    if mac <= 5:
        return str(int(mac))
    if mac <= 10:
        return "6-10"
    if mac <= 20:
        return "11-20"
    if mac <= 50:
        return "21-50"
    if mac <= 100:
        return "51-100"
    return ">100"


def retain_sites_by_p_excess(sites: pd.DataFrame, alpha: float = RETAIN_ALPHA) -> pd.DataFrame:
    """Keep MAC>0 sites with ``p_excess >= alpha`` (fail to reject the selfing model)."""
    mac = pd.to_numeric(sites["MAC"], errors="coerce")
    p = pd.to_numeric(sites["p_excess"], errors="coerce")
    return sites[(mac > 0) & (p >= float(alpha))].copy()


def drop_sites_by_p_excess(sites: pd.DataFrame, alpha: float = RETAIN_ALPHA) -> pd.DataFrame:
    """Keep MAC>0 sites with ``p_excess < alpha`` (reject the selfing model)."""
    mac = pd.to_numeric(sites["MAC"], errors="coerce")
    p = pd.to_numeric(sites["p_excess"], errors="coerce")
    return sites[(mac > 0) & (p < float(alpha))].copy()


def _heatmap_x_tick_step(mac_max: int) -> int:
    if mac_max <= 100:
        return 10
    if mac_max <= 200:
        return 20
    if mac_max <= 500:
        return 50
    return 100


def _plot_mac_het_excess_heatmaps(
    sites: pd.DataFrame,
    retained: pd.DataFrame,
    dropped: pd.DataFrame,
    output_prefix: str,
    alpha: float,
    f_data: float,
    mac_max: int,
) -> None:
    tick = _heatmap_x_tick_step(mac_max)
    span = f"0-{mac_max}"
    stem_span = f"0_{mac_max}"
    _plot_mac_hetfrac_heatmap(
        sites, output_prefix, mac_max=mac_max, x_tick_step=tick, stem=f"heatmap_mac_hetfrac.{stem_span}"
    )
    _plot_mac_h_heatmap(
        sites, output_prefix, mac_max=mac_max, x_tick_step=tick, stem=f"heatmap_mac_h.{stem_span}"
    )
    _plot_mac_hetfrac_heatmap(
        retained,
        output_prefix,
        mac_max=mac_max,
        x_tick_step=tick,
        title=(
            f"Observed heterozygosity per MAC ({span})\n"
            f"retained at P_excess>={alpha:g} (F={f_data:.4g})"
        ),
        stem=f"heatmap_mac_hetfrac.retained.{stem_span}",
    )
    _plot_mac_h_heatmap(
        retained,
        output_prefix,
        mac_max=mac_max,
        x_tick_step=tick,
        title=(
            f"Heterozygous count (H) per MAC ({span})\n"
            f"retained at P_excess>={alpha:g} (F={f_data:.4g})"
        ),
        stem=f"heatmap_mac_h.retained.{stem_span}",
    )
    _plot_mac_hetfrac_heatmap(
        dropped,
        output_prefix,
        mac_max=mac_max,
        x_tick_step=tick,
        title=(
            f"Observed heterozygosity per MAC ({span})\n"
            f"dropped at P_excess<{alpha:g} (F={f_data:.4g})"
        ),
        stem=f"heatmap_mac_hetfrac.dropped.{stem_span}",
    )
    _plot_mac_h_heatmap(
        dropped,
        output_prefix,
        mac_max=mac_max,
        x_tick_step=tick,
        title=(
            f"Heterozygous count (H) per MAC ({span})\n"
            f"dropped at P_excess<{alpha:g} (F={f_data:.4g})"
        ),
        stem=f"heatmap_mac_h.dropped.{stem_span}",
    )
    if mac_max in HEATMAP_MAC_FILTER_OVERLAYS:
        n_diploid = None
        if "Total_Samples" in sites.columns:
            n_med = pd.to_numeric(sites["Total_Samples"], errors="coerce").median()
            if pd.notna(n_med) and float(n_med) > 0:
                n_diploid = float(n_med)
        curve_x, curve_y = excess_het_cutoff_heatmap_coords(
            mac_max, f_data, alpha=alpha, n_diploid=n_diploid
        )
        _plot_mac_h_heatmap(
            sites,
            output_prefix,
            mac_max=mac_max,
            x_tick_step=tick,
            title=(
                f"Heterozygous count (H) per MAC ({span})\n"
                f"excess-het cutoff and MAC<{HEATMAP_MAC_MIN_KEEP} "
                f"(P_excess<{alpha:g}, F={f_data:.4g})"
            ),
            stem=f"heatmap_mac_h.filter.{stem_span}",
            cmap=SEQUENTIAL_WHITE_RED,
            gamma=HEATMAP_FRACTION_GAMMA,
            vmax=1.0,
            overlays=[
                {
                    "x": curve_x,
                    "y": curve_y,
                    "color": "#111111",
                    "linestyle": "-",
                    "linewidth": 2.0,
                    "label": f"Excess-het cutoff (P < {alpha:g})",
                }
            ],
            vlines=[
                {
                    "x": float(HEATMAP_MAC_MIN_KEEP),
                    "color": TWO_CATEGORY_RED,
                    "linestyle": "--",
                    "linewidth": 1.8,
                    "label": f"MAC < {HEATMAP_MAC_MIN_KEEP}",
                }
            ],
        )


def ana_mac_hetfrac_heatmap_retained(
    input_file: str,
    output_prefix: str = "variant_mac",
    het_path: str | None = None,
    f_panel: float | None = None,
    alpha: float = RETAIN_ALPHA,
) -> None:
    """MAC vs H and MAC vs H/MAC heatmaps at ``p_excess >= alpha``.

    Uses median sample F from PLINK2 ``--het`` when ``het_path`` is set.
    Writes unfiltered / retained / dropped MAC-vs-het-fraction and MAC-vs-H
    heatmaps at MAC windows 0-100, 0-200, 0-500, and 0-1000
    (``heatmap_mac_hetfrac*.0_{100,200,500,1000}.png``, ``heatmap_mac_h*.0_{100,200,500,1000}.png``)
    plus MAC-vs-H overlays (``heatmap_mac_h.filter.0_{100,200,500}.png``) with a
    stretched white-to-red scale, the excess-het H cutoff, and a MAC=5 line.
    Also writes ``{output_prefix}.het_retain.th.tsv``. Does not emit a filtered pfile.
    """
    df_raw = load_df_from_plink_gcount(input_file)
    if het_path:
        f_data = load_median_inbreeding_from_het(het_path)
        f_source = "het_median"
    elif f_panel is not None:
        f_data = float(f_panel)
        f_source = "supplied"
    else:
        fail("Provide het_path (PLINK2 --het) or f_panel for the retained MAC heatmap")

    sites = site_excess_het_table(df_raw, f_data)
    retained = retain_sites_by_p_excess(sites, alpha=alpha)
    n_tested = int((pd.to_numeric(sites["MAC"], errors="coerce") > 0).sum())
    n_retain = int(len(retained))
    stats = {
        "F_panel": f_data,
        "F_source": f_source,
        "p_threshold": float(alpha),
        "N_tested_MAC_gt_0": n_tested,
        "N_retain": n_retain,
        "Frac_retain": (n_retain / n_tested) if n_tested else 0.0,
    }
    save_thresholds(stats, f"{output_prefix}.het_retain.th.tsv")
    print(stats)
    dropped = drop_sites_by_p_excess(sites, alpha=alpha)
    for mac_max in HEATMAP_MAC_MAXES:
        _plot_mac_het_excess_heatmaps(
            sites, retained, dropped, output_prefix, alpha, f_data, mac_max
        )


HET_FILTER_PASS = 0
HET_FILTER_FAIL = 1
HET_FILTER_LABELS = {
    HET_FILTER_PASS: "pass het filter",
    HET_FILTER_FAIL: "fail het filter",
}
HET_FILTER_COLORS = {
    HET_FILTER_PASS: "#4c78a8",
    HET_FILTER_FAIL: "#c44e52",
}
HET_FILTER_HUE_ORDER = (
    HET_FILTER_LABELS[HET_FILTER_PASS],
    HET_FILTER_LABELS[HET_FILTER_FAIL],
)
HET_FILTER_COLORS_BY_LABEL = {
    HET_FILTER_LABELS[HET_FILTER_PASS]: HET_FILTER_COLORS[HET_FILTER_PASS],
    HET_FILTER_LABELS[HET_FILTER_FAIL]: HET_FILTER_COLORS[HET_FILTER_FAIL],
}
POPDEP_HET_BINS_X = 250
POPDEP_HET_BINS_Y = 150
_POPDEP_HET_SCATTERS = (
    ("Depth_Mean", "Depth_SD", "Mean Depth", "Depth SD"),
    ("RelativeDepth_Mean", "RelativeDepth_SD", "Relative Mean Depth", "Relative Depth SD"),
    ("Depth_Mean", "Depth_CV", "Mean Depth", "Depth CV"),
    ("RelativeDepth_Mean", "RelativeDepth_CV", "Relative Mean Depth", "Relative Depth CV"),
)


def _sites_with_het_filter_category(sites: pd.DataFrame, alpha: float) -> pd.DataFrame:
    """MAC>0 sites with integer het-filter flag (0=pass drawn first, 1=fail on top)."""
    mac = pd.to_numeric(sites["MAC"], errors="coerce")
    tested = sites[mac > 0].copy()
    p = pd.to_numeric(tested["p_excess"], errors="coerce")
    tested["het_filter"] = np.where(p >= float(alpha), HET_FILTER_PASS, HET_FILTER_FAIL).astype(int)
    return tested


def _join_popdep_het_filter(
    sites: pd.DataFrame,
    popdep: pd.DataFrame,
    alpha: float,
) -> pd.DataFrame:
    tested = _sites_with_het_filter_category(sites, alpha)
    if "ID" not in tested.columns or "ID" not in popdep.columns:
        fail("gcount and popdep tables must both have an ID column")
    depth_cols = [
        "ID",
        "Depth_Mean",
        "Depth_SD",
        "Depth_CV",
        "RelativeDepth_Mean",
        "RelativeDepth_SD",
        "RelativeDepth_CV",
    ]
    missing = [c for c in depth_cols if c not in popdep.columns]
    if missing:
        fail(f"popdep table missing columns: {missing}")
    merged = tested.merge(popdep[depth_cols], on="ID", how="inner")
    for col in depth_cols[1:]:
        merged[col] = pd.to_numeric(merged[col], errors="coerce")
    n_depth = int(merged["Depth_Mean"].notna().sum())
    if n_depth == 0:
        fail(
            "popdep Depth_Mean has 0 numeric values; re-annotate variants from a "
            "tabix popdep grid before plotting het-filter scatters"
        )
    return merged


def ana_popdep_scatter_by_het_filter(
    gcount_file: str,
    popdep_file: str,
    output_prefix: str = "variant.popdep.het_filter",
    het_path: str | None = None,
    f_panel: float | None = None,
    alpha: float = RETAIN_ALPHA,
) -> None:
    """Popdep mean-vs-SD and mean-vs-CV scatters coloured by the wheat excess-het
    retain/drop call.

    Writes linear and log-log JointGrid pairs (scatter plus 100% fill x/y
    marginals: within-bin pass/fail counts). Fail is drawn last; pass and fail
    share the same point alpha. Marginal bins are 250 (x) by 150 (y).
    ``{output_prefix}.scatter_Depth_SD_vs_Depth_Mean.png`` /
    ``{output_prefix}.scatter_Depth_CV_vs_Depth_Mean.png`` (and ``*.log.png``)
    plus the relative-depth equivalents. Pass sites (steelblue) are drawn first;
    fail sites (crimson) on top. Does not emit a filtered pfile.
    """
    df_raw = load_df_from_plink_gcount(gcount_file)
    if het_path:
        f_data = load_median_inbreeding_from_het(het_path)
    elif f_panel is not None:
        f_data = float(f_panel)
    else:
        fail("Provide het_path (PLINK2 --het) or f_panel for popdep het-filter scatters")

    sites = site_excess_het_table(df_raw, f_data)
    popdep = load_df_from_tsv_popdep_info(popdep_file)
    if popdep is None:
        fail(f"Failed to load popdep info: {popdep_file}")
    merged = _join_popdep_het_filter(sites, popdep, alpha)
    n_pass = int((merged["het_filter"] == HET_FILTER_PASS).sum())
    n_fail = int((merged["het_filter"] == HET_FILTER_FAIL).sum())
    print(
        {
            "F_panel": f_data,
            "p_threshold": float(alpha),
            "N_joined": int(len(merged)),
            "N_pass": n_pass,
            "N_fail": n_fail,
        }
    )

    wrote = 0
    for x_col, y_col, x_label, y_label in _POPDEP_HET_SCATTERS:
        plot_df = merged.dropna(subset=[x_col, y_col])
        if plot_df.empty:
            logger.warning("Skipping %s vs %s scatter; no numeric rows", y_col, x_col)
            continue
        for log_axes in (False, True):
            frame = plot_df
            if log_axes:
                frame = plot_df[(plot_df[x_col] > 0) & (plot_df[y_col] > 0)]
                if frame.empty:
                    logger.warning(
                        "Skipping log-log %s vs %s scatter; no positive rows", y_col, x_col
                    )
                    continue
            scale = ", log-log" if log_axes else ""
            stem = f"{output_prefix}.scatter_{y_col}_vs_{x_col}"
            filename = f"{stem}.log.png" if log_axes else f"{stem}.png"
            plot_frame = frame.copy()
            plot_frame["het_filter_label"] = plot_frame["het_filter"].map(HET_FILTER_LABELS)
            plot_joint_regression(
                plot_frame,
                x_col=x_col,
                y_col=y_col,
                group_col="het_filter_label",
                x_label=x_label,
                y_label=y_label,
                filename=filename,
                title=(
                    f"{y_label} vs {x_label} by het filter{scale} "
                    f"(P_excess>={alpha:g}, F={f_data:.4g})"
                ),
                figure_size=(10, 8),
                palette=HET_FILTER_COLORS_BY_LABEL,
                hue_order=HET_FILTER_HUE_ORDER,
                log_x=log_axes,
                log_y=log_axes,
                draw_regression=False,
                alpha=0.25,
                s=4,
                rasterized=True,
                legend_title=None,
                bins_x=POPDEP_HET_BINS_X,
                bins_y=POPDEP_HET_BINS_Y,
                min_points=1,
                marginal_multiple="fill",
                stack_last_on_top=True,
            )
            wrote += 1
    if wrote == 0:
        fail("No popdep het-filter scatters written; depth columns were empty")


def site_excess_het_table(gcount_df: pd.DataFrame, f_panel: float) -> pd.DataFrame:
    """Per-site MAC, H, expected H, excess-het tail probability, and all-het flag."""
    return _annotate_excess_het(_compute_mac_table(gcount_df), f_panel)


def _annotate_excess_het(mac_df: pd.DataFrame, f_panel: float) -> pd.DataFrame:
    out = mac_df.copy()
    mac = pd.to_numeric(out["MAC"], errors="coerce").fillna(0).to_numpy(dtype=float)
    h = pd.to_numeric(out["HET_REF_ALT_CTS"], errors="coerce").fillna(0).to_numpy(dtype=float)
    n = pd.to_numeric(out["Total_Samples"], errors="coerce").to_numpy(dtype=float)
    f = _clip_f(f_panel)
    out["H"] = h
    out["F_panel"] = f
    out["E_H"] = mac * (1.0 - mac / np.clip(2.0 * n, 1.0, None)) * (1.0 - f)
    out["p_excess"] = p_excess_het_binomial(h, mac, f, n_diploid=n)
    out["all_het_low_mac"] = (mac > 0) & (mac <= ALL_HET_MAC_MAX) & np.isclose(h, mac)
    return out


def sweep_excess_het_thresholds(
    sites: pd.DataFrame,
    f_panel: float,
    thresholds: Sequence[float] = DEFAULT_P_THRESHOLDS,
    f_source: str = "supplied",
) -> pd.DataFrame:
    """Retain/drop counts at each tail-probability threshold (tested sites: MAC>0)."""
    tested = sites[pd.to_numeric(sites["MAC"], errors="coerce") > 0].copy()
    n_tested = int(len(tested))
    n_all_het = int(tested["all_het_low_mac"].sum()) if "all_het_low_mac" in tested.columns else 0
    p = pd.to_numeric(tested["p_excess"], errors="coerce") if n_tested else pd.Series(dtype=float)
    rows = []
    for alpha in thresholds:
        alpha_f = float(alpha)
        n_retain = int((p >= alpha_f).sum()) if n_tested else 0
        n_drop = n_tested - n_retain
        rows.append(
            {
                "f_source": f_source,
                "f_panel": float(f_panel),
                "p_threshold": alpha_f,
                "n_tested": n_tested,
                "n_retain": n_retain,
                "n_drop": n_drop,
                "frac_retain": (n_retain / n_tested) if n_tested else np.nan,
                "n_all_het_mac_le_10": n_all_het,
            }
        )
    return pd.DataFrame(rows)


def sweep_excess_het_by_mac(
    sites: pd.DataFrame,
    f_panel: float,
    thresholds: Sequence[float] = DEFAULT_P_THRESHOLDS,
    f_source: str = "supplied",
) -> pd.DataFrame:
    """Retain counts stratified by MAC bin at each threshold."""
    tested = sites[pd.to_numeric(sites["MAC"], errors="coerce") > 0].copy()
    if tested.empty:
        return pd.DataFrame(
            columns=[
                "f_source",
                "f_panel",
                "p_threshold",
                "mac_bin",
                "n_tested",
                "n_retain",
                "n_drop",
                "frac_retain",
            ]
        )
    tested = tested.copy()
    tested["mac_bin"] = tested["MAC"].astype(int).map(_mac_stratum)
    rows = []
    for alpha in thresholds:
        alpha_f = float(alpha)
        keep = pd.to_numeric(tested["p_excess"], errors="coerce") >= alpha_f
        for mac_bin, grp in tested.groupby("mac_bin", sort=False):
            n_tested = int(len(grp))
            n_retain = int(keep.loc[grp.index].sum())
            rows.append(
                {
                    "f_source": f_source,
                    "f_panel": float(f_panel),
                    "p_threshold": alpha_f,
                    "mac_bin": mac_bin,
                    "n_tested": n_tested,
                    "n_retain": n_retain,
                    "n_drop": n_tested - n_retain,
                    "frac_retain": (n_retain / n_tested) if n_tested else np.nan,
                }
            )
    order = [str(i) for i in range(1, 6)] + ["6-10", "11-20", "21-50", "51-100", ">100"]
    out = pd.DataFrame(rows)
    out["mac_bin"] = pd.Categorical(out["mac_bin"], categories=order, ordered=True)
    return out.sort_values(["p_threshold", "mac_bin"]).reset_index(drop=True)


def max_allowed_het_table(
    f_panel: float,
    alpha: float = 1e-6,
    mac_max: int = MAC_MAX_ALLOWED_TABLE,
    n_diploid: float | None = None,
) -> pd.DataFrame:
    rows = []
    for mac in range(1, int(mac_max) + 1):
        allowed = max_allowed_het(mac, f_panel, alpha, n_diploid=n_diploid)
        e_h = mac * (1.0 - f_panel) if n_diploid is None else mac * (1.0 - mac / (2.0 * n_diploid)) * (1.0 - f_panel)
        rows.append(
            {
                "MAC": mac,
                "f_panel": float(f_panel),
                "p_threshold": float(alpha),
                "E_H": e_h,
                "max_allowed_H": allowed,
            }
        )
    return pd.DataFrame(rows)


def _parse_float_list(raw: str | Sequence[float] | None, default: Sequence[float]) -> tuple[float, ...]:
    if raw is None:
        return tuple(float(x) for x in default)
    if isinstance(raw, str):
        parts = [p.strip() for p in raw.split(",") if p.strip()]
        return tuple(float(p) for p in parts)
    return tuple(float(x) for x in raw)


def _all_het_flag_summary(sites: pd.DataFrame) -> pd.DataFrame:
    tested = sites[pd.to_numeric(sites["MAC"], errors="coerce") > 0]
    low = tested[(tested["MAC"] >= 1) & (tested["MAC"] <= ALL_HET_MAC_MAX)].copy()
    if low.empty:
        return pd.DataFrame({"MAC": [], "n_sites": [], "n_all_het": [], "frac_all_het": []})
    low["MAC"] = low["MAC"].astype(int)
    grouped = (
        low.groupby("MAC", sort=True)
        .agg(n_sites=("MAC", "size"), n_all_het=("all_het_low_mac", "sum"))
        .reset_index()
    )
    grouped["frac_all_het"] = grouped["n_all_het"] / grouped["n_sites"]
    return grouped


def _plot_retain_sweep(sweep: pd.DataFrame, output_prefix: str) -> None:
    work = sweep.copy()
    work["neg_log10_p"] = -np.log10(pd.to_numeric(work["p_threshold"], errors="coerce"))
    work["f_label"] = work.apply(
        lambda r: f"F={r['f_panel']:.4g} ({r['f_source']})",
        axis=1,
    )
    wide_frac = work.pivot_table(
        index="neg_log10_p",
        columns="f_label",
        values="frac_retain",
        aggfunc="first",
    ).sort_index().reset_index()
    y_specs = [
        {"y_col": col, "label": col, "linewidth": 1.5, "marker": "o"}
        for col in wide_frac.columns
        if col != "neg_log10_p"
    ]
    if y_specs:
        plot_multi_line_series(
            wide_frac,
            x_col="neg_log10_p",
            y_specs=y_specs,
            title="Sites retained vs excess-heterozygosity tail threshold",
            filename=f"{output_prefix}.sweep.frac_retain.png",
            x_label="-log10(P_excess threshold)",
            y_label="Fraction of MAC>0 sites retained",
        )
    wide_n = work.pivot_table(
        index="neg_log10_p",
        columns="f_label",
        values="n_retain",
        aggfunc="first",
    ).sort_index().reset_index()
    y_specs_n = [
        {"y_col": col, "label": col, "linewidth": 1.5, "marker": "o"}
        for col in wide_n.columns
        if col != "neg_log10_p"
    ]
    if y_specs_n:
        plot_multi_line_series(
            wide_n,
            x_col="neg_log10_p",
            y_specs=y_specs_n,
            title="Number of sites retained vs excess-heterozygosity tail threshold",
            filename=f"{output_prefix}.sweep.n_retain.png",
            x_label="-log10(P_excess threshold)",
            y_label="Sites retained (MAC>0)",
        )


def ana_mac_het_excess_sweep(
    input_file: str,
    output_prefix: str = "variant_mac_het",
    het_path: str | None = None,
    f_panel: float | None = None,
    f_grid: str | Sequence[float] | None = DEFAULT_F_GRID,
    thresholds: str | Sequence[float] | None = DEFAULT_P_THRESHOLDS,
) -> None:
    """Sweep wheat excess-het tail thresholds on a PLINK2 gcount table.

    Writes retain/drop tables and plots. Does not emit a filtered pfile.
    """
    logger.info("Excess-het probability sweep: %s", input_file)
    df_raw = load_df_from_plink_gcount(input_file)
    f_values = list(_parse_float_list(f_grid, DEFAULT_F_GRID))
    p_values = _parse_float_list(thresholds, DEFAULT_P_THRESHOLDS)

    f_data: float | None = None
    f_source_data = "het_median"
    if het_path:
        f_data = load_median_inbreeding_from_het(het_path)
        logger.info("F_panel median from %s: %s", het_path, f_data)
    elif f_panel is not None:
        f_data = float(f_panel)
        f_source_data = "supplied"
    if f_data is None:
        fail("Provide het_path (PLINK2 --het) or f_panel to run the excess-het sweep")

    preview = _compute_mac_table(df_raw)
    n_tested_preview = int((pd.to_numeric(preview["MAC"], errors="coerce") > 0).sum())
    p_list = list(p_values)
    if n_tested_preview:
        bonf = 0.05 / n_tested_preview
        if all(abs(bonf - t) > 0 for t in p_list):
            p_list.append(bonf)
    p_values = tuple(sorted(set(p_list)))

    models: list[tuple[str, float]] = [(f_source_data, f_data)]
    for f_g in f_values:
        if abs(f_g - f_data) > 1e-9:
            models.append(("grid", float(f_g)))

    sweep_frames = []
    mac_frames = []
    sites_data = None
    for source, f_val in models:
        sites = _annotate_excess_het(preview, f_val)
        if source == f_source_data:
            sites_data = sites
        sweep_frames.append(sweep_excess_het_thresholds(sites, f_val, p_values, f_source=source))
        mac_frames.append(sweep_excess_het_by_mac(sites, f_val, p_values, f_source=source))

    sweep = pd.concat(sweep_frames, ignore_index=True)
    by_mac = pd.concat(mac_frames, ignore_index=True)
    save_df_to_tsv(sweep, f"{output_prefix}.sweep.tsv", float_format=_FLOAT_TSV)
    save_df_to_tsv(by_mac, f"{output_prefix}.sweep_by_mac.tsv", float_format=_FLOAT_TSV)

    n_median = float(np.nanmedian(pd.to_numeric(sites_data["Total_Samples"], errors="coerce")))
    allowed = max_allowed_het_table(f_data, alpha=1e-6, n_diploid=n_median if np.isfinite(n_median) else None)
    save_df_to_tsv(allowed, f"{output_prefix}.max_allowed_het.1e-6.tsv", float_format=_FLOAT_TSV)

    flag_summary = _all_het_flag_summary(sites_data)
    save_df_to_tsv(flag_summary, f"{output_prefix}.all_het_flag.info.tsv", float_format=_FLOAT_TSV)

    n_tested = int((pd.to_numeric(sites_data["MAC"], errors="coerce") > 0).sum())
    n_all_het = int(sites_data["all_het_low_mac"].sum())
    stats = {
        "F_panel": f_data,
        "F_source": f_source_data,
        "N_median": n_median,
        "N_tested_MAC_gt_0": n_tested,
        "N_all_het_MAC_le_10": n_all_het,
        "Frac_all_het_MAC_le_10": (n_all_het / n_tested) if n_tested else 0.0,
        "Bonferroni_0.05": (0.05 / n_tested) if n_tested else np.nan,
    }
    save_thresholds(stats, f"{output_prefix}.th.tsv")
    print(stats)
    print(sweep.to_string(index=False))

    _plot_retain_sweep(sweep, output_prefix)
    if not flag_summary.empty:
        plot_bar_chart(
            names=[str(x) for x in flag_summary["MAC"].tolist()],
            values=flag_summary["n_all_het"].astype(float).tolist(),
            title="Sites with all minor alleles in heterozygotes (MAC 1-10; flagged, not deleted)",
            ylabel="Site count",
            filename=f"{output_prefix}.all_het_flag.png",
            ylim=None,
        )
    plot_bar_chart(
        names=allowed["MAC"].astype(str).tolist()[:: max(1, len(allowed) // 20)],
        values=allowed["max_allowed_H"].astype(float).tolist()[:: max(1, len(allowed) // 20)],
        title=f"Max allowed heterozygotes at P_excess>=1e-6 (F={f_data:.4g})",
        ylabel="Max H retained",
        filename=f"{output_prefix}.max_allowed_het.1e-6.png",
        ylim=None,
        rotate_xlabels=45,
    )
    logger.info("Wrote excess-het sweep under %s (no filtered pfile)", output_prefix)


def combine_mac_het_excess_sweeps(
    sweep_files: Sequence[str],
    output_prefix: str = "all.variant.mac_het",
) -> None:
    """Sum per-subgenome retain/drop sweeps into a panel-wide table and plot."""
    if not sweep_files:
        fail("No sweep TSV files to combine")
    frames = [pd.read_csv(path, sep="\t") for path in sweep_files]
    df = pd.concat(frames, ignore_index=True)
    grouped = (
        df.groupby(["f_source", "f_panel", "p_threshold"], as_index=False)
        .agg(
            n_tested=("n_tested", "sum"),
            n_retain=("n_retain", "sum"),
            n_drop=("n_drop", "sum"),
            n_all_het_mac_le_10=("n_all_het_mac_le_10", "sum"),
        )
        .sort_values(["f_source", "f_panel", "p_threshold"])
        .reset_index(drop=True)
    )
    grouped["frac_retain"] = grouped["n_retain"] / grouped["n_tested"]
    save_df_to_tsv(grouped, f"{output_prefix}.sweep.tsv", float_format=_FLOAT_TSV)
    _plot_retain_sweep(grouped, output_prefix)
    print(grouped.to_string(index=False))
