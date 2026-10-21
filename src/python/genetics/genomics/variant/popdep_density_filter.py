"""RelativeDepth mean-vs-SD grid density filter for PopDep trusted sites.

Stream TIGER popdepth grids into a 2D histogram of RelativeDepth_Mean vs
RelativeDepth_SD. Linear-range filters keep the densest 60% / 80% of
non-empty cells. Log-range filters keep densest non-empty cells until
60% / 80% / 90% of counted sites remain, and report the implied cell-top
fraction. A separate linear site-coverage filter uses RelativeDepth mean
and SD on [0.02, 6] (no log) with site-coverage keep at 60 / 80 / 90 /
95 / 97 / 98 / 99%. Mean-vs-CV filters RelMean on [0, 6] with observed
CV min/max (no CV clip) and keep 60 / 70 / 80 / 85 / 90 / 95 / 97 / 98 / 99%.
The n500 reliable-site libraries use RelMean [0, 6] × RelCV [0, 25] on a
500×500 grid (250_000 cells) and keep 60% / 80% of in-range sites.
Points outside the histogram range are not trusted.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from genetics.genomics.variant.popdep import _parse_popdep_grid_row, _popdep_grid_column_map
from genetics.genomics.variant.popdep_region_class import (
    REGION_CLASS_ALPHA,
    REGION_CLASS_COLORS,
    REGION_CLASS_HUE_ORDER,
    REGION_CLASS_LEGEND_LABELS,
)
from infra.utils.graph import (
    SEQUENTIAL_WHITE_RED,
    TWO_CATEGORY_BLUE,
    plot_heatmap_custom,
    plot_joint_regression,
    plot_scatter_with_thresholds,
)
from infra.utils.io import load_df_from_tsv, save_df_to_tsv, save_thresholds

KEEP_FRAC_STRICT = 0.6
KEEP_FRAC_LOOSE = 0.8
KEEP_FRAC_WIDE = 0.9
LOG_SITE_KEEP_FRACS = (KEEP_FRAC_STRICT, KEEP_FRAC_LOOSE, KEEP_FRAC_WIDE)
DEFAULT_BINS_X = 250
DEFAULT_BINS_Y = 150
DEFAULT_X_RANGE = (0.0, 5.0)
DEFAULT_Y_RANGE = (0.0, 5.0)
LINEAR_RELDEPTH_PLOT_RANGE = (0.0, 10.0)
LINEAR_SITE_RANGE = (0.02, 6.0)
LINEAR_SITE_KEEP_FRACS = (
    KEEP_FRAC_STRICT,
    KEEP_FRAC_LOOSE,
    KEEP_FRAC_WIDE,
    0.95,
    0.97,
    0.98,
    0.99,
)
MEAN_RANGE_0_6 = (0.0, 6.0)
CV_RANGE_0_25 = (0.0, 25.0)
CV_RELIABLE_BINS = 500
CV_RELIABLE_N_CELLS = CV_RELIABLE_BINS * CV_RELIABLE_BINS
CV_RELIABLE_KEEP_FRACS = (KEEP_FRAC_STRICT, KEEP_FRAC_LOOSE)
CV_SITE_KEEP_FRACS = (
    KEEP_FRAC_STRICT,
    0.7,
    KEEP_FRAC_LOOSE,
    0.85,
    KEEP_FRAC_WIDE,
    0.95,
    0.97,
    0.98,
    0.99,
)
HIST_BATCH = 1_000_000
LOG_N_CELLS = 90_000
LOG_COL_X = "_log_rel_mean"
LOG_COL_Y = "_log_rel_sd"
SCATTER_CAP = 2_000_000
X_COL = "RelativeDepth_Mean"
Y_COL = "RelativeDepth_SD"
CV_COL = "RelativeDepth_CV"
DENSITY_KEEP_LABEL = "Kept"
DENSITY_DROP_LABEL = "Dropped"
DENSITY_DROP_COLOR = "#bdbdbd"
DENSITY_FILTER_HUE_ORDER = (DENSITY_DROP_LABEL, DENSITY_KEEP_LABEL)
DENSITY_FILTER_PALETTE = {
    DENSITY_DROP_LABEL: DENSITY_DROP_COLOR,
    DENSITY_KEEP_LABEL: TWO_CATEGORY_BLUE,
}
DENSITY_FILTER_HUE_ALPHA = {
    DENSITY_DROP_LABEL: 0.2,
    DENSITY_KEEP_LABEL: 0.2,
}
_MEAN_SD_SCATTERS = (
    ("Depth_Mean", "Depth_SD", "Mean depth", "Depth SD"),
    ("RelativeDepth_Mean", "RelativeDepth_SD", "Relative mean depth", "Relative depth SD"),
)
_RELDEPTH_SCATTERS = (_MEAN_SD_SCATTERS[1],)
_RELDEPTH_CV_SCATTERS = (
    (X_COL, CV_COL, "Relative mean depth", "Relative depth CV"),
)


def keep_frac_pct(frac: float) -> int:
    """Percent integer for a keep fraction (0.95 → 95)."""
    return int(round(100.0 * float(frac)))


def _bin_edges(limit: tuple[float, float], n_bins: int) -> np.ndarray:
    lo, hi = float(limit[0]), float(limit[1])
    if hi <= lo:
        raise ValueError(f"range must be increasing, got {limit}")
    return np.linspace(lo, hi, int(n_bins) + 1)


def histogram_rel_mean_sd(
    x,
    y,
    *,
    x_range: tuple[float, float] = DEFAULT_X_RANGE,
    y_range: tuple[float, float] = DEFAULT_Y_RANGE,
    bins_x: int = DEFAULT_BINS_X,
    bins_y: int = DEFAULT_BINS_Y,
) -> tuple[np.ndarray, int, int]:
    """Return ``(counts[ix, iy], n_valid, n_oob)``; NaN rows are ignored."""
    xs = np.asarray(x, dtype=float)
    ys = np.asarray(y, dtype=float)
    finite = np.isfinite(xs) & np.isfinite(ys)
    xs = xs[finite]
    ys = ys[finite]
    x_edges = _bin_edges(x_range, bins_x)
    y_edges = _bin_edges(y_range, bins_y)
    inside = (
        (xs >= x_edges[0])
        & (xs < x_edges[-1])
        & (ys >= y_edges[0])
        & (ys < y_edges[-1])
    )
    n_oob = int((~inside).sum())
    n_valid = int(inside.sum())
    if n_valid == 0:
        return np.zeros((int(bins_x), int(bins_y)), dtype=np.int64), 0, n_oob
    counts, _, _ = np.histogram2d(xs[inside], ys[inside], bins=[x_edges, y_edges])
    return counts.astype(np.int64, copy=False), n_valid, n_oob


def assign_bin_indices(
    x,
    y,
    *,
    x_range: tuple[float, float] = DEFAULT_X_RANGE,
    y_range: tuple[float, float] = DEFAULT_Y_RANGE,
    bins_x: int | None = None,
    bins_y: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return 0-based ``(ix, iy)``; out-of-range or non-finite is ``-1``."""
    xs = np.asarray(x, dtype=float)
    ys = np.asarray(y, dtype=float)
    n_x = int(bins_x) if bins_x is not None else DEFAULT_BINS_X
    n_y = int(bins_y) if bins_y is not None else DEFAULT_BINS_Y
    x_edges = _bin_edges(x_range, n_x)
    y_edges = _bin_edges(y_range, n_y)
    ix = np.digitize(xs, x_edges) - 1
    iy = np.digitize(ys, y_edges) - 1
    ok = (
        np.isfinite(xs)
        & np.isfinite(ys)
        & (ix >= 0)
        & (ix < n_x)
        & (iy >= 0)
        & (iy < n_y)
        & (xs < x_edges[-1])
        & (ys < y_edges[-1])
    )
    ix = ix.astype(np.int64, copy=False)
    iy = iy.astype(np.int64, copy=False)
    ix[~ok] = -1
    iy[~ok] = -1
    return ix, iy


def trusted_mask_from_counts(counts: np.ndarray, keep_frac: float) -> np.ndarray:
    """Keep the densest ``keep_frac`` of non-empty cells (empty cells never pass)."""
    frac = float(keep_frac)
    if not 0.0 < frac <= 1.0:
        raise ValueError(f"keep_frac must be in (0, 1], got {keep_frac}")
    flat = np.asarray(counts, dtype=np.int64).ravel()
    nonempty = np.flatnonzero(flat > 0)
    mask = np.zeros(counts.shape, dtype=bool)
    if nonempty.size == 0:
        return mask
    n_keep = max(1, int(np.ceil(frac * nonempty.size)))
    order = nonempty[np.argsort(flat[nonempty], kind="mergesort")[::-1]]
    mask.ravel()[order[:n_keep]] = True
    return mask


def trusted_mask_from_site_coverage(
    counts: np.ndarray,
    site_frac: float,
) -> tuple[np.ndarray, dict[str, float | int]]:
    """Keep densest non-empty cells until they hold ``site_frac`` of counted sites.

    Empty cells are ignored. Returns the mask and the implied cell-top fraction
    ``n_cells_kept / n_nonempty``.
    """
    frac = float(site_frac)
    if not 0.0 < frac <= 1.0:
        raise ValueError(f"site_frac must be in (0, 1], got {site_frac}")
    flat = np.asarray(counts, dtype=np.int64).ravel()
    nonempty = np.flatnonzero(flat > 0)
    mask = np.zeros(counts.shape, dtype=bool)
    empty_info: dict[str, float | int] = {
        "n_cells_kept": 0,
        "n_nonempty": 0,
        "n_sites_kept": 0,
        "n_sites_total": int(flat.sum()),
        "cell_top_frac": 0.0,
        "site_frac_target": frac,
        "site_frac_kept": 0.0,
    }
    if nonempty.size == 0:
        return mask, empty_info
    order = nonempty[np.argsort(flat[nonempty], kind="mergesort")[::-1]]
    ranked = flat[order]
    total = int(ranked.sum())
    target = max(1, int(np.ceil(frac * total)))
    cum = np.cumsum(ranked)
    n_keep = int(np.searchsorted(cum, target, side="left") + 1)
    n_keep = min(max(n_keep, 1), int(order.size))
    mask.ravel()[order[:n_keep]] = True
    n_sites = int(cum[n_keep - 1])
    info: dict[str, float | int] = {
        "n_cells_kept": n_keep,
        "n_nonempty": int(nonempty.size),
        "n_sites_kept": n_sites,
        "n_sites_total": total,
        "cell_top_frac": n_keep / float(nonempty.size),
        "site_frac_target": frac,
        "site_frac_kept": n_sites / float(total) if total else 0.0,
    }
    return mask, info


def site_coverage_curve(counts: np.ndarray) -> pd.DataFrame:
    """Cumulative site fraction vs nonempty-cell top fraction, densest cells first.

    Row 0 is keep-nothing (0%, 0%). Each later row adds the next densest
    non-empty cell. Empty cells are omitted.
    """
    columns = (
        "n_cells_kept",
        "n_nonempty",
        "cell_top_frac",
        "n_sites_kept",
        "n_sites_total",
        "site_frac_kept",
        "site_pct",
        "cell_top_pct",
    )
    flat = np.asarray(counts, dtype=np.int64).ravel()
    nonempty = np.flatnonzero(flat > 0)
    n_nonempty = int(nonempty.size)
    total = int(flat.sum())
    origin = {
        "n_cells_kept": 0,
        "n_nonempty": n_nonempty,
        "cell_top_frac": 0.0,
        "n_sites_kept": 0,
        "n_sites_total": total,
        "site_frac_kept": 0.0,
        "site_pct": 0.0,
        "cell_top_pct": 0.0,
    }
    if n_nonempty == 0:
        return pd.DataFrame([origin], columns=list(columns))
    order = nonempty[np.argsort(flat[nonempty], kind="mergesort")[::-1]]
    ranked = flat[order]
    n_cells = np.arange(1, n_nonempty + 1, dtype=np.int64)
    n_sites = np.cumsum(ranked)
    site_frac = n_sites / float(total) if total else np.zeros(n_nonempty, dtype=float)
    cell_frac = n_cells / float(n_nonempty)
    body = pd.DataFrame(
        {
            "n_cells_kept": n_cells,
            "n_nonempty": n_nonempty,
            "cell_top_frac": cell_frac,
            "n_sites_kept": n_sites.astype(np.int64, copy=False),
            "n_sites_total": total,
            "site_frac_kept": site_frac,
            "site_pct": 100.0 * site_frac,
            "cell_top_pct": 100.0 * cell_frac,
        }
    )
    return pd.concat([pd.DataFrame([origin]), body], ignore_index=True)


_SITE_MARK_COLORS = {
    KEEP_FRAC_STRICT: "#c44e52",
    KEEP_FRAC_LOOSE: "#4c72b0",
    0.7: "#64b5cd",
    0.85: "#cc78bc",
    KEEP_FRAC_WIDE: "#55a868",
    0.95: "#dd8452",
    0.97: "#8172b3",
    0.98: "#937860",
    0.99: "#4c4c4c",
}


def plot_site_coverage_curve(
    curve: pd.DataFrame,
    filename: str | Path,
    *,
    log_y: bool = False,
    site_marks: tuple[float, ...] = LOG_SITE_KEEP_FRACS,
) -> None:
    """Line of retained-site % vs the implied nonempty-cell top %."""
    frame = curve.copy()
    if log_y:
        frame = frame.loc[frame["cell_top_pct"] > 0].copy()
    vlines: list[dict[str, float | str]] = []
    hlines: list[dict[str, float | str]] = []
    for frac in site_marks:
        site_pct = 100.0 * float(frac)
        color = _SITE_MARK_COLORS.get(float(frac), "#7f7f7f")
        label = f"sites {site_pct:.0f}%"
        vlines.append(
            {"value": site_pct, "color": color, "linestyle": "--", "label": label}
        )
        hit = frame.loc[frame["site_pct"] >= site_pct]
        if hit.empty:
            continue
        cell_pct = float(hit.iloc[0]["cell_top_pct"])
        hlines.append(
            {
                "value": cell_pct,
                "color": color,
                "linestyle": ":",
                "label": f"cells top {cell_pct:.2f}%",
            }
        )
    plot_scatter_with_thresholds(
        frame,
        x_col="site_pct",
        y_col="cell_top_pct",
        title="Nonempty-cell top % vs retained PopDep sites",
        filename=str(filename),
        xlabel="Retained sites (%)",
        ylabel="Corresponding nonempty-cell top (%)",
        color=TWO_CATEGORY_BLUE,
        connect=True,
        s=1,
        alpha=0.0,
        xlim=(0.0, 100.0),
        ylim=None if log_y else (0.0, 100.0),
        thresholds_v=vlines or None,
        thresholds_h=hlines or None,
        log_y=log_y,
        figure_size=(10, 6),
    )


def sites_in_trusted_cells(
    sites: pd.DataFrame,
    mask: np.ndarray,
    *,
    x_col: str = X_COL,
    y_col: str = Y_COL,
    x_range: tuple[float, float] = DEFAULT_X_RANGE,
    y_range: tuple[float, float] = DEFAULT_Y_RANGE,
) -> pd.DataFrame:
    """Return rows whose RelativeDepth cell is True in ``mask``."""
    keep = flag_sites_in_trusted_cells(
        sites,
        mask,
        x_col=x_col,
        y_col=y_col,
        x_range=x_range,
        y_range=y_range,
    )
    return sites.loc[keep].copy()


def _open_popdep_text(path: str | Path):
    path_s = str(path)
    if path_s.endswith((".gz", ".bgz")):
        return gzip.open(path_s, "rt", encoding="utf-8")
    return open(path_s, "rt", encoding="utf-8")


def accumulate_rel_mean_sd_histogram(
    path: str | Path,
    *,
    x_range: tuple[float, float] = DEFAULT_X_RANGE,
    y_range: tuple[float, float] = DEFAULT_Y_RANGE,
    bins_x: int = DEFAULT_BINS_X,
    bins_y: int = DEFAULT_BINS_Y,
    batch_size: int = HIST_BATCH,
) -> tuple[np.ndarray, dict[str, int]]:
    """Stream one TIGER popdepth table into a RelativeDepth mean-vs-SD histogram."""
    counts = np.zeros((int(bins_x), int(bins_y)), dtype=np.int64)
    n_rows = 0
    n_valid = 0
    n_oob = 0
    n_nan = 0
    xs: list[float] = []
    ys: list[float] = []

    def _flush() -> None:
        nonlocal n_valid, n_oob, xs, ys
        if not xs:
            return
        part, n_ok, n_out = histogram_rel_mean_sd(
            xs, ys, x_range=x_range, y_range=y_range, bins_x=bins_x, bins_y=bins_y
        )
        counts.__iadd__(part)
        n_valid += n_ok
        n_oob += n_out
        xs = []
        ys = []

    with _open_popdep_text(path) as handle:
        header = handle.readline()
        if not header:
            return counts, {"n_rows": 0, "n_valid": 0, "n_oob": 0, "n_nan": 0}
        col_map = _popdep_grid_column_map(header.rstrip("\n").split("\t"))
        for line in handle:
            n_rows += 1
            parsed = _parse_popdep_grid_row(line.rstrip("\n").split("\t"), col_map)
            if parsed is None:
                n_nan += 1
                continue
            _pos, _mean, _sd, rel_mean, rel_sd = parsed
            if not (np.isfinite(rel_mean) and np.isfinite(rel_sd)):
                n_nan += 1
                continue
            xs.append(rel_mean)
            ys.append(rel_sd)
            if len(xs) >= int(batch_size):
                _flush()
        _flush()
    return counts, {
        "n_rows": int(n_rows),
        "n_valid": int(n_valid),
        "n_oob": int(n_oob),
        "n_nan": int(n_nan),
    }


def merge_histograms(parts: list[np.ndarray]) -> np.ndarray:
    if not parts:
        raise ValueError("merge_histograms: no arrays to add")
    out = np.zeros_like(parts[0], dtype=np.int64)
    for part in parts:
        if part.shape != out.shape:
            raise ValueError(f"histogram shape mismatch: {part.shape} vs {out.shape}")
        out += np.asarray(part, dtype=np.int64)
    return out


def save_histogram_npz(
    path: str | Path,
    counts: np.ndarray,
    stats: dict[str, int],
    *,
    x_range: tuple[float, float] = DEFAULT_X_RANGE,
    y_range: tuple[float, float] = DEFAULT_Y_RANGE,
) -> None:
    np.savez_compressed(
        path,
        counts=np.asarray(counts, dtype=np.int64),
        n_rows=np.int64(stats.get("n_rows", 0)),
        n_valid=np.int64(stats.get("n_valid", 0)),
        n_oob=np.int64(stats.get("n_oob", 0)),
        n_nan=np.int64(stats.get("n_nan", 0)),
        x_range=np.asarray(x_range, dtype=float),
        y_range=np.asarray(y_range, dtype=float),
    )


def load_histogram_npz(path: str | Path) -> tuple[np.ndarray, dict[str, int], tuple[float, float], tuple[float, float]]:
    with np.load(path) as handle:
        counts = np.asarray(handle["counts"], dtype=np.int64)
        stats = {
            "n_rows": int(handle["n_rows"]),
            "n_valid": int(handle["n_valid"]),
            "n_oob": int(handle["n_oob"]),
            "n_nan": int(handle["n_nan"]),
        }
        x_range = (float(handle["x_range"][0]), float(handle["x_range"][1]))
        y_range = (float(handle["y_range"][0]), float(handle["y_range"][1]))
    return counts, stats, x_range, y_range


def write_chr_histogram(
    popdep_path: str | Path,
    output_npz: str | Path,
    *,
    x_range: tuple[float, float] = DEFAULT_X_RANGE,
    y_range: tuple[float, float] = DEFAULT_Y_RANGE,
    bins_x: int = DEFAULT_BINS_X,
    bins_y: int = DEFAULT_BINS_Y,
) -> dict[str, int]:
    counts, stats = accumulate_rel_mean_sd_histogram(
        popdep_path, x_range=x_range, y_range=y_range, bins_x=bins_x, bins_y=bins_y
    )
    save_histogram_npz(output_npz, counts, stats, x_range=x_range, y_range=y_range)
    print(
        f"[Info] {popdep_path}: rows={stats['n_rows']} valid={stats['n_valid']} "
        f"oob={stats['n_oob']} nan={stats['n_nan']} -> {output_npz}"
    )
    return stats


def reliable_cv_hist_ranges() -> tuple[tuple[float, float], tuple[float, float]]:
    """Half-open RelMean [0, 6] and RelCV [0, 25] that still contain the endpoints."""
    return (
        expand_range_to_include_max(*MEAN_RANGE_0_6),
        expand_range_to_include_max(*CV_RANGE_0_25),
    )


def accumulate_rel_mean_cv_histogram_grid(
    path: str | Path,
    *,
    x_range: tuple[float, float] | None = None,
    y_range: tuple[float, float] | None = None,
    bins_x: int = CV_RELIABLE_BINS,
    bins_y: int = CV_RELIABLE_BINS,
    batch_size: int = HIST_BATCH,
) -> tuple[np.ndarray, dict[str, int]]:
    """Stream one TIGER popdepth table into a RelMean vs RelCV histogram."""
    hist_x, hist_y = reliable_cv_hist_ranges()
    if x_range is not None:
        hist_x = (float(x_range[0]), float(x_range[1]))
    if y_range is not None:
        hist_y = (float(y_range[0]), float(y_range[1]))
    counts = np.zeros((int(bins_x), int(bins_y)), dtype=np.int64)
    n_rows = 0
    n_valid = 0
    n_oob = 0
    n_nan = 0
    xs: list[float] = []
    ys: list[float] = []

    def _flush() -> None:
        nonlocal n_valid, n_oob, xs, ys
        if not xs:
            return
        part, n_ok, n_out = histogram_rel_mean_sd(
            xs, ys, x_range=hist_x, y_range=hist_y, bins_x=bins_x, bins_y=bins_y
        )
        counts.__iadd__(part)
        n_valid += n_ok
        n_oob += n_out
        xs = []
        ys = []

    with _open_popdep_text(path) as handle:
        header = handle.readline()
        if not header:
            return counts, {"n_rows": 0, "n_valid": 0, "n_oob": 0, "n_nan": 0}
        col_map = _popdep_grid_column_map(header.rstrip("\n").split("\t"))
        for line in handle:
            n_rows += 1
            parsed = _parse_popdep_grid_row(line.rstrip("\n").split("\t"), col_map)
            if parsed is None:
                n_nan += 1
                continue
            _pos, _mean, _sd, rel_mean, rel_sd = parsed
            if not (
                np.isfinite(rel_mean)
                and np.isfinite(rel_sd)
                and float(rel_mean) > 0.0
            ):
                n_nan += 1
                continue
            xs.append(float(rel_mean))
            ys.append(float(rel_sd) / float(rel_mean))
            if len(xs) >= int(batch_size):
                _flush()
        _flush()
    return counts, {
        "n_rows": int(n_rows),
        "n_valid": int(n_valid),
        "n_oob": int(n_oob),
        "n_nan": int(n_nan),
    }


def write_chr_mean_cv_histogram(
    popdep_path: str | Path,
    output_npz: str | Path,
    *,
    x_range: tuple[float, float] | None = None,
    y_range: tuple[float, float] | None = None,
    bins_x: int = CV_RELIABLE_BINS,
    bins_y: int = CV_RELIABLE_BINS,
) -> dict[str, int]:
    """Write one chromosome RelMean vs RelCV histogram npz (n500 500×500 default)."""
    hist_x, hist_y = reliable_cv_hist_ranges()
    if x_range is not None:
        hist_x = (float(x_range[0]), float(x_range[1]))
    if y_range is not None:
        hist_y = (float(y_range[0]), float(y_range[1]))
    counts, stats = accumulate_rel_mean_cv_histogram_grid(
        popdep_path, x_range=hist_x, y_range=hist_y, bins_x=bins_x, bins_y=bins_y
    )
    save_histogram_npz(output_npz, counts, stats, x_range=hist_x, y_range=hist_y)
    print(
        f"[Info] {popdep_path}: RelMean-CV rows={stats['n_rows']} valid={stats['n_valid']} "
        f"oob={stats['n_oob']} nan={stats['n_nan']} -> {output_npz}"
    )
    return stats


def write_cv_reliable_density_outputs(
    counts: np.ndarray,
    stats: dict[str, int],
    output_prefix: str | Path,
    *,
    x_range: tuple[float, float] | None = None,
    y_range: tuple[float, float] | None = None,
) -> dict[str, np.ndarray | dict[str, float | int]]:
    """Site-coverage keep 60/80 on a RelMean vs RelCV histogram."""
    hist_x, hist_y = reliable_cv_hist_ranges()
    if x_range is not None:
        hist_x = (float(x_range[0]), float(x_range[1]))
    if y_range is not None:
        hist_y = (float(y_range[0]), float(y_range[1]))
    return write_linear_site_density_outputs(
        counts,
        stats,
        output_prefix,
        x_range=hist_x,
        y_range=hist_y,
        keep_fracs=CV_RELIABLE_KEEP_FRACS,
        heatmap_stem="rel_mean_cv",
        heatmap_ylabel="RelativeDepth_CV",
        heatmap_title="PopDep RelMean vs RelCV cell density (n500, 500×500)",
        keep_scale="linear_cv_reliable",
    )


class _BedRunWriter:
    """Merge consecutive 1-based positions into 0-based half-open BED rows."""

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        if str(self.path).endswith(".gz"):
            self.handle = gzip.open(self.path, "wt", encoding="utf-8")
        else:
            self.handle = self.path.open("w", encoding="utf-8")
        self.chrom: str | None = None
        self.start: int | None = None
        self.end: int | None = None
        self.n_bp = 0
        self.n_intervals = 0

    def add(self, chrom: str, pos_1based: int) -> None:
        start = int(pos_1based) - 1
        end = int(pos_1based)
        self.n_bp += 1
        if self.chrom == chrom and self.end is not None and start <= self.end:
            if end > self.end:
                self.end = end
            return
        self.flush()
        self.chrom = str(chrom)
        self.start = start
        self.end = end

    def flush(self) -> None:
        if self.chrom is None or self.start is None or self.end is None:
            return
        self.handle.write(f"{self.chrom}\t{self.start}\t{self.end}\n")
        self.n_intervals += 1
        self.chrom = None
        self.start = None
        self.end = None

    def close(self) -> dict[str, int]:
        self.flush()
        self.handle.close()
        return {"n_bp": int(self.n_bp), "n_intervals": int(self.n_intervals)}


def write_chr_cv_reliable_beds(
    popdep_path: str | Path,
    masks_npz: str | Path,
    output_prefix: str | Path,
    *,
    keep_pcts: tuple[int, ...] = (60, 80),
    batch_size: int = HIST_BATCH,
) -> dict[str, int]:
    """Write per-chromosome BED libraries of PopDep sites in keep-60/80 cells."""
    packed = np.load(masks_npz)
    x_range = (float(packed["x_range"][0]), float(packed["x_range"][1]))
    y_range = (float(packed["y_range"][0]), float(packed["y_range"][1]))
    masks = {int(pct): np.asarray(packed[f"mask_{int(pct)}"], dtype=bool) for pct in keep_pcts}
    bins_x, bins_y = next(iter(masks.values())).shape
    prefix = Path(output_prefix)
    writers = {
        int(pct): _BedRunWriter(Path(f"{prefix}.keep{int(pct)}.bed.gz")) for pct in keep_pcts
    }
    n_rows = 0
    n_valid = 0
    n_oob = 0
    n_nan = 0
    chroms: list[str] = []
    positions: list[int] = []
    means: list[float] = []
    cvs: list[float] = []

    def _flush_batch() -> None:
        nonlocal n_valid, n_oob, chroms, positions, means, cvs
        if not means:
            return
        ix, iy = assign_bin_indices(
            means,
            cvs,
            x_range=x_range,
            y_range=y_range,
            bins_x=bins_x,
            bins_y=bins_y,
        )
        for chrom, pos, x_i, y_i in zip(chroms, positions, ix, iy, strict=True):
            if int(x_i) < 0 or int(y_i) < 0:
                n_oob += 1
                continue
            n_valid += 1
            for pct, mask in masks.items():
                if bool(mask[int(x_i), int(y_i)]):
                    writers[pct].add(chrom, pos)
        chroms = []
        positions = []
        means = []
        cvs = []

    with _open_popdep_text(popdep_path) as handle:
        header = handle.readline()
        if not header:
            for writer in writers.values():
                writer.close()
            return {"n_rows": 0, "n_valid": 0, "n_oob": 0, "n_nan": 0}
        col_map = _popdep_grid_column_map(header.rstrip("\n").split("\t"))
        for line in handle:
            n_rows += 1
            parts = line.rstrip("\n").split("\t")
            parsed = _parse_popdep_grid_row(parts, col_map)
            if parsed is None:
                n_nan += 1
                continue
            pos, _mean, _sd, rel_mean, rel_sd = parsed
            if not (
                np.isfinite(rel_mean)
                and np.isfinite(rel_sd)
                and float(rel_mean) > 0.0
            ):
                n_nan += 1
                continue
            chroms.append(parts[0])
            positions.append(int(pos))
            means.append(float(rel_mean))
            cvs.append(float(rel_sd) / float(rel_mean))
            if len(means) >= int(batch_size):
                _flush_batch()
        _flush_batch()
    stats = {
        "n_rows": int(n_rows),
        "n_valid": int(n_valid),
        "n_oob": int(n_oob),
        "n_nan": int(n_nan),
    }
    for pct, writer in writers.items():
        closed = writer.close()
        stats[f"n_kept{pct}"] = int(closed["n_bp"])
        stats[f"n_intervals{pct}"] = int(closed["n_intervals"])
    save_df_to_tsv(pd.DataFrame([stats]), f"{prefix}.lib.info.tsv")
    print(
        f"[Info] {popdep_path}: reliable beds rows={stats['n_rows']} "
        f"kept60={stats.get('n_kept60', 0)} kept80={stats.get('n_kept80', 0)} -> {prefix}",
        flush=True,
    )
    return stats


def concat_gzip_beds(bed_paths, output_bed: str | Path) -> int:
    """Concatenate chromosome BED libraries; writes gzip when the path ends with .gz."""
    out_path = Path(output_bed)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    paths = sorted((Path(p) for p in bed_paths), key=lambda path: path.name)
    n_lines = 0
    opener = gzip.open if str(out_path).endswith(".gz") else open
    with opener(out_path, "wt", encoding="utf-8") as handle:
        for path in paths:
            reader = gzip.open if str(path).endswith(".gz") else open
            with reader(path, "rt", encoding="utf-8") as src:
                for line in src:
                    handle.write(line)
                    n_lines += 1
    print(f"[Info] concatenated {len(paths)} BED files ({n_lines} intervals) -> {out_path}")
    return n_lines


def merge_lib_info_tables(info_paths, output_tsv: str | Path) -> pd.DataFrame:
    """Stack per-chromosome reliable-library count tables."""
    frames = [load_df_from_tsv(str(path)) for path in sorted(info_paths, key=lambda p: Path(p).name)]
    if not frames:
        raise ValueError("merge_lib_info_tables: no lib.info tables")
    out = pd.concat(frames, ignore_index=True)
    save_df_to_tsv(out, str(output_tsv))
    return out


def merge_histogram_files(npz_paths, output_npz: str | Path) -> tuple[np.ndarray, dict[str, int]]:
    parts: list[np.ndarray] = []
    totals = {"n_rows": 0, "n_valid": 0, "n_oob": 0, "n_nan": 0}
    x_range = DEFAULT_X_RANGE
    y_range = DEFAULT_Y_RANGE
    for path in npz_paths:
        counts, stats, x_range, y_range = load_histogram_npz(path)
        parts.append(counts)
        for key in totals:
            totals[key] += int(stats[key])
    merged = merge_histograms(parts)
    save_histogram_npz(output_npz, merged, totals, x_range=x_range, y_range=y_range)
    return merged, totals


def summarize_trusted_cells(
    counts: np.ndarray,
    stats: dict[str, int],
    *,
    keep_fracs: tuple[float, ...] = (KEEP_FRAC_STRICT, KEEP_FRAC_LOOSE),
) -> pd.DataFrame:
    rows = [
        {
            "set": "all_grid_rows",
            "n_cells": int(counts.size),
            "n_sites": int(stats["n_rows"]),
            "frac_sites": 1.0,
        }
    ]
    n_valid = max(int(stats["n_valid"]), 1)
    for frac in keep_fracs:
        mask = trusted_mask_from_counts(counts, frac)
        n_sites = int(counts[mask].sum())
        rows.append(
            {
                "set": f"top_{int(round(frac * 100))}pct_cells",
                "n_cells": int(mask.sum()),
                "n_sites": n_sites,
                "frac_sites": n_sites / n_valid,
            }
        )
    return pd.DataFrame(rows)


def plot_density_grid(
    counts: np.ndarray,
    filename: str,
    *,
    x_range: tuple[float, float] = DEFAULT_X_RANGE,
    y_range: tuple[float, float] = DEFAULT_Y_RANGE,
    mask_60: np.ndarray | None = None,
    mask_80: np.ndarray | None = None,
) -> None:
    """Heatmap of log1p cell counts; y increases upward in the matrix top-to-bottom flip."""
    nx, ny = counts.shape
    x_edges = _bin_edges(x_range, nx)
    y_edges = _bin_edges(y_range, ny)
    matrix = np.log1p(counts.astype(float)).T[::-1]
    x_labels = [f"{0.5 * (x_edges[i] + x_edges[i + 1]):.2f}" for i in range(nx)]
    y_labels = [f"{0.5 * (y_edges[j] + y_edges[j + 1]):.2f}" for j in range(ny - 1, -1, -1)]
    overlays = []
    if mask_80 is not None:
        ys, xs = _mask_boundary_heatmap_coords(mask_80.T[::-1])
        if xs:
            overlays.append(
                {
                    "x": xs,
                    "y": ys,
                    "color": "#4c72b0",
                    "linestyle": "none",
                    "marker": "s",
                    "markersize": 1.2,
                    "label": "80% cells",
                }
            )
    if mask_60 is not None:
        ys, xs = _mask_boundary_heatmap_coords(mask_60.T[::-1])
        if xs:
            overlays.append(
                {
                    "x": xs,
                    "y": ys,
                    "color": "#c44e52",
                    "linestyle": "none",
                    "marker": "s",
                    "markersize": 1.2,
                    "label": "60% cells",
                }
            )
    plot_heatmap_custom(
        matrix,
        x_labels,
        y_labels,
        title="PopDep RelativeDepth mean vs SD cell density (log1p count)",
        filename=filename,
        xlabel="RelativeDepth_Mean",
        ylabel="RelativeDepth_SD",
        cmap=SEQUENTIAL_WHITE_RED,
        cbar_label="log1p(sites / cell)",
        overlays=overlays or None,
    )


def _mask_boundary_heatmap_coords(mask: np.ndarray) -> tuple[list[float], list[float]]:
    """Cheap outline: one point at the centre of every kept heatmap cell."""
    rows, cols = np.nonzero(mask)
    return (rows + 0.5).tolist(), (cols + 0.5).tolist()


def write_density_outputs(
    counts: np.ndarray,
    stats: dict[str, int],
    output_prefix: str | Path,
    *,
    x_range: tuple[float, float] = DEFAULT_X_RANGE,
    y_range: tuple[float, float] = DEFAULT_Y_RANGE,
) -> dict[str, np.ndarray]:
    output_prefix = Path(output_prefix)
    mask_60 = trusted_mask_from_counts(counts, KEEP_FRAC_STRICT)
    mask_80 = trusted_mask_from_counts(counts, KEEP_FRAC_LOOSE)
    summary = summarize_trusted_cells(counts, stats)
    save_df_to_tsv(summary, f"{output_prefix}.trusted_cells.tsv")
    save_histogram_npz(f"{output_prefix}.counts.npz", counts, stats, x_range=x_range, y_range=y_range)
    np.savez_compressed(f"{output_prefix}.masks.npz", mask_60=mask_60, mask_80=mask_80)
    plot_density_grid(
        counts,
        f"{output_prefix}.rel_mean_sd.density.png",
        x_range=x_range,
        y_range=y_range,
        mask_60=mask_60,
        mask_80=mask_80,
    )
    save_thresholds(
        {
            "n_rows": int(stats["n_rows"]),
            "n_valid": int(stats["n_valid"]),
            "n_oob": int(stats["n_oob"]),
            "n_nan": int(stats["n_nan"]),
            "n_sites_top60": int(counts[mask_60].sum()),
            "n_sites_top80": int(counts[mask_80].sum()),
            "n_cells_top60": int(mask_60.sum()),
            "n_cells_top80": int(mask_80.sum()),
            "keep_frac_strict": KEEP_FRAC_STRICT,
            "keep_frac_loose": KEEP_FRAC_LOOSE,
            "bins_x": int(counts.shape[0]),
            "bins_y": int(counts.shape[1]),
            "x_min": float(x_range[0]),
            "x_max": float(x_range[1]),
            "y_min": float(y_range[0]),
            "y_max": float(y_range[1]),
        },
        f"{output_prefix}.th.tsv",
    )
    print(f"[Info] Wrote density outputs under {output_prefix}")
    return {"mask_60": mask_60, "mask_80": mask_80}


def flag_sites_in_trusted_cells(
    sites: pd.DataFrame,
    mask: np.ndarray,
    *,
    x_col: str = X_COL,
    y_col: str = Y_COL,
    x_range: tuple[float, float] = DEFAULT_X_RANGE,
    y_range: tuple[float, float] = DEFAULT_Y_RANGE,
) -> np.ndarray:
    """Boolean keep flags; out-of-range or non-finite RelativeDepth is False."""
    bins_x, bins_y = mask.shape
    ix, iy = assign_bin_indices(
        pd.to_numeric(sites[x_col], errors="coerce"),
        pd.to_numeric(sites[y_col], errors="coerce"),
        x_range=x_range,
        y_range=y_range,
        bins_x=bins_x,
        bins_y=bins_y,
    )
    ok = (ix >= 0) & (iy >= 0)
    keep = np.zeros(len(sites), dtype=bool)
    keep[ok] = mask[ix[ok], iy[ok]]
    return keep


def label_density_filter(
    sites: pd.DataFrame,
    mask: np.ndarray,
    *,
    col: str = "density_filter",
    x_col: str = X_COL,
    y_col: str = Y_COL,
    x_range: tuple[float, float] = DEFAULT_X_RANGE,
    y_range: tuple[float, float] = DEFAULT_Y_RANGE,
) -> pd.DataFrame:
    """Add Kept/Dropped labels from a RelativeDepth trusted-cell mask."""
    out = sites.copy()
    keep = flag_sites_in_trusted_cells(
        out, mask, x_col=x_col, y_col=y_col, x_range=x_range, y_range=y_range
    )
    out[col] = np.where(keep, DENSITY_KEEP_LABEL, DENSITY_DROP_LABEL)
    return out


def _counts_row(labeled: pd.DataFrame, keep_pct: int) -> dict[str, float | int]:
    n_keep = int((labeled["density_filter"] == DENSITY_KEEP_LABEL).sum())
    n_drop = int((labeled["density_filter"] == DENSITY_DROP_LABEL).sum())
    n = n_keep + n_drop
    return {
        "keep_pct": int(keep_pct),
        "n_sites": n,
        "n_kept": n_keep,
        "n_dropped": n_drop,
        "frac_kept": (n_keep / n) if n else 0.0,
    }


def plot_density_keep_drop_scatters(
    sites: pd.DataFrame,
    mask: np.ndarray,
    output_prefix: str | Path,
    keep_pct: int,
    *,
    x_range: tuple[float, float] = DEFAULT_X_RANGE,
    y_range: tuple[float, float] = DEFAULT_Y_RANGE,
    x_col: str = X_COL,
    y_col: str = Y_COL,
    scatter_pairs: tuple[tuple[str, str, str, str], ...] | None = None,
    title_note: str | None = None,
    extra_counts: dict[str, float | int] | None = None,
    linear_x_lim: tuple[float, float] | None = None,
    linear_y_lim: tuple[float, float] | None = None,
    plot_log_axes: bool = True,
) -> pd.DataFrame:
    """Mean-vs-SD JointGrids (linear + optional log) coloured Kept vs Dropped.

    ``linear_x_lim`` / ``linear_y_lim`` fix the linear-axis view (and drop
    points outside that window). Log JointGrids keep auto range.
    ``plot_log_axes`` false writes only the linear JointGrids.
    """
    output_prefix = Path(output_prefix)
    labeled = label_density_filter(
        sites, mask, x_col=x_col, y_col=y_col, x_range=x_range, y_range=y_range
    )
    row = _counts_row(labeled, keep_pct)
    if extra_counts:
        row.update(extra_counts)
    counts = pd.DataFrame([row])
    save_df_to_tsv(counts, f"{output_prefix}.keep{int(keep_pct)}.counts.tsv")
    stem = f"{output_prefix}.keep{int(keep_pct)}"
    title_frac = title_note or f"top {int(keep_pct)}% density cells"
    pairs = scatter_pairs if scatter_pairs is not None else _MEAN_SD_SCATTERS
    for plot_x, plot_y, x_label, y_label in pairs:
        plot_df = labeled.dropna(subset=[plot_x, plot_y, "density_filter"])
        if plot_df.empty:
            continue
        for log_axes in ((False, True) if plot_log_axes else (False,)):
            frame = plot_df
            x_lim = None
            y_lim = None
            if log_axes:
                frame = plot_df[(plot_df[plot_x] > 0) & (plot_df[plot_y] > 0)]
                if frame.empty:
                    continue
            else:
                if linear_x_lim is not None:
                    x_lim = linear_x_lim
                    frame = frame[
                        (frame[plot_x] >= float(x_lim[0]))
                        & (frame[plot_x] <= float(x_lim[1]))
                    ]
                if linear_y_lim is not None:
                    y_lim = linear_y_lim
                    frame = frame[
                        (frame[plot_y] >= float(y_lim[0]))
                        & (frame[plot_y] <= float(y_lim[1]))
                    ]
                if frame.empty:
                    continue
            if log_axes:
                scale = ", log"
            elif x_lim == LINEAR_RELDEPTH_PLOT_RANGE and y_lim == LINEAR_RELDEPTH_PLOT_RANGE:
                scale = ", 0-10"
            elif x_lim == LINEAR_SITE_RANGE and y_lim == LINEAR_SITE_RANGE:
                scale = ", 0.02-6"
            elif x_lim == MEAN_RANGE_0_6:
                scale = ", mean 0-6"
            else:
                scale = ""
            pair = f"{plot_y}_vs_{plot_x}"
            filename = f"{stem}.scatter_{pair}.log.png" if log_axes else f"{stem}.scatter_{pair}.png"
            plot_joint_regression(
                frame,
                x_col=plot_x,
                y_col=plot_y,
                group_col="density_filter",
                x_label=x_label,
                y_label=y_label,
                filename=filename,
                title=f"{y_label} vs {x_label.lower()} ({title_frac}{scale})",
                figure_size=(10, 8),
                palette=DENSITY_FILTER_PALETTE,
                hue_order=list(DENSITY_FILTER_HUE_ORDER),
                log_x=log_axes,
                log_y=log_axes,
                x_lim=x_lim,
                y_lim=y_lim,
                draw_regression=False,
                alpha=0.2,
                s=8,
                rasterized=True,
                legend_title=None,
                bins_x=DEFAULT_BINS_X,
                bins_y=DEFAULT_BINS_Y,
                min_points=1,
                marginal_multiple="stack",
                stack_last_on_top=True,
                hue_alpha=DENSITY_FILTER_HUE_ALPHA,
                also_stack_counts=True,
            )
    return labeled


def ana_test_popdep_density_keep_plots(
    sites: pd.DataFrame | str | Path,
    masks_npz: str | Path,
    output_prefix: str | Path,
    *,
    x_range: tuple[float, float] = DEFAULT_X_RANGE,
    y_range: tuple[float, float] = DEFAULT_Y_RANGE,
) -> None:
    """Colour test_popdep points by genome 60%/80% trusted cells; plot mean vs SD."""
    table = sites if isinstance(sites, pd.DataFrame) else load_df_from_tsv(str(sites))
    packed = np.load(masks_npz)
    for keep_pct, key in ((60, "mask_60"), (80, "mask_80")):
        if key not in packed.files:
            raise KeyError(f"{masks_npz} missing {key}")
        plot_density_keep_drop_scatters(
            table,
            packed[key],
            output_prefix,
            keep_pct,
            x_range=x_range,
            y_range=y_range,
        )


def square_bins_for_n_cells(n_cells: int) -> tuple[int, int]:
    """Near-square factorisation; 90_000 → (300, 300)."""
    n = int(n_cells)
    if n < 1:
        raise ValueError(f"n_cells must be positive, got {n_cells}")
    side = int(round(np.sqrt(n)))
    other = int(round(n / side)) if side else n
    if side * other == n:
        return int(side), int(other)
    side = int(np.floor(np.sqrt(n)))
    while side > 1 and n % side != 0:
        side -= 1
    return int(side), int(n // side)


def log10_positive(values) -> np.ndarray:
    """log10 for finite positive values; others NaN."""
    xs = np.asarray(values, dtype=float)
    out = np.full(xs.shape, np.nan, dtype=float)
    ok = np.isfinite(xs) & (xs > 0)
    out[ok] = np.log10(xs[ok])
    return out


def expand_range_to_include_max(lo: float, hi: float) -> tuple[float, float]:
    """Half-open histogram range that still contains ``hi``."""
    left = float(lo)
    right = float(hi)
    if right <= left:
        right = left + 1e-6
    span = right - left
    return left, right + max(span * 1e-9, 1e-12)


def histogram_log_log(
    x,
    y,
    *,
    x_range: tuple[float, float] | None = None,
    y_range: tuple[float, float] | None = None,
    bins_x: int = 300,
    bins_y: int = 300,
) -> tuple[np.ndarray, int, int]:
    """2D histogram of log10(x) vs log10(y); non-positive values are ignored."""
    lx = log10_positive(x)
    ly = log10_positive(y)
    finite = np.isfinite(lx) & np.isfinite(ly)
    if not bool(finite.any()):
        return np.zeros((int(bins_x), int(bins_y)), dtype=np.int64), 0, 0
    if x_range is None:
        x_range = expand_range_to_include_max(float(lx[finite].min()), float(lx[finite].max()))
    if y_range is None:
        y_range = expand_range_to_include_max(float(ly[finite].min()), float(ly[finite].max()))
    return histogram_rel_mean_sd(
        lx, ly, x_range=x_range, y_range=y_range, bins_x=bins_x, bins_y=bins_y
    )


def add_log_reldepth_columns(sites: pd.DataFrame) -> pd.DataFrame:
    out = sites.copy()
    out[LOG_COL_X] = log10_positive(out[X_COL])
    out[LOG_COL_Y] = log10_positive(out[Y_COL])
    return out


def _site_table_paths(panel) -> list[Path]:
    if isinstance(panel, pd.DataFrame):
        return []
    if isinstance(panel, (str, Path)):
        root = Path(panel)
        if root.is_dir():
            return sorted(root.glob("*.popdep.sites.tsv"))
        return [root]
    paths: list[Path] = []
    for item in panel:
        path = Path(item)
        if path.is_dir():
            paths.extend(sorted(path.glob("*.popdep.sites.tsv")))
        else:
            paths.append(path)
    return paths


def _iter_reldepth_chunks(panel, *, chunksize: int = HIST_BATCH):
    if isinstance(panel, pd.DataFrame):
        yield panel
        return
    cols = [X_COL, Y_COL]
    for path in _site_table_paths(panel):
        for chunk in pd.read_csv(path, sep="\t", usecols=cols, chunksize=int(chunksize)):
            yield chunk


def scan_log_reldepth_limits(panel) -> tuple[tuple[float, float], tuple[float, float], dict[str, int]]:
    """Pass 1: log10 min/max of positive RelativeDepth mean and SD."""
    xmin, xmax = np.inf, -np.inf
    ymin, ymax = np.inf, -np.inf
    n_rows = 0
    n_pos = 0
    n_skip = 0
    for chunk in _iter_reldepth_chunks(panel):
        n_rows += int(len(chunk))
        lx = log10_positive(chunk[X_COL])
        ly = log10_positive(chunk[Y_COL])
        ok = np.isfinite(lx) & np.isfinite(ly)
        n_pos += int(ok.sum())
        n_skip += int((~ok).sum())
        if not bool(ok.any()):
            continue
        xmin = min(xmin, float(lx[ok].min()))
        xmax = max(xmax, float(lx[ok].max()))
        ymin = min(ymin, float(ly[ok].min()))
        ymax = max(ymax, float(ly[ok].max()))
    if not np.isfinite(xmin):
        raise ValueError("scan_log_reldepth_limits: no positive RelativeDepth rows")
    x_range = expand_range_to_include_max(xmin, xmax)
    y_range = expand_range_to_include_max(ymin, ymax)
    stats = {"n_rows": n_rows, "n_pos": n_pos, "n_skip": n_skip}
    return x_range, y_range, stats


def accumulate_log_reldepth_histogram(
    panel,
    *,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
    bins_x: int,
    bins_y: int,
) -> tuple[np.ndarray, dict[str, int]]:
    """Pass 2: stream RelativeDepth into a log10 mean-vs-SD histogram."""
    counts = np.zeros((int(bins_x), int(bins_y)), dtype=np.int64)
    n_rows = 0
    n_valid = 0
    n_oob = 0
    n_nan = 0
    for chunk in _iter_reldepth_chunks(panel):
        n_rows += int(len(chunk))
        part, n_ok, n_out = histogram_log_log(
            chunk[X_COL],
            chunk[Y_COL],
            x_range=x_range,
            y_range=y_range,
            bins_x=bins_x,
            bins_y=bins_y,
        )
        n_skip = int(len(chunk) - n_ok - n_out)
        n_nan += n_skip
        counts += part
        n_valid += n_ok
        n_oob += n_out
    return counts, {
        "n_rows": n_rows,
        "n_valid": n_valid,
        "n_oob": n_oob,
        "n_nan": n_nan,
    }


def accumulate_linear_reldepth_histogram(
    panel,
    *,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
    bins_x: int,
    bins_y: int,
) -> tuple[np.ndarray, dict[str, int]]:
    """Stream RelativeDepth mean vs SD into a linear histogram."""
    counts = np.zeros((int(bins_x), int(bins_y)), dtype=np.int64)
    n_rows = 0
    n_valid = 0
    n_oob = 0
    n_nan = 0
    for chunk in _iter_reldepth_chunks(panel):
        n_rows += int(len(chunk))
        part, n_ok, n_out = histogram_rel_mean_sd(
            chunk[X_COL],
            chunk[Y_COL],
            x_range=x_range,
            y_range=y_range,
            bins_x=bins_x,
            bins_y=bins_y,
        )
        n_skip = int(len(chunk) - n_ok - n_out)
        n_nan += n_skip
        counts += part
        n_valid += n_ok
        n_oob += n_out
    return counts, {
        "n_rows": n_rows,
        "n_valid": n_valid,
        "n_oob": n_oob,
        "n_nan": n_nan,
    }


def relative_depth_cv(mean, sd) -> np.ndarray:
    """RelSD / RelMean; non-positive or non-finite mean is NaN."""
    m = np.asarray(mean, dtype=float)
    s = np.asarray(sd, dtype=float)
    out = np.full(m.shape, np.nan, dtype=float)
    ok = np.isfinite(m) & np.isfinite(s) & (m > 0)
    out[ok] = s[ok] / m[ok]
    return out


def add_reldepth_cv_column(sites: pd.DataFrame) -> pd.DataFrame:
    out = sites.copy()
    if CV_COL not in out.columns:
        out[CV_COL] = np.nan
    miss = pd.to_numeric(out[CV_COL], errors="coerce").isna()
    if bool(miss.any()):
        out.loc[miss, CV_COL] = relative_depth_cv(out.loc[miss, X_COL], out.loc[miss, Y_COL])
    return out


def scan_mean_cv_limits(
    panel,
    *,
    x_range: tuple[float, float],
) -> tuple[tuple[float, float], dict[str, int]]:
    """Pass 1: observed CV min/max for sites with RelMean in ``x_range``."""
    ymin, ymax = np.inf, -np.inf
    n_rows = 0
    n_cv = 0
    n_skip = 0
    x_lo, x_hi = float(x_range[0]), float(x_range[1])
    for chunk in _iter_reldepth_chunks(panel):
        n_rows += int(len(chunk))
        cv = relative_depth_cv(chunk[X_COL], chunk[Y_COL])
        mean = np.asarray(chunk[X_COL], dtype=float)
        in_x = np.isfinite(mean) & (mean >= x_lo) & (mean < x_hi)
        ok = in_x & np.isfinite(cv)
        n_cv += int(ok.sum())
        n_skip += int((~ok).sum())
        if not bool(ok.any()):
            continue
        ymin = min(ymin, float(cv[ok].min()))
        ymax = max(ymax, float(cv[ok].max()))
    if not np.isfinite(ymin):
        raise ValueError("scan_mean_cv_limits: no finite CV rows in mean range")
    return expand_range_to_include_max(ymin, ymax), {
        "n_rows": n_rows,
        "n_cv": n_cv,
        "n_skip": n_skip,
    }


def accumulate_mean_cv_histogram(
    panel,
    *,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
    bins_x: int,
    bins_y: int,
) -> tuple[np.ndarray, dict[str, int]]:
    """Stream RelMean vs RelCV into a linear histogram."""
    counts = np.zeros((int(bins_x), int(bins_y)), dtype=np.int64)
    n_rows = 0
    n_valid = 0
    n_oob = 0
    n_nan = 0
    for chunk in _iter_reldepth_chunks(panel):
        n_rows += int(len(chunk))
        cv = relative_depth_cv(chunk[X_COL], chunk[Y_COL])
        part, n_ok, n_out = histogram_rel_mean_sd(
            chunk[X_COL],
            cv,
            x_range=x_range,
            y_range=y_range,
            bins_x=bins_x,
            bins_y=bins_y,
        )
        n_skip = int(len(chunk) - n_ok - n_out)
        n_nan += n_skip
        counts += part
        n_valid += n_ok
        n_oob += n_out
    return counts, {
        "n_rows": n_rows,
        "n_valid": n_valid,
        "n_oob": n_oob,
        "n_nan": n_nan,
    }


def _log10_cell_labels(edges: np.ndarray, n_show: int = 8) -> list[str]:
    n = len(edges) - 1
    labels = [""] * n
    if n <= 0:
        return labels
    step = max(1, n // max(int(n_show) - 1, 1))
    for i in list(range(0, n, step)) + [n - 1]:
        mid = 0.5 * (edges[i] + edges[i + 1])
        labels[i] = f"{10.0 ** mid:.3g}"
    return labels


def _linear_cell_labels(edges: np.ndarray, n_show: int = 8) -> list[str]:
    n = len(edges) - 1
    labels = [""] * n
    if n <= 0:
        return labels
    step = max(1, n // max(int(n_show) - 1, 1))
    for i in list(range(0, n, step)) + [n - 1]:
        mid = 0.5 * (edges[i] + edges[i + 1])
        labels[i] = f"{mid:.3g}"
    return labels


def plot_log_density_grid(
    counts: np.ndarray,
    filename: str,
    *,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
    mask_60: np.ndarray | None = None,
    mask_80: np.ndarray | None = None,
    mask_90: np.ndarray | None = None,
) -> None:
    """Heatmap of all log10 RelMean/RelSD cell counts (every input site)."""
    nx, ny = counts.shape
    x_edges = _bin_edges(x_range, nx)
    y_edges = _bin_edges(y_range, ny)
    matrix = np.log1p(counts.astype(float)).T[::-1]
    x_labels = _log10_cell_labels(x_edges)
    y_labels = _log10_cell_labels(y_edges)[::-1]
    overlays = []
    for mask, frac, label in (
        (mask_90, KEEP_FRAC_WIDE, "sites 90%"),
        (mask_80, KEEP_FRAC_LOOSE, "sites 80%"),
        (mask_60, KEEP_FRAC_STRICT, "sites 60%"),
    ):
        if mask is None:
            continue
        ys, xs = _mask_boundary_heatmap_coords(mask.T[::-1])
        if not xs:
            continue
        overlays.append(
            {
                "x": xs,
                "y": ys,
                "color": _SITE_MARK_COLORS[frac],
                "linestyle": "none",
                "marker": "s",
                "markersize": 0.8,
                "label": label,
            }
        )
    plot_heatmap_custom(
        matrix,
        x_labels,
        y_labels,
        title="PopDep log10 RelativeDepth mean vs SD cell density (all sites)",
        filename=filename,
        xlabel="RelativeDepth_Mean",
        ylabel="RelativeDepth_SD",
        cmap=SEQUENTIAL_WHITE_RED,
        cbar_label="log1p(sites / cell)",
        overlays=overlays or None,
        x_tick_step=max(1, nx // 8),
    )


def plot_linear_density_grid(
    counts: np.ndarray,
    filename: str,
    *,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
    masks: dict[int, np.ndarray] | None = None,
    xlabel: str = "RelativeDepth_Mean",
    ylabel: str = "RelativeDepth_SD",
    title: str | None = None,
) -> None:
    """Heatmap of RelMean vs SD or CV cell counts on a linear axis range."""
    nx, ny = counts.shape
    x_edges = _bin_edges(x_range, nx)
    y_edges = _bin_edges(y_range, ny)
    matrix = np.log1p(counts.astype(float)).T[::-1]
    x_labels = _linear_cell_labels(x_edges)
    y_labels = _linear_cell_labels(y_edges)[::-1]
    overlays = []
    for pct in sorted(masks or {}, reverse=True):
        mask = (masks or {})[pct]
        if mask is None:
            continue
        ys, xs = _mask_boundary_heatmap_coords(mask.T[::-1])
        if not xs:
            continue
        frac = pct / 100.0
        overlays.append(
            {
                "x": xs,
                "y": ys,
                "color": _SITE_MARK_COLORS.get(frac, "#7f7f7f"),
                "linestyle": "none",
                "marker": "s",
                "markersize": 0.8,
                "label": f"sites {pct}%",
            }
        )
    plot_heatmap_custom(
        matrix,
        x_labels,
        y_labels,
        title=title or "PopDep RelativeDepth mean vs SD cell density (all sites in range)",
        filename=filename,
        xlabel=xlabel,
        ylabel=ylabel,
        cmap=SEQUENTIAL_WHITE_RED,
        cbar_label="log1p(sites / cell)",
        overlays=overlays or None,
        x_tick_step=max(1, nx // 8),
    )


def write_log_density_outputs(
    counts: np.ndarray,
    stats: dict[str, int],
    output_prefix: str | Path,
    *,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
) -> dict[str, np.ndarray | dict[str, float | int]]:
    output_prefix = Path(output_prefix)
    masks: dict[str, np.ndarray | dict[str, float | int]] = {}
    infos: dict[int, dict[str, float | int]] = {}
    for frac in LOG_SITE_KEEP_FRACS:
        pct = int(round(100.0 * float(frac)))
        mask, info = trusted_mask_from_site_coverage(counts, frac)
        masks[f"mask_{pct}"] = mask
        infos[pct] = info
        print(
            f"[Info] site coverage {pct}% -> cells top {100.0 * float(info['cell_top_frac']):.4f}% "
            f"({info['n_cells_kept']}/{info['n_nonempty']} nonempty; "
            f"{info['n_sites_kept']}/{info['n_sites_total']} sites)",
            flush=True,
        )
    info_60, info_80, info_90 = infos[60], infos[80], infos[90]
    mask_60 = masks["mask_60"]
    mask_80 = masks["mask_80"]
    mask_90 = masks["mask_90"]
    summary_rows = [
        {
            "set": "all_nonempty_cells",
            "site_frac_target": 1.0,
            "n_cells": int(info_80["n_nonempty"]),
            "n_nonempty": int(info_80["n_nonempty"]),
            "cell_top_frac": 1.0,
            "n_sites": int(info_80["n_sites_total"]),
            "frac_sites": 1.0,
        }
    ]
    for pct, info in infos.items():
        summary_rows.append(
            {
                "set": f"sites_{pct}pct",
                "site_frac_target": float(info["site_frac_target"]),
                "n_cells": int(info["n_cells_kept"]),
                "n_nonempty": int(info["n_nonempty"]),
                "cell_top_frac": float(info["cell_top_frac"]),
                "n_sites": int(info["n_sites_kept"]),
                "frac_sites": float(info["site_frac_kept"]),
            }
        )
    summary = pd.DataFrame(summary_rows)
    save_df_to_tsv(summary, f"{output_prefix}.trusted_cells.tsv")
    save_histogram_npz(
        f"{output_prefix}.counts.npz", counts, stats, x_range=x_range, y_range=y_range
    )
    np.savez_compressed(
        f"{output_prefix}.masks.npz",
        mask_60=mask_60,
        mask_80=mask_80,
        mask_90=mask_90,
        x_range=np.asarray(x_range, dtype=float),
        y_range=np.asarray(y_range, dtype=float),
    )
    plot_log_density_grid(
        counts,
        f"{output_prefix}.rel_mean_sd.log.density.png",
        x_range=x_range,
        y_range=y_range,
        mask_60=mask_60,
        mask_80=mask_80,
        mask_90=mask_90,
    )
    curve = site_coverage_curve(counts)
    save_df_to_tsv(curve, f"{output_prefix}.site_coverage_curve.tsv")
    plot_site_coverage_curve(
        curve, f"{output_prefix}.site_coverage_curve.png"
    )
    plot_site_coverage_curve(
        curve, f"{output_prefix}.site_coverage_curve.log.png", log_y=True
    )
    save_thresholds(
        {
            "n_rows": int(stats["n_rows"]),
            "n_valid": int(stats["n_valid"]),
            "n_oob": int(stats["n_oob"]),
            "n_nan": int(stats["n_nan"]),
            "n_sites_top60": int(info_60["n_sites_kept"]),
            "n_sites_top80": int(info_80["n_sites_kept"]),
            "n_sites_top90": int(info_90["n_sites_kept"]),
            "n_cells_top60": int(info_60["n_cells_kept"]),
            "n_cells_top80": int(info_80["n_cells_kept"]),
            "n_cells_top90": int(info_90["n_cells_kept"]),
            "n_nonempty_cells": int(info_80["n_nonempty"]),
            "cell_top_frac_60": float(info_60["cell_top_frac"]),
            "cell_top_frac_80": float(info_80["cell_top_frac"]),
            "cell_top_frac_90": float(info_90["cell_top_frac"]),
            "n_cells": int(counts.size),
            "keep_frac_strict": KEEP_FRAC_STRICT,
            "keep_frac_loose": KEEP_FRAC_LOOSE,
            "keep_frac_wide": KEEP_FRAC_WIDE,
            "keep_rule": "site_coverage",
            "bins_x": int(counts.shape[0]),
            "bins_y": int(counts.shape[1]),
            "log10_x_min": float(x_range[0]),
            "log10_x_max": float(x_range[1]),
            "log10_y_min": float(y_range[0]),
            "log10_y_max": float(y_range[1]),
        },
        f"{output_prefix}.th.tsv",
    )
    print(f"[Info] Wrote log-density outputs under {output_prefix}")
    out: dict[str, np.ndarray | dict[str, float | int]] = {
        "mask_60": mask_60,
        "mask_80": mask_80,
        "mask_90": mask_90,
        "info_60": info_60,
        "info_80": info_80,
        "info_90": info_90,
    }
    return out


def write_linear_site_density_outputs(
    counts: np.ndarray,
    stats: dict[str, int],
    output_prefix: str | Path,
    *,
    x_range: tuple[float, float],
    y_range: tuple[float, float],
    keep_fracs: tuple[float, ...] = LINEAR_SITE_KEEP_FRACS,
    heatmap_stem: str = "rel_mean_sd",
    heatmap_ylabel: str = "RelativeDepth_SD",
    heatmap_title: str | None = None,
    keep_scale: str = "linear",
) -> dict[str, np.ndarray | dict[str, float | int]]:
    output_prefix = Path(output_prefix)
    masks: dict[str, np.ndarray | dict[str, float | int]] = {}
    infos: dict[int, dict[str, float | int]] = {}
    overlay_masks: dict[int, np.ndarray] = {}
    for frac in keep_fracs:
        pct = keep_frac_pct(frac)
        mask, info = trusted_mask_from_site_coverage(counts, frac)
        masks[f"mask_{pct}"] = mask
        infos[pct] = info
        overlay_masks[pct] = mask
        print(
            f"[Info] site coverage {pct}% -> cells top {100.0 * float(info['cell_top_frac']):.4f}% "
            f"({info['n_cells_kept']}/{info['n_nonempty']} nonempty; "
            f"{info['n_sites_kept']}/{info['n_sites_total']} sites)",
            flush=True,
        )
    n_nonempty = int(next(iter(infos.values()))["n_nonempty"])
    n_sites_total = int(next(iter(infos.values()))["n_sites_total"])
    summary_rows = [
        {
            "set": "all_nonempty_cells",
            "site_frac_target": 1.0,
            "n_cells": n_nonempty,
            "n_nonempty": n_nonempty,
            "cell_top_frac": 1.0,
            "n_sites": n_sites_total,
            "frac_sites": 1.0,
        }
    ]
    for pct, info in infos.items():
        summary_rows.append(
            {
                "set": f"sites_{pct}pct",
                "site_frac_target": float(info["site_frac_target"]),
                "n_cells": int(info["n_cells_kept"]),
                "n_nonempty": int(info["n_nonempty"]),
                "cell_top_frac": float(info["cell_top_frac"]),
                "n_sites": int(info["n_sites_kept"]),
                "frac_sites": float(info["site_frac_kept"]),
            }
        )
    summary = pd.DataFrame(summary_rows)
    save_df_to_tsv(summary, f"{output_prefix}.trusted_cells.tsv")
    save_histogram_npz(
        f"{output_prefix}.counts.npz", counts, stats, x_range=x_range, y_range=y_range
    )
    np.savez_compressed(
        f"{output_prefix}.masks.npz",
        x_range=np.asarray(x_range, dtype=float),
        y_range=np.asarray(y_range, dtype=float),
        **{f"mask_{pct}": overlay_masks[pct] for pct in overlay_masks},
    )
    plot_linear_density_grid(
        counts,
        f"{output_prefix}.{heatmap_stem}.density.png",
        x_range=x_range,
        y_range=y_range,
        masks=overlay_masks,
        ylabel=heatmap_ylabel,
        title=heatmap_title,
    )
    curve = site_coverage_curve(counts)
    save_df_to_tsv(curve, f"{output_prefix}.site_coverage_curve.tsv")
    plot_site_coverage_curve(
        curve,
        f"{output_prefix}.site_coverage_curve.png",
        site_marks=keep_fracs,
    )
    th: dict[str, float | int | str] = {
        "n_rows": int(stats["n_rows"]),
        "n_valid": int(stats["n_valid"]),
        "n_oob": int(stats["n_oob"]),
        "n_nan": int(stats["n_nan"]),
        "n_nonempty_cells": n_nonempty,
        "n_cells": int(counts.size),
        "keep_frac_strict": KEEP_FRAC_STRICT,
        "keep_frac_loose": KEEP_FRAC_LOOSE,
        "keep_frac_wide": KEEP_FRAC_WIDE,
        "keep_rule": "site_coverage",
        "keep_scale": keep_scale,
        "bins_x": int(counts.shape[0]),
        "bins_y": int(counts.shape[1]),
        "x_min": float(x_range[0]),
        "x_max": float(x_range[1]),
        "y_min": float(y_range[0]),
        "y_max": float(y_range[1]),
    }
    for pct, info in infos.items():
        th[f"n_sites_top{pct}"] = int(info["n_sites_kept"])
        th[f"n_cells_top{pct}"] = int(info["n_cells_kept"])
        th[f"cell_top_frac_{pct}"] = float(info["cell_top_frac"])
    save_thresholds(th, f"{output_prefix}.th.tsv")
    print(f"[Info] Wrote linear site-density outputs under {output_prefix}")
    out: dict[str, np.ndarray | dict[str, float | int]] = dict(masks)
    for pct, info in infos.items():
        out[f"info_{pct}"] = info
    return out


def _load_scatter_frame(panel, scatter, scatter_cap: int) -> pd.DataFrame:
    if scatter is not None:
        return scatter if isinstance(scatter, pd.DataFrame) else load_df_from_tsv(str(scatter))
    if isinstance(panel, pd.DataFrame):
        return panel
    paths = _site_table_paths(panel)
    if not paths:
        raise ValueError("no site tables to plot")
    if len(paths) == 1:
        frame = load_df_from_tsv(str(paths[0]))
        if len(frame) <= int(scatter_cap):
            return frame
        return frame.sample(n=int(scatter_cap), random_state=42)
    parts: list[pd.DataFrame] = []
    remain = int(scatter_cap)
    for path in paths:
        if remain <= 0:
            break
        chunk = pd.read_csv(path, sep="\t", nrows=remain)
        parts.append(chunk)
        remain -= int(len(chunk))
    return pd.concat(parts, ignore_index=True)


def ana_test_popdep_log_density(
    panel,
    output_prefix: str | Path,
    *,
    n_cells: int = LOG_N_CELLS,
    scatter=None,
    scatter_cap: int = SCATTER_CAP,
    counts_npz: str | Path | None = None,
    plot_keep_drop: bool = True,
) -> dict[str, np.ndarray | dict[str, float | int]]:
    """log10 RelMean vs RelSD histogram on all panel sites; keep/drop JointGrids.

    Axis range is the observed min/max of log10 values (every positive site).
    ``n_cells`` is split into a near-square grid (90_000 → 300×300). Trusted
    cells are the densest non-empty cells whose cumulative sites reach 60%,
    80% and 90% of counted sites (not those fractions of cells). ``counts_npz`` skips the
    two-pass histogram. Scatter plots may use a subsample. Linear RelDepth
    JointGrids use axes [0, 10]; log JointGrids keep the observed log10
    range. ``plot_keep_drop`` false skips the JointGrids (curve and heatmap
    still write).
    """
    if counts_npz:
        counts, stats, x_range, y_range = load_histogram_npz(counts_npz)
        print(
            f"[Info] reuse log histogram {counts_npz} bins={counts.shape[0]}x{counts.shape[1]}",
            flush=True,
        )
    else:
        bins_x, bins_y = square_bins_for_n_cells(n_cells)
        print(
            f"[Info] log10 RelDepth histogram n_cells={n_cells} bins={bins_x}x{bins_y}",
            flush=True,
        )
        x_range, y_range, _scan = scan_log_reldepth_limits(panel)
        print(
            f"[Info] log10 range mean={x_range} sd={y_range} rows={_scan['n_rows']} "
            f"positive={_scan['n_pos']}",
            flush=True,
        )
        counts, stats = accumulate_log_reldepth_histogram(
            panel, x_range=x_range, y_range=y_range, bins_x=bins_x, bins_y=bins_y
        )
    masks = write_log_density_outputs(
        counts, stats, output_prefix, x_range=x_range, y_range=y_range
    )
    if not plot_keep_drop:
        return masks
    scatter_df = add_log_reldepth_columns(
        _load_scatter_frame(panel, scatter, scatter_cap)
    )
    for keep_pct, key in ((60, "mask_60"), (80, "mask_80"), (90, "mask_90")):
        info = masks[f"info_{keep_pct}"]
        cell_pct = 100.0 * float(info["cell_top_frac"])
        print(
            f"[Info] sites {keep_pct}% -> cells top {cell_pct:.4f}% "
            f"({info['n_cells_kept']}/{info['n_nonempty']} nonempty)",
            flush=True,
        )
        plot_density_keep_drop_scatters(
            scatter_df,
            masks[key],
            output_prefix,
            keep_pct,
            x_col=LOG_COL_X,
            y_col=LOG_COL_Y,
            x_range=x_range,
            y_range=y_range,
            scatter_pairs=_RELDEPTH_SCATTERS,
            linear_x_lim=LINEAR_RELDEPTH_PLOT_RANGE,
            linear_y_lim=LINEAR_RELDEPTH_PLOT_RANGE,
            title_note=(
                f"sites {keep_pct}%; cells top {cell_pct:.2f}%"
            ),
            extra_counts={
                "cell_top_frac": float(info["cell_top_frac"]),
                "n_cells_kept": int(info["n_cells_kept"]),
                "n_nonempty": int(info["n_nonempty"]),
                "n_sites_kept_hist": int(info["n_sites_kept"]),
                "n_sites_total_hist": int(info["n_sites_total"]),
                "site_frac_kept_hist": float(info["site_frac_kept"]),
            },
        )
    return masks


def ana_test_popdep_linear_density(
    panel,
    output_prefix: str | Path,
    *,
    n_cells: int = LOG_N_CELLS,
    scatter=None,
    scatter_cap: int = SCATTER_CAP,
    counts_npz: str | Path | None = None,
    plot_keep_drop: bool = True,
    x_range: tuple[float, float] = LINEAR_SITE_RANGE,
    y_range: tuple[float, float] = LINEAR_SITE_RANGE,
) -> dict[str, np.ndarray | dict[str, float | int]]:
    """Linear RelMean vs RelSD histogram; keep/drop JointGrids, no log axes.

    Only sites inside ``x_range`` × ``y_range`` (default [0.02, 6]) enter the
    60/80/90/95/97/98/99% site-coverage filter. Points outside that window
    are dropped. ``n_cells`` is split into a near-square grid (90_000 →
    300×300). Scatter JointGrids use the same linear RelDepth axes.
    """
    hist_x = expand_range_to_include_max(float(x_range[0]), float(x_range[1]))
    hist_y = expand_range_to_include_max(float(y_range[0]), float(y_range[1]))
    if counts_npz:
        counts, stats, hist_x, hist_y = load_histogram_npz(counts_npz)
        print(
            f"[Info] reuse linear histogram {counts_npz} bins={counts.shape[0]}x{counts.shape[1]}",
            flush=True,
        )
    else:
        bins_x, bins_y = square_bins_for_n_cells(n_cells)
        print(
            f"[Info] linear RelDepth histogram n_cells={n_cells} bins={bins_x}x{bins_y} "
            f"mean={hist_x} sd={hist_y}",
            flush=True,
        )
        counts, stats = accumulate_linear_reldepth_histogram(
            panel, x_range=hist_x, y_range=hist_y, bins_x=bins_x, bins_y=bins_y
        )
    masks = write_linear_site_density_outputs(
        counts, stats, output_prefix, x_range=hist_x, y_range=hist_y
    )
    if not plot_keep_drop:
        return masks
    scatter_df = _load_scatter_frame(panel, scatter, scatter_cap)
    for frac in LINEAR_SITE_KEEP_FRACS:
        keep_pct = keep_frac_pct(frac)
        key = f"mask_{keep_pct}"
        info = masks[f"info_{keep_pct}"]
        cell_pct = 100.0 * float(info["cell_top_frac"])
        print(
            f"[Info] sites {keep_pct}% -> cells top {cell_pct:.4f}% "
            f"({info['n_cells_kept']}/{info['n_nonempty']} nonempty)",
            flush=True,
        )
        plot_density_keep_drop_scatters(
            scatter_df,
            masks[key],
            output_prefix,
            keep_pct,
            x_col=X_COL,
            y_col=Y_COL,
            x_range=hist_x,
            y_range=hist_y,
            scatter_pairs=_RELDEPTH_SCATTERS,
            linear_x_lim=x_range,
            linear_y_lim=y_range,
            plot_log_axes=False,
            title_note=(
                f"sites {keep_pct}%; cells top {cell_pct:.2f}%"
            ),
            extra_counts={
                "cell_top_frac": float(info["cell_top_frac"]),
                "n_cells_kept": int(info["n_cells_kept"]),
                "n_nonempty": int(info["n_nonempty"]),
                "n_sites_kept_hist": int(info["n_sites_kept"]),
                "n_sites_total_hist": int(info["n_sites_total"]),
                "site_frac_kept_hist": float(info["site_frac_kept"]),
            },
        )
    return masks


def plot_cv_vs_mean_by_class(
    sites: pd.DataFrame,
    filename: str | Path,
    *,
    x_lim: tuple[float, float] = MEAN_RANGE_0_6,
) -> None:
    """RelDepth CV vs mean JointGrid coloured by genomic class; mean clipped to ``x_lim``."""
    if "region_class" not in sites.columns or CV_COL not in sites.columns:
        print("[Info] skip cv_vs_mean_by_class: missing region_class or CV", flush=True)
        return
    frame = sites.dropna(subset=[X_COL, CV_COL, "region_class"]).copy()
    frame["region_class"] = frame["region_class"].astype(str)
    frame = frame[
        (frame[X_COL] >= float(x_lim[0])) & (frame[X_COL] <= float(x_lim[1]))
    ]
    if frame.empty:
        print("[Info] skip cv_vs_mean_by_class: no rows in mean window", flush=True)
        return
    plot_joint_regression(
        frame,
        x_col=X_COL,
        y_col=CV_COL,
        group_col="region_class",
        x_label="Relative mean depth",
        y_label="Relative depth CV",
        filename=str(filename),
        title="Relative depth CV vs relative mean depth by genomic class (mean 0-6)",
        figure_size=(10, 8),
        palette=REGION_CLASS_COLORS,
        hue_order=list(REGION_CLASS_HUE_ORDER),
        x_lim=x_lim,
        draw_regression=False,
        alpha=0.45,
        s=12,
        rasterized=True,
        legend_title="Genomic class",
        bins_x=DEFAULT_BINS_X,
        bins_y=DEFAULT_BINS_Y,
        min_points=1,
        marginal_multiple="fill",
        stack_last_on_top=True,
        hue_alpha=REGION_CLASS_ALPHA,
        also_stack_counts=True,
        hue_labels=REGION_CLASS_LEGEND_LABELS,
    )


def ana_test_popdep_cv_density(
    panel,
    output_prefix: str | Path,
    *,
    n_cells: int = LOG_N_CELLS,
    scatter=None,
    scatter_cap: int = SCATTER_CAP,
    counts_npz: str | Path | None = None,
    plot_keep_drop: bool = True,
    x_range: tuple[float, float] = MEAN_RANGE_0_6,
) -> dict[str, np.ndarray | dict[str, float | int]]:
    """RelMean vs RelCV histogram; keep/drop JointGrids and class CV plot.

    Mean axis is ``x_range`` (default [0, 6]). CV uses the observed min/max
    among in-range positive-mean sites (no CV clip). Site-coverage keep
    60/70/80/85/90/95/97/98/99%. Class JointGrid uses mean 0–6 and free CV.
    """
    hist_x = expand_range_to_include_max(float(x_range[0]), float(x_range[1]))
    if counts_npz:
        counts, stats, hist_x, hist_y = load_histogram_npz(counts_npz)
        print(
            f"[Info] reuse mean-CV histogram {counts_npz} bins={counts.shape[0]}x{counts.shape[1]}",
            flush=True,
        )
    else:
        bins_x, bins_y = square_bins_for_n_cells(n_cells)
        hist_y, scan = scan_mean_cv_limits(panel, x_range=hist_x)
        print(
            f"[Info] RelMean vs CV histogram n_cells={n_cells} bins={bins_x}x{bins_y} "
            f"mean={hist_x} cv={hist_y} rows={scan['n_rows']} cv_ok={scan['n_cv']}",
            flush=True,
        )
        counts, stats = accumulate_mean_cv_histogram(
            panel, x_range=hist_x, y_range=hist_y, bins_x=bins_x, bins_y=bins_y
        )
    masks = write_linear_site_density_outputs(
        counts,
        stats,
        output_prefix,
        x_range=hist_x,
        y_range=hist_y,
        keep_fracs=CV_SITE_KEEP_FRACS,
        heatmap_stem="rel_mean_cv",
        heatmap_ylabel="RelativeDepth_CV",
        heatmap_title="PopDep RelativeDepth mean vs CV cell density (all sites in mean range)",
        keep_scale="mean_cv",
    )
    scatter_df = add_reldepth_cv_column(_load_scatter_frame(panel, scatter, scatter_cap))
    plot_cv_vs_mean_by_class(
        scatter_df,
        f"{output_prefix}.reldepth.cv_vs_mean_by_class.png",
        x_lim=x_range,
    )
    if not plot_keep_drop:
        return masks
    for frac in CV_SITE_KEEP_FRACS:
        keep_pct = keep_frac_pct(frac)
        key = f"mask_{keep_pct}"
        info = masks[f"info_{keep_pct}"]
        cell_pct = 100.0 * float(info["cell_top_frac"])
        print(
            f"[Info] sites {keep_pct}% -> cells top {cell_pct:.4f}% "
            f"({info['n_cells_kept']}/{info['n_nonempty']} nonempty)",
            flush=True,
        )
        plot_density_keep_drop_scatters(
            scatter_df,
            masks[key],
            output_prefix,
            keep_pct,
            x_col=X_COL,
            y_col=CV_COL,
            x_range=hist_x,
            y_range=hist_y,
            scatter_pairs=_RELDEPTH_CV_SCATTERS,
            linear_x_lim=x_range,
            linear_y_lim=None,
            plot_log_axes=False,
            title_note=(f"sites {keep_pct}%; cells top {cell_pct:.2f}%"),
            extra_counts={
                "cell_top_frac": float(info["cell_top_frac"]),
                "n_cells_kept": int(info["n_cells_kept"]),
                "n_nonempty": int(info["n_nonempty"]),
                "n_sites_kept_hist": int(info["n_sites_kept"]),
                "n_sites_total_hist": int(info["n_sites_total"]),
                "site_frac_kept_hist": float(info["site_frac_kept"]),
            },
        )
    return masks

