"""Feature engineering from raw MAPQ-depth bins."""

from __future__ import annotations

import numpy as np
import pandas as pd

from breachr.config import CallConfig
from breachr.constants import MQ_LOW_LABELS, SUBGENOME_AUTOSOMES


def subgenome_from_logical_chr(logical_chr: str) -> str:
    """Return A/B/D suffix from a logical chromosome name."""
    if len(logical_chr) < 2:
        return ""
    return logical_chr[-1]


def add_mapq_fractions(df: pd.DataFrame) -> pd.DataFrame:
    """Add low/mid/high MAPQ fractions and depth support flags."""
    out = df.copy()
    reads = pd.to_numeric(out["total_reads"], errors="coerce").fillna(0).astype(float)
    reads_safe = reads.replace(0, np.nan)
    out["mq_low_frac"] = out[list(MQ_LOW_LABELS)].sum(axis=1) / reads_safe
    out["mq_low_frac"] = out["mq_low_frac"].fillna(0.0)
    out["subgenome"] = out["logical_chr"].map(subgenome_from_logical_chr)
    return out


def compute_subgenome_baselines(df: pd.DataFrame, config: CallConfig) -> dict[str, float]:
    """Median ``bin_mean_depth`` per subgenome across autosomal bins."""
    baselines: dict[str, float] = {}
    for sg, chroms in SUBGENOME_AUTOSOMES.items():
        sub = df[df["logical_chr"].isin(chroms)]
        depths = pd.to_numeric(sub["bin_mean_depth"], errors="coerce")
        depths = depths[(depths > 0) & np.isfinite(depths)]
        if depths.empty:
            baselines[sg] = float("nan")
        else:
            baselines[sg] = float(depths.median())
    return baselines


def add_depth_relatives(df: pd.DataFrame, baselines: dict[str, float]) -> pd.DataFrame:
    """Normalize depth by subgenome baseline."""
    out = df.copy()
    base = out["subgenome"].map(baselines).astype(float)
    depth = pd.to_numeric(out["bin_mean_depth"], errors="coerce").astype(float)
    out["depth_rel"] = depth / base.replace(0, np.nan)
    out["depth_rel"] = out["depth_rel"].replace([np.inf, -np.inf], np.nan)
    return out


def add_smoothed_depth(df: pd.DataFrame, config: CallConfig) -> pd.DataFrame:
    """Rolling median of ``depth_rel`` along each logical chromosome."""
    out = df.copy()
    smoothed = []
    for _chr, grp in out.groupby("logical_chr", sort=False):
        series = pd.to_numeric(grp["depth_rel"], errors="coerce")
        roll = series.rolling(
            config.smooth_window_bins,
            center=True,
            min_periods=1,
        ).median()
        smoothed.append(roll)
    out["depth_rel_smooth"] = pd.concat(smoothed).sort_index()
    return out


def sample_mq_low_frac_reference(df: pd.DataFrame, config: CallConfig) -> float:
    """Robust genome-wide reference for low-MAPQ fraction on supported bins."""
    reads = pd.to_numeric(df["total_reads"], errors="coerce").fillna(0)
    supported = df.loc[reads >= config.min_reads_per_bin, "mq_low_frac"]
    supported = pd.to_numeric(supported, errors="coerce").dropna()
    if supported.empty:
        return 0.0
    return float(supported.median())
