"""Chromosome-level dosage profiling for aneuploidy vs focal CNV."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from breachr.config import CallConfig


@dataclass(frozen=True)
class ChrDosageProfile:
    """Summary of per-bin relative depth on one logical chromosome."""

    n_bins: int
    median_rel: float
    std_rel: float
    iqr_rel: float
    disomic_frac: float
    loss_frac: float
    gain_frac: float
    longest_disomic_run_bins: int
    longest_loss_run_bins: int


def _longest_run(mask: np.ndarray) -> int:
    best = cur = 0
    for value in mask.astype(bool):
        if value:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best


def chr_dosage_profile(chr_df: pd.DataFrame, config: CallConfig) -> ChrDosageProfile | None:
    """Compute dosage shape metrics for one chromosome."""
    depth = pd.to_numeric(chr_df["depth_rel"], errors="coerce").dropna()
    if depth.empty:
        return None
    values = depth.to_numpy(dtype=float)
    q25, q75 = np.quantile(values, [0.25, 0.75])
    disomic = (values >= config.disomic_depth_min) & (values <= config.disomic_depth_max)
    loss = values <= config.deletion_depth_rel_max
    gain = values >= config.duplication_depth_rel_min
    return ChrDosageProfile(
        n_bins=int(len(values)),
        median_rel=float(np.median(values)),
        std_rel=float(np.std(values)),
        iqr_rel=float(q75 - q25),
        disomic_frac=float(disomic.mean()),
        loss_frac=float(loss.mean()),
        gain_frac=float(gain.mean()),
        longest_disomic_run_bins=_longest_run(disomic),
        longest_loss_run_bins=_longest_run(loss),
    )


def is_uniform_monosomic(profile: ChrDosageProfile, config: CallConfig) -> bool:
    """Whole-chromosome ~0.5x dose with little disomic signal."""
    if profile.n_bins < config.min_chr_bins_for_aneuploidy:
        return False
    if not (config.uniform_monosomic_median_min <= profile.median_rel <= config.uniform_monosomic_median_max):
        return False
    if profile.disomic_frac > config.chr_disomic_frac_max_for_aneuploidy:
        return False
    if profile.iqr_rel > config.uniform_aneuploidy_iqr_max:
        return False
    cv = profile.std_rel / max(profile.median_rel, 1e-6)
    return profile.loss_frac >= config.uniform_monosomic_loss_frac_min or cv <= config.uniform_aneuploidy_cv_max


def is_uniform_trisomic(profile: ChrDosageProfile, config: CallConfig) -> bool:
    """Whole-chromosome ~1.5x dose with little disomic signal."""
    if profile.n_bins < config.min_chr_bins_for_aneuploidy:
        return False
    if not (config.uniform_trisomic_median_min <= profile.median_rel <= config.uniform_trisomic_median_max):
        return False
    if profile.disomic_frac > config.chr_disomic_frac_max_for_aneuploidy:
        return False
    if profile.iqr_rel > config.uniform_aneuploidy_iqr_max:
        return False
    cv = profile.std_rel / max(profile.median_rel, 1e-6)
    return profile.gain_frac >= config.uniform_trisomic_gain_frac_min or cv <= config.uniform_aneuploidy_cv_max


def is_partial_cnv_pattern(profile: ChrDosageProfile, config: CallConfig) -> bool:
    """Mix of disomic and loss/gain bins — focal CNV, not whole-chr aneuploidy."""
    if profile.disomic_frac >= config.chr_partial_cnv_disomic_frac_min:
        if profile.loss_frac >= config.chr_partial_cnv_altered_frac_min:
            return True
        if profile.gain_frac >= config.chr_partial_cnv_altered_frac_min:
            return True
    min_run = max(2, int(profile.n_bins * config.chr_partial_cnv_disomic_run_frac_min))
    if profile.longest_disomic_run_bins >= min_run and profile.loss_frac >= 0.20:
        return True
    if profile.longest_disomic_run_bins >= min_run and profile.gain_frac >= 0.20:
        return True
    return False
