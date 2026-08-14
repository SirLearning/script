"""Subgenome-level absent or reference-ploidy mismatch detection."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from breachr.config import CallConfig
from breachr.constants import SUBGENOME_AUTOSOMES


@dataclass(frozen=True)
class SubgenomeState:
    """Detected state for one A/B/D subgenome."""

    subgenome: str
    state: str
    evidence: str
    raw_depth_median: float
    depth_rel_median: float
    mq_low_frac_median: float


def detect_subgenome_states(df: pd.DataFrame, config: CallConfig) -> dict[str, SubgenomeState]:
    """Flag absent subgenomes or tetraploid-vs-hexaploid reference mismatch."""
    states: dict[str, SubgenomeState] = {}
    raw_by_sg: dict[str, float] = {}
    rel_by_sg: dict[str, float] = {}
    mq_by_sg: dict[str, float] = {}

    for sg, chroms in SUBGENOME_AUTOSOMES.items():
        sub = df[df["logical_chr"].isin(chroms)]
        if sub.empty:
            continue
        raw = pd.to_numeric(sub["bin_mean_depth"], errors="coerce").dropna()
        rel = pd.to_numeric(sub["depth_rel"], errors="coerce").dropna()
        mq = pd.to_numeric(sub["mq_low_frac"], errors="coerce").dropna()
        reads = pd.to_numeric(sub["total_reads"], errors="coerce").fillna(0)
        raw_by_sg[sg] = float(raw.median()) if not raw.empty else 0.0
        rel_by_sg[sg] = float(rel.median()) if not rel.empty else float("nan")
        mq_by_sg[sg] = float(mq.median()) if not mq.empty else 0.0

        if raw.empty or float(reads.median()) <= config.subgenome_absent_reads_median_max:
            states[sg] = SubgenomeState(
                subgenome=sg,
                state="subgenome_absent",
                evidence="near_zero_reads_or_depth_on_subgenome",
                raw_depth_median=raw_by_sg[sg],
                depth_rel_median=rel_by_sg[sg],
                mq_low_frac_median=mq_by_sg[sg],
            )

    ref_sg = "A"
    if ref_sg not in raw_by_sg or raw_by_sg[ref_sg] <= 0:
        ref_sg = max(raw_by_sg, key=raw_by_sg.get, default="A")
    ref_depth = raw_by_sg.get(ref_sg, 0.0)

    for sg in ("A", "B", "D"):
        if sg in states or sg not in raw_by_sg or ref_depth <= 0:
            continue
        depth_ratio = raw_by_sg[sg] / ref_depth
        if (
            depth_ratio <= config.subgenome_depth_ratio_mismatch_max
            and mq_by_sg.get(sg, 0.0) >= config.subgenome_mq_low_mismatch_min
        ):
            states[sg] = SubgenomeState(
                subgenome=sg,
                state="subgenome_reference_mismatch",
                evidence=(
                    f"raw_depth_ratio_to_{ref_sg}={depth_ratio:.2f};"
                    f"elevated_low_mq={mq_by_sg[sg]:.3f};"
                    "tetraploid_or_lower_in_hexaploid_reference"
                ),
                raw_depth_median=raw_by_sg[sg],
                depth_rel_median=rel_by_sg.get(sg, float("nan")),
                mq_low_frac_median=mq_by_sg.get(sg, 0.0),
            )
    return states
