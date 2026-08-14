"""Load and validate ``*.mapq_depth_bin.tsv`` inputs."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from breachr.constants import MAPQ_DEPTH_BIN_COLUMNS, MQ_BIN_LABELS
from infra.utils.io import load_df_generic


def validate_mapq_depth_bin_df(df: pd.DataFrame) -> None:
    """Ensure a dataframe matches the published MAPQ-depth bin schema."""
    if df is None or df.empty:
        raise ValueError("MAPQ-depth bin table is empty or missing.")
    missing = [c for c in MAPQ_DEPTH_BIN_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"MAPQ-depth bin table missing columns: {missing}")
    for col in ("bin_idx", "bin_start", "bin_end", "total_reads", "bin_mean_depth", *MQ_BIN_LABELS):
        if not pd.api.types.is_numeric_dtype(df[col]):
            df[col] = pd.to_numeric(df[col], errors="coerce")


def load_mapq_depth_bin_tsv(path: str | Path) -> pd.DataFrame:
    """Read one sample ``*.mapq_depth_bin.tsv`` file."""
    path = Path(path)
    df = load_df_generic(str(path))
    validate_mapq_depth_bin_df(df)
    out = df[list(MAPQ_DEPTH_BIN_COLUMNS)].copy()
    out["logical_chr"] = out["logical_chr"].astype(str)
    return out.sort_values(["logical_chr", "bin_idx"]).reset_index(drop=True)
