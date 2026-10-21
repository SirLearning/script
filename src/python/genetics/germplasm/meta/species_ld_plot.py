"""Plot comparative LD decay distance vs genome size from a species metadata table."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from infra.utils.graph import (
    TWO_CATEGORY_BLUE,
    TWO_CATEGORY_RED,
    plot_scatter_labeled,
)
from infra.utils.io import load_df_generic, save_df_to_tsv

REQUIRED_COLUMNS = (
    "species_scientific",
    "ld_kb",
    "genome_size_mb",
    "mating_system",
)

MATING_SYSTEM_MAP = {
    "自交": "Selfing",
    "异交": "Outcrossing",
    "selfing": "Selfing",
    "outcrossing": "Outcrossing",
    "Selfing": "Selfing",
    "Outcrossing": "Outcrossing",
}

MATING_SYSTEM_COLORS = {
    "Selfing": TWO_CATEGORY_BLUE,
    "Outcrossing": TWO_CATEGORY_RED,
}

MATING_SYSTEM_ORDER = ("Selfing", "Outcrossing")


def _as_plot_included(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    text = series.astype(str).str.strip().str.upper()
    return text.isin({"TRUE", "1", "YES", "T"})


def _normalize_mating_system(value) -> str:
    key = str(value).strip()
    if key not in MATING_SYSTEM_MAP:
        raise ValueError(f"Unknown mating_system value: {value!r}")
    return MATING_SYSTEM_MAP[key]


def load_species_ld_table(csv_path: str | Path) -> pd.DataFrame:
    """
    Load a species LD-vs-genome table and keep rows marked for plotting.

    Mating-system labels are mapped to English (Selfing / Outcrossing).
    """
    raw = load_df_generic(str(csv_path))
    missing = [col for col in REQUIRED_COLUMNS if col not in raw.columns]
    if missing:
        raise ValueError(f"Missing expected columns in species LD CSV: {missing}")

    out = raw.copy()
    if "plot_included" in out.columns:
        out = out.loc[_as_plot_included(out["plot_included"])].copy()

    out["ld_kb"] = pd.to_numeric(out["ld_kb"], errors="coerce")
    out["genome_size_mb"] = pd.to_numeric(out["genome_size_mb"], errors="coerce")
    out["mating_system"] = out["mating_system"].map(_normalize_mating_system)
    out["species"] = out["species_scientific"].astype(str).str.strip()
    keep = ["species", "ld_kb", "genome_size_mb", "mating_system"]
    if "species_cn" in out.columns:
        keep.append("species_cn")
    if "genome_assembly" in out.columns:
        keep.append("genome_assembly")
    if "ncbi_assembly_accession" in out.columns:
        keep.append("ncbi_assembly_accession")
    return out.dropna(subset=["species", "ld_kb", "genome_size_mb", "mating_system"]).reset_index(drop=True)[keep]


def plot_species_ld_vs_genome(
    species_df: pd.DataFrame,
    output_png: str | Path,
    *,
    labeled: bool = True,
    log_x: bool = True,
) -> Path:
    """Scatter of LD decay distance vs genome size, colored by mating system.

    The canvas is 2:1 (width:height). Genome size is log-scaled. The LD axis
    is log by default; set ``log_x=False`` for a linear kb scale.
    """
    plot_scatter_labeled(
        species_df,
        x_col="ld_kb",
        y_col="genome_size_mb",
        hue_col="mating_system",
        label_col="species" if labeled else None,
        title="LD decay distance vs genome size",
        filename=str(output_png),
        hue_order=MATING_SYSTEM_ORDER,
        hue_colors=MATING_SYSTEM_COLORS,
        xlabel="LD decay distance (kb)",
        ylabel="Genome size (Mb)",
        log_x=log_x,
        log_y=True,
        italic_labels=True,
        s=160,
        label_fontsize=11,
        figure_size=(16, 8),
    )
    return Path(output_png)


def run_species_ld_plot(csv_path: str | Path, output_dir: str | Path) -> dict[str, Path]:
    """Parse the species LD table and emit log-x and linear-x scatter PNGs."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    species_df = load_species_ld_table(csv_path)
    table_tsv = out_dir / "species_ld_vs_genome.tsv"
    save_df_to_tsv(species_df, str(table_tsv))
    unlabeled = plot_species_ld_vs_genome(
        species_df, out_dir / "species_ld_vs_genome.png", labeled=False, log_x=True
    )
    labeled = plot_species_ld_vs_genome(
        species_df, out_dir / "species_ld_vs_genome.labeled.png", labeled=True, log_x=True
    )
    unlabeled_lin = plot_species_ld_vs_genome(
        species_df, out_dir / "species_ld_vs_genome.linear_x.png", labeled=False, log_x=False
    )
    labeled_lin = plot_species_ld_vs_genome(
        species_df, out_dir / "species_ld_vs_genome.linear_x.labeled.png", labeled=True, log_x=False
    )
    return {
        "table_tsv": table_tsv,
        "scatter_plot": unlabeled,
        "scatter_plot_labeled": labeled,
        "scatter_plot_linear_x": unlabeled_lin,
        "scatter_plot_linear_x_labeled": labeled_lin,
    }
