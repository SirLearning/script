"""Plot confirmed whole-chromosome aneuploidy counts from BreachR reconciliation tables."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from breachr.constants import LOGICAL_AUTOSOMES
from breachr.reconcile import RECONCILE_OUTPUT_COLUMNS
from genetics.germplasm.sample.mosdepth_aneuploidy import CN_CLASS_BG_COLORS
from infra.utils.graph import plot_stacked_bar_chart
from infra.utils.io import save_df_to_tsv

CONFIRMED_RECONCILIATION = "consistent_whole_chr_aneuploidy"

MOSDEPTH_TO_PLOT_LABEL = {
    "null": "nullisomy",
    "monosomic": "monosomy",
    "trisomic": "trisomy",
    "tetrasomic+": "tetrasomy",
}

STACK_ORDER = ("nullisomy", "monosomy", "trisomy", "tetrasomy")

STACK_COLORS = {
    "nullisomy": CN_CLASS_BG_COLORS["null"],
    "monosomy": CN_CLASS_BG_COLORS["monosomic"],
    "trisomy": CN_CLASS_BG_COLORS["trisomic"],
    "tetrasomy": CN_CLASS_BG_COLORS["tetrasomic+"],
}

RECONCILE_ALL_FILENAME = "breachr_reconcile_all.tsv"


def _resolve_reconcile_table(breachr_input_dir: Path) -> Path:
    """Return the BreachR reconciliation table under ``breachr_input_dir``."""
    direct = breachr_input_dir / RECONCILE_ALL_FILENAME
    if direct.is_file():
        return direct
    candidates = sorted(breachr_input_dir.glob(f"*{RECONCILE_ALL_FILENAME}"))
    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        raise FileNotFoundError(
            f"Multiple reconciliation tables under {breachr_input_dir}: "
            f"{', '.join(p.name for p in candidates)}",
        )
    raise FileNotFoundError(
        f"Missing {RECONCILE_ALL_FILENAME} under {breachr_input_dir}",
    )


def _load_reconcile_table(reconcile_path: Path) -> pd.DataFrame:
    """Load reconciliation TSV without treating mosdepth label ``null`` as NA."""
    return pd.read_csv(reconcile_path, sep="\t", keep_default_na=False, na_values=[""])


def load_confirmed_aneuploidy_events(reconcile_path: Path | str) -> pd.DataFrame:
    """Load confirmed non-disomic whole-chromosome aneuploidy events."""
    path = Path(reconcile_path)
    df = _load_reconcile_table(path)
    missing_cols = [c for c in RECONCILE_OUTPUT_COLUMNS if c not in df.columns]
    if missing_cols:
        raise ValueError(
            f"{path} missing columns: {', '.join(missing_cols)}",
        )

    confirmed = df[df["reconciliation"] == CONFIRMED_RECONCILIATION].copy()
    confirmed = confirmed[~confirmed["mosdepth_label"].astype(str).str.lower().eq("disomic")]
    confirmed["aneuploidy_class"] = confirmed["mosdepth_label"].map(MOSDEPTH_TO_PLOT_LABEL)
    unknown = confirmed["aneuploidy_class"].isna()
    if unknown.any():
        labels = sorted(confirmed.loc[unknown, "mosdepth_label"].astype(str).unique())
        print(
            f"[Warning] Dropping {int(unknown.sum())} confirmed rows with "
            f"unmapped mosdepth labels: {labels}",
            file=sys.stderr,
        )
        confirmed = confirmed.loc[~unknown].copy()
    return confirmed.reset_index(drop=True)


def summarize_confirmed_aneuploidy_by_chromosome(events: pd.DataFrame) -> pd.DataFrame:
    """Count confirmed aneuploidy events per chromosome and copy-number class."""
    if events.empty:
        counts = pd.DataFrame(0, index=list(LOGICAL_AUTOSOMES), columns=list(STACK_ORDER))
        counts.index.name = "logical_chr"
        return counts.reset_index()

    counts = (
        events.groupby(["logical_chr", "aneuploidy_class"], observed=True)
        .size()
        .unstack(fill_value=0)
        .reindex(index=list(LOGICAL_AUTOSOMES), columns=list(STACK_ORDER), fill_value=0)
    )
    counts.index.name = "logical_chr"
    return counts.reset_index()


def plot_confirmed_aneuploidy_chromosome_distribution(
    breachr_input_dir: Path | str,
    output_dir: Path | str,
    *,
    plot_id: str = "benchmark",
    plot_title: str | None = None,
) -> dict[str, Path]:
    """
    Read BreachR reconciliation output, filter confirmed non-disomic aneuploidy,
    and publish a stacked chromosome distribution plot plus event/count tables.
    """
    breachr_dir = Path(breachr_input_dir)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    reconcile_path = _resolve_reconcile_table(breachr_dir)
    events = load_confirmed_aneuploidy_events(reconcile_path)
    counts = summarize_confirmed_aneuploidy_by_chromosome(events)

    events_path = out_dir / f"{plot_id}.confirmed_aneuploidy.events.tsv"
    counts_path = out_dir / f"{plot_id}.confirmed_aneuploidy.counts.tsv"
    plot_path = out_dir / f"{plot_id}.confirmed_aneuploidy.chromosome_distribution.png"

    save_df_to_tsv(events, str(events_path))
    save_df_to_tsv(counts, str(counts_path))

    title = plot_title or (
        f"Confirmed whole-chromosome aneuploidy by chromosome ({plot_id}; n={len(events)})"
    )
    count_df = counts.set_index("logical_chr")[list(STACK_ORDER)]
    plot_stacked_bar_chart(
        count_df,
        title=title,
        filename=str(plot_path),
        x_label="Chromosome",
        y_label="Confirmed aneuploidy events",
        stack_col_title="Copy-number class",
        stack_colors=STACK_COLORS,
        figure_size=(16, 6),
        rotate_xlabels=90,
    )

    print(f"Read {reconcile_path} -> {len(events)} confirmed events")
    return {
        "plot": plot_path,
        "events": events_path,
        "counts": counts_path,
        "reconcile": reconcile_path,
    }
