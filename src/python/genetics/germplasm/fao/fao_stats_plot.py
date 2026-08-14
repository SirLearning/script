"""Parse FAOSTAT CSV exports and plot global wheat production and supply share trends."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from infra.utils.graph import plot_dual_axis_line_series, plot_multi_line_series
from infra.utils.io import load_df_generic, save_df_to_tsv

ELEMENT_PRODUCTION = "Production"
ELEMENT_AREA_HARVESTED = "Area harvested"
ELEMENT_FOOD_SUPPLY = "Food supply (kcal/capita/day)"
ELEMENT_PROTEIN_SUPPLY = "Protein supply quantity (g/capita/day)"

ITEM_WHEAT_PRODUCTS = "Wheat and products"
ITEM_CEREALS = "Cereals - Excluding Beer"
ITEM_GRAND_TOTAL = "Grand Total"


def _numeric_value(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _pivot_element_item(df: pd.DataFrame, elements: tuple[str, ...], items: tuple[str, ...]) -> pd.DataFrame:
    subset = df[df["Element"].isin(elements) & df["Item"].isin(items)].copy()
    subset["Value"] = _numeric_value(subset["Value"])
    subset["Year"] = _numeric_value(subset["Year"]).astype("Int64")
    return subset.pivot_table(
        index="Year",
        columns=["Element", "Item"],
        values="Value",
        aggfunc="first",
    ).sort_index()


def load_fao_production_series(production_csv: str | Path) -> pd.DataFrame:
    """Return yearly wheat production (t) and area harvested (ha) from a FAOSTAT QCL export."""
    raw = load_df_generic(str(production_csv))
    wide = _pivot_element_item(raw, (ELEMENT_PRODUCTION, ELEMENT_AREA_HARVESTED), ("Wheat",))
    out = pd.DataFrame({"year": wide.index.astype(int)})
    out["production_t"] = wide[(ELEMENT_PRODUCTION, "Wheat")].values
    out["area_harvested_ha"] = wide[(ELEMENT_AREA_HARVESTED, "Wheat")].values
    return out.dropna(subset=["production_t", "area_harvested_ha"]).reset_index(drop=True)


def load_fao_supply_percentages(supply_csv: str | Path) -> pd.DataFrame:
    """
    Return yearly wheat share of cereal and total food/protein supply (%).

    Percentages are Wheat and products divided by the matching element for
    Cereals - Excluding Beer or Grand Total.
    """
    raw = load_df_generic(str(supply_csv))
    elements = (ELEMENT_FOOD_SUPPLY, ELEMENT_PROTEIN_SUPPLY)
    items = (ITEM_WHEAT_PRODUCTS, ITEM_CEREALS, ITEM_GRAND_TOTAL)
    wide = _pivot_element_item(raw, elements, items)

    def _pct(numerator_col, denominator_col) -> pd.Series:
        num = wide[numerator_col]
        den = wide[denominator_col]
        return (num / den) * 100.0

    out = pd.DataFrame({"year": wide.index.astype(int)})
    out["cereal_food_pct"] = _pct(
        (ELEMENT_FOOD_SUPPLY, ITEM_WHEAT_PRODUCTS),
        (ELEMENT_FOOD_SUPPLY, ITEM_CEREALS),
    ).values
    out["cereal_protein_pct"] = _pct(
        (ELEMENT_PROTEIN_SUPPLY, ITEM_WHEAT_PRODUCTS),
        (ELEMENT_PROTEIN_SUPPLY, ITEM_CEREALS),
    ).values
    out["total_food_pct"] = _pct(
        (ELEMENT_FOOD_SUPPLY, ITEM_WHEAT_PRODUCTS),
        (ELEMENT_FOOD_SUPPLY, ITEM_GRAND_TOTAL),
    ).values
    out["total_protein_pct"] = _pct(
        (ELEMENT_PROTEIN_SUPPLY, ITEM_WHEAT_PRODUCTS),
        (ELEMENT_PROTEIN_SUPPLY, ITEM_GRAND_TOTAL),
    ).values
    return out.dropna().reset_index(drop=True)


def plot_fao_wheat_production_area(production_df: pd.DataFrame, output_png: str | Path) -> Path:
    """Dual-axis line chart: production (left, t) vs area harvested (right, ha)."""
    plot_dual_axis_line_series(
        production_df,
        x_col="year",
        left_spec={
            "y_col": "production_t",
            "label": "Production",
            "color": "#1f77b4",
            "y_label": "Production (t)",
        },
        right_spec={
            "y_col": "area_harvested_ha",
            "label": "Area harvested",
            "color": "#ff7f0e",
            "y_label": "Area harvested (ha)",
        },
        title="Global wheat production and area harvested",
        filename=str(output_png),
        x_label="Year",
        rotate_xlabels=45,
    )
    return Path(output_png)


def plot_fao_wheat_cereal_percentage(supply_pct_df: pd.DataFrame, output_png: str | Path) -> Path:
    """Wheat share of cereal food and protein supply (%)."""
    plot_data = supply_pct_df.rename(
        columns={
            "cereal_food_pct": "food_supply_pct",
            "cereal_protein_pct": "protein_supply_pct",
        }
    )
    plot_multi_line_series(
        plot_data,
        x_col="year",
        y_specs=[
            {
                "y_col": "food_supply_pct",
                "label": "Food supply (kcal/capita/day)",
                "color": "#1f77b4",
                "marker": "o",
            },
            {
                "y_col": "protein_supply_pct",
                "label": "Protein supply (g/capita/day)",
                "color": "#d62728",
                "marker": "s",
            },
        ],
        title="Wheat share of cereal supply (world)",
        filename=str(output_png),
        x_label="Year",
        y_label="Share of cereal supply (%)",
        rotate_xlabels=45,
    )
    return Path(output_png)


def plot_fao_wheat_total_percentage(supply_pct_df: pd.DataFrame, output_png: str | Path) -> Path:
    """Wheat share of total food and protein supply (%)."""
    plot_data = supply_pct_df.rename(
        columns={
            "total_food_pct": "food_supply_pct",
            "total_protein_pct": "protein_supply_pct",
        }
    )
    plot_multi_line_series(
        plot_data,
        x_col="year",
        y_specs=[
            {
                "y_col": "food_supply_pct",
                "label": "Food supply (kcal/capita/day)",
                "color": "#2ca02c",
                "marker": "o",
            },
            {
                "y_col": "protein_supply_pct",
                "label": "Protein supply (g/capita/day)",
                "color": "#9467bd",
                "marker": "s",
            },
        ],
        title="Wheat share of total food supply (world)",
        filename=str(output_png),
        x_label="Year",
        y_label="Share of total supply (%)",
        rotate_xlabels=45,
    )
    return Path(output_png)


def run_fao_stats_plots(
    production_csv: str | Path,
    supply_csv: str | Path,
    output_dir: str | Path,
) -> dict[str, Path]:
    """
    Parse FAOSTAT exports, write summary TSVs, and emit three PNG figures.

    Returns a mapping of artefact stem to output path.
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    production_df = load_fao_production_series(production_csv)
    supply_pct_df = load_fao_supply_percentages(supply_csv)

    prod_tsv = out_dir / "wheat_production_area.tsv"
    supply_tsv = out_dir / "wheat_supply_percentage.tsv"
    save_df_to_tsv(production_df, str(prod_tsv))
    save_df_to_tsv(supply_pct_df, str(supply_tsv))

    plots = {
        "wheat_production_area": plot_fao_wheat_production_area(
            production_df, out_dir / "wheat_production_area.png"
        ),
        "wheat_cereal_percentage": plot_fao_wheat_cereal_percentage(
            supply_pct_df, out_dir / "wheat_cereal_percentage.png"
        ),
        "wheat_total_percentage": plot_fao_wheat_total_percentage(
            supply_pct_df, out_dir / "wheat_total_percentage.png"
        ),
    }
    return {"production_tsv": prod_tsv, "supply_tsv": supply_tsv, **plots}
