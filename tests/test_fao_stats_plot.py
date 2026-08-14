"""Unit tests for FAOSTAT wheat stats parsing and plotting."""

from pathlib import Path

import pandas as pd

from genetics.germplasm.fao.fao_stats_plot import (
    load_fao_production_series,
    load_fao_supply_percentages,
    run_fao_stats_plots,
)

PRODUCTION_CSV = """Domain Code,Domain,Area Code (M49),Area,Element Code,Element,Item Code (CPC),Item,Year Code,Year,Unit,Value,Flag,Flag Description,Note
"QCL","Crops and livestock products","001","World","5312","Area harvested","0111","Wheat","1961","1961","ha","100","A","Official value",""
"QCL","Crops and livestock products","001","World","5510","Production","0111","Wheat","1961","1961","t","200","A","Official value",""
"QCL","Crops and livestock products","001","World","5312","Area harvested","0111","Wheat","1962","1962","ha","110","A","Official value",""
"QCL","Crops and livestock products","001","World","5510","Production","0111","Wheat","1962","1962","t","220","A","Official value",""
"""

SUPPLY_CSV = """Domain Code,Domain,Area Code (M49),Area,Element Code,Element,Item Code (FBS),Item,Year Code,Year,Unit,Value,Flag,Flag Description,Note
"FBS","Food Balances (2010-)","001","World","664","Food supply (kcal/capita/day)","S2901","Grand Total","2010","2010","kcal/cap/d","1000","E","Estimated value",""
"FBS","Food Balances (2010-)","001","World","674","Protein supply quantity (g/capita/day)","S2901","Grand Total","2010","2010","g/cap/d","100","E","Estimated value",""
"FBS","Food Balances (2010-)","001","World","664","Food supply (kcal/capita/day)","S2905","Cereals - Excluding Beer","2010","2010","kcal/cap/d","500","E","Estimated value",""
"FBS","Food Balances (2010-)","001","World","674","Protein supply quantity (g/capita/day)","S2905","Cereals - Excluding Beer","2010","2010","g/cap/d","50","E","Estimated value",""
"FBS","Food Balances (2010-)","001","World","664","Food supply (kcal/capita/day)","S2511","Wheat and products","2010","2010","kcal/cap/d","250","E","Estimated value",""
"FBS","Food Balances (2010-)","001","World","674","Protein supply quantity (g/capita/day)","S2511","Wheat and products","2010","2010","g/cap/d","25","E","Estimated value",""
"FBS","Food Balances (2010-)","001","World","664","Food supply (kcal/capita/day)","S2901","Grand Total","2011","2011","kcal/cap/d","1100","E","Estimated value",""
"FBS","Food Balances (2010-)","001","World","674","Protein supply quantity (g/capita/day)","S2901","Grand Total","2011","2011","g/cap/d","110","E","Estimated value",""
"FBS","Food Balances (2010-)","001","World","664","Food supply (kcal/capita/day)","S2905","Cereals - Excluding Beer","2011","2011","kcal/cap/d","550","E","Estimated value",""
"FBS","Food Balances (2010-)","001","World","674","Protein supply quantity (g/capita/day)","S2905","Cereals - Excluding Beer","2011","2011","g/cap/d","55","E","Estimated value",""
"FBS","Food Balances (2010-)","001","World","664","Food supply (kcal/capita/day)","S2511","Wheat and products","2011","2011","kcal/cap/d","275","E","Estimated value",""
"FBS","Food Balances (2010-)","001","World","674","Protein supply quantity (g/capita/day)","S2511","Wheat and products","2011","2011","g/cap/d","27.5","E","Estimated value",""
"""


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_load_fao_production_series(tmp_path):
    csv_path = _write(tmp_path / "production.csv", PRODUCTION_CSV)
    df = load_fao_production_series(csv_path)
    assert list(df.columns) == ["year", "production_t", "area_harvested_ha"]
    assert len(df) == 2
    assert df.loc[0, "production_t"] == 200.0
    assert df.loc[1, "area_harvested_ha"] == 110.0


def test_load_fao_supply_percentages(tmp_path):
    csv_path = _write(tmp_path / "supply.csv", SUPPLY_CSV)
    df = load_fao_supply_percentages(csv_path)
    assert len(df) == 2
    row_2010 = df[df["year"] == 2010].iloc[0]
    assert row_2010["cereal_food_pct"] == 50.0
    assert row_2010["cereal_protein_pct"] == 50.0
    assert row_2010["total_food_pct"] == 25.0
    assert row_2010["total_protein_pct"] == 25.0


def test_run_fao_stats_plots_writes_outputs(tmp_path):
    production_csv = _write(tmp_path / "production.csv", PRODUCTION_CSV)
    supply_csv = _write(tmp_path / "supply.csv", SUPPLY_CSV)
    out_dir = tmp_path / "out"
    outputs = run_fao_stats_plots(production_csv, supply_csv, out_dir)

    for key in (
        "production_tsv",
        "supply_tsv",
        "wheat_production_area",
        "wheat_cereal_percentage",
        "wheat_total_percentage",
    ):
        assert outputs[key].exists()

    prod = pd.read_csv(outputs["production_tsv"], sep="\t")
    assert "production_t" in prod.columns
