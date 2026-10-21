"""Unit tests for species LD-vs-genome scatter parsing and labeled plotting."""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from genetics.germplasm.meta.species_ld_plot import (
    load_species_ld_table,
    run_species_ld_plot,
)
from infra.utils.graph import place_nonoverlapping_labels

SPECIES_LD_CSV = """species_cn,species_scientific,ld_kb,ld_reported,ld_basis,ld_scope,ld_source_url,ld_validation,genome_size_mb,ncbi_assembly_accession,genome_size_source_url,genome_size_basis,mating_system,plot_included
玉米,Zea mays,7.07,x,x,x,x,x,2182.076,x,x,x,异交,TRUE
水稻,Oryza sativa,150,x,x,x,x,x,385.711,x,x,x,自交,TRUE
白云杉,Picea glauca,0.065,x,x,x,x,x,26571.973,x,x,x,异交,FALSE
大豆,Glycine max,326,x,x,x,x,x,978.387,x,x,x,自交,TRUE
人类,Homo sapiens,17.3,x,x,x,x,x,3099.7,x,x,x,异交,TRUE
"""


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_load_species_ld_table_keeps_plot_rows_and_maps_mating(tmp_path):
    csv_path = _write(tmp_path / "species_ld.csv", SPECIES_LD_CSV)
    df = load_species_ld_table(csv_path)
    assert list(df["species"]) == ["Zea mays", "Oryza sativa", "Glycine max", "Homo sapiens"]
    assert list(df["mating_system"]) == ["Outcrossing", "Selfing", "Selfing", "Outcrossing"]
    assert df.loc[0, "ld_kb"] == 7.07
    assert df.loc[1, "genome_size_mb"] == 385.711
    assert df.loc[df["species"] == "Homo sapiens", "ld_kb"].iloc[0] == 17.3


def test_place_nonoverlapping_labels_moves_close_labels_apart():
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter([1.0, 1.01], [1.0, 1.01])
    ax.set_xlim(0.5, 1.5)
    ax.set_ylim(0.5, 1.5)
    texts = place_nonoverlapping_labels(
        ax,
        [1.0, 1.01],
        [1.0, 1.01],
        ["Alpha species", "Beta species"],
        fontsize=12,
    )
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    b0 = texts[0].get_window_extent(renderer)
    b1 = texts[1].get_window_extent(renderer)
    plt.close(fig)
    c0 = ((b0.x0 + b0.x1) * 0.5, (b0.y0 + b0.y1) * 0.5)
    c1 = ((b1.x0 + b1.x1) * 0.5, (b1.y0 + b1.y1) * 0.5)
    dist = ((c0[0] - c1[0]) ** 2 + (c0[1] - c1[1]) ** 2) ** 0.5
    assert dist > 20


def test_run_species_ld_plot_writes_outputs(tmp_path):
    csv_path = _write(tmp_path / "species_ld.csv", SPECIES_LD_CSV)
    out_dir = tmp_path / "out"
    outputs = run_species_ld_plot(csv_path, out_dir)
    assert outputs["table_tsv"].exists()
    assert outputs["scatter_plot"].exists()
    assert outputs["scatter_plot_labeled"].exists()
    assert outputs["scatter_plot_linear_x"].exists()
    assert outputs["scatter_plot_linear_x_labeled"].exists()
    for key in (
        "scatter_plot",
        "scatter_plot_labeled",
        "scatter_plot_linear_x",
        "scatter_plot_linear_x_labeled",
    ):
        img = mpimg.imread(str(outputs[key]))
        height, width = img.shape[0], img.shape[1]
        assert height > 100
        assert width > 100
        # 2:1 canvas; tight bbox plus the legend still leaves a wide figure.
        assert width / height >= 1.6
