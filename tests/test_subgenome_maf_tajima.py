"""Stacked A/B/D MAF and Tajima's D distribution plots."""

from pathlib import Path

import pandas as pd
import pytest

from genetics.genomics.variant.subgenome_maf_tajima import (
    load_maf_from_afreq,
    load_tajima_d_table,
    run_subgenome_maf_plots,
    run_subgenome_tajima_plots,
    summarize_numeric_by_subgenome,
)
from infra.utils.graph import plot_faceted_histograms


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


AFREQ_A = """#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT
1\t1-1\tA\tT\t0.02\t100
1\t1-2\tC\tG\t0.80\t100
"""

AFREQ_B = """#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT
7\t7-1\tA\tT\t0.10\t100
7\t7-2\tC\tG\t0.40\t100
"""

AFREQ_D = """#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT
13\t13-1\tA\tT\t0.25\t100
13\t13-2\tC\tG\t0.49\t100
"""

TAJIMA_A = """CHROM\tBIN_START\tN_SNPS\tTajimaD
1\t0\t12\t-1.20
1\t100000\t8\tnan
1\t200000\t15\t0.40
"""

TAJIMA_B = """CHROM\tBIN_START\tN_SNPS\tTajimaD
7\t0\t10\t-0.50
7\t100000\t11\t-0.10
"""

TAJIMA_D = """CHROM\tBIN_START\tN_SNPS\tTajimaD
13\t0\t9\t0.80
13\t100000\t14\t1.10
"""


def test_load_maf_from_afreq_uses_minor_allele(tmp_path):
    path = _write(tmp_path / "A.info.afreq", AFREQ_A)
    df = load_maf_from_afreq(path, "A")
    assert df["MAF"].tolist() == pytest.approx([0.02, 0.20])
    assert set(df["subgenome"]) == {"A"}


def test_load_tajima_d_table_drops_nan(tmp_path):
    path = _write(tmp_path / "A.variant.tajima_d.Tajima.D", TAJIMA_A)
    df = load_tajima_d_table(path, "A")
    assert len(df) == 2
    assert df["TajimaD"].tolist() == [-1.20, 0.40]


def test_summarize_numeric_by_subgenome_order():
    df = pd.DataFrame(
        {
            "subgenome": ["D", "A", "B", "A"],
            "MAF": [0.1, 0.2, 0.3, 0.4],
        }
    )
    summary = summarize_numeric_by_subgenome(df, "MAF")
    assert list(summary["subgenome"]) == ["A", "B", "D"]
    assert summary.loc[0, "n"] == 2
    assert summary.loc[0, "mean"] == pytest.approx(0.3)


def test_plot_faceted_histograms_writes_png(tmp_path):
    df = pd.DataFrame(
        {
            "MAF": [0.01, 0.02, 0.10, 0.20, 0.30, 0.40],
            "subgenome": ["A", "A", "B", "B", "D", "D"],
        }
    )
    out = tmp_path / "stacked.png"
    plot_faceted_histograms(
        data=df,
        x_col="MAF",
        facet_col="subgenome",
        facet_order=["A", "B", "D"],
        title="MAF",
        filename=str(out),
        xlim=(0.0, 0.5),
        log_scale=True,
    )
    assert out.exists()
    assert out.stat().st_size > 0


def test_plot_faceted_histograms_shared_ylim_plain_titles(tmp_path, monkeypatch):
    import matplotlib.pyplot as plt

    recorded = {}
    orig_savefig = plt.Figure.savefig

    def wrapped(self, *args, **kwargs):
        recorded["titles"] = [ax.get_title() for ax in self.axes]
        recorded["ylims"] = [tuple(ax.get_ylim()) for ax in self.axes]
        return orig_savefig(self, *args, **kwargs)

    monkeypatch.setattr(plt.Figure, "savefig", wrapped)
    df = pd.DataFrame(
        {
            "MAF": [0.01] * 8 + [0.20, 0.30, 0.40, 0.45],
            "subgenome": ["A"] * 8 + ["B", "B", "D", "D"],
        }
    )
    plot_faceted_histograms(
        data=df,
        x_col="MAF",
        facet_col="subgenome",
        facet_order=["A", "B", "D"],
        title="MAF",
        filename=str(tmp_path / "shared.png"),
        facet_labels={"A": "Subgenome A", "B": "Subgenome B", "D": "Subgenome D"},
        xlim=(0.0, 0.5),
        sharey=True,
    )
    assert recorded["titles"] == ["Subgenome A", "Subgenome B", "Subgenome D"]
    assert all("n=" not in title for title in recorded["titles"])
    assert recorded["ylims"][0] == recorded["ylims"][1] == recorded["ylims"][2]


def test_run_subgenome_maf_plots_writes_outputs(tmp_path):
    paths = [
        _write(tmp_path / "A.info.afreq", AFREQ_A),
        _write(tmp_path / "B.info.afreq", AFREQ_B),
        _write(tmp_path / "D.info.afreq", AFREQ_D),
    ]
    outputs = run_subgenome_maf_plots(paths, tmp_path / "all.variant.maf")
    for key in ("summary", "dist", "dist_log"):
        assert outputs[key].exists()
    summary = pd.read_csv(outputs["summary"], sep="\t")
    assert list(summary["subgenome"]) == ["A", "B", "D"]
    assert int(summary.loc[0, "n"]) == 2


def test_run_subgenome_tajima_plots_writes_outputs(tmp_path):
    paths = [
        _write(tmp_path / "A.variant.tajima_d.Tajima.D", TAJIMA_A),
        _write(tmp_path / "B.variant.tajima_d.Tajima.D", TAJIMA_B),
        _write(tmp_path / "D.variant.tajima_d.Tajima.D", TAJIMA_D),
    ]
    outputs = run_subgenome_tajima_plots(paths, tmp_path / "all.variant.tajima_d")
    for key in ("summary", "dist", "info"):
        assert outputs[key].exists()
    info = pd.read_csv(outputs["info"], sep="\t")
    assert set(info["subgenome"]) == {"A", "B", "D"}
    assert len(info) == 6
