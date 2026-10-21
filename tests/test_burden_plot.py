"""Gene-level burden Manhattan and QQ."""

from pathlib import Path

import pandas as pd

from genetics.gwas.association_plot import plot_burden_results


def test_plot_burden_results_writes_qq_and_manhattan(tmp_path: Path):
    gff = tmp_path / "genes.gff3"
    gff.write_text(
        "\n".join([
            "chr1A\tIWGSC\tgene\t100\t300\t.\t+\t.\tID=TraesCS1A02G000100",
            "chr2A\tIWGSC\tgene\t1000\t2000\t.\t+\t.\tID=TraesCS2A02G000100",
            "chr2A\tIWGSC\tgene\t5000\t6000\t.\t+\t.\tID=TraesCS2A02G000200",
        ]) + "\n",
        encoding="utf-8",
    )
    burden = tmp_path / "A.PH.burden.tsv"
    pd.DataFrame({
        "gene": ["TraesCS1A02G000100", "TraesCS2A02G000100", "TraesCS2A02G000200", "NO_SUCH_GENE"],
        "n_variants": [1, 2, 1, 1],
        "n_carriers": [4, 21, 3, 2],
        "p_value": [0.4, 1e-12, 0.2, 0.01],
    }).to_csv(burden, sep="\t", index=False)

    prefix = tmp_path / "A.PH.plof"
    placed = plot_burden_results(str(burden), str(prefix), str(gff))

    assert (tmp_path / "A.PH.plof.qq.png").is_file()
    assert (tmp_path / "A.PH.plof.manhattan.png").is_file()
    info = pd.read_csv(tmp_path / "A.PH.plof.burden.plot.tsv", sep="\t")
    assert set(info["gene"]) == {
        "TraesCS1A02G000100", "TraesCS2A02G000100", "TraesCS2A02G000200",
    }
    assert "NO_SUCH_GENE" not in set(placed["gene"])
    assert info.loc[info["gene"] == "TraesCS2A02G000100", "chrom"].item() == "chr2A"
