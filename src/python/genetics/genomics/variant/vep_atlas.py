"""Subgenome atlas of unique-site VEP consequences, frequencies, and window counts.

Genome-wide tables come from ``test_wap_refchr``; exon-interval tables from
``test_wap_exon`` (GFF ``exon``, including UTRs). One variant keeps its most
severe VEP consequence. Density plots that lack callable bases are labelled as
counts per reference window, not per callable Mb.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import numpy as np
import pandas as pd

from genetics.genomics.variant.coding_mask import (
    MISSENSE_TERMS,
    PLOF_HIGH_TERMS,
    SEVERITY_ORDER,
    VEP_META_PREFIX,
    most_severe_consequence,
)
from infra.utils.graph import (
    LEGEND_FONT_SIZE,
    SUBGENOME_COLORS,
    SUBGENOME_ORDER,
    SUBGENOME_PURPLE,
    TICK_FONT_SIZE,
    TITLE_FONT_SIZE,
    TWO_CATEGORY_BLUE,
    TWO_CATEGORY_RED,
    X_LABEL_FONT_SIZE,
    Y_LABEL_FONT_SIZE,
    plot_circos_chromosome_tracks,
    plot_heatmap_custom,
    plot_stacked_proportion_bar,
)
from infra.utils.io import load_df_from_tsv, save_df_to_tsv, save_thresholds
from infra.wheat.ref_v1 import get_chromosomes_of_subgenome, get_ref_chromosome_length

CONSEQUENCE_GROUPS = (
    "predicted LoF",
    "splice-related",
    "missense",
    "synonymous",
    "UTR",
    "intron",
    "upstream/downstream",
    "intergenic",
    "other",
)
IMPACT_ORDER = ("HIGH", "MODERATE", "LOW", "MODIFIER")
FREQ_CLASSES = (
    "singleton",
    "doubleton",
    "MAC 3-10",
    "MAF <0.1%",
    "0.1%-0.5%",
    "0.5%-1%",
    "1%-5%",
    ">=5%",
)
HEATMAP_FUNCTION_ROWS = (
    "intergenic",
    "upstream/downstream",
    "intron",
    "UTR",
    "synonymous",
    "missense",
    "splice-related",
    "predicted LoF",
)
MAF_RARE = 0.01
MAF_COMMON = 0.05
DEFAULT_BIN_BP = 10_000_000
PERICENTROMERE_FLANK_BP = 50_000_000
COMMON_COLOR = TWO_CATEGORY_BLUE
RARE_COLOR = "#e69f00"
FUNCTIONAL_COLOR = TWO_CATEGORY_RED
GENE_COLOR = SUBGENOME_PURPLE
REGION_COLORS = {
    "Arm": "#d9d9d9",
    "Pericentromere": "#56B4E9",
    "Centromere": "#1a1a1a",
}
_SPLICE_REGION_TERMS = frozenset(
    {
        "splice_region_variant",
        "splice_donor_5th_base_variant",
        "splice_donor_region_variant",
        "splice_polypyrimidine_tract_variant",
    }
)
_HIGH_TERMS = PLOF_HIGH_TERMS | frozenset({"stop_lost", "start_lost", "transcript_ablation"})
_MODERATE_TERMS = MISSENSE_TERMS | frozenset(
    {"inframe_insertion", "inframe_deletion", "protein_altering_variant"}
)
_LOW_TERMS = frozenset(
    {
        "synonymous_variant",
        "start_retained_variant",
        "stop_retained_variant",
        "incomplete_terminal_codon_variant",
    }
) | _SPLICE_REGION_TERMS
_VEP_COLUMNS = ("Uploaded_variation", "Location", "Gene", "Consequence")


def consequence_group(term: str | float) -> str:
    token = "" if pd.isna(term) else str(term).replace("&", ",").split(",")[0].strip()
    if token in PLOF_HIGH_TERMS:
        return "predicted LoF"
    if token in _SPLICE_REGION_TERMS:
        return "splice-related"
    if token in MISSENSE_TERMS:
        return "missense"
    if token == "synonymous_variant":
        return "synonymous"
    if "UTR" in token:
        return "UTR"
    if "intron" in token:
        return "intron"
    if token in {"upstream_gene_variant", "downstream_gene_variant"}:
        return "upstream/downstream"
    if "intergenic" in token:
        return "intergenic"
    return "other"


def impact_from_consequence(term: str | float) -> str:
    token = "" if pd.isna(term) else str(term).replace("&", ",").split(",")[0].strip()
    if token in _HIGH_TERMS:
        return "HIGH"
    if token in _MODERATE_TERMS:
        return "MODERATE"
    if token in _LOW_TERMS:
        return "LOW"
    return "MODIFIER"


def frequency_class(*, mac: float, maf: float) -> str:
    if pd.isna(mac) or mac < 0:
        return ""
    mac_i = int(mac)
    if mac_i == 1:
        return "singleton"
    if mac_i == 2:
        return "doubleton"
    if 3 <= mac_i <= 10:
        return "MAC 3-10"
    if pd.isna(maf):
        return ""
    if maf < 0.001:
        return "MAF <0.1%"
    if maf < 0.005:
        return "0.1%-0.5%"
    if maf < 0.01:
        return "0.5%-1%"
    if maf < 0.05:
        return "1%-5%"
    return ">=5%"


def _read_vep_columns(path: str | Path) -> pd.DataFrame:
    meta_lines = 0
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(VEP_META_PREFIX):
                break
            meta_lines += 1
    wanted = set(_VEP_COLUMNS)

    def _keep(name: str) -> bool:
        return str(name).lstrip("#") in wanted

    frame = pd.read_csv(path, sep="\t", skiprows=meta_lines, usecols=_keep, dtype=str)
    frame.columns = [str(column).lstrip("#") for column in frame.columns]
    missing = [name for name in _VEP_COLUMNS if name not in frame.columns]
    if missing:
        raise ValueError(f"{path} missing VEP columns: {missing}")
    return frame


def collapse_vep_sites(vep_path: str | Path) -> pd.DataFrame:
    """One row per ``Uploaded_variation``, keeping the most severe consequence."""

    raw = _read_vep_columns(vep_path)
    work = pd.DataFrame(
        {
            "id": raw["Uploaded_variation"].astype(str),
            "location": raw["Location"].astype(str),
            "gene": raw["Gene"].astype(str),
            "consequence": raw["Consequence"].map(
                lambda text: most_severe_consequence(str(text).replace("&", ","))
            ),
        }
    )
    rank = {term: i for i, term in enumerate(SEVERITY_ORDER)}
    work["_rank"] = work["consequence"].map(lambda term: rank.get(term, 1000))
    work = work.sort_values(["id", "_rank"]).drop_duplicates("id", keep="first")
    work["group"] = work["consequence"].map(consequence_group)
    work["impact"] = work["consequence"].map(impact_from_consequence)
    return work.drop(columns="_rank").reset_index(drop=True)


def subgenome_ref_bp(subgenome: str) -> int:
    return int(
        sum(get_ref_chromosome_length(chrom) for chrom in get_chromosomes_of_subgenome(subgenome))
    )


def subgenome_count_table(sites_by_sg: Mapping[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for sg in SUBGENOME_ORDER:
        if sg not in sites_by_sg:
            continue
        n = int(len(sites_by_sg[sg]))
        ref_bp = subgenome_ref_bp(sg)
        rows.append(
            {
                "subgenome": sg,
                "n_variants": n,
                "ref_bp": ref_bp,
                "n_per_ref_gb": n / (ref_bp / 1e9) if ref_bp else np.nan,
            }
        )
    out = pd.DataFrame(rows).set_index("subgenome")
    return out


def _proportion_table(sites_by_sg: Mapping[str, pd.DataFrame], column: str, categories: tuple[str, ...]) -> pd.DataFrame:
    rows = []
    for sg in SUBGENOME_ORDER:
        if sg not in sites_by_sg:
            continue
        counts = sites_by_sg[sg][column].value_counts()
        total = max(int(counts.sum()), 1)
        rows.append({name: float(counts.get(name, 0)) / total for name in categories})
    return pd.DataFrame(rows, index=[sg for sg in SUBGENOME_ORDER if sg in sites_by_sg])


def attach_frequency(sites: pd.DataFrame, afreq_path: str | Path, acount_path: str | Path) -> pd.DataFrame:
    afreq = pd.read_csv(afreq_path, sep="\t", dtype={"ID": str})
    afreq.columns = [str(column).lstrip("#") for column in afreq.columns]
    acount = pd.read_csv(acount_path, sep="\t", dtype={"ID": str})
    acount.columns = [str(column).lstrip("#") for column in acount.columns]
    if "ID" not in afreq.columns or "ALT_FREQS" not in afreq.columns:
        raise ValueError(f"{afreq_path} must contain ID and ALT_FREQS")
    if "ID" not in acount.columns or "ALT_CTS" not in acount.columns:
        raise ValueError(f"{acount_path} must contain ID and ALT_CTS")
    freq = pd.to_numeric(afreq["ALT_FREQS"], errors="coerce")
    maf = pd.Series(np.minimum(freq.to_numpy(dtype=float), 1.0 - freq.to_numpy(dtype=float)), index=afreq["ID"].astype(str))
    mac = pd.Series(pd.to_numeric(acount["ALT_CTS"], errors="coerce").to_numpy(dtype=float), index=acount["ID"].astype(str))
    out = sites.copy()
    out["maf"] = out["id"].map(maf)
    out["mac"] = out["id"].map(mac)
    out["freq_class"] = [
        frequency_class(mac=mac_i, maf=maf_i) for mac_i, maf_i in zip(out["mac"], out["maf"])
    ]
    out["is_common"] = out["maf"] >= MAF_COMMON
    out["is_rare"] = (out["mac"] >= 2) & (out["maf"] < MAF_RARE)
    if "group" in out.columns:
        out["is_rare_functional"] = out["is_rare"] & (out["group"] == "predicted LoF")
    else:
        out["is_rare_functional"] = False
    return out


def _stacked_bars(ax, prop_df: pd.DataFrame, colors: list[str]) -> None:
    work = prop_df.fillna(0.0)
    bottom = np.zeros(len(work))
    x_pos = np.arange(len(work))
    for column, color in zip(work.columns, colors):
        vals = work[column].to_numpy(dtype=float)
        ax.bar(x_pos, vals, bottom=bottom, label=str(column), color=color, width=0.65, edgecolor="white", linewidth=0.4)
        bottom = bottom + vals
    ax.set_ylim(0.0, 1.0)
    ax.set_xticks(x_pos)
    ax.set_xticklabels([str(x) for x in work.index], fontsize=TICK_FONT_SIZE)
    ax.set_ylabel("Fraction of variants", fontsize=Y_LABEL_FONT_SIZE)
    ax.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_vep_atlas_composition(
    counts: pd.DataFrame,
    genome_groups: pd.DataFrame,
    exon_groups: pd.DataFrame,
    genome_impact: pd.DataFrame,
    filename: str,
) -> None:
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_style("white")
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    sg = list(counts.index)
    colors = [SUBGENOME_COLORS[name] for name in sg]
    axes[0, 0].bar(sg, counts["n_per_ref_gb"].to_numpy(dtype=float), color=colors, width=0.65)
    axes[0, 0].set_title("Variants per reference-genome Gb", fontsize=TITLE_FONT_SIZE)
    axes[0, 0].set_ylabel("Sites / Gb", fontsize=Y_LABEL_FONT_SIZE)
    axes[0, 0].set_xlabel("Subgenome", fontsize=X_LABEL_FONT_SIZE)
    axes[0, 0].tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    for i, name in enumerate(sg):
        axes[0, 0].text(
            i,
            counts.loc[name, "n_per_ref_gb"],
            f"n={int(counts.loc[name, 'n_variants']):,}",
            ha="center",
            va="bottom",
            fontsize=TICK_FONT_SIZE,
        )
    axes[0, 0].spines["top"].set_visible(False)
    axes[0, 0].spines["right"].set_visible(False)

    group_colors = sns.color_palette("deep", len(CONSEQUENCE_GROUPS))
    _stacked_bars(axes[0, 1], genome_groups[list(CONSEQUENCE_GROUPS)], group_colors)
    axes[0, 1].set_title("Genome-wide VEP consequence", fontsize=TITLE_FONT_SIZE)
    axes[0, 1].legend(fontsize=LEGEND_FONT_SIZE, bbox_to_anchor=(1.02, 1.0), loc="upper left", frameon=False)

    _stacked_bars(axes[1, 0], exon_groups[list(CONSEQUENCE_GROUPS)], group_colors)
    axes[1, 0].set_title("Exon-interval VEP consequence (includes UTR)", fontsize=TITLE_FONT_SIZE)
    axes[1, 0].set_xlabel("Subgenome", fontsize=X_LABEL_FONT_SIZE)

    impact_colors = ["#c44e52", "#e69f00", "#4c72b0", "#bdbdbd"]
    _stacked_bars(axes[1, 1], genome_impact[list(IMPACT_ORDER)], impact_colors)
    axes[1, 1].set_title("VEP IMPACT (from consequence)", fontsize=TITLE_FONT_SIZE)
    axes[1, 1].legend(fontsize=LEGEND_FONT_SIZE, loc="upper right", frameon=False)

    fig.suptitle(
        "Distribution of genetic variants and predicted consequences in the A, B and D subgenomes",
        fontsize=TITLE_FONT_SIZE,
    )
    fig.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _freq_function_matrix(sites: pd.DataFrame) -> pd.DataFrame:
    work = sites.loc[sites["freq_class"].isin(FREQ_CLASSES) & sites["group"].isin(HEATMAP_FUNCTION_ROWS)]
    if work.empty:
        return pd.DataFrame(0.0, index=HEATMAP_FUNCTION_ROWS, columns=FREQ_CLASSES)
    counts = work.groupby(["group", "freq_class"]).size().unstack(fill_value=0)
    counts = counts.reindex(index=HEATMAP_FUNCTION_ROWS, columns=FREQ_CLASSES, fill_value=0)
    totals = counts.sum(axis=0).replace(0, np.nan)
    return counts.div(totals, axis=1).fillna(0.0)


def add_rare_enrichment(windows: pd.DataFrame) -> pd.DataFrame:
    """log2 fold of the window rare fraction versus the genome-wide rare fraction."""

    out = windows.copy()
    n_all = float(out["n_variants"].sum())
    n_rare = float(out["n_rare"].sum())
    global_frac = n_rare / n_all if n_all else np.nan
    local = out["n_rare"].to_numpy(dtype=float) / np.where(
        out["n_variants"].to_numpy(dtype=float) > 0,
        out["n_variants"].to_numpy(dtype=float),
        np.nan,
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        out["rare_enrichment"] = np.log2(local / global_frac)
    return out


def _ordered_windows(windows: pd.DataFrame, subgenome: str) -> pd.DataFrame:
    from genetics.genomics.variant.circos_ld_map import WHEAT_CIRCOS_CHR_ORDER

    chroms = [name for name in WHEAT_CIRCOS_CHR_ORDER if str(name).endswith(subgenome)]
    rank = {name: i for i, name in enumerate(chroms)}
    part = windows.loc[windows["subgenome"].astype(str) == subgenome].copy()
    part["_ord"] = part["chrom"].astype(str).map(rank)
    return part.sort_values(["_ord", "bin_start"], kind="mergesort").drop(columns="_ord")


def _draw_window_bars(ax, part: pd.DataFrame, values, colors, ylabel: str) -> None:
    x_pos = np.arange(len(part))
    ax.bar(x_pos, np.nan_to_num(values, nan=0.0), color=colors, width=1.0, linewidth=0)
    ax.axhline(0.0, color="black", linewidth=0.6)
    ax.set_ylabel(ylabel, fontsize=Y_LABEL_FONT_SIZE)
    ax.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xticks([])
    if len(part):
        chrom_arr = part["chrom"].astype(str).to_numpy()
        changes = np.concatenate(([True], chrom_arr[1:] != chrom_arr[:-1]))
        y_top = ax.get_ylim()[1]
        for idx, chrom in zip(np.flatnonzero(changes), chrom_arr[changes]):
            ax.axvline(idx - 0.5, color="#bdbdbd", linewidth=0.4)
            ax.text(idx, y_top, chrom, fontsize=TICK_FONT_SIZE, ha="left", va="top")


def plot_common_rare_mirror(windows: pd.DataFrame, filename: str) -> None:
    """Common counts above the axis and rare counts below, one panel per subgenome."""

    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(16, 10), sharex=False)
    for ax, sg in zip(axes, SUBGENOME_ORDER):
        part = _ordered_windows(windows, sg)
        common = part["n_common"].to_numpy(dtype=float) if len(part) else np.array([])
        rare = part["n_rare"].to_numpy(dtype=float) if len(part) else np.array([])
        x_pos = np.arange(len(part))
        ax.bar(x_pos, common, color=COMMON_COLOR, width=1.0, linewidth=0, label="Common")
        ax.bar(x_pos, -rare, color=RARE_COLOR, width=1.0, linewidth=0, label="Rare")
        ax.axhline(0.0, color="black", linewidth=0.6)
        ax.set_ylabel(f"{sg} count / window", fontsize=Y_LABEL_FONT_SIZE)
        ax.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xticks([])
        if len(part):
            chrom_arr = part["chrom"].astype(str).to_numpy()
            changes = np.concatenate(([True], chrom_arr[1:] != chrom_arr[:-1]))
            y_top = ax.get_ylim()[1]
            for idx, chrom in zip(np.flatnonzero(changes), chrom_arr[changes]):
                ax.axvline(idx - 0.5, color="#bdbdbd", linewidth=0.4)
                ax.text(idx, y_top, chrom, fontsize=TICK_FONT_SIZE, ha="left", va="top")
        if sg == "A":
            ax.legend(fontsize=LEGEND_FONT_SIZE, loc="upper right", frameon=False)
    axes[-1].set_xlabel("10 Mb windows along each subgenome", fontsize=X_LABEL_FONT_SIZE)
    fig.suptitle(
        "Common versus rare variant counts per 10 Mb window (not callable-normalised)",
        fontsize=TITLE_FONT_SIZE,
    )
    fig.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_rare_enrichment(windows: pd.DataFrame, filename: str) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(16, 10), sharex=False)
    for ax, sg in zip(axes, SUBGENOME_ORDER):
        part = _ordered_windows(windows, sg)
        values = part["rare_enrichment"].to_numpy(dtype=float) if len(part) else np.array([])
        colors = [RARE_COLOR if (pd.notna(v) and v >= 0) else COMMON_COLOR for v in values]
        _draw_window_bars(ax, part, values, colors, f"{sg} rare E")
    axes[-1].set_xlabel("10 Mb windows along each subgenome", fontsize=X_LABEL_FONT_SIZE)
    fig.suptitle(
        r"Rare-variant enrichment $\log_2$ (window rare fraction / genome-wide rare fraction)",
        fontsize=TITLE_FONT_SIZE,
    )
    fig.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_window_subgenome_violins(windows: pd.DataFrame, filename: str) -> None:
    """A/B/D window-count distributions (not callable-normalised)."""

    import matplotlib.pyplot as plt
    import seaborn as sns

    metrics = (
        ("n_common", "Common / 10 Mb window"),
        ("n_rare", "Rare / 10 Mb window"),
        ("n_rare_functional", "Exon rare predicted LoF / 10 Mb window"),
    )
    sns.set_style("white")
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharex=True)
    work = windows.loc[windows["subgenome"].isin(SUBGENOME_ORDER)].copy()
    for ax, (column, title) in zip(axes, metrics):
        sns.violinplot(
            data=work,
            x="subgenome",
            y=column,
            hue="subgenome",
            order=list(SUBGENOME_ORDER),
            palette=SUBGENOME_COLORS,
            inner=None,
            cut=0,
            legend=False,
            ax=ax,
        )
        sns.boxplot(
            data=work,
            x="subgenome",
            y=column,
            order=list(SUBGENOME_ORDER),
            width=0.2,
            color="white",
            fliersize=2,
            ax=ax,
        )
        ax.set_title(title, fontsize=TITLE_FONT_SIZE)
        ax.set_xlabel("Subgenome", fontsize=X_LABEL_FONT_SIZE)
        ax.set_ylabel("Count / window", fontsize=Y_LABEL_FONT_SIZE)
        ax.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle(
        "Subgenome window counts are not corrected for callable bases",
        fontsize=TITLE_FONT_SIZE,
    )
    fig.tight_layout()
    fig.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _window_table(sites: pd.DataFrame, bin_bp: int) -> pd.DataFrame:
    work = sites.dropna(subset=["chrom", "pos"]).copy()
    work["bin_start"] = (work["pos"].astype(int) // bin_bp) * bin_bp
    grouped = work.groupby(["chrom", "bin_start"], sort=True)
    rows = []
    for (chrom, start), part in grouped:
        rows.append(
            {
                "chrom": chrom,
                "bin_start": int(start),
                "bin_end": int(start) + bin_bp,
                "subgenome": str(chrom)[-1] if str(chrom) else "",
                "n_variants": int(len(part)),
                "n_common": int(part["is_common"].sum()) if "is_common" in part else 0,
                "n_rare": int(part["is_rare"].sum()) if "is_rare" in part else 0,
                "n_rare_functional": int(part["is_rare_functional"].sum()) if "is_rare_functional" in part else 0,
            }
        )
    return pd.DataFrame(rows)


def _attach_pvar_coordinates(sites: pd.DataFrame, pvar_path: str | Path) -> pd.DataFrame:
    from genetics.genomics.variant.circos_ld_map import load_block_snps_from_pvar

    coords = load_block_snps_from_pvar(str(pvar_path)).set_index("id")
    out = sites.copy()
    out["chrom"] = out["id"].map(coords["chrom"])
    out["pos"] = out["id"].map(coords["pos"])
    return out


def chromosome_region_intervals(
    *,
    flank_bp: int = PERICENTROMERE_FLANK_BP,
    chrom_lengths: Mapping[str, int] | None = None,
) -> pd.DataFrame:
    """Arm / pericentromere / centromere intervals on IWGSC v1 chromosomes."""

    from genetics.genomics.variant.circos_ld_map import WHEAT_CIRCOS_CHR_ORDER, wheat_circos_chrom_lengths
    from infra.wheat.ref_v1 import get_centromere_end, get_centromere_start

    lengths = dict(chrom_lengths) if chrom_lengths is not None else wheat_circos_chrom_lengths()
    flank = int(flank_bp)
    rows: list[dict] = []
    for chrom in WHEAT_CIRCOS_CHR_ORDER:
        if chrom not in lengths:
            continue
        length = int(lengths[chrom])
        c0 = get_centromere_start(chrom)
        c1 = get_centromere_end(chrom)
        if c0 is None or c1 is None:
            rows.append(
                {
                    "chrom": chrom,
                    "start_pos": 0,
                    "end_pos": length,
                    "label": "Arm",
                    "color": REGION_COLORS["Arm"],
                }
            )
            continue
        c0_i = int(c0)
        c1_i = int(c1)
        peri0 = max(0, c0_i - flank)
        peri1 = min(length, c1_i + flank)
        for label, start, end in (
            ("Arm", 0, peri0),
            ("Pericentromere", peri0, c0_i),
            ("Centromere", c0_i, c1_i),
            ("Pericentromere", c1_i, peri1),
            ("Arm", peri1, length),
        ):
            if end <= start:
                continue
            rows.append(
                {
                    "chrom": chrom,
                    "start_pos": start,
                    "end_pos": end,
                    "label": label,
                    "color": REGION_COLORS[label],
                }
            )
    return pd.DataFrame(rows)


def _histogram_values(windows: pd.DataFrame, column: str) -> pd.DataFrame:
    out = windows.copy()
    out["value"] = pd.to_numeric(out[column], errors="coerce").fillna(0.0)
    return out


def rare_common_share_table(counts: pd.DataFrame) -> pd.DataFrame:
    """Within each subgenome, the share of sites that are common or rare.

    The denominator is ``n_common + n_rare`` on that subgenome. Sites in neither
    class are left out of the bar.
    """
    work = counts.copy()
    if "subgenome" in work.columns:
        work = work.set_index("subgenome")
    order = [sg for sg in ("A", "B", "D") if sg in work.index]
    if not order:
        raise ValueError("counts has no A, B, or D subgenome rows")
    rare = pd.to_numeric(work.loc[order, "n_rare"], errors="coerce")
    common = pd.to_numeric(work.loc[order, "n_common"], errors="coerce")
    total = (rare + common).replace(0, np.nan)
    return pd.DataFrame(
        {
            "Common": (common / total).to_numpy(dtype=float),
            "Rare": (rare / total).to_numpy(dtype=float),
        },
        index=order,
    )


def plot_subgenome_rare_common_share(
    counts: pd.DataFrame,
    filename: str,
    *,
    title: str = "Rare and common share by subgenome",
) -> pd.DataFrame:
    """One 100% stacked bar per subgenome: common (red) under rare (blue)."""
    prop = rare_common_share_table(counts)
    plot_stacked_proportion_bar(
        prop,
        title,
        filename,
        x_label="Subgenome",
        y_label="Share of rare and common sites",
        stack_col_title="Class",
        palette_map={"Common": TWO_CATEGORY_RED, "Rare": TWO_CATEGORY_BLUE},
        annotate=True,
        percent_axis=True,
        figure_size=(7, 5),
    )
    return prop


def plot_vep_atlas_circos(
    windows: str | Path | pd.DataFrame,
    filename: str,
    bin_bp: int = DEFAULT_BIN_BP,
    *,
    ld_blocks: str | Path | pd.DataFrame | None = None,
    gene_gff: str | Path | None = None,
    peri_flank_bp: int = PERICENTROMERE_FLANK_BP,
) -> None:
    """Circos: thin cytobands, LD blocks, then genes / LoF / common / rare histograms."""

    from genetics.genomics.gene.gff_bed import bin_gff_feature_counts
    from genetics.genomics.variant.circos_ld_map import (
        LD_BLOCK_COLOR,
        WHEAT_CIRCOS_CHR_ORDER,
        wheat_circos_chrom_lengths,
    )

    if isinstance(windows, (str, Path)):
        windows = load_df_from_tsv(str(windows))
    if isinstance(ld_blocks, (str, Path)):
        ld_blocks = load_df_from_tsv(str(ld_blocks))

    lengths = wheat_circos_chrom_lengths()
    regions = chromosome_region_intervals(flank_bp=peri_flank_bp, chrom_lengths=lengths)
    tracks = [
        {
            "kind": "karyotype",
            "r_inner": 0.945,
            "r_outer": 0.972,
            "intervals": regions,
            "base_color": REGION_COLORS["Arm"],
            "show_ticks": True,
            "show_legend": False,
        }
    ]
    if ld_blocks is not None and not ld_blocks.empty:
        tracks.append(
            {
                "kind": "karyotype",
                "label": "LD block",
                "r_inner": 0.875,
                "r_outer": 0.930,
                "intervals": ld_blocks,
                "color": LD_BLOCK_COLOR,
                "interval_color": LD_BLOCK_COLOR,
                "base_color": "#f0f0f0",
                "show_ticks": False,
            }
        )
    if gene_gff is not None:
        genes = bin_gff_feature_counts(gene_gff, bin_bp=int(bin_bp), feature_type="gene")
        tracks.append(
            {
                "kind": "histogram",
                "r_inner": 0.73,
                "r_outer": 0.86,
                "color": "#e69f00",
                "label": "Genes",
                "data": genes,
            }
        )
    tracks.extend(
        [
            {
                "kind": "histogram",
                "r_inner": 0.58,
                "r_outer": 0.71,
                "color": SUBGENOME_PURPLE,
                "label": "Predicted LoF",
                "data": _histogram_values(windows, "n_rare_functional"),
            },
            {
                "kind": "histogram",
                "r_inner": 0.43,
                "r_outer": 0.56,
                "color": TWO_CATEGORY_RED,
                "label": "Common",
                "data": _histogram_values(windows, "n_common"),
            },
            {
                "kind": "histogram",
                "r_inner": 0.28,
                "r_outer": 0.41,
                "color": TWO_CATEGORY_BLUE,
                "label": "Rare",
                "data": _histogram_values(windows, "n_rare"),
            },
        ]
    )
    plot_circos_chromosome_tracks(
        lengths,
        tracks,
        filename,
        chrom_order=WHEAT_CIRCOS_CHR_ORDER,
        legend_loc="center",
        wrap_gap_frac=0.0,
        chrom_label_pad=0.06,
        figure_size=(13, 13),
    )


def _plot_window_circos(
    windows: pd.DataFrame,
    filename: str,
    bin_bp: int,
    *,
    ld_blocks: str | Path | pd.DataFrame | None = None,
    gene_gff: str | Path | None = None,
) -> None:
    plot_vep_atlas_circos(
        windows,
        filename,
        bin_bp=bin_bp,
        ld_blocks=ld_blocks,
        gene_gff=gene_gff,
    )


def merge_genome_exon_windows(genome_windows: pd.DataFrame, exon_windows: pd.DataFrame) -> pd.DataFrame:
    """Keep genome common/rare counts; replace LoF counts with the exon-interval table."""

    keys = ["chrom", "bin_start"]
    genome = genome_windows.drop(columns=["n_rare_functional"], errors="ignore").copy()
    lof = (
        exon_windows.reindex(columns=["chrom", "bin_start", "n_rare_functional"])
        .dropna(subset=keys)
        .groupby(keys, as_index=False)["n_rare_functional"]
        .sum()
    )
    out = genome.merge(lof, on=keys, how="left")
    out["n_rare_functional"] = out["n_rare_functional"].fillna(0).astype(int)
    return out


def _subset_keep_ids(
    sites_by_sg: dict[str, pd.DataFrame],
    keep_ids: set[str] | None,
) -> dict[str, pd.DataFrame]:
    if keep_ids is None:
        return sites_by_sg
    keep = {str(variant_id) for variant_id in keep_ids}
    return {
        sg: frame.loc[frame["id"].astype(str).isin(keep)].reset_index(drop=True)
        for sg, frame in sites_by_sg.items()
    }


def run_vep_atlas(
    genome_vep: Mapping[str, str],
    exon_vep: Mapping[str, str],
    output_prefix: str,
    *,
    genome_afreq: Mapping[str, str] | None = None,
    genome_acount: Mapping[str, str] | None = None,
    genome_pvar: Mapping[str, str] | None = None,
    exon_afreq: Mapping[str, str] | None = None,
    exon_acount: Mapping[str, str] | None = None,
    exon_pvar: Mapping[str, str] | None = None,
    ld_blocks: str | Path | pd.DataFrame | None = None,
    gene_gff: str | Path | None = None,
    bin_bp: int = DEFAULT_BIN_BP,
    genome_keep_ids: set[str] | None = None,
    exon_keep_ids: set[str] | None = None,
) -> dict:
    genome_sites = _subset_keep_ids(
        {sg: collapse_vep_sites(path) for sg, path in genome_vep.items()},
        genome_keep_ids,
    )
    exon_sites = _subset_keep_ids(
        {sg: collapse_vep_sites(path) for sg, path in exon_vep.items()},
        exon_keep_ids,
    )
    counts = subgenome_count_table(genome_sites)
    genome_groups = _proportion_table(genome_sites, "group", CONSEQUENCE_GROUPS)
    exon_groups = _proportion_table(exon_sites, "group", CONSEQUENCE_GROUPS)
    genome_impact = _proportion_table(genome_sites, "impact", IMPACT_ORDER)

    if genome_afreq and genome_acount:
        for sg in list(genome_sites):
            if sg in genome_afreq and sg in genome_acount:
                genome_sites[sg] = attach_frequency(genome_sites[sg], genome_afreq[sg], genome_acount[sg])
                counts.loc[sg, "n_rare"] = int(genome_sites[sg]["is_rare"].sum())
                counts.loc[sg, "n_common"] = int(genome_sites[sg]["is_common"].sum())

    if exon_afreq and exon_acount:
        for sg in list(exon_sites):
            if sg in exon_afreq and sg in exon_acount:
                exon_sites[sg] = attach_frequency(exon_sites[sg], exon_afreq[sg], exon_acount[sg])
                counts.loc[sg, "n_rare_functional"] = int(exon_sites[sg]["is_rare_functional"].sum())

    plot_vep_atlas_composition(
        counts,
        genome_groups,
        exon_groups,
        genome_impact,
        f"{output_prefix}.composition.png",
    )
    save_df_to_tsv(counts.reset_index(), f"{output_prefix}.counts.tsv")
    save_df_to_tsv(genome_groups.reset_index(), f"{output_prefix}.consequence.genome.tsv")
    save_df_to_tsv(exon_groups.reset_index(), f"{output_prefix}.consequence.exon.tsv")
    save_df_to_tsv(genome_impact.reset_index(), f"{output_prefix}.impact.genome.tsv")

    freq_sites = pd.concat(
        [frame.assign(subgenome=sg) for sg, frame in genome_sites.items() if "freq_class" in frame.columns],
        ignore_index=True,
    ) if any("freq_class" in frame.columns for frame in genome_sites.values()) else pd.DataFrame()
    if not freq_sites.empty:
        matrix = _freq_function_matrix(freq_sites)
        plot_heatmap_custom(
            matrix.to_numpy(dtype=float),
            list(matrix.columns),
            list(matrix.index),
            title="Consequence group within frequency class (genome-wide unique sites)",
            filename=f"{output_prefix}.freq_function.png",
            xlabel="Frequency class",
            ylabel="Most-severe VEP group",
            cmap="YlOrRd",
            cbar_label="Fraction of class",
            vmin=0.0,
            vmax=1.0,
        )
        save_df_to_tsv(matrix.reset_index(), f"{output_prefix}.freq_function.tsv")

    if genome_pvar:
        window_frames = []
        for sg, frame in genome_sites.items():
            if sg not in genome_pvar:
                continue
            located = _attach_pvar_coordinates(frame, genome_pvar[sg])
            window_frames.append(_window_table(located, bin_bp))
        if window_frames:
            windows = pd.concat(window_frames, ignore_index=True)
            if exon_pvar:
                exon_frames = []
                for sg, frame in exon_sites.items():
                    if sg not in exon_pvar:
                        continue
                    located = _attach_pvar_coordinates(frame, exon_pvar[sg])
                    exon_frames.append(_window_table(located, bin_bp))
                if exon_frames:
                    windows = merge_genome_exon_windows(windows, pd.concat(exon_frames, ignore_index=True))
            windows = add_rare_enrichment(windows)
            save_df_to_tsv(windows, f"{output_prefix}.windows.tsv")
            _plot_window_circos(
                windows,
                f"{output_prefix}.circos_windows.png",
                bin_bp,
                ld_blocks=ld_blocks,
                gene_gff=gene_gff,
            )
            plot_common_rare_mirror(windows, f"{output_prefix}.mirror.png")
            plot_rare_enrichment(windows, f"{output_prefix}.enrichment.png")
            plot_window_subgenome_violins(windows, f"{output_prefix}.windows.violin.png")

    stats = {
        "n_genome_sites": int(sum(len(frame) for frame in genome_sites.values())),
        "n_exon_sites": int(sum(len(frame) for frame in exon_sites.values())),
        "bin_bp": int(bin_bp),
        "note": "n_per_ref_gb uses IWGSC v1 chromosome length, not callable bases; n_rare_functional is rare predicted LoF from exon-interval VEP",
    }
    save_thresholds(stats, f"{output_prefix}.th.tsv")
    print(
        f"vep_atlas: genome={stats['n_genome_sites']} exon={stats['n_exon_sites']} "
        f"wrote {output_prefix}.composition.png"
    )
    return stats
