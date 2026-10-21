"""Per-variant site annotation from the wheat featHighRes conservation/epigenome BED.

The ``00data/11epi_anno/chr*_featHighRes.bed.gz`` files carry 30 per-base
features on IWGSC RefSeq v1.0 whole-chromosome (linear) coordinates, which is
the same system as the ``*_refchr`` pfiles. A Nextflow process streams one BED
against one pfile chromosome and writes the intersected rows; this module gives
the column contract, derives analysis classes from the raw features, and
summarizes the result.

Degeneracy encoding
-------------------
The upstream README lists ``degeneration`` as taking 0, 1 or 2 without saying
what the codes mean. Cross-tabulating 271k chr1A CDS bases against fold
degeneracy computed from the reference codons gives an unambiguous mapping::

    true fold     degen=0   degen=1   degen=2
    0-fold           2563         0    172426
    2-fold          49437        10        11
    3-fold           4722         1         0
    4-fold            629     41512         0

so ``2`` is non-degenerate and ``1`` — not ``2`` — marks the 4-fold sites used
as the neutral reference. See ``DEGENERACY_LABELS``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from infra.utils.graph import plot_categorical_bar_chart, plot_grouped_density_lines
from infra.utils.io import load_df_generic, save_df_to_tsv, save_thresholds

# Feature columns in BED order, following section #02 of the 11epi_anno README.
EPI_FEATURE_COLUMNS: tuple[str, ...] = (
    "heptamer_USscore",
    "CDS_HC",
    "CDS_LC",
    "ncRNA",
    "degeneration",
    "CGI",
    "TE",
    "triticeae_phyloP",
    "poaceae_phyloP",
    "monocots_phyloP",
    "triticeae_phastCons",
    "poaceae_phastCons",
    "monocots_phastCons",
    "GERP29",
    "CE_triticeae",
    "CE_BOP",
    "CE_poaceae",
    "CE_monocots",
    "PPH2",
    "LIST2",
    "SIFT",
    "DAP",
    "DHS",
    "ATAC",
    "MSF",
    "RNAPII",
    "pNET",
    "GRO",
    "epiRE",
    "histoneState_level",
)

KEY_COLUMNS: tuple[str, ...] = ("ID", "chrom", "pos")
SITE_ANNO_COLUMNS: tuple[str, ...] = KEY_COLUMNS + EPI_FEATURE_COLUMNS

BINARY_FEATURES = frozenset(
    {
        "CDS_HC",
        "CDS_LC",
        "ncRNA",
        "CGI",
        "TE",
        "CE_triticeae",
        "CE_BOP",
        "CE_poaceae",
        "CE_monocots",
        "DAP",
        "DHS",
        "ATAC",
        "MSF",
        "RNAPII",
        "pNET",
        "GRO",
        "epiRE",
    }
)

# Chromatin / binding tracks that together define an accessible-or-bound site.
REGULATORY_FEATURES: tuple[str, ...] = (
    "DAP",
    "DHS",
    "ATAC",
    "MSF",
    "RNAPII",
    "pNET",
    "GRO",
    "epiRE",
)

CONSERVATION_FEATURES: tuple[str, ...] = (
    "triticeae_phyloP",
    "poaceae_phyloP",
    "monocots_phyloP",
    "triticeae_phastCons",
    "poaceae_phastCons",
    "monocots_phastCons",
    "GERP29",
)

DEGENERACY_LABELS = {0: "2fold", 1: "4fold", 2: "0fold"}
DEGENERACY_ORDER = ("0fold", "2fold", "4fold", "noncoding")

ELEMENT_ORDER = ("CDS_HC", "CDS_LC", "ncRNA", "TE", "CGI", "intergenic")

DEFAULT_PHYLOP_MIN = 2.0
DEFAULT_PHASTCONS_MIN = 0.8
DEFAULT_GERP_MIN = 2.0

# Column added by add_site_classes.
COL_DEGENERACY = "degeneracy"
COL_ELEMENT = "element"
COL_CODING = "coding"
COL_CONSERVED = "conserved"
COL_REGULATORY = "regulatory"
COL_NEUTRAL = "neutral_4fold"

SITE_CLASS_COLUMNS: tuple[str, ...] = (
    COL_DEGENERACY,
    COL_ELEMENT,
    COL_CODING,
    COL_CONSERVED,
    COL_REGULATORY,
    COL_NEUTRAL,
)


def site_anno_dtypes() -> dict[str, str]:
    """Compact dtypes so a whole-subgenome table stays in memory."""

    dtypes: dict[str, str] = {"ID": "string", "chrom": "string", "pos": "int64"}
    for name in EPI_FEATURE_COLUMNS:
        dtypes[name] = "int8" if name in BINARY_FEATURES else "float32"
    dtypes["degeneration"] = "float32"
    dtypes["histoneState_level"] = "float32"
    return dtypes


def load_site_anno_table(path: str | Path) -> pd.DataFrame:
    """Read a site-annotation table that already carries a header row."""

    table = load_df_generic(str(path))
    missing = [name for name in SITE_ANNO_COLUMNS if name not in table.columns]
    if missing:
        raise ValueError(f"site_anno table missing columns: {', '.join(missing)}")
    return table


def read_site_anno_chunk(path: str | Path) -> pd.DataFrame:
    """Read one headerless per-chromosome chunk written by ``site_anno_extract``."""

    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=list(SITE_ANNO_COLUMNS),
        dtype=site_anno_dtypes(),
        na_values=["NA", "."],
    )


def load_site_anno_chunks(paths) -> pd.DataFrame:
    """Concatenate per-chromosome chunks into one subgenome table."""

    frames = [read_site_anno_chunk(path) for path in paths]
    frames = [frame for frame in frames if not frame.empty]
    if not frames:
        raise ValueError("site_anno: every chunk was empty")
    return pd.concat(frames, ignore_index=True)


def classify_degeneracy(table: pd.DataFrame) -> pd.Series:
    """Map the raw ``degeneration`` code to a fold label, ``noncoding`` off CDS."""

    codes = pd.to_numeric(table["degeneration"], errors="coerce")
    coding = pd.to_numeric(table["CDS_HC"], errors="coerce").fillna(0).astype(int) == 1
    labels = codes.map(DEGENERACY_LABELS)
    return labels.where(coding, "noncoding").fillna("noncoding").astype("string")


def classify_element(table: pd.DataFrame) -> pd.Series:
    """Assign each site to one sequence element, most specific first."""

    element = pd.Series("intergenic", index=table.index, dtype="object")
    for name in reversed(ELEMENT_ORDER[:-1]):
        flag = pd.to_numeric(table[name], errors="coerce").fillna(0).astype(int) == 1
        element = element.mask(flag, name)
    return element.astype("string")


def flag_conserved(
    table: pd.DataFrame,
    *,
    phylop_min: float = DEFAULT_PHYLOP_MIN,
    phastcons_min: float = DEFAULT_PHASTCONS_MIN,
    gerp_min: float = DEFAULT_GERP_MIN,
) -> pd.Series:
    """Conserved across monocots by phyloP and phastCons, or by GERP alone."""

    phylop = pd.to_numeric(table["monocots_phyloP"], errors="coerce")
    phastcons = pd.to_numeric(table["monocots_phastCons"], errors="coerce")
    gerp = pd.to_numeric(table["GERP29"], errors="coerce")
    by_alignment = phylop.ge(phylop_min) & phastcons.ge(phastcons_min)
    return (by_alignment | gerp.ge(gerp_min)).fillna(False)


def flag_regulatory(table: pd.DataFrame) -> pd.Series:
    """True when any chromatin accessibility or binding track is set."""

    flags = [
        pd.to_numeric(table[name], errors="coerce").fillna(0).astype(int) == 1
        for name in REGULATORY_FEATURES
        if name in table.columns
    ]
    if not flags:
        return pd.Series(False, index=table.index)
    combined = flags[0]
    for flag in flags[1:]:
        combined = combined | flag
    return combined


def add_site_classes(
    table: pd.DataFrame,
    *,
    phylop_min: float = DEFAULT_PHYLOP_MIN,
    phastcons_min: float = DEFAULT_PHASTCONS_MIN,
    gerp_min: float = DEFAULT_GERP_MIN,
) -> pd.DataFrame:
    """Return a copy with degeneracy, element, conservation and regulatory classes."""

    out = table.copy()
    out[COL_DEGENERACY] = classify_degeneracy(out)
    out[COL_ELEMENT] = classify_element(out)
    out[COL_CODING] = out[COL_ELEMENT].isin(("CDS_HC", "CDS_LC"))
    out[COL_CONSERVED] = flag_conserved(
        out,
        phylop_min=phylop_min,
        phastcons_min=phastcons_min,
        gerp_min=gerp_min,
    )
    out[COL_REGULATORY] = flag_regulatory(out)
    out[COL_NEUTRAL] = out[COL_DEGENERACY].eq("4fold") & ~out[COL_CONSERVED]
    return out


def summarize_site_anno(table: pd.DataFrame) -> dict[str, int | float]:
    """Counts per class plus median conservation, for the thresholds file."""

    stats: dict[str, int | float] = {"n_sites": int(len(table))}
    for label in DEGENERACY_ORDER:
        stats[f"n_{label}"] = int((table[COL_DEGENERACY] == label).sum())
    for label in ELEMENT_ORDER:
        stats[f"n_element_{label}"] = int((table[COL_ELEMENT] == label).sum())
    stats["n_conserved"] = int(table[COL_CONSERVED].sum())
    stats["n_regulatory"] = int(table[COL_REGULATORY].sum())
    stats["n_neutral_4fold"] = int(table[COL_NEUTRAL].sum())
    for name in ("monocots_phyloP", "monocots_phastCons", "GERP29"):
        values = pd.to_numeric(table[name], errors="coerce")
        stats[f"median_{name}"] = float(np.nanmedian(values)) if len(values) else float("nan")
    return stats


def plot_site_anno(table: pd.DataFrame, output_prefix: str | Path) -> None:
    """Class-count bars plus a conservation density split by degeneracy."""

    counts = [int((table[COL_ELEMENT] == label).sum()) for label in ELEMENT_ORDER]
    counts += [
        int(table[COL_CONSERVED].sum()),
        int(table[COL_REGULATORY].sum()),
        int(table[COL_NEUTRAL].sum()),
    ]
    plot_categorical_bar_chart(
        categories=[*ELEMENT_ORDER, "conserved", "regulatory", "neutral_4fold"],
        values=[float(value) for value in counts],
        title="Site annotation classes",
        ylabel="Variants",
        filename=f"{output_prefix}.site_anno.classes.png",
        rotate_xlabels=30,
    )

    coding = table[table[COL_DEGENERACY].isin(("0fold", "2fold", "4fold"))]
    if coding.empty:
        return
    palette = {"0fold": "#c0392b", "2fold": "#e6a23c", "4fold": "#4c72b0"}
    plot_grouped_density_lines(
        data=coding,
        x_col="monocots_phyloP",
        hue_col=COL_DEGENERACY,
        hue_order=["0fold", "2fold", "4fold"],
        palette_map=palette,
        title="Monocot phyloP by codon degeneracy",
        filename=f"{output_prefix}.site_anno.phylop.png",
        xlabel="monocots phyloP",
    )


def run_site_anno(
    chunk_paths,
    output_prefix: str | Path,
    *,
    phylop_min: float = DEFAULT_PHYLOP_MIN,
    phastcons_min: float = DEFAULT_PHASTCONS_MIN,
    gerp_min: float = DEFAULT_GERP_MIN,
) -> pd.DataFrame:
    """Merge per-chromosome chunks, derive classes, write table, thresholds, plots."""

    classed = add_site_classes(
        load_site_anno_chunks(chunk_paths),
        phylop_min=phylop_min,
        phastcons_min=phastcons_min,
        gerp_min=gerp_min,
    )
    keep = [*KEY_COLUMNS, *SITE_CLASS_COLUMNS, *CONSERVATION_FEATURES]
    save_df_to_tsv(classed[keep], f"{output_prefix}.site_anno.tsv")
    stats = summarize_site_anno(classed)
    save_thresholds(stats, f"{output_prefix}.site_anno.th.tsv")
    plot_site_anno(classed, output_prefix)
    print(
        f"site_anno: sites={stats['n_sites']} 4fold={stats['n_4fold']} "
        f"0fold={stats['n_0fold']} conserved={stats['n_conserved']} "
        f"regulatory={stats['n_regulatory']} neutral={stats['n_neutral_4fold']}"
    )
    return classed
