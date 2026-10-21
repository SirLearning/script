"""Project coding-consequence masks (pLoF / missense) from a VEP-like table."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from infra.utils.graph import plot_bar_chart
from infra.utils.io import save_df_to_tsv, save_thresholds

VEP_META_PREFIX = "##"

MASK_PLOF_HIGH = "plof_high"
MASK_MIS_SEVERE = "mis_severe"
MASK_MIS_CONSENSUS = "mis_consensus"
MASK_PROTEIN_ALTERING = "protein_altering"
MASK_OTHER = "other"

PLOF_HIGH_TERMS = frozenset(
    {
        "stop_gained",
        "frameshift_variant",
        "splice_acceptor_variant",
        "splice_donor_variant",
    }
)
MISSENSE_TERMS = frozenset({"missense_variant"})
PROTEIN_ALTERING_TERMS = PLOF_HIGH_TERMS | MISSENSE_TERMS | frozenset(
    {
        "stop_lost",
        "start_lost",
        "inframe_insertion",
        "inframe_deletion",
        "protein_altering_variant",
    }
)
SEVERITY_ORDER = (
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "missense_variant",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
    "synonymous_variant",
)

WHEAT_GENE_RE = re.compile(r"TraesCS(\d)([ABD])", re.IGNORECASE)
PROTEIN_POS_RE = re.compile(r"(?P<pos>\d+)(?:-\d+)?/(?P<length>\d+)")
TERMINAL_FRAC = 0.95
DELETERIOUS_SIFT = frozenset({"deleterious", "deleterious_low_confidence"})
DAMAGING_POLYPHEN = frozenset({"probably_damaging", "possibly_damaging"})


def most_severe_consequence(consequence: str | float) -> str:
    if pd.isna(consequence):
        return ""
    parts = [part.strip() for part in str(consequence).split(",") if part.strip()]
    if not parts:
        return ""
    for term in SEVERITY_ORDER:
        if term in parts:
            return term
    return parts[0]


def protein_position_frac(protein_position: str | float) -> float | None:
    if pd.isna(protein_position):
        return None
    match = PROTEIN_POS_RE.search(str(protein_position))
    if match is None:
        return None
    length = int(match.group("length"))
    if length <= 0:
        return None
    return int(match.group("pos")) / length


def subgenome_from_gene_id(gene_id: str | float) -> str:
    if pd.isna(gene_id):
        return ""
    match = WHEAT_GENE_RE.search(str(gene_id))
    if match is None:
        return ""
    return match.group(2).upper()


def _is_deleterious_sift(value: str | float) -> bool:
    if pd.isna(value):
        return False
    token = str(value).split("(")[0].strip().lower()
    return token in DELETERIOUS_SIFT


def _is_damaging_polyphen(value: str | float) -> bool:
    if pd.isna(value):
        return False
    token = str(value).split("(")[0].strip().lower()
    return token in DAMAGING_POLYPHEN


def assign_coding_mask(row: pd.Series, *, terminal_frac: float = TERMINAL_FRAC) -> str:
    term = most_severe_consequence(row.get("Consequence", ""))
    frac = protein_position_frac(row.get("Protein_position", pd.NA))
    terminal = frac is not None and frac >= terminal_frac

    if term in PLOF_HIGH_TERMS:
        if terminal and term in {"stop_gained", "frameshift_variant"}:
            return MASK_PROTEIN_ALTERING
        return MASK_PLOF_HIGH

    if term in MISSENSE_TERMS:
        sift_hit = _is_deleterious_sift(row.get("SIFT", pd.NA))
        pph_hit = _is_damaging_polyphen(row.get("PolyPhen", pd.NA))
        if sift_hit and pph_hit:
            return MASK_MIS_CONSENSUS
        if sift_hit or pph_hit:
            return MASK_MIS_SEVERE
        return MASK_PROTEIN_ALTERING

    if term in PROTEIN_ALTERING_TERMS:
        return MASK_PROTEIN_ALTERING
    return MASK_OTHER


def classify_vep_table(
    vep_df: pd.DataFrame,
    *,
    terminal_frac: float = TERMINAL_FRAC,
) -> pd.DataFrame:
    if vep_df is None or vep_df.empty:
        raise ValueError("Empty VEP table")
    work = vep_df.copy()
    if "Consequence" not in work.columns:
        raise ValueError("VEP table missing Consequence column")
    gene_col = "Gene" if "Gene" in work.columns else ("SYMBOL" if "SYMBOL" in work.columns else None)
    work["consequence_severe"] = work["Consequence"].map(most_severe_consequence)
    work["coding_mask"] = work.apply(lambda row: assign_coding_mask(row, terminal_frac=terminal_frac), axis=1)
    work["terminal_truncation"] = [
        (protein_position_frac(val) or 0.0) >= terminal_frac
        and most_severe_consequence(cons) in {"stop_gained", "frameshift_variant"}
        for val, cons in zip(work.get("Protein_position", pd.Series(index=work.index)), work["Consequence"])
    ]
    if gene_col:
        work["subgenome"] = work[gene_col].map(subgenome_from_gene_id)
    else:
        work["subgenome"] = ""
    return work


def summarize_coding_masks(classified: pd.DataFrame) -> dict[str, int]:
    counts = classified["coding_mask"].value_counts().to_dict()
    stats = {
        "n_rows": int(len(classified)),
        "n_plof_high": int(counts.get(MASK_PLOF_HIGH, 0)),
        "n_mis_severe": int(counts.get(MASK_MIS_SEVERE, 0)),
        "n_mis_consensus": int(counts.get(MASK_MIS_CONSENSUS, 0)),
        "n_protein_altering": int(counts.get(MASK_PROTEIN_ALTERING, 0)),
        "n_other": int(counts.get(MASK_OTHER, 0)),
        "n_terminal_truncation": int(classified["terminal_truncation"].sum()),
    }
    for sg in ("A", "B", "D"):
        stats[f"n_subgenome_{sg}"] = int((classified["subgenome"] == sg).sum())
    return stats


def plot_coding_masks(classified: pd.DataFrame, output_prefix: str | Path) -> None:
    order = [MASK_PLOF_HIGH, MASK_MIS_CONSENSUS, MASK_MIS_SEVERE, MASK_PROTEIN_ALTERING, MASK_OTHER]
    counts = classified["coding_mask"].value_counts()
    names = [name for name in order if name in counts.index]
    values = [float(counts[name]) for name in names]
    plot_bar_chart(
        names=names,
        values=values,
        title="Coding masks (VEP consequences)",
        ylabel="Annotation rows",
        filename=f"{output_prefix}.mask.counts.png",
        ylim=(0.0, (max(values) * 1.15) if values else 1.05),
        rotate_xlabels=30,
    )


def load_vep_tab_table(vep_table: str | Path) -> pd.DataFrame:
    """Read a VEP ``--tab`` file.

    VEP writes a block of ``##`` metadata lines before the column header, which
    itself starts with a single ``#``. A generic reader takes the first metadata
    line as the header and then fails on the real columns.
    """

    path = Path(vep_table)
    meta_lines = 0
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(VEP_META_PREFIX):
                break
            meta_lines += 1
    frame = pd.read_csv(path, sep="\t", skiprows=meta_lines, low_memory=False)
    frame.columns = [str(column).lstrip("#") for column in frame.columns]
    return frame


def run_coding_mask(
    vep_table: str | Path,
    output_prefix: str | Path,
    *,
    terminal_frac: float = TERMINAL_FRAC,
) -> pd.DataFrame:
    raw = load_vep_tab_table(vep_table)
    classified = classify_vep_table(raw, terminal_frac=terminal_frac)
    save_df_to_tsv(classified, f"{output_prefix}.coding_mask.tsv")
    stats = summarize_coding_masks(classified)
    save_thresholds(stats, f"{output_prefix}.coding_mask.th.tsv")
    plot_coding_masks(classified, output_prefix)
    print(
        f"coding_mask: n={stats['n_rows']} plof_high={stats['n_plof_high']} "
        f"mis_consensus={stats['n_mis_consensus']}"
    )
    return classified
