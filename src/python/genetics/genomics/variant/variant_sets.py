"""Analysis variant sets: coding masks intersected with cohort allele frequency.

Joins a ``coding_mask`` table (from a VEP consequence table) with a ``dual_af``
table (WGS-core vs WAP frequencies) and emits one named set per mask/frequency
combination, plus PLINK2 ``--extract`` lists for downstream burden or
association runs.

The join key is the variant ID. VEP reports it as ``Uploaded_variation`` and
PLINK2 as ``ID``; both sides must use the same variant naming, otherwise the
intersection is empty.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from genetics.genomics.variant.coding_mask import (
    MASK_MIS_CONSENSUS,
    MASK_MIS_SEVERE,
    MASK_PLOF_HIGH,
    MASK_PROTEIN_ALTERING,
)
from genetics.genomics.variant.site_anno import (
    COL_CODING,
    COL_CONSERVED,
    COL_DEGENERACY,
    COL_REGULATORY,
    SITE_CLASS_COLUMNS,
)
from infra.utils.graph import plot_bar_chart
from infra.utils.io import load_df_generic, save_df_to_tsv, save_thresholds

VARIANT_ID_CANDIDATES = ("ID", "Uploaded_variation", "SNP", "Variant")

FREQ_ANY = "any"
FREQ_RARE_WGS = "rare_wgs"
FREQ_RARE_WAP = "rare_wap"
FREQ_RARE_BOTH = "rare_both"
FREQ_COMMON_WAP = "common_wap"
FREQUENCY_CLASSES = (FREQ_ANY, FREQ_RARE_WGS, FREQ_RARE_WAP, FREQ_RARE_BOTH, FREQ_COMMON_WAP)

DEFAULT_COMMON_MIN = 0.05
SET_COLUMN_PREFIX = "set_"

MIS_SEVERE_MASKS = (MASK_MIS_CONSENSUS, MASK_MIS_SEVERE)
PROTEIN_ALTERING_MASKS = (
    MASK_PLOF_HIGH,
    MASK_MIS_CONSENSUS,
    MASK_MIS_SEVERE,
    MASK_PROTEIN_ALTERING,
)


@dataclass(frozen=True)
class VariantSetSpec:
    """One named analysis set.

    A set is the intersection of an optional coding-mask class, a cohort
    frequency class, and — when a site-annotation or gene-expression table is
    joined — site classes from ``genetics.genomics.variant.site_anno``. An empty
    ``masks`` means the coding mask is not constrained, which is what the
    non-coding and neutral sets need.
    """

    name: str
    masks: tuple[str, ...]
    frequency: str
    site_flags: tuple[str, ...] = ()
    site_not_flags: tuple[str, ...] = ()
    degeneracy: str | None = None
    require_expressed: bool = False

    @property
    def needs_site_anno(self) -> bool:
        return bool(self.site_flags or self.site_not_flags or self.degeneracy)


DEFAULT_SET_SPECS: tuple[VariantSetSpec, ...] = (
    VariantSetSpec("plof_rare_wgs", (MASK_PLOF_HIGH,), FREQ_RARE_WGS),
    VariantSetSpec("plof_rare_wap", (MASK_PLOF_HIGH,), FREQ_RARE_WAP),
    VariantSetSpec("plof_common_wap", (MASK_PLOF_HIGH,), FREQ_COMMON_WAP),
    VariantSetSpec("mis_consensus_rare_wgs", (MASK_MIS_CONSENSUS,), FREQ_RARE_WGS),
    VariantSetSpec("mis_severe_rare_wgs", MIS_SEVERE_MASKS, FREQ_RARE_WGS),
    VariantSetSpec("protein_altering_rare_both", PROTEIN_ALTERING_MASKS, FREQ_RARE_BOTH),
)

# Activated only when a site_anno (and, for the last one, gene_expr) table is joined.
SITE_SET_SPECS: tuple[VariantSetSpec, ...] = (
    # 4-fold degenerate and unconstrained: the neutral reference for MAF spectra
    # and Tajima's D, which otherwise have no within-genome baseline.
    VariantSetSpec(
        "neutral_4fold",
        (),
        FREQ_ANY,
        site_not_flags=(COL_CONSERVED,),
        degeneracy="4fold",
    ),
    VariantSetSpec(
        "conserved_noncoding_rare_wap",
        (),
        FREQ_RARE_WAP,
        site_flags=(COL_CONSERVED,),
        site_not_flags=(COL_CODING,),
    ),
    VariantSetSpec("regulatory_rare_wap", (), FREQ_RARE_WAP, site_flags=(COL_REGULATORY,)),
    VariantSetSpec(
        "plof_conserved_rare_wap",
        (MASK_PLOF_HIGH,),
        FREQ_RARE_WAP,
        site_flags=(COL_CONSERVED,),
    ),
    VariantSetSpec("plof_expressed_rare_wap", (MASK_PLOF_HIGH,), FREQ_RARE_WAP, require_expressed=True),
)

ALL_SET_SPECS: tuple[VariantSetSpec, ...] = DEFAULT_SET_SPECS + SITE_SET_SPECS

GENE_ID_CANDIDATES = ("Gene", "gene_id")
COL_EXPRESSED = "expressed"


def normalize_variant_id(df: pd.DataFrame, *, source: str) -> pd.DataFrame:
    """Return a copy whose variant identifier column is named ``ID``."""

    if df is None or df.empty:
        raise ValueError(f"Empty {source} table")
    for candidate in VARIANT_ID_CANDIDATES:
        if candidate in df.columns:
            out = df.copy()
            if candidate != "ID":
                out = out.rename(columns={candidate: "ID"})
            out["ID"] = out["ID"].astype(str)
            return out.drop_duplicates(subset=["ID"], keep="first")
    raise ValueError(
        f"{source} table has no variant ID column (looked for {', '.join(VARIANT_ID_CANDIDATES)})"
    )


def select_set_specs(names: str | None) -> tuple[VariantSetSpec, ...]:
    """Pick specs by comma-separated name; ``None`` keeps every known set."""

    if names is None or not str(names).strip():
        return ALL_SET_SPECS
    wanted = [part.strip() for part in str(names).split(",") if part.strip()]
    by_name = {spec.name: spec for spec in ALL_SET_SPECS}
    unknown = [name for name in wanted if name not in by_name]
    if unknown:
        raise ValueError(f"Unknown variant set(s): {', '.join(unknown)}")
    return tuple(by_name[name] for name in wanted)


def applicable_set_specs(
    specs: tuple[VariantSetSpec, ...],
    merged: pd.DataFrame,
) -> tuple[VariantSetSpec, ...]:
    """Drop specs whose site or expression inputs were not joined."""

    has_site = COL_DEGENERACY in merged.columns
    has_expression = COL_EXPRESSED in merged.columns
    kept = []
    for spec in specs:
        if spec.needs_site_anno and not has_site:
            print(f"variant_sets: skipping {spec.name} (no site_anno table joined)")
            continue
        if spec.require_expressed and not has_expression:
            print(f"variant_sets: skipping {spec.name} (no gene_expr table joined)")
            continue
        kept.append(spec)
    if not kept:
        raise ValueError("variant_sets: no set spec is applicable to the joined inputs")
    return tuple(kept)


def frequency_mask(
    merged: pd.DataFrame,
    frequency: str,
    *,
    common_min: float = DEFAULT_COMMON_MIN,
) -> pd.Series:
    """Boolean row mask for one frequency class."""

    if frequency == FREQ_ANY:
        return pd.Series(True, index=merged.index)
    if frequency == FREQ_COMMON_WAP:
        return pd.to_numeric(merged["MAF_wap"], errors="coerce").ge(common_min).fillna(False)
    if frequency in {FREQ_RARE_WGS, FREQ_RARE_WAP, FREQ_RARE_BOTH}:
        if frequency not in merged.columns:
            raise ValueError(f"dual_af table missing {frequency} column")
        return merged[frequency].astype("boolean").fillna(False).astype(bool)
    raise ValueError(f"Unknown frequency class: {frequency}")


def join_mask_and_af(coding_mask_df: pd.DataFrame, dual_af_df: pd.DataFrame) -> pd.DataFrame:
    """Inner-join coding masks to cohort frequencies on variant ID."""

    masks = normalize_variant_id(coding_mask_df, source="coding_mask")
    if "coding_mask" not in masks.columns:
        raise ValueError("coding_mask table missing coding_mask column")
    freqs = normalize_variant_id(dual_af_df, source="dual_af")
    keep = ["ID", "coding_mask"]
    for optional in ("consequence_severe", "Gene", "SYMBOL", "Feature", "subgenome"):
        if optional in masks.columns:
            keep.append(optional)
    return masks[keep].merge(freqs, on="ID", how="inner")


def join_site_anno(merged: pd.DataFrame, site_anno_df: pd.DataFrame) -> pd.DataFrame:
    """Left-join site classes on variant ID; unannotated variants keep NA."""

    sites = normalize_variant_id(site_anno_df, source="site_anno")
    available = [column for column in SITE_CLASS_COLUMNS if column in sites.columns]
    if not available:
        raise ValueError(
            "site_anno table has no class columns; run genetics.genomics.variant.site_anno first"
        )
    return merged.merge(sites[["ID", *available]], on="ID", how="left")


def join_gene_expression(merged: pd.DataFrame, gene_expr_df: pd.DataFrame) -> pd.DataFrame:
    """Left-join per-gene expression onto the VEP gene column."""

    gene_column = next((name for name in GENE_ID_CANDIDATES if name in merged.columns), None)
    if gene_column is None:
        raise ValueError(
            "coding_mask table has no gene column, so gene expression cannot be joined"
        )
    if "gene_id" not in gene_expr_df.columns:
        raise ValueError("gene_expr table missing gene_id column")
    keep = ["gene_id", COL_EXPRESSED]
    for optional in ("max_tpm", "tau", "breadth"):
        if optional in gene_expr_df.columns:
            keep.append(optional)
    expression = gene_expr_df[keep].drop_duplicates(subset=["gene_id"], keep="first")
    out = merged.merge(expression, left_on=gene_column, right_on="gene_id", how="left")
    out[COL_EXPRESSED] = _as_bool(out[COL_EXPRESSED])
    return out.drop(columns=["gene_id"])


def _as_bool(column: pd.Series) -> pd.Series:
    """Coerce a possibly-NA flag column to plain bool, treating NA as False."""

    return column.astype("boolean").fillna(False).astype(bool)


def site_mask(merged: pd.DataFrame, spec: VariantSetSpec) -> pd.Series:
    """Boolean row mask for a spec's site-annotation and expression requirements."""

    keep = pd.Series(True, index=merged.index)
    for column in spec.site_flags:
        keep &= _as_bool(merged[column])
    for column in spec.site_not_flags:
        keep &= ~_as_bool(merged[column])
    if spec.degeneracy is not None:
        keep &= _as_bool(merged[COL_DEGENERACY].eq(spec.degeneracy))
    if spec.require_expressed:
        keep &= _as_bool(merged[COL_EXPRESSED])
    return keep


def build_variant_sets(
    merged: pd.DataFrame,
    specs: tuple[VariantSetSpec, ...] = DEFAULT_SET_SPECS,
    *,
    common_min: float = DEFAULT_COMMON_MIN,
) -> pd.DataFrame:
    """Add one boolean ``set_<name>`` column per spec."""

    out = merged.copy()
    for spec in specs:
        in_set = frequency_mask(out, spec.frequency, common_min=common_min)
        if spec.masks:
            in_set &= out["coding_mask"].isin(spec.masks)
        in_set &= site_mask(out, spec)
        out[f"{SET_COLUMN_PREFIX}{spec.name}"] = in_set
    return out


def summarize_variant_sets(
    sets_df: pd.DataFrame,
    specs: tuple[VariantSetSpec, ...] = DEFAULT_SET_SPECS,
    *,
    common_min: float = DEFAULT_COMMON_MIN,
) -> dict[str, int | float]:
    stats: dict[str, int | float] = {
        "n_annotated_variants": int(len(sets_df)),
        "common_min": common_min,
    }
    for spec in specs:
        stats[f"n_{spec.name}"] = int(sets_df[f"{SET_COLUMN_PREFIX}{spec.name}"].sum())
    return stats


def write_variant_set_extracts(
    sets_df: pd.DataFrame,
    specs: tuple[VariantSetSpec, ...],
    output_prefix: str | Path,
) -> dict[str, int]:
    """Write one PLINK2 ``--extract`` ID list per set."""

    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    counts: dict[str, int] = {}
    for spec in specs:
        ids = sets_df.loc[sets_df[f"{SET_COLUMN_PREFIX}{spec.name}"], "ID"].astype(str).tolist()
        path = Path(f"{prefix}.{spec.name}.extract")
        with path.open("w", encoding="utf-8") as handle:
            for variant_id in ids:
                handle.write(f"{variant_id}\n")
        counts[spec.name] = len(ids)
    return counts


def plot_variant_sets(
    sets_df: pd.DataFrame,
    specs: tuple[VariantSetSpec, ...],
    output_prefix: str | Path,
) -> None:
    names = [spec.name for spec in specs]
    values = [float(sets_df[f"{SET_COLUMN_PREFIX}{spec.name}"].sum()) for spec in specs]
    plot_bar_chart(
        names=names,
        values=values,
        title="Analysis variant sets (coding mask x cohort frequency)",
        ylabel="Variants",
        filename=f"{output_prefix}.variant_sets.counts.png",
        ylim=(0.0, (max(values) * 1.15) if values and max(values) > 0 else 1.05),
        rotate_xlabels=30,
    )


def run_variant_sets(
    coding_mask_table: str | Path,
    dual_af_table: str | Path,
    output_prefix: str | Path,
    *,
    set_names: str | None = None,
    common_min: float = DEFAULT_COMMON_MIN,
    site_anno_table: str | Path | None = None,
    gene_expr_table: str | Path | None = None,
) -> pd.DataFrame:
    """Build analysis sets and write table, extract lists, thresholds, and plot."""

    merged = join_mask_and_af(
        load_df_generic(str(coding_mask_table)),
        load_df_generic(str(dual_af_table)),
    )
    if site_anno_table:
        merged = join_site_anno(merged, load_df_generic(str(site_anno_table)))
    if gene_expr_table:
        merged = join_gene_expression(merged, load_df_generic(str(gene_expr_table)))
    specs = applicable_set_specs(select_set_specs(set_names), merged)
    sets_df = build_variant_sets(merged, specs, common_min=common_min)
    save_df_to_tsv(sets_df, f"{output_prefix}.variant_sets.tsv")
    stats = summarize_variant_sets(sets_df, specs, common_min=common_min)
    save_thresholds(stats, f"{output_prefix}.variant_sets.th.tsv")
    write_variant_set_extracts(sets_df, specs, output_prefix)
    plot_variant_sets(sets_df, specs, output_prefix)
    print(
        f"variant_sets: annotated={stats['n_annotated_variants']} "
        + " ".join(f"{spec.name}={stats[f'n_{spec.name}']}" for spec in specs)
    )
    return sets_df
