"""Classify MAC vs H into selfing (bottom) vs all-het (diagonal) bands.

Diploid counting identity: ``MAC = H + 2A``.  Then ``R_H = H/MAC`` is 0 when
every minor allele is in homozygotes and 1 when ``A=0`` (all minor alleles in
heterozygotes).  The two-band MAC–H heatmap is treated as a mixture of those
mechanisms rather than extra-binomial overdispersion.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from genetics.genomics.variant.mac import _compute_mac_table
from genetics.genomics.variant.mac_het_excess import load_median_inbreeding_from_het
from genetics.genomics.variant.mq import load_variant_mq_info
from genetics.genomics.variant.popdep import load_df_from_tsv_popdep_info
from genetics.genomics.variant.variant_utils import (
    load_df_from_plink_gcount,
    load_df_from_plink_variant,
)
from infra.utils.errors import fail
from infra.utils.graph import (
    SEQUENTIAL_WHITE_RED,
    SUBGENOME_PURPLE,
    TWO_CATEGORY_BLUE,
    TWO_CATEGORY_RED,
    plot_binned_mean_heatmap,
    plot_distribution_with_stats,
    plot_grouped_density_lines,
    plot_multi_line_series,
)
from infra.utils.io import save_df_to_tsv, save_thresholds

logger = logging.getLogger(__name__)

RH_BOTTOM = "bottom"
RH_MIDDLE = "middle"
RH_DIAGONAL = "diagonal"
RH_MAC_LT_MIN = "mac_lt_min"
RH_CLASS_ORDER = (RH_BOTTOM, RH_MIDDLE, RH_DIAGONAL)
RH_CLASS_COLORS = {
    RH_BOTTOM: TWO_CATEGORY_BLUE,
    RH_MIDDLE: SUBGENOME_PURPLE,
    RH_DIAGONAL: TWO_CATEGORY_RED,
}
RH_MAC_MIN = 5
RH_HEATMAP_MAC_MAX = 100
_FLOAT_TSV = "%.8g"
_POPDEP_KEEP = (
    "ID",
    "Depth_Mean",
    "Depth_SD",
    "Depth_CV",
    "RelativeDepth_Mean",
    "RelativeDepth_SD",
    "RelativeDepth_CV",
)


def annotate_rh_class(mac_df: pd.DataFrame, mac_min: int = RH_MAC_MIN) -> pd.DataFrame:
    """Add ``R_H = H/MAC``, ``rh_class``, and ``all_het`` (MAC>=mac_min and H=MAC)."""
    out = mac_df.copy()
    if "H" not in out.columns:
        if "HET_REF_ALT_CTS" not in out.columns:
            fail("annotate_rh_class needs H or HET_REF_ALT_CTS")
        out["H"] = pd.to_numeric(out["HET_REF_ALT_CTS"], errors="coerce").fillna(0)
    else:
        out["H"] = pd.to_numeric(out["H"], errors="coerce").fillna(0)
    out["MAC"] = pd.to_numeric(out["MAC"], errors="coerce").fillna(0)
    mac = out["MAC"].to_numpy(dtype=float)
    h = out["H"].to_numpy(dtype=float)
    rh = np.full(mac.shape, np.nan)
    ok = mac > 0
    rh[ok] = h[ok] / mac[ok]
    out["R_H"] = rh
    mac_i = np.rint(mac)
    h_i = np.rint(h)
    classes = np.full(mac.shape, RH_MAC_LT_MIN, dtype=object)
    band = mac >= float(mac_min)
    classes[band & (h_i == 0)] = RH_BOTTOM
    classes[band & (h_i == mac_i) & (mac_i > 0)] = RH_DIAGONAL
    classes[band & (h_i > 0) & (h_i < mac_i)] = RH_MIDDLE
    out["rh_class"] = classes
    out["all_het"] = band & (h_i == mac_i) & (mac_i > 0)
    return out


def inbreeding_from_het_fraction(
    sites: pd.DataFrame,
    mac_min: int = RH_MAC_MIN,
    exclude_diagonal: bool = False,
    rh_classes: Sequence[str] | None = None,
) -> dict[str, float]:
    """Moment estimator of F from ``H ~ Binomial(MAC, (1-q)(1-F))``.

    ``F_mom = 1 - sum(H) / sum(MAC * (1-q))`` with ``q = MAC/(2N)``.
    """
    work = sites[pd.to_numeric(sites["MAC"], errors="coerce") >= float(mac_min)].copy()
    if rh_classes is not None:
        work = work[work["rh_class"].isin(list(rh_classes))]
    elif exclude_diagonal:
        work = work[work["rh_class"] != RH_DIAGONAL]
    empty = {
        "n_sites": 0,
        "sum_H": 0.0,
        "sum_MAC": 0.0,
        "mean_R_H": np.nan,
        "F_mom": np.nan,
    }
    if work.empty:
        return empty
    mac = pd.to_numeric(work["MAC"], errors="coerce").to_numpy(dtype=float)
    h = pd.to_numeric(work["H"], errors="coerce").to_numpy(dtype=float)
    n = pd.to_numeric(work["Total_Samples"], errors="coerce").to_numpy(dtype=float)
    q = np.divide(mac, np.clip(2.0 * n, 1.0, None), out=np.zeros_like(mac), where=np.isfinite(n))
    weight = mac * np.clip(1.0 - q, 0.0, 1.0)
    wsum = float(np.nansum(weight))
    hsum = float(np.nansum(h))
    p_het = hsum / wsum if wsum > 0 else np.nan
    mean_rh = float(np.nanmean(np.divide(h, mac, out=np.full_like(h, np.nan), where=mac > 0)))
    return {
        "n_sites": int(len(work)),
        "sum_H": hsum,
        "sum_MAC": float(np.nansum(mac)),
        "mean_R_H": mean_rh,
        "F_mom": float(1.0 - p_het) if np.isfinite(p_het) else np.nan,
    }


def allhet_fraction_by_mac(
    sites: pd.DataFrame,
    mac_max: int = RH_HEATMAP_MAC_MAX,
    f_panel: float | None = None,
    n_diploid: float | None = None,
) -> pd.DataFrame:
    """Observed P(H=MAC) per integer MAC, plus binomial all-het probability if F given."""
    tested = sites[pd.to_numeric(sites["MAC"], errors="coerce") > 0].copy()
    tested["MAC"] = pd.to_numeric(tested["MAC"], errors="coerce").round().astype(int)
    rows = []
    for mac in range(1, int(mac_max) + 1):
        sub = tested[tested["MAC"] == mac]
        n_sites = int(len(sub))
        n_all = int((pd.to_numeric(sub["R_H"], errors="coerce") >= 1.0 - 1e-12).sum()) if n_sites else 0
        expected = np.nan
        if f_panel is not None and n_diploid is not None and float(n_diploid) > 0:
            q = mac / (2.0 * float(n_diploid))
            p = max(0.0, min(1.0, (1.0 - q) * (1.0 - float(f_panel))))
            expected = float(p**mac)
        elif f_panel is not None:
            p = max(0.0, min(1.0, 1.0 - float(f_panel)))
            expected = float(p**mac)
        rows.append(
            {
                "MAC": mac,
                "n_sites": n_sites,
                "n_all_het": n_all,
                "frac_all_het": (n_all / n_sites) if n_sites else np.nan,
                "expected_frac_all_het": expected,
            }
        )
    return pd.DataFrame(rows)


def summarize_rh_classes(sites: pd.DataFrame, mac_min: int = RH_MAC_MIN) -> pd.DataFrame:
    work = sites.copy()
    work["MAC"] = pd.to_numeric(work["MAC"], errors="coerce")
    rows = []
    for cls in (RH_MAC_LT_MIN,) + RH_CLASS_ORDER:
        sub = work[work["rh_class"] == cls]
        rows.append(
            {
                "rh_class": cls,
                "n_sites": int(len(sub)),
                "mac_min": int(mac_min) if cls != RH_MAC_LT_MIN else 0,
                "mean_R_H": float(pd.to_numeric(sub["R_H"], errors="coerce").mean()) if len(sub) else np.nan,
                "median_MAC": float(sub["MAC"].median()) if len(sub) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def _het_ohet_series(work: pd.DataFrame) -> pd.Series:
    """Heterozygous genotype count per sample from a PLINK2 ``--het`` table."""
    cols = {str(c).replace("#", ""): c for c in work.columns}
    if "O_HET" in work.columns:
        return pd.to_numeric(work["O_HET"], errors="coerce")
    hom_col = cols.get("O_HOM", cols.get("O(HOM)"))
    obs_col = cols.get("OBS_CT")
    if hom_col is not None and obs_col is not None:
        return pd.to_numeric(work[obs_col], errors="coerce") - pd.to_numeric(
            work[hom_col], errors="coerce"
        )
    fail(f"PLINK2 --het table missing O_HET / O(HOM)+OBS_CT; got {list(work.columns)}")


def sample_het_concentration(het_df: pd.DataFrame) -> dict[str, float]:
    """How concentrated all-het heterozygous calls are across samples."""
    work = het_df.copy()
    work.columns = [str(c).replace("#", "") for c in work.columns]
    o = _het_ohet_series(work).fillna(0.0)
    n = int(len(o))
    total = float(o.sum())
    if n == 0 or total <= 0:
        return {
            "n_samples": n,
            "sum_O_HET": total,
            "mean_O_HET": np.nan,
            "median_O_HET": np.nan,
            "cv_O_HET": np.nan,
            "max_frac": np.nan,
            "n_for_half_hets": 0,
        }
    ordered = o.sort_values(ascending=False).to_numpy(dtype=float)
    cum = np.cumsum(ordered) / total
    n50 = int(np.searchsorted(cum, 0.5, side="left") + 1)
    mean = float(o.mean())
    std = float(o.std(ddof=0))
    return {
        "n_samples": n,
        "sum_O_HET": total,
        "mean_O_HET": mean,
        "median_O_HET": float(o.median()),
        "cv_O_HET": (std / mean) if mean > 0 else np.nan,
        "max_frac": float(ordered[0] / total),
        "n_for_half_hets": n50,
    }


def _rename_popdep(df: pd.DataFrame, suffix: str) -> pd.DataFrame:
    keep = [c for c in _POPDEP_KEEP if c in df.columns]
    out = df[keep].copy()
    rename = {c: f"{c}{suffix}" for c in keep if c != "ID"}
    return out.rename(columns=rename)


def _join_covariates(
    sites: pd.DataFrame,
    popdep: pd.DataFrame | None = None,
    popdep_mapq30: pd.DataFrame | None = None,
    mq: pd.DataFrame | None = None,
    vmiss: pd.DataFrame | None = None,
) -> pd.DataFrame:
    merged = sites.copy()
    if "ID" not in merged.columns:
        fail("gcount-derived sites need an ID column to join covariates")
    if popdep is not None:
        merged = merged.merge(_rename_popdep(popdep, ""), on="ID", how="left")
    if popdep_mapq30 is not None:
        merged = merged.merge(_rename_popdep(popdep_mapq30, "_mapq30"), on="ID", how="left")
    if mq is not None:
        mq_work = mq.copy()
        mq_work.columns = [str(c).replace("#", "") for c in mq_work.columns]
        if "ID" in mq_work.columns and "MQ" in mq_work.columns:
            merged = merged.merge(mq_work[["ID", "MQ"]], on="ID", how="left")
    if vmiss is not None:
        miss = vmiss.copy()
        miss.columns = [str(c).replace("#", "") for c in miss.columns]
        cols = [c for c in ("ID", "F_MISS", "MISSING_CT") if c in miss.columns]
        if "ID" in cols:
            merged = merged.merge(miss[cols], on="ID", how="left")
    for col in (
        "Depth_Mean",
        "Depth_Mean_mapq30",
        "RelativeDepth_Mean",
        "RelativeDepth_Mean_mapq30",
        "MQ",
        "F_MISS",
    ):
        if col in merged.columns:
            merged[col] = pd.to_numeric(merged[col], errors="coerce")
    return merged


def covariate_enrichment_by_rh_class(merged: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "Depth_Mean",
        "Depth_Mean_mapq30",
        "RelativeDepth_Mean",
        "RelativeDepth_Mean_mapq30",
        "MQ",
        "F_MISS",
        "R_H",
        "MAC",
    ]
    rows = []
    for cls in RH_CLASS_ORDER:
        sub = merged[merged["rh_class"] == cls]
        row: dict[str, float | str | int] = {"rh_class": cls, "n_sites": int(len(sub))}
        for col in metrics:
            if col not in merged.columns:
                continue
            vals = pd.to_numeric(sub[col], errors="coerce")
            row[f"median_{col}"] = float(vals.median()) if vals.notna().any() else np.nan
            row[f"mean_{col}"] = float(vals.mean()) if vals.notna().any() else np.nan
        rows.append(row)
    out = pd.DataFrame(rows)
    if RH_BOTTOM in set(out["rh_class"]) and RH_DIAGONAL in set(out["rh_class"]):
        bottom = out.set_index("rh_class").loc[RH_BOTTOM]
        for col in list(out.columns):
            if not col.startswith("median_") or col == "median_R_H":
                continue
            ratio_col = col.replace("median_", "ratio_diag_over_bottom_")
            b = bottom.get(col, np.nan)
            d = out.set_index("rh_class").loc[RH_DIAGONAL].get(col, np.nan)
            out[ratio_col] = np.nan
            if np.isfinite(b) and float(b) != 0:
                out.loc[out["rh_class"] == RH_DIAGONAL, ratio_col] = float(d) / float(b)
    return out


def _integer_edges(lo: int, hi: int) -> np.ndarray:
    return np.arange(float(lo) - 0.5, float(hi) + 1.5, 1.0)


def _write_extract(sites: pd.DataFrame, path: str) -> int:
    ids = sites.loc[sites["all_het"].astype(bool), "ID"].astype(str)
    Path(path).write_text("\n".join(ids.tolist()) + ("\n" if len(ids) else ""))
    return int(len(ids))


def _plot_rh_diagnostics(
    sites: pd.DataFrame,
    by_mac: pd.DataFrame,
    merged: pd.DataFrame,
    output_prefix: str,
    f_het: float | None,
    f_ex: float | None,
) -> None:
    band = sites[sites["rh_class"].isin(RH_CLASS_ORDER)].copy()
    if not band.empty and band["R_H"].notna().any():
        rh_vals = pd.to_numeric(band["R_H"], errors="coerce").dropna()
        plot_distribution_with_stats(
            band.dropna(subset=["R_H"]),
            col="R_H",
            title="H/MAC among MAC>=5 sites (bottom ≈0, diagonal =1)",
            filename=f"{output_prefix}.rh_hist.png",
            mean_val=float(rh_vals.mean()),
            median_val=float(rh_vals.median()),
            x_label="R_H = H / MAC",
            y_label="Site count",
            bins=50,
            xlim=(0.0, 1.0),
        )
    y_specs = [
        {"y_col": "frac_all_het", "label": "Observed P(H=MAC)", "color": TWO_CATEGORY_RED, "linewidth": 1.8},
    ]
    if "expected_frac_all_het" in by_mac.columns and by_mac["expected_frac_all_het"].notna().any():
        y_specs.append(
            {
                "y_col": "expected_frac_all_het",
                "label": f"Selfing P(H=MAC) at F={f_het:.4g}" if f_het is not None else "Selfing P(H=MAC)",
                "color": TWO_CATEGORY_BLUE,
                "linestyle": "--",
                "linewidth": 1.5,
            }
        )
    plot_multi_line_series(
        by_mac,
        x_col="MAC",
        y_specs=y_specs,
        title="All-heterozygote site fraction vs MAC",
        filename=f"{output_prefix}.frac_allhet_by_mac.png",
        x_label="Minor allele count (MAC)",
        y_label="Fraction of sites with H = MAC",
        ylim=(0.0, 1.05),
    )
    class_df = merged[merged["rh_class"].isin(RH_CLASS_ORDER)].copy()
    if class_df.empty:
        return
    for col, xlabel, stem in (
        ("Depth_Mean", "Mean depth (unfiltered)", "depth"),
        ("Depth_Mean_mapq30", "Mean depth (MAPQ>=30)", "depth_mapq30"),
        ("MQ", "Site mapping quality", "mq"),
        ("F_MISS", "Missing rate", "miss"),
    ):
        if col not in class_df.columns:
            continue
        plot_df = class_df.dropna(subset=[col])
        if plot_df.empty or plot_df[col].nunique() < 2:
            continue
        present = [c for c in RH_CLASS_ORDER if c in set(plot_df["rh_class"])]
        if len(present) < 2:
            continue
        plot_grouped_density_lines(
            plot_df,
            x_col=col,
            hue_col="rh_class",
            hue_order=present,
            palette_map={k: RH_CLASS_COLORS[k] for k in present},
            title=f"{xlabel} by R_H class",
            filename=f"{output_prefix}.{stem}_by_rh_class.png",
            xlabel=xlabel,
            ylabel="Density",
        )
    window = class_df[
        (pd.to_numeric(class_df["MAC"], errors="coerce") >= 1)
        & (pd.to_numeric(class_df["MAC"], errors="coerce") <= RH_HEATMAP_MAC_MAX)
        & (pd.to_numeric(class_df["H"], errors="coerce") <= RH_HEATMAP_MAC_MAX)
    ].copy()
    edges = _integer_edges(0, RH_HEATMAP_MAC_MAX)
    f_note = f"F_het={f_het:.4g}" if f_het is not None else ""
    f_ex_note = f"F_excl_diag={f_ex:.4g}" if f_ex is not None else ""
    subtitle = ", ".join(x for x in (f_note, f_ex_note) if x)
    for col, cbar, stem in (
        ("Depth_Mean", "Mean depth", "depth"),
        ("Depth_Mean_mapq30", "MAPQ>=30 mean depth", "depth_mapq30"),
        ("F_MISS", "Missing rate", "miss"),
        ("MQ", "Mapping quality", "mq"),
    ):
        if col not in window.columns:
            continue
        heat = window.dropna(subset=["MAC", "H", col])
        if len(heat) < 10:
            continue
        plot_binned_mean_heatmap(
            heat,
            x_col="MAC",
            y_col="H",
            value_col=col,
            x_label="Minor allele count (MAC)",
            y_label="Heterozygous sample count (H)",
            cbar_label=cbar,
            title=f"Mean {cbar.lower()} in MAC–H bins" + (f" ({subtitle})" if subtitle else ""),
            filename=f"{output_prefix}.heatmap_mac_h.{stem}.png",
            x_edges=edges,
            y_edges=edges,
            cmap=SEQUENTIAL_WHITE_RED,
            x_label_fmt="{:.0f}",
            y_label_fmt="{:.0f}",
            x_tick_step=10,
            y_tick_step=10,
        )


def ana_mac_het_rh_diagnostics(
    gcount_file: str,
    output_prefix: str = "variant.mac_het.rh",
    het_path: str | None = None,
    f_panel: float | None = None,
    popdep_file: str | None = None,
    popdep_mapq30_file: str | None = None,
    mq_file: str | None = None,
    vmiss_file: str | None = None,
    mac_min: int = RH_MAC_MIN,
    mac_max: int = RH_HEATMAP_MAC_MAX,
) -> None:
    """Classify rare-site MAC–H bands and write F / enrichment diagnostics.

    Does not emit a filtered pfile. Optional MAPQ>=30 popdep is joined as a
    technical covariate (same site IDs; genotypes stay on the unfiltered callset).
    """
    df_raw = load_df_from_plink_gcount(gcount_file)
    sites = annotate_rh_class(_compute_mac_table(df_raw), mac_min=mac_min)
    f_het = None
    if het_path:
        f_het = load_median_inbreeding_from_het(het_path)
    elif f_panel is not None:
        f_het = float(f_panel)

    n_med = float(np.nanmedian(pd.to_numeric(sites["Total_Samples"], errors="coerce")))
    n_diploid = n_med if np.isfinite(n_med) and n_med > 0 else None
    f_all = inbreeding_from_het_fraction(sites, mac_min=mac_min, exclude_diagonal=False)
    f_ex = inbreeding_from_het_fraction(sites, mac_min=mac_min, exclude_diagonal=True)
    f_bottom = inbreeding_from_het_fraction(sites, mac_min=mac_min, rh_classes=(RH_BOTTOM,))
    classes = summarize_rh_classes(sites, mac_min=mac_min)
    by_mac = allhet_fraction_by_mac(sites, mac_max=mac_max, f_panel=f_het, n_diploid=n_diploid)

    popdep = load_df_from_tsv_popdep_info(popdep_file) if popdep_file else None
    popdep30 = load_df_from_tsv_popdep_info(popdep_mapq30_file) if popdep_mapq30_file else None
    mq = load_variant_mq_info(mq_file) if mq_file else None
    vmiss = load_df_from_plink_variant(vmiss_file) if vmiss_file else None
    merged = _join_covariates(sites, popdep=popdep, popdep_mapq30=popdep30, mq=mq, vmiss=vmiss)
    enrich = covariate_enrichment_by_rh_class(merged)

    n_extract = _write_extract(sites, f"{output_prefix}.all_het.extract")
    save_df_to_tsv(classes, f"{output_prefix}.rh_class.info.tsv", float_format=_FLOAT_TSV)
    save_df_to_tsv(by_mac, f"{output_prefix}.frac_allhet_by_mac.info.tsv", float_format=_FLOAT_TSV)
    save_df_to_tsv(enrich, f"{output_prefix}.rh_enrich.info.tsv", float_format=_FLOAT_TSV)
    n_band = int((sites["rh_class"].isin(RH_CLASS_ORDER)).sum())
    n_diag = int((sites["rh_class"] == RH_DIAGONAL).sum())
    stats = {
        "F_het_median": f_het if f_het is not None else np.nan,
        "F_mom_all": f_all["F_mom"],
        "F_mom_exclude_diagonal": f_ex["F_mom"],
        "F_mom_bottom": f_bottom["F_mom"],
        "N_median": n_diploid if n_diploid is not None else np.nan,
        "mac_min": int(mac_min),
        "n_mac_ge_min": n_band,
        "n_bottom": int((sites["rh_class"] == RH_BOTTOM).sum()),
        "n_middle": int((sites["rh_class"] == RH_MIDDLE).sum()),
        "n_diagonal": n_diag,
        "frac_diagonal_mac_ge_min": (n_diag / n_band) if n_band else 0.0,
        "n_all_het_extract": n_extract,
        "mean_R_H_all": f_all["mean_R_H"],
        "mean_R_H_exclude_diagonal": f_ex["mean_R_H"],
    }
    save_thresholds(stats, f"{output_prefix}.th.tsv")
    print(stats)
    _plot_rh_diagnostics(sites, by_mac, merged, output_prefix, f_het, f_ex["F_mom"])


def ana_allhet_sample_burden(het_file: str, output_prefix: str = "variant.mac_het.rh_sample") -> None:
    """Per-sample heterozygous burden on ``all_het`` sites from PLINK2 ``--het``."""
    het = pd.read_csv(het_file, sep=r"\s+")
    stats = sample_het_concentration(het)
    save_thresholds(stats, f"{output_prefix}.th.tsv")
    print(stats)
    work = het.copy()
    work.columns = [str(c).replace("#", "") for c in work.columns]
    work["O_HET"] = _het_ohet_series(work)
    o = pd.to_numeric(work["O_HET"], errors="coerce").dropna()
    plot_distribution_with_stats(
        work,
        col="O_HET",
        title="Per-sample heterozygous count on all-het (H=MAC, MAC>=5) sites",
        filename=f"{output_prefix}.ohet_hist.png",
        mean_val=float(o.mean()) if len(o) else None,
        median_val=float(o.median()) if len(o) else None,
        x_label="O(HET) among all-het sites",
        y_label="Sample count",
        bins=40,
    )


def combine_mac_het_rh_thresholds(
    th_files: Sequence[str],
    output_prefix: str = "all.variant.mac_het.rh",
) -> None:
    """Stack per-subgenome R_H F / class-count threshold tables."""
    if not th_files:
        fail("No R_H threshold TSV files to combine")
    frames = []
    for path in th_files:
        df = pd.read_csv(path, sep="\t")
        stem = Path(path).name
        sg = stem.split(".", 1)[0]
        df.insert(0, "subgenome", sg)
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    save_df_to_tsv(out, f"{output_prefix}.th.info.tsv", float_format=_FLOAT_TSV)
    print(out.to_string(index=False))
