"""Variant filter cascade: optional PopDep density cells → het → MAC≥5 → WAP MAC≥2."""

from __future__ import annotations

import gzip
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from genetics.genomics.variant.mac import _compute_mac_table
from genetics.genomics.variant.mac_het_excess import (
    RETAIN_ALPHA,
    load_median_inbreeding_from_het,
    site_excess_het_table,
)
from genetics.genomics.variant.maf_miss_group import _maf_from_afreq_df
from genetics.genomics.variant.variant_utils import (
    load_df_from_plink_gcount,
    load_df_from_plink_variant,
)
from infra.utils.graph import plot_bar_chart, plot_stacked_histogram, plot_stacked_proportion_barh
from infra.utils.io import save_df_to_tsv, save_thresholds

# Genome-scale reference totals (post-popdep 80% → het → MAC≥5 → WAP+MAC≥2).
GENOME_PIPELINE_TOTALS = {
    "after_popdep": 2.24e9,
    "het_pass": 1.71e9,
    "mac_ge5": 1.26e9,
    "wap_mac_ge2": 777e6,
}

STEP_ORDER = ("after_popdep", "het_pass", "mac_ge5", "wap_mac_ge2")
MAF_CLASS_SHORT = ("rare", "common")

def _out_path(output_prefix: str | Path, name: str) -> Path:
    return Path(f"{output_prefix}.{name}")


RARE_COMMON_LABELS = ("Rare (MAF < 0.05)", "Common (MAF ≥ 0.05)")
RARE_COMMON_COLORS = {
    RARE_COMMON_LABELS[0]: "#4c78a8",
    RARE_COMMON_LABELS[1]: "#e45756",
}


def assign_rare_common_class(maf: pd.Series, *, threshold: float = 0.05) -> pd.Series:
    m = pd.to_numeric(maf, errors="coerce")
    return pd.Series(
        np.where(m < threshold, RARE_COMMON_LABELS[0], RARE_COMMON_LABELS[1]),
        index=maf.index,
    )


def _write_id_list(ids: pd.Series, path: Path) -> None:
    path.write_text("\n".join(ids.astype(str)) + "\n", encoding="utf-8")


def _write_keep_file(sample_ids: list[str], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for iid in sample_ids:
            handle.write(f"0\t{iid}\n")


def _load_wap_keep_samples(wap_keep_path: str | Path) -> list[str]:
    keep = pd.read_csv(wap_keep_path, sep=r"\s+", header=None, names=["FID", "IID"])
    return keep["IID"].astype(str).tolist()


def _wap_mac_from_plink(
    pfile_prefix: str | Path,
    variant_ids: pd.Series,
    wap_keep_path: str | Path,
    *,
    plink2: str = "plink2",
    min_mac: int = 2,
) -> pd.DataFrame:
    """Return ID, MAC_WAP, MAF_WAP for variants passing MAC≥min_mac in WAP panel."""
    prefix = Path(pfile_prefix)
    if variant_ids.empty:
        return pd.DataFrame(columns=["ID", "MAC_WAP", "MAF_WAP"])

    wap_samples = _load_wap_keep_samples(wap_keep_path)
    with tempfile.TemporaryDirectory(prefix="filter_cascade_wap_") as tmp:
        tmp_dir = Path(tmp)
        extract_path = tmp_dir / "extract.ids"
        keep_path = tmp_dir / "wap.keep"
        out_prefix = tmp_dir / "wap_gcount"
        _write_id_list(variant_ids, extract_path)
        _write_keep_file(wap_samples, keep_path)

        cmd = [
            plink2,
            "--pfile",
            str(prefix),
            "--keep",
            str(keep_path),
            "--extract",
            str(extract_path),
            "--geno-counts",
            "--out",
            str(out_prefix),
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        gcount_path = Path(f"{out_prefix}.gcount")
        if not gcount_path.exists():
            raise FileNotFoundError(f"plink2 geno-counts missing output: {gcount_path}")

        mac_df = _compute_mac_table(load_df_from_plink_gcount(str(gcount_path)))
        mac_df = mac_df[pd.to_numeric(mac_df["MAC"], errors="coerce") >= int(min_mac)].copy()
        mac_df = mac_df.rename(columns={"MAC": "MAC_WAP", "MAF": "MAF_WAP"})
        return mac_df[["ID", "MAC_WAP", "MAF_WAP"]]


def _load_subgenome_frame(
    process_dir: Path,
    subgenome: str,
    *,
    f_panel: float,
    alpha: float,
) -> pd.DataFrame:
    variant_dir = process_dir / "variant"
    gcount = load_df_from_plink_gcount(str(variant_dir / f"{subgenome}.info.gcount"))
    sites = site_excess_het_table(gcount, f_panel)
    afreq_path = variant_dir / f"{subgenome}.info.afreq"
    if afreq_path.exists():
        afreq = _maf_from_afreq_df(load_df_from_plink_variant(str(afreq_path)))
        if "MAF" in sites.columns:
            sites = sites.drop(columns=["MAF"])
        sites = sites.merge(afreq[["ID", "MAF"]], on="ID", how="left")
    else:
        sites["MAF"] = pd.to_numeric(sites["MAF"], errors="coerce")
    sites["subgenome"] = subgenome
    sites["het_pass"] = (pd.to_numeric(sites["MAC"], errors="coerce") > 0) & (
        pd.to_numeric(sites["p_excess"], errors="coerce") >= float(alpha)
    )
    sites["mac_ge5"] = pd.to_numeric(sites["MAC"], errors="coerce") >= 5
    return sites


def filter_cascade_test_thin(
    process_dir: str | Path,
    *,
    het_paths: dict[str, str | Path],
    wap_keep_path: str | Path | None,
    subgenomes: tuple[str, ...] = ("A", "B", "D", "Others"),
    alpha: float = RETAIN_ALPHA,
    plink2: str = "plink2",
    popdep_keep_ids: set[str] | None = None,
    popdep_keep_ids_60: set[str] | None = None,
    skip_wap: bool = False,
    min_mac: int = 5,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Apply optional PopDep trusted cells → het → MAC≥5 → WAP MAC≥2.

    ``popdep_keep_ids`` is the 80% density set used in the cascade.
    ``popdep_keep_ids_60`` is counted only. Returns (final_sites, by_subgenome, totals).
    """
    process_dir = Path(process_dir)
    step_rows: list[dict] = []
    parts: list[pd.DataFrame] = []

    raw_total = 0
    popdep60_total = 0
    popdep80_total = 0
    het_total = 0
    mac5_total = 0
    use_popdep = popdep_keep_ids is not None

    for sg in subgenomes:
        het_path = het_paths.get(sg)
        if het_path is None:
            raise ValueError(f"Missing het path for subgenome {sg}")
        f_panel = load_median_inbreeding_from_het(str(het_path))
        sites = _load_subgenome_frame(process_dir, sg, f_panel=f_panel, alpha=alpha)
        n_raw = int(len(sites))
        raw_total += n_raw

        ids = sites["ID"].astype(str)
        n_popdep_60 = int(ids.isin(popdep_keep_ids_60).sum()) if popdep_keep_ids_60 is not None else n_raw
        popdep60_total += n_popdep_60

        if use_popdep:
            sites = sites[ids.isin(popdep_keep_ids)].copy()
        n_popdep_80 = int(len(sites))
        popdep80_total += n_popdep_80

        het_sites = sites[sites["het_pass"]].copy()
        n_het = int(len(het_sites))
        het_total += n_het

        mac5_sites = het_sites[het_sites["mac_ge5"]].copy()
        n_mac5 = int(len(mac5_sites))
        mac5_total += n_mac5

        if skip_wap or wap_keep_path is None:
            mac_num = pd.to_numeric(het_sites["MAC"], errors="coerce")
            final = het_sites[mac_num >= int(min_mac)].copy()
            final["MAC_WAP"] = pd.to_numeric(final["MAC"], errors="coerce")
            final["MAF_WAP"] = pd.to_numeric(final["MAF"], errors="coerce")
        else:
            pfile_prefix = process_dir / f"{sg}_test.plink2"
            wap_pass = _wap_mac_from_plink(
                pfile_prefix,
                mac5_sites["ID"],
                wap_keep_path,
                plink2=plink2,
                min_mac=2,
            )
            final = mac5_sites.merge(wap_pass, on="ID", how="inner")
            final["MAF"] = pd.to_numeric(final["MAF_WAP"], errors="coerce")
        parts.append(final)

        step_rows.append(
            {
                "subgenome": sg,
                "n_raw": n_raw,
                "n_popdep_60": n_popdep_60,
                "n_popdep_80": n_popdep_80,
                "n_het_pass": n_het,
                "n_mac_ge5": n_mac5,
                "n_wap_mac_ge2": int(len(final)),
            }
        )

    final_all = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    summary = pd.DataFrame(step_rows)
    if use_popdep:
        n_final = int(len(final_all))
        totals_df = pd.DataFrame(
            {
                "step": ["raw", "popdep_60", "popdep_80", "het_pass", "mac_ge5", "wap_mac_ge2"],
                "n_sites": [
                    raw_total,
                    popdep60_total,
                    popdep80_total,
                    het_total,
                    mac5_total,
                    n_final,
                ],
            }
        )
        totals_df["frac_retain_of_raw"] = totals_df["n_sites"] / raw_total if raw_total else np.nan
        totals_df["frac_drop_in_step"] = [
            np.nan,
            np.nan,
            1 - popdep80_total / raw_total if raw_total else np.nan,
            1 - het_total / popdep80_total if popdep80_total else np.nan,
            1 - mac5_total / het_total if het_total else np.nan,
            1 - n_final / mac5_total if mac5_total else np.nan,
        ]
    else:
        totals_df = pd.DataFrame(
            {
                "step": ["raw", "het_pass", "mac_ge5", "wap_mac_ge2"],
                "n_sites": [raw_total, het_total, mac5_total, int(len(final_all))],
            }
        )
        totals_df["frac_retain_of_raw"] = totals_df["n_sites"] / raw_total if raw_total else np.nan
        totals_df["frac_drop_in_step"] = [
            np.nan,
            1 - het_total / raw_total if raw_total else np.nan,
            1 - mac5_total / het_total if het_total else np.nan,
            1 - len(final_all) / mac5_total if mac5_total else np.nan,
        ]
    return final_all, summary, totals_df


def _maf_class_series(maf: pd.Series, *, threshold: float = 0.05) -> pd.Series:
    m = pd.to_numeric(maf, errors="coerce")
    return pd.Series(np.where(m < threshold, "rare", "common"), index=maf.index)


def summarize_stepwise_maf_class_counts(
    process_dir: str | Path,
    *,
    het_paths: dict[str, str | Path],
    final_sites: pd.DataFrame,
    subgenomes: tuple[str, ...] = ("A", "B", "D"),
    alpha: float = RETAIN_ALPHA,
    maf_threshold: float = 0.05,
) -> pd.DataFrame:
    """Count rare/common at each cascade step (full-panel MAF until WAP; WAP MAF at final)."""
    rows: list[dict] = []
    process_dir = Path(process_dir)

    for sg in subgenomes:
        het_path = het_paths.get(sg)
        if het_path is None:
            raise ValueError(f"Missing het path for subgenome {sg}")
        f_panel = load_median_inbreeding_from_het(str(het_path))
        sites = _load_subgenome_frame(process_dir, sg, f_panel=f_panel, alpha=alpha)
        sites["maf_class"] = _maf_class_series(sites["MAF"], threshold=maf_threshold)

        step_masks = {
            "after_popdep": pd.Series(True, index=sites.index),
            "het_pass": sites["het_pass"],
            "mac_ge5": sites["het_pass"] & sites["mac_ge5"],
        }
        for step, mask in step_masks.items():
            sub = sites.loc[mask]
            for cls in MAF_CLASS_SHORT:
                n_cls = int((sub["maf_class"] == cls).sum())
                rows.append({"subgenome": sg, "step": step, "maf_class": cls, "n_sites": n_cls})
            rows.append({"subgenome": sg, "step": step, "maf_class": "total", "n_sites": int(len(sub))})

        fin = final_sites[final_sites["subgenome"] == sg].copy()
        fin["maf_class"] = _maf_class_series(fin["MAF_WAP"], threshold=maf_threshold)
        for cls in MAF_CLASS_SHORT:
            rows.append(
                {
                    "subgenome": sg,
                    "step": "wap_mac_ge2",
                    "maf_class": cls,
                    "n_sites": int((fin["maf_class"] == cls).sum()),
                }
            )
        rows.append(
            {"subgenome": sg, "step": "wap_mac_ge2", "maf_class": "total", "n_sites": int(len(fin))}
        )

    out = pd.DataFrame(rows)
    totals = out[out["maf_class"] == "total"][["subgenome", "step", "n_sites"]].rename(
        columns={"n_sites": "n_step_total"}
    )
    out = out.merge(totals, on=["subgenome", "step"], how="left")
    out["frac_within_step"] = np.where(
        out["maf_class"].isin(MAF_CLASS_SHORT),
        out["n_sites"] / out["n_step_total"],
        np.nan,
    )
    return out.drop(columns=["n_step_total"])


def subgenome_final_maf_proportions(
    final_sites: pd.DataFrame,
    *,
    subgenomes: tuple[str, ...] = ("A", "B", "D"),
    maf_col: str = "MAF_WAP",
    maf_threshold: float = 0.05,
) -> pd.DataFrame:
    """Return fraction rare/common per subgenome (rows sum to 1)."""
    prop_rows: list[dict] = []
    for sg in subgenomes:
        sub = final_sites[final_sites["subgenome"] == sg]
        maf = pd.to_numeric(sub[maf_col], errors="coerce")
        n = int(maf.notna().sum())
        if n == 0:
            continue
        prop_rows.append(
            {
                "subgenome": sg,
                "frac_rare": float((maf < maf_threshold).sum()) / n,
                "frac_common": float((maf >= maf_threshold).sum()) / n,
            }
        )
    prop = pd.DataFrame(prop_rows).set_index("subgenome")
    return prop[["frac_rare", "frac_common"]].rename(
        columns={"frac_rare": RARE_COMMON_LABELS[0], "frac_common": RARE_COMMON_LABELS[1]}
    )


def _step_transition_retention(class_counts: pd.DataFrame) -> pd.DataFrame:
    """Per-step rare/common retention vs previous step (ABD pooled)."""
    cls = class_counts[class_counts["maf_class"].isin(MAF_CLASS_SHORT)].copy()
    pooled = cls.groupby(["step", "maf_class"], as_index=False)["n_sites"].sum()
    wide = pooled.pivot(index="step", columns="maf_class", values="n_sites").reindex(STEP_ORDER)
    rows: list[dict] = []
    prev_step = None
    for step in STEP_ORDER:
        if prev_step is None:
            prev_step = step
            continue
        for maf_class in MAF_CLASS_SHORT:
            prev_n = float(wide.loc[prev_step, maf_class])
            cur_n = float(wide.loc[step, maf_class])
            rows.append(
                {
                    "transition": f"{prev_step}_to_{step}",
                    "maf_class": maf_class,
                    "retain_frac": cur_n / prev_n if prev_n else np.nan,
                }
            )
        prev_step = step
    return pd.DataFrame(rows)


def project_genome_pipeline_by_maf_class(
    class_counts: pd.DataFrame,
    *,
    pipeline_totals: dict[str, float] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Forward-project genome totals using test_thin rare/common retention rates.

    Returns (step_projection, comparison_vs_unstratified).
    """
    pipeline_totals = pipeline_totals or GENOME_PIPELINE_TOTALS
    cls = class_counts[class_counts["maf_class"].isin(MAF_CLASS_SHORT)].copy()
    pooled = cls.groupby(["step", "maf_class"], as_index=False)["n_sites"].sum()
    wide = pooled.pivot(index="step", columns="maf_class", values="n_sites").reindex(STEP_ORDER).fillna(0)

    start = STEP_ORDER[0]
    n_rare = float(wide.loc[start, "rare"])
    n_common = float(wide.loc[start, "common"])
    scale = pipeline_totals[start] / (n_rare + n_common) if (n_rare + n_common) else 1.0
    state = {"rare": n_rare * scale, "common": n_common * scale}

    proj_rows: list[dict] = []
    for step in STEP_ORDER:
        n_total = state["rare"] + state["common"]
        proj_rows.append(
            {
                "step": step,
                "maf_class": "rare",
                "n_sites": state["rare"],
                "frac_rare": state["rare"] / n_total if n_total else np.nan,
            }
        )
        proj_rows.append(
            {
                "step": step,
                "maf_class": "common",
                "n_sites": state["common"],
                "frac_rare": state["rare"] / n_total if n_total else np.nan,
            }
        )
        proj_rows.append(
            {
                "step": step,
                "maf_class": "total",
                "n_sites": n_total,
                "frac_rare": state["rare"] / n_total if n_total else np.nan,
            }
        )
        if step == STEP_ORDER[-1]:
            break
        next_step = STEP_ORDER[STEP_ORDER.index(step) + 1]
        for maf_class in MAF_CLASS_SHORT:
            prev_n = float(wide.loc[step, maf_class])
            next_n = float(wide.loc[next_step, maf_class])
            retain = next_n / prev_n if prev_n else 1.0
            state[maf_class] *= retain

    projection = pd.DataFrame(proj_rows)
    compare = pd.DataFrame(
        [
            {
                "step": step,
                "n_unstratified": pipeline_totals[step],
                "n_stratified": float(projection.loc[(projection["step"] == step) & (projection["maf_class"] == "total"), "n_sites"].iloc[0]),
            }
            for step in STEP_ORDER
        ]
    )
    compare["n_stratified_rare"] = [
        float(projection.loc[(projection["step"] == s) & (projection["maf_class"] == "rare"), "n_sites"].iloc[0])
        for s in STEP_ORDER
    ]
    compare["n_stratified_common"] = [
        float(projection.loc[(projection["step"] == s) & (projection["maf_class"] == "common"), "n_sites"].iloc[0])
        for s in STEP_ORDER
    ]
    compare["delta_total"] = compare["n_stratified"] - compare["n_unstratified"]
    return projection, compare


def plot_subgenome_rare_common_proportion_bar(
    final_sites: pd.DataFrame,
    output_prefix: str | Path,
    *,
    subgenomes: tuple[str, ...] = ("A", "B", "D"),
    maf_col: str = "MAF_WAP",
    maf_threshold: float = 0.05,
) -> pd.DataFrame:
    """Horizontal 100% stacked bars of rare vs common fractions for each subgenome."""
    prop_df = subgenome_final_maf_proportions(
        final_sites,
        subgenomes=subgenomes,
        maf_col=maf_col,
        maf_threshold=maf_threshold,
    )
    if prop_df.empty:
        print(f"[Warning] No subgenome proportions to plot for {output_prefix}")
        return prop_df

    plot_order = [sg for sg in reversed(subgenomes) if sg in prop_df.index]
    plot_stacked_proportion_barh(
        prop_df,
        title="Rare vs common fraction after het + MAC≥5 + WAP MAC≥2 (WAP MAF)",
        filename=str(_out_path(output_prefix, "subgenome_maf_rare_common_prop.png")),
        y_label="Subgenome",
        x_label="Fraction of variants",
        stack_col_title="MAF class",
        palette_map=RARE_COMMON_COLORS,
        figure_size=(8, 2.8),
        row_order=plot_order,
    )
    save_df_to_tsv(
        prop_df.reindex([sg for sg in subgenomes if sg in prop_df.index]).reset_index(),
        str(_out_path(output_prefix, "subgenome_maf_rare_common_prop.tsv")),
    )
    return prop_df


def plot_maf_rare_common_distribution(
    sites: pd.DataFrame,
    output_prefix: str | Path,
    *,
    maf_col: str = "MAF",
    maf_threshold: float = 0.05,
    bins: int = 100,
    title: str | None = None,
) -> None:
    """Stacked MAF histogram: rare vs common split at ``maf_threshold``."""
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    plot_df = sites[[maf_col]].copy()
    plot_df[maf_col] = pd.to_numeric(plot_df[maf_col], errors="coerce")
    plot_df = plot_df.dropna()
    plot_df = plot_df[(plot_df[maf_col] >= 0) & (plot_df[maf_col] <= 0.5)]
    if plot_df.empty:
        print(f"[Warning] No MAF values to plot for {output_prefix}")
        return

    plot_df["maf_class"] = assign_rare_common_class(plot_df[maf_col], threshold=maf_threshold)
    n_rare = int((plot_df[maf_col] < maf_threshold).sum())
    n_common = int((plot_df[maf_col] >= maf_threshold).sum())
    n_total = len(plot_df)

    plot_stacked_histogram(
        plot_df,
        x_col=maf_col,
        hue_col="maf_class",
        hue_order=list(RARE_COMMON_LABELS),
        palette_map=RARE_COMMON_COLORS,
        title=(
            title
            if title is not None
            else (
                f"MAF distribution after het + MAC≥5 + WAP MAC≥2 "
                f"(n={n_total:,}; rare={n_rare:,}, common={n_common:,})"
            )
        ),
        filename=str(_out_path(output_prefix, "maf_rare_common.png")),
        xlabel="Minor allele frequency (WAP panel)",
        ylabel="Variant count",
        xlim=(0.0, 0.5),
        legend_title="MAF class",
    )

    save_df_to_tsv(
        pd.DataFrame(
            [
                {
                    "n_total": n_total,
                    "n_rare_maf_lt_0p05": n_rare,
                    "n_common_maf_ge_0p05": n_common,
                    "frac_rare": n_rare / n_total if n_total else np.nan,
                    "frac_common": n_common / n_total if n_total else np.nan,
                    "maf_threshold": maf_threshold,
                }
            ]
        ),
        str(_out_path(output_prefix, "maf_rare_common.th.tsv")),
    )


def run_filter_cascade_maf_plot(
    process_dir: str | Path,
    output_prefix: str | Path,
    *,
    het_paths: dict[str, str | Path],
    wap_keep_path: str | Path,
    alpha: float = RETAIN_ALPHA,
    plink2: str = "plink2",
) -> pd.DataFrame:
    """Run cascade filters, write step table + MAF plot."""
    output_prefix = Path(output_prefix)
    final, by_sg, totals = filter_cascade_test_thin(
        process_dir,
        het_paths=het_paths,
        wap_keep_path=wap_keep_path,
        alpha=alpha,
        plink2=plink2,
    )
    save_df_to_tsv(by_sg, str(_out_path(output_prefix, "filter_steps_by_subgenome.tsv")))
    save_df_to_tsv(totals, str(_out_path(output_prefix, "filter_steps_total.tsv")))
    save_thresholds(
        {
            "N_raw": int(totals.loc[totals["step"] == "raw", "n_sites"].iloc[0]),
            "N_het_pass": int(totals.loc[totals["step"] == "het_pass", "n_sites"].iloc[0]),
            "N_mac_ge5": int(totals.loc[totals["step"] == "mac_ge5", "n_sites"].iloc[0]),
            "N_wap_mac_ge2": int(totals.loc[totals["step"] == "wap_mac_ge2", "n_sites"].iloc[0]),
            "alpha_het": float(alpha),
        },
        str(_out_path(output_prefix, "filter.th.tsv")),
    )
    if not final.empty:
        save_df_to_tsv(
            final[["subgenome", "ID", "MAC", "MAF", "MAC_WAP", "MAF_WAP", "p_excess"]],
            str(_out_path(output_prefix, "filtered_sites.tsv")),
        )
    plot_maf_rare_common_distribution(final, output_prefix)
    class_counts = summarize_stepwise_maf_class_counts(
        process_dir,
        het_paths=het_paths,
        final_sites=final,
        subgenomes=("A", "B", "D"),
        alpha=alpha,
    )
    save_df_to_tsv(class_counts, str(_out_path(output_prefix, "stepwise_maf_class_counts.tsv")))
    retention = _step_transition_retention(class_counts)
    save_df_to_tsv(retention, str(_out_path(output_prefix, "stepwise_maf_class_retention.tsv")))
    projection, compare = project_genome_pipeline_by_maf_class(class_counts)
    save_df_to_tsv(projection, str(_out_path(output_prefix, "genome_pipeline_maf_class_projection.tsv")))
    save_df_to_tsv(compare, str(_out_path(output_prefix, "genome_pipeline_maf_class_compare.tsv")))
    plot_subgenome_rare_common_proportion_bar(final, output_prefix, subgenomes=("A", "B", "D"))
    return final


def estimate_full_library_from_thin(
    totals: pd.DataFrame,
    *,
    thin_rate: float,
) -> pd.DataFrame:
    """Scale test-library step counts by ``1 / thin_rate`` (test_thin default 0.001)."""
    rate = float(thin_rate)
    if rate <= 0:
        raise ValueError(f"thin_rate must be positive, got {thin_rate}")
    out = totals.copy()
    out["n_estimate_full"] = out["n_sites"] / rate
    return out


def _keep_ids_from_popdep_table(
    popdep_path: str | Path,
    mask: np.ndarray,
) -> set[str]:
    from genetics.genomics.variant.popdep_density_filter import sites_in_trusted_cells
    from infra.utils.io import load_df_from_tsv

    table = load_df_from_tsv(str(popdep_path))
    table.columns = [c.replace("#", "") if isinstance(c, str) else c for c in table.columns]
    kept = sites_in_trusted_cells(table, mask)
    if kept.empty or "ID" not in kept.columns:
        return set()
    return set(kept["ID"].astype(str))


def run_filter_cascade_popdep_density(
    process_dir: str | Path,
    output_prefix: str | Path,
    *,
    het_paths: dict[str, str | Path],
    wap_keep_path: str | Path,
    mask_60: np.ndarray,
    mask_80: np.ndarray,
    popdep_tables: dict[str, str | Path],
    genome_density_stats: dict[str, int] | None = None,
    thin_rate: float = 0.001,
    alpha: float = RETAIN_ALPHA,
    plink2: str = "plink2",
    skip_wap: bool = False,
) -> pd.DataFrame:
    """Filter test pfiles by PopDep 80% cells, then het / MAC≥5 / WAP MAC≥2."""
    output_prefix = Path(output_prefix)
    keep_60: set[str] = set()
    keep_80: set[str] = set()
    for _sg, path in popdep_tables.items():
        keep_60 |= _keep_ids_from_popdep_table(path, mask_60)
        keep_80 |= _keep_ids_from_popdep_table(path, mask_80)

    final, by_sg, totals = filter_cascade_test_thin(
        process_dir,
        het_paths=het_paths,
        wap_keep_path=wap_keep_path,
        alpha=alpha,
        plink2=plink2,
        popdep_keep_ids=keep_80,
        popdep_keep_ids_60=keep_60,
        skip_wap=skip_wap,
    )
    estimate = estimate_full_library_from_thin(totals, thin_rate=thin_rate)
    save_df_to_tsv(by_sg, str(_out_path(output_prefix, "filter_steps_by_subgenome.tsv")))
    save_df_to_tsv(totals, str(_out_path(output_prefix, "filter_steps_total.tsv")))
    save_df_to_tsv(estimate, str(_out_path(output_prefix, "filter_steps_genome_estimate.tsv")))
    th = {
        "N_raw": int(totals.loc[totals["step"] == "raw", "n_sites"].iloc[0]),
        "N_popdep_60": int(totals.loc[totals["step"] == "popdep_60", "n_sites"].iloc[0]),
        "N_popdep_80": int(totals.loc[totals["step"] == "popdep_80", "n_sites"].iloc[0]),
        "N_het_pass": int(totals.loc[totals["step"] == "het_pass", "n_sites"].iloc[0]),
        "N_mac_ge5": int(totals.loc[totals["step"] == "mac_ge5", "n_sites"].iloc[0]),
        "N_wap_mac_ge2": int(totals.loc[totals["step"] == "wap_mac_ge2", "n_sites"].iloc[0]),
        "alpha_het": float(alpha),
        "thin_rate": float(thin_rate),
        "N_estimate_full_final": float(
            estimate.loc[estimate["step"] == "wap_mac_ge2", "n_estimate_full"].iloc[0]
        ),
    }
    if genome_density_stats:
        th["genome_popdep_n_rows"] = int(genome_density_stats.get("n_rows", 0))
        th["genome_popdep_n_valid"] = int(genome_density_stats.get("n_valid", 0))
        th["genome_popdep_n_top60"] = int(genome_density_stats.get("n_sites_top60", 0))
        th["genome_popdep_n_top80"] = int(genome_density_stats.get("n_sites_top80", 0))
    save_thresholds(th, str(_out_path(output_prefix, "filter.th.tsv")))
    if not final.empty:
        cols = [c for c in ("subgenome", "ID", "MAC", "MAF", "MAC_WAP", "MAF_WAP", "p_excess") if c in final.columns]
        save_df_to_tsv(final[cols], str(_out_path(output_prefix, "filtered_sites.tsv")))
    if not final.empty and "MAF" in final.columns:
        plot_maf_rare_common_distribution(final, output_prefix)
    return final


def _chrom_key_from_bed_name(path: Path) -> str:
    """``chr001.popdep.cv_reliable.keep80.bed.gz`` → ``\"1\"`` (PLINK #CHROM)."""
    token = path.name.split(".", 1)[0]
    if not (token.startswith("chr") and token[3:].isdigit()):
        raise ValueError(f"BED name must start with chrNNN: {path.name}")
    return str(int(token[3:]))


def load_bed_intervals(path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """0-based half-open starts and ends from one chromosome BED."""
    starts: list[int] = []
    ends: list[int] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith(("#", "track", "browser")):
                continue
            _chrom, start_s, end_s = line.rstrip("\n").split("\t")[:3]
            starts.append(int(start_s))
            ends.append(int(end_s))
    return np.asarray(starts, dtype=np.int64), np.asarray(ends, dtype=np.int64)


def variant_ids_in_reliable_beds(
    pvar_paths,
    bed_paths,
) -> set[str]:
    """Variant IDs whose 1-based POS lies in a PopDep reliable-site BED.

    BED files are named ``chrNNN.*`` with NNN equal to the PLINK ``#CHROM``.
    Coordinates are 0-based half-open.
    """
    beds: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for path in bed_paths:
        beds[_chrom_key_from_bed_name(Path(path))] = load_bed_intervals(path)
    grouped: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for pvar in pvar_paths:
        with open(pvar, encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                chrom, pos_s, variant_id = line.split("\t", 3)[:3]
                grouped[chrom].append((int(pos_s), variant_id))
    kept: set[str] = set()
    for chrom, items in grouped.items():
        intervals = beds.get(chrom)
        if intervals is None or intervals[0].size == 0:
            continue
        starts, ends = intervals
        coord = np.fromiter((pos - 1 for pos, _vid in items), dtype=np.int64, count=len(items))
        idx = np.searchsorted(ends, coord, side="right")
        inside = idx < starts.size
        clipped = np.minimum(idx, starts.size - 1)
        inside &= starts[clipped] <= coord
        for hit, (_pos, variant_id) in zip(inside, items, strict=True):
            if bool(hit):
                kept.add(variant_id)
    return kept


def _single_popdep_step_tables(
    by_sg: pd.DataFrame,
    totals: pd.DataFrame,
    keep_pct: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Drop the unused 60% count-only row and name the filter step."""
    step = f"popdep_{int(keep_pct)}"
    totals_out = totals.loc[totals["step"] != "popdep_60"].copy()
    totals_out["step"] = totals_out["step"].replace({"popdep_80": step})
    by_out = by_sg.drop(columns=["n_popdep_60"], errors="ignore").rename(
        columns={"n_popdep_80": f"n_{step}"}
    )
    return by_out, totals_out


def run_filter_cascade_cv_reliable(
    process_dir: str | Path,
    output_prefix: str | Path,
    *,
    bed_paths,
    het_paths: dict[str, str | Path],
    wap_keep_path: str | Path,
    keep_pct: int,
    subgenomes: tuple[str, ...] = ("A", "B", "D", "Others"),
    thin_rate: float = 0.001,
    alpha: float = RETAIN_ALPHA,
    plink2: str = "plink2",
    skip_wap: bool = False,
) -> pd.DataFrame:
    """PopDep reliable-BED keep, then het → MAC≥5 → WAP MAC≥2, scaled by thin_rate."""
    process_dir = Path(process_dir)
    output_prefix = Path(output_prefix)
    pvars = [process_dir / f"{sg}_test.plink2.pvar" for sg in subgenomes]
    missing = [str(path) for path in pvars if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing pvar: {', '.join(missing)}")
    keep_ids = variant_ids_in_reliable_beds(pvars, bed_paths)
    print(
        f"[Info] popdep keep{int(keep_pct)}: {len(keep_ids)} test_thin variants in reliable BED",
        flush=True,
    )
    final, by_sg, totals = filter_cascade_test_thin(
        process_dir,
        het_paths=het_paths,
        wap_keep_path=None if skip_wap else wap_keep_path,
        subgenomes=subgenomes,
        alpha=alpha,
        plink2=plink2,
        popdep_keep_ids=keep_ids,
        skip_wap=skip_wap,
    )
    by_sg, totals = _single_popdep_step_tables(by_sg, totals, keep_pct)
    estimate = estimate_full_library_from_thin(totals, thin_rate=thin_rate)
    save_df_to_tsv(by_sg, str(_out_path(output_prefix, "filter_steps_by_subgenome.tsv")))
    save_df_to_tsv(totals, str(_out_path(output_prefix, "filter_steps_total.tsv")))
    save_df_to_tsv(estimate, str(_out_path(output_prefix, "filter_steps_genome_estimate.tsv")))
    step = f"popdep_{int(keep_pct)}"
    save_thresholds(
        {
            "keep_pct": int(keep_pct),
            "N_raw": int(totals.loc[totals["step"] == "raw", "n_sites"].iloc[0]),
            "N_popdep": int(totals.loc[totals["step"] == step, "n_sites"].iloc[0]),
            "N_het_pass": int(totals.loc[totals["step"] == "het_pass", "n_sites"].iloc[0]),
            "N_mac_ge5": int(totals.loc[totals["step"] == "mac_ge5", "n_sites"].iloc[0]),
            "N_wap_mac_ge2": int(totals.loc[totals["step"] == "wap_mac_ge2", "n_sites"].iloc[0]),
            "thin_rate": float(thin_rate),
            "N_estimate_full_final": float(
                estimate.loc[estimate["step"] == "wap_mac_ge2", "n_estimate_full"].iloc[0]
            ),
            "alpha_het": float(alpha),
        },
        str(_out_path(output_prefix, "filter.th.tsv")),
    )
    if not final.empty:
        cols = [c for c in ("subgenome", "ID", "MAC", "MAF", "MAC_WAP", "MAF_WAP", "p_excess") if c in final.columns]
        save_df_to_tsv(final[cols], str(_out_path(output_prefix, "filtered_sites.tsv")))
    if not final.empty and "MAF" in final.columns and not skip_wap:
        plot_maf_rare_common_distribution(final, output_prefix)
    return final


def parse_segment_variant_id(variant_id: str) -> tuple[str, int] | None:
    """``{chr_id}-{local_pos}`` from a pvar ID, or None when the ID is not that form."""
    chrom, sep, pos = str(variant_id).partition("-")
    if not sep or not chrom.isdigit() or not pos.isdigit():
        return None
    return chrom, int(pos)


def variant_ids_in_reliable_beds_from_ids(pvar_paths, bed_paths) -> set[str]:
    """Variant IDs whose arm-local position falls in a PopDep reliable-site BED.

    The ID is ``{chr_id}-{local_pos}`` (1-based on the arm). ``*_refchr`` pvars
    store a linear ``#CHROM`` and ``POS``, so membership uses the ID rather than
    those columns. BED files are named ``chrNNN.*`` and are 0-based half-open.
    """
    beds: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for path in bed_paths:
        beds[_chrom_key_from_bed_name(Path(path))] = load_bed_intervals(path)
    grouped: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for pvar in pvar_paths:
        with open(pvar, encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                variant_id = line.split("\t", 3)[2]
                parsed = parse_segment_variant_id(variant_id)
                if parsed is None:
                    continue
                chrom, pos = parsed
                grouped[chrom].append((pos, variant_id))
    kept: set[str] = set()
    for chrom, items in grouped.items():
        intervals = beds.get(chrom)
        if intervals is None or intervals[0].size == 0:
            continue
        starts, ends = intervals
        coord = np.fromiter((pos - 1 for pos, _vid in items), dtype=np.int64, count=len(items))
        idx = np.searchsorted(ends, coord, side="right")
        inside = idx < starts.size
        clipped = np.minimum(idx, starts.size - 1)
        inside &= starts[clipped] <= coord
        for hit, (_pos, variant_id) in zip(inside, items, strict=True):
            if bool(hit):
                kept.add(variant_id)
    return kept


def write_plink2_het(pfile_prefix: str | Path, out_prefix: str | Path, plink2: str = "plink2") -> Path:
    """Sample inbreeding table from an existing pfile. Writes only ``out_prefix.het``."""
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [plink2, "--pfile", str(pfile_prefix), "--het", "--out", str(out_prefix)],
        check=True,
    )
    het_path = Path(f"{out_prefix}.het")
    if not het_path.exists():
        raise FileNotFoundError(f"plink2 --het did not write {het_path}")
    return het_path


def discover_pfile_subgenomes(process_dir: str | Path) -> tuple[str, ...]:
    """Subgenomes that have both a pvar and a geno-counts table."""
    process_dir = Path(process_dir)
    found = []
    for sg in ("A", "B", "D", "Others"):
        pvar = process_dir / f"{sg}_test.plink2.pvar"
        gcount = process_dir / "variant" / f"{sg}.info.gcount"
        if pvar.exists() and gcount.exists():
            found.append(sg)
    if not found:
        raise FileNotFoundError(f"No subgenome pvar+gcount under {process_dir}")
    return tuple(found)


def het_paths_from_pfiles(
    process_dir: str | Path,
    subgenomes: tuple[str, ...],
    het_dir: str | Path,
    plink2: str = "plink2",
) -> dict[str, Path]:
    """Median-F inputs for each subgenome, written under ``het_dir`` (not the pfile dir)."""
    process_dir = Path(process_dir)
    het_dir = Path(het_dir)
    paths: dict[str, Path] = {}
    for sg in subgenomes:
        print(f"[Info] plink2 --het {sg}", flush=True)
        paths[sg] = write_plink2_het(process_dir / f"{sg}_test.plink2", het_dir / sg, plink2=plink2)
    return paths


def _plot_filter_step_counts(totals: pd.DataFrame, output_prefix: Path) -> None:
    steps = totals["step"].astype(str).tolist()
    counts = [int(value) for value in totals["n_sites"].tolist()]
    ymax = max(counts) if counts else 1
    plot_bar_chart(
        steps,
        counts,
        title="Sites retained after each filter",
        ylabel="Variant count",
        filename=str(_out_path(output_prefix, "filter_steps.png")),
        ylim=(0.0, ymax * 1.2),
        rotate_xlabels=20,
        figure_size=(8, 5),
    )


def run_cv_reliable_three_step(
    process_dir: str | Path,
    output_prefix: str | Path,
    *,
    bed_paths,
    keep_pct: int,
    het_paths: dict[str, str | Path] | None = None,
    subgenomes: tuple[str, ...] | None = None,
    alpha: float = RETAIN_ALPHA,
    plink2: str = "plink2",
    het_dir: str | Path | None = None,
    min_mac: int = 5,
) -> pd.DataFrame:
    """PopDep reliable BED, excess het, then MAC≥``min_mac``. Does not rewrite pfiles.

    Inbreeding F comes from ``het_paths`` or from ``plink2 --het`` on each
    subgenome pfile (written under ``het_dir``). The separate WAP MAC≥2 step is
    omitted because these panels are already the WAP sample set.
    """
    process_dir = Path(process_dir)
    output_prefix = Path(output_prefix)
    if subgenomes is None:
        subgenomes = discover_pfile_subgenomes(process_dir)
    if het_paths is None:
        het_paths = het_paths_from_pfiles(
            process_dir,
            subgenomes,
            Path(het_dir) if het_dir is not None else output_prefix.parent / f"het_keep{int(keep_pct)}",
            plink2=plink2,
        )
    pvars = [process_dir / f"{sg}_test.plink2.pvar" for sg in subgenomes]
    keep_ids = variant_ids_in_reliable_beds_from_ids(pvars, bed_paths)
    print(
        f"[Info] popdep keep{int(keep_pct)}: {len(keep_ids)} variants in reliable BED",
        flush=True,
    )
    final, by_sg, totals = filter_cascade_test_thin(
        process_dir,
        het_paths=het_paths,
        wap_keep_path=None,
        subgenomes=subgenomes,
        alpha=alpha,
        plink2=plink2,
        popdep_keep_ids=keep_ids,
        skip_wap=True,
        min_mac=int(min_mac),
    )
    by_sg, totals = _single_popdep_step_tables(by_sg, totals, keep_pct)
    mac_step = f"mac_ge{int(min_mac)}"
    if int(min_mac) == 5:
        by_sg = by_sg.drop(columns=["n_wap_mac_ge2"], errors="ignore")
        totals = totals.loc[totals["step"] != "wap_mac_ge2"].copy()
    else:
        by_sg = by_sg.copy()
        by_sg[f"n_{mac_step}"] = by_sg["n_wap_mac_ge2"]
        by_sg = by_sg.drop(columns=["n_mac_ge5", "n_wap_mac_ge2"], errors="ignore")
        totals = totals.copy()
        totals.loc[totals["step"] == "mac_ge5", "step"] = mac_step
        totals.loc[totals["step"] == mac_step, "n_sites"] = int(len(final))
        totals = totals.loc[totals["step"] != "wap_mac_ge2"].copy()
        n_sites = totals["n_sites"].to_numpy(dtype=float)
        raw_n = float(n_sites[0]) if len(n_sites) else np.nan
        totals["frac_retain_of_raw"] = n_sites / raw_n if raw_n else np.nan
        drops = [np.nan]
        for prev, cur in zip(n_sites[:-1], n_sites[1:]):
            drops.append(1.0 - cur / prev if prev else np.nan)
        totals["frac_drop_in_step"] = drops
    save_df_to_tsv(by_sg, str(_out_path(output_prefix, "filter_steps_by_subgenome.tsv")))
    save_df_to_tsv(totals, str(_out_path(output_prefix, "filter_steps_total.tsv")))
    step = f"popdep_{int(keep_pct)}"
    save_thresholds(
        {
            "keep_pct": int(keep_pct),
            "N_raw": int(totals.loc[totals["step"] == "raw", "n_sites"].iloc[0]),
            "N_popdep": int(totals.loc[totals["step"] == step, "n_sites"].iloc[0]),
            "N_het_pass": int(totals.loc[totals["step"] == "het_pass", "n_sites"].iloc[0]),
            f"N_{mac_step}": int(totals.loc[totals["step"] == mac_step, "n_sites"].iloc[0]),
            "alpha_het": float(alpha),
            "min_mac": int(min_mac),
            "steps": f"popdep,het,{mac_step}",
        },
        str(_out_path(output_prefix, "filter.th.tsv")),
    )
    cols = [c for c in ("subgenome", "ID", "MAC", "MAF", "p_excess") if c in final.columns or c == "ID"]
    if final.empty:
        save_df_to_tsv(pd.DataFrame(columns=cols), str(_out_path(output_prefix, "filtered_sites.tsv")))
    else:
        save_df_to_tsv(final[cols], str(_out_path(output_prefix, "filtered_sites.tsv")))
    _plot_filter_step_counts(totals, output_prefix)
    if not final.empty and "MAF" in final.columns:
        plot_maf_rare_common_distribution(
            final,
            output_prefix,
            title=(
                f"MAF after PopDep keep{int(keep_pct)} + het + MAC≥{int(min_mac)} "
                f"(n={len(final):,})"
            ),
        )
    return final

