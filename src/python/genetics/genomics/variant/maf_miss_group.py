"""Variant MAF / missing plots stratified by sample_group (dominant missing group)."""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from genetics.genomics.variant.variant_utils import load_df_from_plink_variant
from infra.utils.graph import (
    plot_joint_regression,
    plot_regression_comparison,
    plot_stacked_distribution,
    plot_stacked_proportion_bar,
)
from infra.utils.io import load_df_from_space_sep_no_header, save_df_to_tsv


def _maf_from_afreq_df(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "ALT_FREQS" in out.columns:
        out["MAF"] = out["ALT_FREQS"].apply(lambda x: x if x <= 0.5 else 1 - x)
    elif "AF" in out.columns:
        out["MAF"] = out["AF"].apply(lambda x: x if x <= 0.5 else 1 - x)
    else:
        raise ValueError("afreq table missing ALT_FREQS / AF column")
    return out


def _load_psam_iids(psam_path: Path) -> set[str]:
    psam = pd.read_csv(psam_path, sep="\t")
    psam.columns = [c.lstrip("#") for c in psam.columns]
    iid_col = "IID" if "IID" in psam.columns else psam.columns[0]
    return set(psam[iid_col].astype(str))


def _write_keep_file(samples: list[str], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for sample in samples:
            handle.write(f"0\t{sample}\n")


def _run_plink_group_vmiss(pfile_prefix: str, keep_path: Path, out_prefix: Path) -> pd.DataFrame:
    cmd = [
        "plink2",
        "--pfile",
        pfile_prefix,
        "--keep",
        str(keep_path),
        "--missing",
        "--out",
        str(out_prefix),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    vmiss_path = Path(f"{out_prefix}.vmiss")
    df = load_df_from_plink_variant(str(vmiss_path))
    if df is None or df.empty:
        return pd.DataFrame(columns=["ID", "F_MISS"])
    return df[["ID", "F_MISS"]].rename(columns={"F_MISS": "F_MISS_group"})


def _group_stratified_vmiss(
    pfile_prefix: str,
    group_file: str,
    *,
    min_samples: int = 5,
) -> pd.DataFrame:
    """Long table: ID, Group, F_MISS_group for each sample_group present in the pfile."""
    psam_path = Path(f"{pfile_prefix}.psam")
    if not psam_path.exists():
        raise FileNotFoundError(f"Missing psam for group-stratified vmiss: {psam_path}")

    panel_samples = _load_psam_iids(psam_path)
    group_df = load_df_from_space_sep_no_header(group_file, ["Sample", "Group"])
    if group_df is None or group_df.empty:
        raise ValueError(f"Could not load group file: {group_file}")

    group_df = group_df.copy()
    group_df["Sample"] = group_df["Sample"].astype(str)
    group_df = group_df[group_df["Sample"].isin(panel_samples)]

    long_parts: list[pd.DataFrame] = []
    with tempfile.TemporaryDirectory(prefix="variant_group_vmiss_") as tmp:
        tmp_dir = Path(tmp)
        for group, members in group_df.groupby("Group"):
            samples = members["Sample"].astype(str).tolist()
            if len(samples) < min_samples:
                print(f"[Warning] Skipping group {group}: only {len(samples)} samples in panel.")
                continue
            keep_path = tmp_dir / f"{group}.keep"
            out_prefix = tmp_dir / f"{group}.vmiss"
            _write_keep_file(samples, keep_path)
            print(f"[Info] Group-stratified vmiss: {group} (n={len(samples)})")
            grp_vmiss = _run_plink_group_vmiss(pfile_prefix, keep_path, out_prefix)
            grp_vmiss["Group"] = str(group)
            long_parts.append(grp_vmiss)

    if not long_parts:
        raise ValueError("No sample groups with enough members for group-stratified vmiss.")
    return pd.concat(long_parts, ignore_index=True)


def _assign_dominant_missing_group(group_long: pd.DataFrame) -> pd.DataFrame:
    wide = group_long.pivot_table(
        index="ID",
        columns="Group",
        values="F_MISS_group",
        aggfunc="first",
    )
    dominant = wide.idxmax(axis=1, skipna=True)
    max_miss = wide.max(axis=1, skipna=True)
    out = dominant.rename("Group").reset_index()
    out["Group_max_F_MISS"] = max_miss.values
    out["Group"] = out["Group"].fillna("Unknown")
    return out


MAF_BIN_LABELS = ("0-0.01", "0.01-0.05", "0.05-0.1", "0.1-0.5")


def _assign_maf_bin(maf: float) -> str | None:
    if pd.isna(maf) or maf < 0 or maf > 0.5:
        return None
    if maf < 0.01:
        return MAF_BIN_LABELS[0]
    if maf < 0.05:
        return MAF_BIN_LABELS[1]
    if maf < 0.1:
        return MAF_BIN_LABELS[2]
    return MAF_BIN_LABELS[3]


def _maf_bin_group_fraction_table(merged: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (counts, row-normalized fractions) indexed by MAF bin, columns = Group."""
    work = merged.copy()
    work["MAF_bin"] = work["MAF"].map(_assign_maf_bin)
    work = work.dropna(subset=["MAF_bin"])
    counts = (
        work.groupby(["MAF_bin", "Group"], observed=True)
        .size()
        .unstack(fill_value=0)
        .reindex(index=list(MAF_BIN_LABELS), fill_value=0)
    )
    fractions = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    return counts, fractions


MAF_MISS_RANGE_SPECS = (
    ("0-0.01", 0.0, 0.01),
    ("0-0.05", 0.0, 0.05),
    ("0-0.5", 0.0, 0.5),
)


def _mask_maf_upper(maf: pd.Series, upper: float) -> pd.Series:
    if upper >= 0.5:
        return (maf >= 0.0) & (maf <= upper)
    return (maf >= 0.0) & (maf < upper)


def _subsample_rows(df: pd.DataFrame, max_points: int, random_seed: int) -> pd.DataFrame:
    if len(df) <= max_points:
        return df
    return df.sample(n=max_points, random_state=np.random.default_rng(random_seed)).reset_index(drop=True)


def _plot_sample_group_maf_miss_regression(
    merged: pd.DataFrame,
    output_prefix: str,
    sample_group: str,
    regression_max_points: int,
    random_seed: int,
) -> None:
    """
    Scatter + OLS/Huber regression: cohort MAF (x) vs group-specific variant F_MISS (y).
    """
    if merged.empty:
        print(f"[Warning] No variants for group {sample_group} MAF vs miss regression.")
        return

    group_slug = sample_group.replace(" ", "_")
    stats_rows: list[dict] = []

    for tag, _lo, hi in MAF_MISS_RANGE_SPECS:
        sub = merged.loc[_mask_maf_upper(merged["MAF"], hi), ["MAF", "F_MISS"]].dropna()
        if len(sub) < 10:
            print(f"[Warning] Too few variants for group {sample_group} MAF {tag} regression.")
            continue

        reg_plot = _subsample_rows(sub, regression_max_points, random_seed)
        x_hi = 0.5 if hi >= 0.5 else hi
        pearson_r = float(sub["MAF"].corr(sub["F_MISS"])) if len(sub) > 1 else np.nan
        stats_rows.append(
            {
                "sample_group": sample_group,
                "maf_range": tag,
                "n_variants": int(len(sub)),
                "regression_points": int(len(reg_plot)),
                "pearson_r": pearson_r,
            }
        )

        plot_regression_comparison(
            reg_plot,
            x_col="MAF",
            y_col="F_MISS",
            x_label="Minor Allele Frequency (MAF, all samples)",
            y_label=f"Variant Missing Rate (sample_group={sample_group})",
            filename=f"{output_prefix}.group_{group_slug}.reg.maf{tag}.png",
            title=(
                f"Variant MAF vs missing (sample_group={sample_group}, MAF {tag}, n={len(sub)})"
            ),
            x_lim=(0.0, x_hi),
            y_lim=(0.0, 1.0),
        )

    if stats_rows:
        save_df_to_tsv(
            pd.DataFrame(stats_rows),
            f"{output_prefix}.group_{group_slug}.reg_by_maf.info.tsv",
        )


def _subsample_for_plot(df: pd.DataFrame, max_points: int, random_seed: int) -> pd.DataFrame:
    if len(df) <= max_points:
        return df
    rng = np.random.default_rng(random_seed)
    parts = []
    for group, sub in df.groupby("Group", sort=False):
        n = max(1, int(round(max_points * len(sub) / len(df))))
        n = min(n, len(sub))
        parts.append(sub.sample(n=n, random_state=rng))
    out = pd.concat(parts, ignore_index=True)
    if len(out) > max_points:
        out = out.sample(n=max_points, random_state=rng).reset_index(drop=True)
    return out


def ana_variant_maf_miss_by_group(
    afreq_path: str,
    vmiss_path: str,
    pfile_prefix: str,
    group_file: str,
    output_prefix: str,
    regression_max_points: int = 50000,
    random_seed: int = 1,
) -> None:
    """
    Variant MAF / F_MISS plots with sample_group coloring.

    Each variant is tagged with the sample_group that has the highest within-group
    missing rate at that site (dominant missing group). Writes stacked histograms
    and a MAF vs F_MISS joint regression (MAF 0-0.01, F_MISS 0-1).
    """
    if not Path(group_file).exists():
        print(f"[Warning] Group file missing; skipping group-stratified variant plots: {group_file}")
        return
    if not Path(f"{pfile_prefix}.pgen").exists():
        print(f"[Warning] Pfile missing; skipping group-stratified variant plots: {pfile_prefix}")
        return

    afreq = _maf_from_afreq_df(load_df_from_plink_variant(afreq_path))
    vmiss = load_df_from_plink_variant(vmiss_path)
    if afreq is None or vmiss is None or afreq.empty or vmiss.empty:
        print("[Warning] Empty afreq or vmiss; skipping group-stratified variant plots.")
        return

    group_long = _group_stratified_vmiss(pfile_prefix, group_file)
    dominant = _assign_dominant_missing_group(group_long)

    merged = (
        afreq[["ID", "MAF"]]
        .merge(vmiss[["ID", "F_MISS"]], on="ID", how="inner")
        .merge(dominant, on="ID", how="inner")
    )
    merged = merged.replace([np.inf, -np.inf], np.nan).dropna(subset=["MAF", "F_MISS", "Group"])
    if merged.empty:
        print("[Warning] No overlapping variants after group merge; skipping plots.")
        return

    save_df_to_tsv(
        merged[["ID", "MAF", "F_MISS", "Group", "Group_max_F_MISS"]],
        f"{output_prefix}.group.info.tsv",
    )

    mean_maf = float(merged["MAF"].mean())
    median_maf = float(merged["MAF"].median())
    mean_miss = float(merged["F_MISS"].mean())
    median_miss = float(merged["F_MISS"].median())

    maf_sub = merged[merged["MAF"] <= 0.01].copy()
    thresh_lines = [
        {"value": 0.001, "label": "Threshold 0.001", "color": "#55a868", "linestyle": "--"},
        {"value": 0.005, "label": "Threshold 0.005", "color": "#c44e52", "linestyle": "--"},
    ]

    plot_stacked_distribution(
        df=maf_sub,
        col="MAF",
        group_col="Group",
        title="Distribution of Minor Allele Frequency (0-0.01) by sample group",
        filename=f"{output_prefix}.maf.dist.0.01.linear.group.png",
        mean_val=mean_maf,
        median_val=median_maf,
        x_label="MAF",
        y_label="Count of Variants",
        bins=100,
        xlim=(0.0, 0.01),
    )
    plot_stacked_distribution(
        df=maf_sub,
        col="MAF",
        group_col="Group",
        title="Distribution of Minor Allele Frequency (0-0.01, log y) by sample group",
        filename=f"{output_prefix}.maf.dist.0.01.log.group.png",
        mean_val=mean_maf,
        median_val=median_maf,
        x_label="MAF",
        y_label="Count of Variants",
        bins=100,
        log_scale=True,
        xlim=(0.0, 0.01),
    )
    plot_stacked_distribution(
        df=merged,
        col="F_MISS",
        group_col="Group",
        title="Distribution of Variant Missing Rate by sample group",
        filename=f"{output_prefix}.miss.dist.group.png",
        mean_val=mean_miss,
        median_val=median_miss,
        x_label="Missing Rate (F_MISS)",
        y_label="Count of Variants",
        bins=50,
        xlim=(0.0, 1.0),
    )

    reg_df = merged[(merged["MAF"] <= 0.01)].copy()
    reg_plot = _subsample_for_plot(reg_df, regression_max_points, random_seed)
    plot_joint_regression(
        df=reg_plot,
        x_col="MAF",
        y_col="F_MISS",
        group_col="Group",
        x_label="Minor Allele Frequency (MAF)",
        y_label="Variant Missing Rate (F_MISS)",
        filename=f"{output_prefix}.reg.0.01.group.png",
        title="Variant MAF vs Missing Rate (MAF 0-0.01) by sample group",
        x_lim=(0.0, 0.01),
        y_lim=(0.0, 1.0),
    )

    counts, fractions = _maf_bin_group_fraction_table(merged)
    long_stats = []
    for maf_bin in MAF_BIN_LABELS:
        n_bin = int(counts.loc[maf_bin].sum()) if maf_bin in counts.index else 0
        for group in counts.columns:
            n = int(counts.loc[maf_bin, group]) if group in counts.columns else 0
            frac = float(fractions.loc[maf_bin, group]) if group in fractions.columns else 0.0
            long_stats.append(
                {
                    "MAF_bin": maf_bin,
                    "Group": group,
                    "n_variants": n,
                    "fraction": frac,
                    "n_variants_in_bin": n_bin,
                }
            )
    save_df_to_tsv(pd.DataFrame(long_stats), f"{output_prefix}.maf_bin.group_fraction.info.tsv")

    plot_stacked_proportion_bar(
        fractions,
        title="Dominant missing sample group fraction by MAF bin",
        filename=f"{output_prefix}.maf_bin.group_fraction.bar.png",
        x_label="Minor allele frequency (MAF) bin",
        y_label="Fraction of variants",
        stack_col_title="Group",
    )

    d_long = group_long.loc[group_long["Group"] == "D", ["ID", "F_MISS_group"]].rename(
        columns={"F_MISS_group": "F_MISS"}
    )
    if d_long.empty:
        print("[Warning] No group D vmiss rows; skipping D-only miss-by-MAF plots.")
    else:
        merged_d = (
            afreq[["ID", "MAF"]]
            .merge(d_long, on="ID", how="inner")
            .replace([np.inf, -np.inf], np.nan)
            .dropna(subset=["MAF", "F_MISS"])
        )
        save_df_to_tsv(merged_d, f"{output_prefix}.group_D.maf_miss.info.tsv")
        _plot_sample_group_maf_miss_regression(
            merged_d,
            output_prefix,
            sample_group="D",
            regression_max_points=regression_max_points,
            random_seed=random_seed,
        )
