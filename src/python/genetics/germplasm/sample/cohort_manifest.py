"""Sample-role manifest: WGS-core (~5000) vs WAP association panel."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from genetics.genomics.variant.group_maf_intersection import write_sample_group_keep_file
from infra.utils.graph import plot_bar_chart, plot_stacked_histogram
from infra.utils.io import load_df_from_space_sep_no_header, save_df_to_tsv, save_thresholds

WGS_CORE_GROUPS = frozenset({"A", "AB", "ABD", "WAP", "HZNU", "Nature"})
ASSOC_GROUPS = frozenset({"WAP"})

ROLE_ASSOC = "assoc_pheno"
ROLE_WGS_CORE = "wgs_core_only"
ROLE_WGS_OTHER = "wgs_other"
ROLE_ORDER = (ROLE_ASSOC, ROLE_WGS_CORE, ROLE_WGS_OTHER)
ROLE_COLORS = {
    ROLE_ASSOC: "#4c78a8",
    ROLE_WGS_CORE: "#54a24b",
    ROLE_WGS_OTHER: "#e45756",
}

KEEP_ASSOC = "assoc.keep"
KEEP_WGS_CORE = "wgs_core.keep"


def cohort_role(sample_group: str) -> str:
    group = str(sample_group)
    if group in ASSOC_GROUPS:
        return ROLE_ASSOC
    if group in WGS_CORE_GROUPS:
        return ROLE_WGS_CORE
    return ROLE_WGS_OTHER


def parse_group_list(raw: str | None, fallback: frozenset[str]) -> frozenset[str]:
    if raw is None or not str(raw).strip():
        return frozenset(fallback)
    return frozenset(part.strip() for part in str(raw).split(",") if part.strip())


def load_sample_group_table(group_file: str | Path) -> pd.DataFrame:
    group_df = load_df_from_space_sep_no_header(str(group_file), ["Sample", "Group"])
    if group_df is None or group_df.empty:
        raise ValueError(f"Could not load sample_group file: {group_file}")
    out = group_df.copy()
    out["Sample"] = out["Sample"].astype(str)
    out["Group"] = out["Group"].astype(str)
    return out


def load_taxa_coverage_table(taxa_bam_map: str | Path) -> pd.DataFrame:
    path = Path(taxa_bam_map)
    if not path.is_file():
        return pd.DataFrame(columns=["Sample", "Coverage"])
    raw = pd.read_csv(path, sep="\t", header=None, dtype=str, comment=None)
    if raw.empty:
        return pd.DataFrame(columns=["Sample", "Coverage"])
    first = str(raw.iloc[0, 0])
    if first in {"Taxa", "Sample"}:
        raw = raw.iloc[1:].reset_index(drop=True)
    out = pd.DataFrame(
        {
            "Sample": raw.iloc[:, 0].astype(str),
            "Coverage": pd.to_numeric(raw.iloc[:, 1], errors="coerce"),
        }
    )
    return out.drop_duplicates(subset=["Sample"], keep="first")


def _collect_taxa_maps(taxa_bam_map_dir: str | Path | None) -> pd.DataFrame:
    if taxa_bam_map_dir is None:
        return pd.DataFrame(columns=["Sample", "Coverage", "taxa_map"])
    root = Path(taxa_bam_map_dir)
    if not root.is_dir():
        return pd.DataFrame(columns=["Sample", "Coverage", "taxa_map"])
    frames: list[pd.DataFrame] = []
    for path in sorted(root.glob("*.taxaBamMap.txt")):
        table = load_taxa_coverage_table(path)
        if table.empty:
            continue
        table = table.copy()
        table["taxa_map"] = path.name
        frames.append(table)
    if not frames:
        return pd.DataFrame(columns=["Sample", "Coverage", "taxa_map"])
    return pd.concat(frames, ignore_index=True)


def build_sample_manifest(
    group_file: str | Path,
    *,
    taxa_bam_map_dir: str | Path | None = None,
    wgs_core_groups: frozenset[str] | None = None,
    assoc_groups: frozenset[str] | None = None,
) -> pd.DataFrame:
    """Join sample_group to optional taxaBamMap coverage and assign cohort roles."""

    wgs_groups = wgs_core_groups or WGS_CORE_GROUPS
    assoc = assoc_groups or ASSOC_GROUPS
    samples = load_sample_group_table(group_file)
    samples["in_assoc"] = samples["Group"].isin(assoc)
    samples["in_wgs_core"] = samples["Group"].isin(wgs_groups)
    samples["cohort_role"] = [
        ROLE_ASSOC if row.in_assoc else (ROLE_WGS_CORE if row.in_wgs_core else ROLE_WGS_OTHER)
        for row in samples.itertuples(index=False)
    ]

    maps = _collect_taxa_maps(taxa_bam_map_dir)
    if maps.empty:
        samples["Coverage"] = pd.NA
        samples["taxa_map"] = pd.NA
        samples["in_taxa_map"] = False
        return samples

    cov = maps.groupby("Sample", as_index=False).agg(
        Coverage=("Coverage", "max"),
        taxa_map=("taxa_map", "first"),
    )
    merged = samples.merge(cov, on="Sample", how="left")
    merged["in_taxa_map"] = merged["taxa_map"].notna()
    return merged


def write_cohort_keep_files(
    manifest: pd.DataFrame,
    output_dir: str | Path,
) -> dict[str, int]:
    """Write PLINK2 FID/IID keep files for assoc and WGS-core panels."""

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    counts: dict[str, int] = {}
    for flag, name in (("in_assoc", KEEP_ASSOC), ("in_wgs_core", KEEP_WGS_CORE)):
        samples = manifest.loc[manifest[flag], "Sample"].astype(str).tolist()
        path = out_dir / name
        with path.open("w", encoding="utf-8") as handle:
            for sample in samples:
                handle.write(f"0\t{sample}\n")
        counts[name] = len(samples)
    return counts


def summarize_cohort_qc(manifest: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for group, sub in manifest.groupby("Group", sort=True):
        rows.append(
            {
                "key": f"group:{group}",
                "n": int(len(sub)),
                "n_assoc": int(sub["in_assoc"].sum()),
                "n_wgs_core": int(sub["in_wgs_core"].sum()),
                "n_in_taxa_map": int(sub["in_taxa_map"].sum()),
                "coverage_median": float(sub["Coverage"].median()) if sub["Coverage"].notna().any() else pd.NA,
            }
        )
    for role in ROLE_ORDER:
        sub = manifest[manifest["cohort_role"] == role]
        rows.append(
            {
                "key": f"role:{role}",
                "n": int(len(sub)),
                "n_assoc": int(sub["in_assoc"].sum()),
                "n_wgs_core": int(sub["in_wgs_core"].sum()),
                "n_in_taxa_map": int(sub["in_taxa_map"].sum()),
                "coverage_median": float(sub["Coverage"].median()) if sub["Coverage"].notna().any() else pd.NA,
            }
        )
    rows.append(
        {
            "key": "total",
            "n": int(len(manifest)),
            "n_assoc": int(manifest["in_assoc"].sum()),
            "n_wgs_core": int(manifest["in_wgs_core"].sum()),
            "n_in_taxa_map": int(manifest["in_taxa_map"].sum()),
            "coverage_median": float(manifest["Coverage"].median()) if manifest["Coverage"].notna().any() else pd.NA,
        }
    )
    return pd.DataFrame(rows)


def write_cohort_qc_thresholds(summary: pd.DataFrame, output_path: str | Path) -> None:
    stats = {str(row.key): int(row.n) for row in summary.itertuples(index=False)}
    save_thresholds(stats, str(output_path))


def plot_cohort_qc(manifest: pd.DataFrame, output_prefix: str | Path) -> None:
    prefix = Path(output_prefix)
    counts = manifest.groupby("Group").size().sort_values(ascending=False)
    plot_bar_chart(
        names=list(counts.index.astype(str)),
        values=[float(v) for v in counts.to_list()],
        title="Samples by sample_group",
        ylabel="Sample count",
        filename=str(prefix) + ".group.counts.png",
        ylim=(0.0, float(counts.max()) * 1.15 if len(counts) else 1.05),
        rotate_xlabels=45,
    )
    if manifest["Coverage"].notna().any():
        plot_df = manifest.loc[manifest["Coverage"].notna(), ["Coverage", "cohort_role"]].copy()
        plot_stacked_histogram(
            data=plot_df,
            x_col="Coverage",
            hue_col="cohort_role",
            hue_order=[role for role in ROLE_ORDER if role in set(plot_df["cohort_role"])],
            palette_map=ROLE_COLORS,
            title="TaxaBamMap coverage by cohort role",
            filename=str(prefix) + ".coverage.role.png",
            xlabel="Coverage",
            ylabel="Sample count",
            legend_title="Cohort role",
        )


def write_sample_group_keeps_from_params(
    group_file: str | Path,
    output_dir: str | Path,
    *,
    wgs_core_groups: str | None = None,
    assoc_groups: str | None = None,
) -> int:
    """Compatibility helper used when only keep files are needed."""

    wgs = parse_group_list(wgs_core_groups, WGS_CORE_GROUPS)
    assoc = parse_group_list(assoc_groups, ASSOC_GROUPS)
    n_assoc = write_sample_group_keep_file(
        group_file, Path(output_dir) / KEEP_ASSOC, include_groups=assoc
    )
    n_wgs = write_sample_group_keep_file(
        group_file, Path(output_dir) / KEEP_WGS_CORE, include_groups=wgs
    )
    return n_assoc + n_wgs


def run_sample_manifest(
    group_file: str | Path,
    output_prefix: str | Path,
    *,
    taxa_bam_map_dir: str | Path | None = None,
    wgs_core_groups: str | None = None,
    assoc_groups: str | None = None,
) -> pd.DataFrame:
    """Build manifest, keep files, QC table, and plots next to ``output_prefix``."""

    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    manifest = build_sample_manifest(
        group_file,
        taxa_bam_map_dir=taxa_bam_map_dir,
        wgs_core_groups=parse_group_list(wgs_core_groups, WGS_CORE_GROUPS),
        assoc_groups=parse_group_list(assoc_groups, ASSOC_GROUPS),
    )
    save_df_to_tsv(manifest, f"{prefix}.manifest.tsv")
    summary = summarize_cohort_qc(manifest)
    save_df_to_tsv(summary, f"{prefix}.qc.tsv")
    write_cohort_qc_thresholds(summary, f"{prefix}.qc.th.tsv")
    write_cohort_keep_files(manifest, prefix.parent)
    plot_cohort_qc(manifest, prefix)
    print(
        f"sample_manifest: n={len(manifest)} assoc={int(manifest['in_assoc'].sum())} "
        f"wgs_core={int(manifest['in_wgs_core'].sum())}"
    )
    return manifest
