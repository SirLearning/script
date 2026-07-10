"""Per-sample-group MAF with a globally aligned minor allele; intersection rare-site counts."""

from __future__ import annotations

import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from genetics.genomics.variant.variant_utils import load_df_from_plink_variant
from infra.utils.io import load_df_from_space_sep_no_header, save_df_to_tsv


def load_sample_group_members(
    group_file: str,
    *,
    panel_samples: set[str] | None = None,
    min_samples: int = 1,
) -> dict[str, list[str]]:
    """Return {group_name: [sample_id, ...]} from a two-column Sample Group file."""

    group_df = load_df_from_space_sep_no_header(group_file, ["Sample", "Group"])
    if group_df is None or group_df.empty:
        raise ValueError(f"Could not load group file: {group_file}")

    group_df = group_df.copy()
    group_df["Sample"] = group_df["Sample"].astype(str)
    group_df["Group"] = group_df["Group"].astype(str)
    if panel_samples is not None:
        group_df = group_df[group_df["Sample"].isin(panel_samples)]

    members: dict[str, list[str]] = {}
    for group, sub in group_df.groupby("Group", sort=True):
        samples = sub["Sample"].astype(str).tolist()
        if len(samples) < min_samples:
            print(f"[Warning] Skipping group {group}: only {len(samples)} samples in panel.")
            continue
        members[str(group)] = samples
    if not members:
        raise ValueError("No sample groups with enough members.")
    return members


def _write_keep_file(samples: list[str], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for sample in samples:
            handle.write(f"0\t{sample}\n")


def run_plink2_freq(
    pfile_prefix: str,
    out_prefix: Path,
    *,
    keep_file: Path | None = None,
    threads: int = 8,
) -> Path:
    """Run ``plink2 --freq`` and return the ``.afreq`` path."""

    cmd = [
        "plink2",
        "--pfile",
        pfile_prefix,
        "--freq",
        "--allow-extra-chr",
        "--threads",
        str(threads),
        "--out",
        str(out_prefix),
    ]
    if keep_file is not None:
        cmd.extend(["--keep", str(keep_file)])
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    afreq_path = Path(f"{out_prefix}.afreq")
    if not afreq_path.exists():
        raise FileNotFoundError(f"plink2 --freq did not write {afreq_path}")
    return afreq_path


def _run_group_plink2_freq_job(
    pfile_prefix: str,
    group: str,
    samples: list[str],
    work_dir: str,
    plink_threads: int,
) -> tuple[str, str]:
    """Worker for parallel per-group ``plink2 --freq`` (must be module-level for pickling)."""

    work = Path(work_dir)
    keep_path = work / f"{group}.keep"
    out_prefix = work / group
    _write_keep_file(samples, keep_path)
    afreq_path = run_plink2_freq(
        pfile_prefix,
        out_prefix,
        keep_file=keep_path,
        threads=plink_threads,
    )
    return group, str(afreq_path)


def _load_afreq_light(path: str | Path) -> pd.DataFrame:
    """Load only ID + ALT_FREQS from a plink2 ``.afreq`` file (no Position mapping)."""

    df = pd.read_csv(path, sep="\t")
    df.columns = [str(c).lstrip("#") for c in df.columns]
    if "ID" not in df.columns or "ALT_FREQS" not in df.columns:
        raise ValueError(f"afreq missing ID/ALT_FREQS columns: {path}")
    out = df[["ID", "ALT_FREQS"]].copy()
    out["ID"] = out["ID"].astype(str)
    if out.empty:
        raise ValueError(f"Empty afreq: {path}")
    out["ALT_FREQS"] = pd.to_numeric(out["ALT_FREQS"], errors="coerce")
    return out.dropna(subset=["ID", "ALT_FREQS"])


def _published_afreq_paths(
    published_afreq_dir: Path,
    chr_id: str,
    groups: dict[str, list[str]],
) -> tuple[Path, dict[str, Path]]:
    """Return expected published afreq paths: ``{chr_id}.global.afreq``, ``{chr_id}.{group}.afreq``."""

    base = published_afreq_dir / chr_id
    group_paths = {group: base / f"{chr_id}.{group}.afreq" for group in groups}
    return base / f"{chr_id}.global.afreq", group_paths


def _afreq_panel_complete(global_path: Path, group_paths: dict[str, Path]) -> bool:
    if not global_path.is_file() or global_path.stat().st_size == 0:
        return False
    return all(path.is_file() and path.stat().st_size > 0 for path in group_paths.values())


def _scratch_afreq_paths(work_dir: Path, groups: dict[str, list[str]]) -> tuple[Path, dict[str, Path]]:
    group_paths = {group: work_dir / f"{group}.afreq" for group in groups}
    return work_dir / "global.afreq", group_paths


def publish_afreq_panel(
    chr_id: str,
    global_path: Path,
    group_paths: dict[str, Path],
    published_afreq_dir: Path,
) -> tuple[Path, dict[str, Path]]:
    """Copy scratch afreq files to ``published_afreq_dir/{chr_id}/{chr_id}.{group}.afreq``."""

    pub_global, pub_groups = _published_afreq_paths(published_afreq_dir, chr_id, group_paths)
    pub_global.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(global_path, pub_global)
    for group, src in group_paths.items():
        shutil.copy2(src, pub_groups[group])
    print(f"[Info] Published afreq panel for {chr_id} -> {pub_global.parent}", flush=True)
    return pub_global, pub_groups


def migrate_plink_freq_scratch_dir(
    scratch_dir: Path,
    chr_id: str,
    published_afreq_dir: Path,
) -> bool:
    """One-off: copy ``plink_freq_scratch/{global,GROUP}.afreq`` into the published layout."""

    if not scratch_dir.is_dir():
        return False
    global_src = scratch_dir / "global.afreq"
    if not global_src.is_file():
        return False
    group_srcs = {
        path.stem: path
        for path in sorted(scratch_dir.glob("*.afreq"))
        if path.stem != "global"
    }
    if not group_srcs:
        return False
    publish_afreq_panel(chr_id, global_src, group_srcs, published_afreq_dir)
    return True


def _load_afreq(path: str | Path) -> pd.DataFrame:
    df = load_df_from_plink_variant(str(path))
    if df is None or df.empty:
        raise ValueError(f"Empty afreq: {path}")
    if "ALT_FREQS" not in df.columns:
        raise ValueError(f"afreq missing ALT_FREQS column: {path}")
    out = df.copy()
    out["ALT_FREQS"] = pd.to_numeric(out["ALT_FREQS"], errors="coerce")
    out["OBS_CT"] = pd.to_numeric(out.get("OBS_CT", pd.NA), errors="coerce")
    return out.dropna(subset=["ID", "ALT_FREQS"])


def global_minor_is_alt(global_afreq: pd.DataFrame) -> pd.Series:
    """
    For each variant, True when ALT is the globally minor allele (ALT_FREQS <= 0.5).
    """

    return global_afreq.set_index("ID")["ALT_FREQS"].le(0.5)


def group_maf_aligned_to_global_minor(
    group_afreq: pd.DataFrame,
    minor_is_alt: pd.Series,
) -> pd.Series:
    """
    MAF of the globally minor allele within this group's allele frequencies.

    When the global minor allele is ALT, use ALT_FREQS; otherwise use 1 - ALT_FREQS.
    """

    aligned = group_afreq.set_index("ID")
    common_ids = aligned.index.intersection(minor_is_alt.index)
    alt_freq = aligned.loc[common_ids, "ALT_FREQS"]
    use_alt = minor_is_alt.loc[common_ids]
    return alt_freq.where(use_alt, 1.0 - alt_freq).rename("group_maf")


def count_intersection_rare_sites(
    global_afreq: pd.DataFrame,
    group_afreqs: dict[str, pd.DataFrame],
    *,
    maf_threshold: float = 0.05,
) -> tuple[int, pd.DataFrame]:
    """
    Count variants whose aligned group MAF is strictly below ``maf_threshold`` in every group.

    Returns (count, per-variant table with one column per group MAF plus ``all_groups_rare``).
    """

    minor_is_alt = global_minor_is_alt(global_afreq)
    tables: list[pd.DataFrame] = []
    for group, afreq in group_afreqs.items():
        maf = group_maf_aligned_to_global_minor(afreq, minor_is_alt)
        tables.append(maf.rename(group))

    if not tables:
        return 0, pd.DataFrame(columns=["ID", "global_minor_is_alt"])

    merged = pd.concat(tables, axis=1)
    merged["global_minor_is_alt"] = minor_is_alt.reindex(merged.index)
    merged = merged.dropna(subset=list(group_afreqs.keys()))
    if merged.empty:
        empty = merged.reset_index(names="ID")
        empty["all_groups_rare"] = pd.Series(dtype=bool)
        return 0, empty

    group_cols = list(group_afreqs.keys())
    merged["all_groups_rare"] = (merged[group_cols] < maf_threshold).all(axis=1)
    count = int(merged["all_groups_rare"].sum())
    return count, merged.reset_index(names="ID")


def _panel_members_from_pfile(
    pfile_prefix: str,
    group_file: str,
    *,
    min_group_samples: int = 1,
) -> dict[str, list[str]]:
    """Resolve sample_group members present in the chromosome pfile panel."""

    prefix_path = Path(f"{pfile_prefix}.pgen")
    if not prefix_path.exists():
        raise FileNotFoundError(f"Missing pgen for chromosome panel: {prefix_path}")

    psam_path = Path(f"{pfile_prefix}.psam")
    psam = pd.read_csv(psam_path, sep="\t")
    psam.columns = [c.lstrip("#") for c in psam.columns]
    iid_col = "IID" if "IID" in psam.columns else psam.columns[0]
    panel_samples = set(psam[iid_col].astype(str))

    return load_sample_group_members(
        group_file,
        panel_samples=panel_samples,
        min_samples=min_group_samples,
    )


def _emit_named_afreq_copies(
    chr_id: str,
    global_path: Path,
    group_paths: dict[str, Path],
    emit_dir: Path,
) -> None:
    """Write ``{chr_id}.global.afreq`` / ``{chr_id}.{group}.afreq`` for NF publishDir."""

    emit_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(global_path, emit_dir / f"{chr_id}.global.afreq")
    for group, path in group_paths.items():
        shutil.copy2(path, emit_dir / f"{chr_id}.{group}.afreq")


def _run_group_plink_freq_jobs(
    pfile_prefix: str,
    group_file: str,
    work_dir: Path,
    *,
    chr_id: str | None = None,
    published_afreq_dir: Path | None = None,
    afreq_emit_dir: Path | None = None,
    min_group_samples: int = 1,
    threads: int = 8,
    group_workers: int = 4,
    plink_threads: int | None = None,
) -> tuple[Path, dict[str, Path]]:
    """Run or reuse global and per-sample-group ``plink2 --freq``; return ``.afreq`` paths."""

    if plink_threads is None:
        plink_threads = max(2, threads // max(1, group_workers))

    members = _panel_members_from_pfile(
        pfile_prefix,
        group_file,
        min_group_samples=min_group_samples,
    )
    work_dir.mkdir(parents=True, exist_ok=True)

    if published_afreq_dir is not None and chr_id is not None:
        pub_global, pub_groups = _published_afreq_paths(published_afreq_dir, chr_id, members)
        if _afreq_panel_complete(pub_global, pub_groups):
            print(
                f"[Info] Reusing published afreq panel for {chr_id} "
                f"({pub_global.parent})",
                flush=True,
            )
            if afreq_emit_dir is not None:
                _emit_named_afreq_copies(chr_id, pub_global, pub_groups, afreq_emit_dir)
            return pub_global, pub_groups

    print(
        f"[Info] Group freq panel: global + {len(members)} groups "
        f"(group_workers={group_workers}, plink_threads={plink_threads})",
        flush=True,
    )
    global_afreq_path = run_plink2_freq(
        pfile_prefix, work_dir / "global", threads=plink_threads
    )

    group_afreq_paths: dict[str, Path] = {}
    worker_count = min(max(1, group_workers), len(members))
    if worker_count == 1:
        for group, samples in members.items():
            grp, afreq_path = _run_group_plink2_freq_job(
                pfile_prefix, group, samples, str(work_dir), plink_threads
            )
            group_afreq_paths[grp] = Path(afreq_path)
            print(f"[Info] Finished group freq: {grp}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=worker_count) as pool:
            futures = {
                pool.submit(
                    _run_group_plink2_freq_job,
                    pfile_prefix,
                    group,
                    samples,
                    str(work_dir),
                    plink_threads,
                ): group
                for group, samples in members.items()
            }
            for fut in as_completed(futures):
                group, afreq_path = fut.result()
                group_afreq_paths[group] = Path(afreq_path)
                print(f"[Info] Finished group freq: {group}", flush=True)

    if published_afreq_dir is not None and chr_id is not None:
        pub_global, pub_groups = publish_afreq_panel(
            chr_id, global_afreq_path, group_afreq_paths, published_afreq_dir
        )
        global_afreq_path = pub_global
        group_afreq_paths = pub_groups

    if afreq_emit_dir is not None and chr_id is not None:
        _emit_named_afreq_copies(chr_id, global_afreq_path, group_afreq_paths, afreq_emit_dir)

    return global_afreq_path, group_afreq_paths


def _load_group_afreqs_sequential(group_afreq_paths: dict[str, Path]) -> dict[str, pd.DataFrame]:
    """Load per-group afreq tables one at a time (lower peak RAM than eager load)."""

    group_afreqs: dict[str, pd.DataFrame] = {}
    for group, path in group_afreq_paths.items():
        group_afreqs[group] = _load_afreq(path)
    return group_afreqs


def _run_group_freq_panel(
    pfile_prefix: str,
    group_file: str,
    work_dir: Path,
    *,
    min_group_samples: int = 1,
    threads: int = 8,
    group_workers: int = 4,
    plink_threads: int | None = None,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    """Run freq jobs and load all afreq tables (for detail TSV export)."""

    global_afreq_path, group_afreq_paths = _run_group_plink_freq_jobs(
        pfile_prefix,
        group_file,
        work_dir,
        min_group_samples=min_group_samples,
        threads=threads,
        group_workers=group_workers,
        plink_threads=plink_threads,
    )
    global_afreq = _load_afreq(global_afreq_path)
    group_afreqs = _load_group_afreqs_sequential(group_afreq_paths)
    return global_afreq, group_afreqs


def rare_inter_mask_from_freq_paths(
    global_afreq_path: Path,
    group_afreq_paths: dict[str, Path],
    *,
    maf_threshold: float = 0.05,
) -> tuple[int, int, pd.Series]:
    """
    Compute the all-group-rare mask loading one group afreq at a time.

    Returns (n_variants, n_all_groups_rare, boolean mask indexed by variant ID).
    """

    global_afreq = _load_afreq_light(global_afreq_path)
    minor_is_alt = global_minor_is_alt(global_afreq)
    n_variants = int(len(global_afreq))
    del global_afreq

    rare_mask = pd.Series(True, index=minor_is_alt.index)
    for group, path in group_afreq_paths.items():
        afreq = _load_afreq_light(path)
        maf = group_maf_aligned_to_global_minor(afreq, minor_is_alt)
        aligned = maf.reindex(rare_mask.index)
        rare_mask &= aligned.notna() & (aligned < maf_threshold)
        del afreq, maf
        print(f"[Info] Applied group MAF filter: {group}", flush=True)

    rare_count = int(rare_mask.sum())
    return n_variants, rare_count, rare_mask


def write_rare_inter_extract_for_chromosome(
    pfile_prefix: str,
    group_file: str,
    extract_path: str,
    *,
    chr_id: str | None = None,
    maf_threshold: float = 0.05,
    summary_path: str | None = None,
    published_afreq_dir: str | Path | None = None,
    published_extract_path: str | Path | None = None,
    published_summary_path: str | Path | None = None,
    afreq_emit_dir: str | Path | None = None,
    min_group_samples: int = 1,
    threads: int = 8,
    group_workers: int = 4,
    plink_threads: int | None = None,
    work_dir: Path | None = None,
) -> dict[str, int | float | str]:
    """
    Write variant IDs rare (MAF < threshold) in every sample_group on this chromosome.

    Uses the globally aligned minor allele (see :func:`count_intersection_rare_sites`).
    Skips plink2 when a complete published afreq panel exists; skips entirely when
    published rare-inter extract + summary already exist.
    """

    if chr_id is None:
        chr_id = Path(extract_path).stem.replace(".rare_inter", "")

    pub_afreq_dir = Path(published_afreq_dir) if published_afreq_dir else None
    emit_dir = Path(afreq_emit_dir) if afreq_emit_dir else None

    if published_extract_path and published_summary_path:
        pub_extract = Path(published_extract_path)
        pub_summary = Path(published_summary_path)
        if pub_extract.is_file() and pub_summary.is_file():
            out_path = Path(extract_path)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(pub_extract, out_path)
            if summary_path:
                shutil.copy2(pub_summary, summary_path)
            if pub_afreq_dir is not None:
                members = _panel_members_from_pfile(
                    pfile_prefix,
                    group_file,
                    min_group_samples=min_group_samples,
                )
                pub_global, pub_groups = _published_afreq_paths(pub_afreq_dir, chr_id, members)
                if _afreq_panel_complete(pub_global, pub_groups) and emit_dir is not None:
                    _emit_named_afreq_copies(chr_id, pub_global, pub_groups, emit_dir)
            summary_df = pd.read_csv(pub_summary, sep="\t")
            summary = summary_df.iloc[0].to_dict()
            print(
                f"[Info] Skipping {chr_id}: published rare-inter extract already exists",
                flush=True,
            )
            return summary

    if work_dir is None:
        work_dir = Path(extract_path).parent / "plink_freq_scratch"
    global_afreq_path, group_afreq_paths = _run_group_plink_freq_jobs(
        pfile_prefix,
        group_file,
        work_dir,
        chr_id=chr_id,
        published_afreq_dir=pub_afreq_dir,
        afreq_emit_dir=emit_dir,
        min_group_samples=min_group_samples,
        threads=threads,
        group_workers=group_workers,
        plink_threads=plink_threads,
    )
    n_variants, rare_count, rare_mask = rare_inter_mask_from_freq_paths(
        global_afreq_path,
        group_afreq_paths,
        maf_threshold=maf_threshold,
    )
    rare_ids = rare_mask[rare_mask].index.astype(str)
    out_path = Path(extract_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(rare_ids) + ("\n" if len(rare_ids) else ""), encoding="utf-8")

    summary = {
        "n_variants": int(n_variants),
        "n_groups": int(len(group_afreq_paths)),
        "groups": ",".join(sorted(group_afreq_paths.keys())),
        "maf_threshold": float(maf_threshold),
        "n_all_groups_rare": int(rare_count),
        "fraction_all_groups_rare": float(rare_count / n_variants) if n_variants else 0.0,
    }
    if summary_path:
        save_df_to_tsv(pd.DataFrame([summary]), summary_path)
    return summary


def analyze_chromosome_group_maf_intersection(
    pfile_prefix: str,
    group_file: str,
    output_prefix: str,
    *,
    maf_threshold: float = 0.05,
    min_group_samples: int = 1,
    threads: int = 8,
    group_workers: int = 4,
    plink_threads: int | None = None,
    work_dir: Path | None = None,
) -> dict[str, int | float]:
    """
    Per-chromosome freq for global + each sample_group; write detail TSV and return summary stats.
    """

    prefix_path = Path(f"{pfile_prefix}.pgen")
    if not prefix_path.exists():
        raise FileNotFoundError(f"Missing pgen for chromosome panel: {prefix_path}")

    if work_dir is None:
        work_dir = Path(output_prefix).parent / "plink_freq_scratch"
    global_afreq, group_afreqs = _run_group_freq_panel(
        pfile_prefix,
        group_file,
        work_dir,
        min_group_samples=min_group_samples,
        threads=threads,
        group_workers=group_workers,
        plink_threads=plink_threads,
    )

    rare_count, detail = count_intersection_rare_sites(
        global_afreq,
        group_afreqs,
        maf_threshold=maf_threshold,
    )
    save_df_to_tsv(detail, f"{output_prefix}.group_maf.detail.tsv")

    summary = {
        "n_variants": int(len(global_afreq)),
        "n_groups": int(len(group_afreqs)),
        "groups": ",".join(sorted(group_afreqs.keys())),
        "maf_threshold": float(maf_threshold),
        "n_all_groups_rare": int(rare_count),
        "fraction_all_groups_rare": float(rare_count / len(global_afreq)) if len(global_afreq) else 0.0,
    }
    save_df_to_tsv(pd.DataFrame([summary]), f"{output_prefix}.group_maf.summary.tsv")
    return summary


def aggregate_group_maf_summaries(summary_paths: list[str | Path]) -> pd.DataFrame:
    """Sum per-chromosome summaries into a genome-wide table."""

    rows: list[dict] = []
    for path in summary_paths:
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            continue
        rows.append(df.iloc[0].to_dict())
    if not rows:
        return pd.DataFrame()

    agg = pd.DataFrame(rows)
    total_variants = int(agg["n_variants"].sum())
    total_rare = int(agg["n_all_groups_rare"].sum())
    out = {
        "n_chromosomes": int(len(agg)),
        "n_groups": int(agg["n_groups"].iloc[0]) if "n_groups" in agg.columns else 0,
        "maf_threshold": float(agg["maf_threshold"].iloc[0]) if "maf_threshold" in agg.columns else 0.05,
        "n_variants": total_variants,
        "n_all_groups_rare": total_rare,
        "fraction_all_groups_rare": float(total_rare / total_variants) if total_variants else 0.0,
        "thin_rate_for_3M": float(3_000_000 / total_rare) if total_rare else float("nan"),
    }
    return pd.DataFrame([out])
