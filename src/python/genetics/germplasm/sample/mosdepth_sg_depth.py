"""Subgenome mean depth from mosdepth summary files (bases sum / length sum)."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

from genetics.germplasm.sample.process import load_df_from_tbm
from infra.utils.io import load_df_generic, save_df_to_tsv
from infra.wheat import ref_v1

DEFAULT_DEPTH_ROOT = "/data/home/tusr1/01projects/vmap4/00data/04depth"
DEFAULT_TBM_DIR = "/data/home/tusr1/01projects/vmap4/00data/05taxaBamMap"
DEFAULT_SG_DEPTH_DIR = f"{DEFAULT_DEPTH_ROOT}/sg_mean_depth"

SUBGENOMES = ("A", "B", "D", "Others")
TBM_TEMPLATE_BY_SUBGENOME = {
    "A": "all.A.taxaBamMap.txt",
    "B": "all.B.taxaBamMap.txt",
    "D": "all.D.taxaBamMap.txt",
    "Others": "all.ALL.taxaBamMap.txt",
}


def _summary_stem_variants(stem: str) -> set[str]:
    """Candidate lookup keys for a mosdepth summary filename stem."""
    variants = {stem}
    if stem.endswith(".bam"):
        variants.add(stem[:-4])
    if stem.endswith(".rmdup.bam"):
        variants.add(stem[: -len(".rmdup.bam")])
        variants.add(stem[: -len(".rmdup")])
    if stem.endswith(".rmdup"):
        variants.add(stem[: -len(".rmdup")])
    if stem.endswith("_deduped.bam"):
        variants.add(stem[: -len("_deduped.bam")])
    if stem.endswith("_deduped"):
        variants.add(stem[: -len("_deduped")])
    return variants


def build_mosdepth_summary_index(depth_root: str | Path) -> dict[str, Path]:
    """Index mosdepth summary paths by filename stem and common taxa aliases."""
    depth_root = Path(depth_root)
    index: dict[str, Path] = {}
    for path in depth_root.rglob("*mosdepth.summary.txt"):
        stem = path.name.replace(".mosdepth.summary.txt", "")
        for key in _summary_stem_variants(stem):
            index.setdefault(key, path)
    return index


def resolve_mosdepth_summary_path(
    sample: str,
    bam_paths: list[str] | None,
    index: dict[str, Path],
) -> Path | None:
    """Resolve a mosdepth summary file for a taxon using sample id and BAM paths."""
    candidates: list[str] = [str(sample)]
    for bam in bam_paths or []:
        base = os.path.basename(str(bam))
        candidates.append(base)
        if base.endswith(".bam"):
            candidates.append(base[:-4])
        if base.endswith(".rmdup.bam"):
            candidates.append(base[: -len(".rmdup.bam")])
        if base.endswith("_deduped.bam"):
            candidates.append(base[: -len("_deduped.bam")])
    seen: set[str] = set()
    for cand in candidates:
        if not cand or cand in seen:
            continue
        seen.add(cand)
        path = index.get(cand)
        if path is not None:
            return path
    return None


def compute_subgenome_mean_depth(summary_path: str | Path, subgenome: str) -> float:
    """
    Mean depth for one subgenome: sum(bases) / sum(length) over segment chrom ids.

    Matches vmap4 ``getRefV1SubChr`` / mosdepth ``chrom`` column (0–44 strings).
    """
    chr_ids = set(ref_v1.get_ref_v1_subgenome_segment_chr_ids(subgenome))
    df = load_df_generic(str(summary_path))
    if df is None or df.empty:
        return float("nan")
    df = df.copy()
    df["chrom"] = df["chrom"].astype(str)
    df = df[df["chrom"] != "total"]
    sub = df[df["chrom"].isin(chr_ids)]
    total_len = float(sub["length"].sum())
    if total_len <= 0:
        return float("nan")
    return float(sub["bases"].sum() / total_len)


def compute_all_subgenome_mean_depths(summary_path: str | Path) -> dict[str, float]:
    return {
        sub: compute_subgenome_mean_depth(summary_path, sub)
        for sub in SUBGENOMES
    }


def build_subgenome_depth_table(
    taxa_bam_path: str | Path,
    depth_root: str | Path = DEFAULT_DEPTH_ROOT,
    index: dict[str, Path] | None = None,
) -> pd.DataFrame:
    """Per-sample subgenome mean depths for all taxa in a taxaBamMap file."""
    if index is None:
        index = build_mosdepth_summary_index(depth_root)
    tbm = load_df_from_tbm(str(taxa_bam_path))
    if tbm is None or tbm.empty:
        return pd.DataFrame()

    rows = []
    for _, row in tbm.iterrows():
        sample = str(row["Sample"])
        bam_paths = row.get("Bam_Path") or []
        summary_path = resolve_mosdepth_summary_path(sample, bam_paths, index)
        entry = {
            "Sample": sample,
            "mosdepth_summary": str(summary_path) if summary_path else "",
        }
        if summary_path is None:
            for sub in SUBGENOMES:
                entry[f"{sub}_depth"] = np.nan
        else:
            depths = compute_all_subgenome_mean_depths(summary_path)
            for sub in SUBGENOMES:
                entry[f"{sub}_depth"] = depths[sub]
        rows.append(entry)
    return pd.DataFrame(rows)


def write_sg_depth_taxa_bam_map(
    source_taxa_bam_path: str | Path,
    depth_table: pd.DataFrame,
    subgenome: str,
    output_path: str | Path,
    depth_col: str | None = None,
) -> int:
    """
    Write taxaBamMap-compatible TSV with Coverage replaced by subgenome mean depth.

    Returns number of taxa written.
    """
    subgenome = str(subgenome)
    depth_col = depth_col or f"{subgenome}_depth"
    depth_by_sample = depth_table.set_index("Sample")[depth_col].to_dict()

    with open(source_taxa_bam_path, "r", encoding="utf-8") as src, open(
        output_path, "w", encoding="utf-8"
    ) as out:
        header = src.readline()
        out.write(header)
        n = 0
        for line in src:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            sample = parts[0]
            depth = depth_by_sample.get(sample, np.nan)
            if pd.isna(depth):
                cov_str = "0"
            else:
                cov_str = f"{float(depth):.4f}".rstrip("0").rstrip(".")
            parts[1] = cov_str
            out.write("\t".join(parts) + "\n")
            n += 1
    return n


def sg_depth_taxa_bam_map_path(
    home_dir: str | Path,
    subgenome: str,
    depth_root: str | Path | None = None,
) -> Path:
    depth_root = Path(depth_root or f"{home_dir}/00data/04depth")
    return depth_root / "sg_mean_depth" / f"all.{subgenome}.sg_depth.taxaBamMap.txt"


def build_all_sg_depth_taxa_bam_maps(
    home_dir: str | Path = "/data/home/tusr1/01projects/vmap4",
    depth_root: str | Path | None = None,
    tbm_dir: str | Path | None = None,
    output_dir: str | Path | None = None,
) -> dict[str, Path]:
    """
    Build per-subgenome taxaBamMap files with mosdepth subgenome mean depth as Coverage.

    Also writes ``all.sg_depth.wide.info.tsv`` (Sample + A/B/D/Others depths).
    """
    home_dir = Path(home_dir)
    depth_root = Path(depth_root or home_dir / "00data/04depth")
    tbm_dir = Path(tbm_dir or home_dir / "00data/05taxaBamMap")
    output_dir = Path(output_dir or depth_root / "sg_mean_depth")
    output_dir.mkdir(parents=True, exist_ok=True)

    index = build_mosdepth_summary_index(depth_root)
    print(f"[Info] mosdepth summary index: {len(index)} lookup keys under {depth_root}")

    all_tbm = tbm_dir / "all.ALL.taxaBamMap.txt"
    depth_table = build_subgenome_depth_table(all_tbm, depth_root=depth_root, index=index)
    save_df_to_tsv(depth_table, str(output_dir / "all.sg_depth.wide.info.tsv"))

    found = depth_table["mosdepth_summary"].astype(str).str.len().gt(0).sum()
    print(f"[Info] mosdepth matched {found}/{len(depth_table)} taxa from {all_tbm}")

    outputs: dict[str, Path] = {}
    for sub in SUBGENOMES:
        template = tbm_dir / TBM_TEMPLATE_BY_SUBGENOME[sub]
        out_path = output_dir / f"all.{sub}.sg_depth.taxaBamMap.txt"
        n = write_sg_depth_taxa_bam_map(template, depth_table, sub, out_path)
        outputs[sub] = out_path
        print(f"[Info] Wrote {n} rows -> {out_path}")
    return outputs


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Build subgenome mean-depth taxaBamMap files")
    parser.add_argument(
        "--home-dir",
        default="/data/home/tusr1/01projects/vmap4",
        help="vmap4 project root (00data lives here)",
    )
    args = parser.parse_args()
    build_all_sg_depth_taxa_bam_maps(home_dir=args.home_dir)
