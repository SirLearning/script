"""RefSeq v1.0 resource inventory (FASTA under home_dir; annotation may be absolute)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from infra.utils.io import save_df_to_tsv, save_thresholds

# Relative to vmap4 home_dir.
REF_V1_RELATIVE_ASSETS = {
    "fasta_A": "00data/03ref/01A/a_iwgscV1.fa.gz",
    "fasta_B": "00data/03ref/05B/b_iwgscV1.fa.gz",
    "fasta_D": "00data/03ref/04D/d_iwgscV1.fa.gz",
    "fasta_AB": "00data/03ref/02AB/ab_iwgscV1.fa.gz",
    "fasta_ABD": "00data/03ref/03ABD/abd_iwgscV1.fa.gz",
}

REQUIRED_ASSETS = frozenset({"fasta_A", "fasta_B", "fasta_D", "fasta_ABD"})


def inventory_ref_v1_assets(
    home_dir: str | Path,
    *,
    extra_assets: dict[str, str | Path] | None = None,
) -> pd.DataFrame:
    root = Path(home_dir)
    rows = []
    for name, rel in REF_V1_RELATIVE_ASSETS.items():
        path = root / rel
        rows.append(
            {
                "asset": name,
                "relative_path": rel,
                "path": str(path),
                "required": name in REQUIRED_ASSETS,
                "present": path.is_file(),
            }
        )
    for name, raw in (extra_assets or {}).items():
        path = Path(raw)
        rows.append(
            {
                "asset": name,
                "relative_path": "",
                "path": str(path),
                "required": False,
                "present": path.is_file(),
            }
        )
    return pd.DataFrame(rows)


def run_ref_v1_inventory(
    home_dir: str | Path,
    output_prefix: str | Path,
    *,
    extra_assets: dict[str, str | Path] | None = None,
) -> pd.DataFrame:
    table = inventory_ref_v1_assets(home_dir, extra_assets=extra_assets)
    save_df_to_tsv(table, f"{output_prefix}.ref_v1.tsv")
    missing_required = table[table["required"] & ~table["present"]]
    annotation = table[table["relative_path"] == ""]
    save_thresholds(
        {
            "n_assets": int(len(table)),
            "n_present": int(table["present"].sum()),
            "n_missing_required": int(len(missing_required)),
            "n_annotation_present": int(annotation["present"].sum()) if not annotation.empty else 0,
            "n_annotation_missing": int((~annotation["present"]).sum()) if not annotation.empty else 0,
        },
        f"{output_prefix}.ref_v1.th.tsv",
    )
    if not missing_required.empty:
        print("[Warning] Missing required RefSeq v1.0 FASTA assets:")
        for row in missing_required.itertuples(index=False):
            print(f"  {row.asset}: {row.path}")
    missing_ann = annotation[~annotation["present"]] if not annotation.empty else pd.DataFrame()
    if not missing_ann.empty:
        print("[Info] Annotation freeze targets not on disk:")
        for row in missing_ann.itertuples(index=False):
            print(f"  {row.asset}: {row.path}")
    return table
