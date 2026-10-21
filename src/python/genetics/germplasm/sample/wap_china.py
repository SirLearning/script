"""WAP passport helpers: country subsets and PLINK2 keep files."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def resolve_r1_id_to_taxa(r1_id: str, available: set[str]) -> str | None:
    """Map a WAP ``r1_id`` (e.g. P0030) onto a genotype IID (P0030 or P0030A)."""
    token = str(r1_id).strip()
    if not token or token.lower() == "nan":
        return None
    if token in available:
        return token
    alias = f"{token}A"
    if alias in available:
        return alias
    return None


def china_wap_taxa(
    germplasm_tsv: str | Path,
    available_iids: set[str] | list[str],
    *,
    country: str = "China",
) -> list[str]:
    """Return genotype IIDs for WAP accessions with the given passport country."""
    table = pd.read_csv(germplasm_tsv, sep="\t", dtype=str)
    if "r1_id" not in table.columns or "country" not in table.columns:
        raise ValueError("germplasm table must have r1_id and country columns")
    want = country.strip().casefold()
    available = {str(sample) for sample in available_iids}
    taxa: list[str] = []
    seen: set[str] = set()
    for _, row in table.iterrows():
        if str(row["country"]).strip().casefold() != want:
            continue
        iid = resolve_r1_id_to_taxa(row["r1_id"], available)
        if iid is None or iid in seen:
            continue
        seen.add(iid)
        taxa.append(iid)
    return taxa


def write_plink2_keep(samples: list[str], path: str | Path) -> None:
    """Write FID/IID keep lines (FID=0) for PLINK2 ``--keep``."""
    if not samples:
        raise ValueError("no samples to write to keep file")
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as handle:
        for sample in samples:
            handle.write(f"0\t{sample}\n")


def load_psam_iids(psam_path: str | Path) -> set[str]:
    iids: set[str] = set()
    with Path(psam_path).open(encoding="utf-8") as handle:
        header = handle.readline()
        iid_col = 0
        if header.startswith("#"):
            cols = header.lstrip("#").split()
            if "IID" in cols:
                iid_col = cols.index("IID")
        for line in handle:
            parts = line.split()
            if len(parts) > iid_col:
                iids.add(parts[iid_col])
    return iids
