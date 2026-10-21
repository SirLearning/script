"""SnpEff ANN VCF → VEP-shaped table so coding_mask can read either annotator."""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

from infra.utils.io import save_df_to_tsv

VEP_COLUMNS = (
    "Uploaded_variation",
    "Location",
    "Allele",
    "Gene",
    "Feature",
    "Feature_type",
    "Consequence",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
)


def _open_text(path: str | Path):
    raw = Path(path)
    if raw.suffix == ".gz" or raw.name.endswith(".vcf.gz"):
        return gzip.open(raw, "rt", encoding="utf-8")
    return raw.open(encoding="utf-8")


def _first_ann(info: str) -> list[str]:
    for field in info.split(";"):
        if field.startswith("ANN="):
            first = field[4:].split(",")[0]
            return first.split("|")
    return []


def snpeff_vcf_to_vep_rows(vcf_path: str | Path) -> pd.DataFrame:
    """One row per variant from the first ANN allele (SnpEff severity order)."""

    rows: list[dict[str, str]] = []
    with _open_text(vcf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            chrom, pos, var_id, _ref, alt, _qual, _filt, info = line.rstrip("\n").split("\t")[:8]
            ann = _first_ann(info)
            allele = ann[0] if len(ann) > 0 else alt.split(",")[0]
            consequence = ann[1] if len(ann) > 1 else ""
            gene = ann[3] if len(ann) > 3 else ""
            feature = ann[6] if len(ann) > 6 else ""
            feature_type = ann[5] if len(ann) > 5 else ""
            hgvs_c = ann[9] if len(ann) > 9 else ""
            hgvs_p = ann[10] if len(ann) > 10 else ""
            rows.append(
                {
                    "Uploaded_variation": var_id if var_id not in {".", ""} else f"{chrom}_{pos}_{allele}",
                    "Location": f"{chrom}:{pos}",
                    "Allele": allele,
                    "Gene": gene,
                    "Feature": feature,
                    "Feature_type": feature_type,
                    "Consequence": consequence,
                    "cDNA_position": "",
                    "CDS_position": hgvs_c,
                    "Protein_position": hgvs_p,
                    "Amino_acids": "",
                    "Codons": "",
                }
            )
    if not rows:
        raise ValueError(f"No SnpEff ANN rows in {vcf_path}")
    return pd.DataFrame(rows, columns=list(VEP_COLUMNS))


def run_snpeff_vcf_to_tab(vcf_path: str | Path, output_prefix: str | Path) -> pd.DataFrame:
    table = snpeff_vcf_to_vep_rows(vcf_path)
    save_df_to_tsv(table, f"{output_prefix}.snpeff.tsv")
    return table
