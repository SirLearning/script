"""Classify putative LoF alleles as derived, ancestral, or ancestry-unknown.

The Triticeae syntenic matrix under ``00data/06maf/syntenic_sites_k2`` stores
one column per Wheat-Aegilops-complex genome. Ancestral state is the consensus
of **lineage outgroups** (einkorn / emmer / tauschii relatives), never Chinese
Spring itself. A LoF allele that matches that consensus is treated as a
suspicious annotation artefact (LOFTEE-style ANC_ALLELE). Missing or
conflicting outgroup bases stay ``unknown`` and are not dropped with the
ancestral class.
"""

from __future__ import annotations

import gzip
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import pandas as pd

from genetics.genomics.alignment.collinear_sites import LINEAGE_SPECIES
from genetics.genomics.gene.gff_bed import gene_lengths_from_gff
from genetics.genomics.variant.coding_mask import (
    PLOF_HIGH_TERMS,
    SEVERITY_ORDER,
    VEP_META_PREFIX,
    most_severe_consequence,
)
from infra.utils.graph import (
    SUBGENOME_COLORS,
    SUBGENOME_ORDER,
    plot_bar_chart,
    plot_correlation_with_regression,
    plot_joint_regression,
)
from infra.utils.io import save_df_to_tsv, save_thresholds
from infra.wheat.ref_v1 import get_chromosome, get_pos_on_chromosome

CLASS_DERIVED = "derived"
CLASS_ANCESTRAL = "ancestral"
CLASS_UNKNOWN = "unknown"
CLASS_ORDER = (CLASS_DERIVED, CLASS_ANCESTRAL, CLASS_UNKNOWN)

_BASES = frozenset("ACGT")
_RANK = {term: i for i, term in enumerate(SEVERITY_ORDER)}


def outgroups_for_subgenome(subgenome: str) -> tuple[str, ...]:
    """Lineage species used to polarize, excluding bread-wheat ``traes*``."""

    key = str(subgenome).upper()
    if key not in LINEAGE_SPECIES:
        raise ValueError(f"Unsupported subgenome: {subgenome}")
    return tuple(name for name in LINEAGE_SPECIES[key] if not name.startswith("traes"))


def infer_ancestral_allele(bases: Mapping[str, str], outgroups: Sequence[str]) -> str:
    """Consensus ACGT among outgroups; empty if none or they conflict."""

    observed: list[str] = []
    for name in outgroups:
        token = str(bases.get(name, "-")).strip().upper()
        if token in _BASES:
            observed.append(token)
    unique = set(observed)
    if len(unique) != 1:
        return ""
    return unique.pop()


def classify_lof_ancestry(lof_allele: str, ancestral_allele: str) -> str:
    """Compare the VEP LoF allele to a polarized ancestral base."""

    lof = str(lof_allele or "").strip().upper()
    anc = str(ancestral_allele or "").strip().upper()
    if lof not in _BASES or anc not in _BASES:
        return CLASS_UNKNOWN
    if lof == anc:
        return CLASS_ANCESTRAL
    return CLASS_DERIVED


def segment_to_matrix_key(segment_id: int, pos_1based: int) -> tuple[str, int]:
    """Map FastCall segment coordinates to syntenic-matrix (chr label, 0-based pos)."""

    chr_label = get_chromosome(int(segment_id))
    if not chr_label or chr_label[-1] not in "ABD":
        raise ValueError(f"No ABD chromosome for segment {segment_id}")
    linear = int(get_pos_on_chromosome(int(segment_id), int(pos_1based)))
    return str(chr_label), linear - 1


def lookup_ancestral_alleles(
    matrix_path: str | Path,
    positions_0based: Iterable[int],
    outgroups: Sequence[str],
) -> dict[int, str]:
    """Stream one ``*.syntenic.matrix.tsv.gz`` and polarize requested positions."""

    wanted = {int(pos) for pos in positions_0based}
    found = {pos: "" for pos in wanted}
    if not wanted:
        return found
    path = Path(matrix_path)
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        try:
            pos_idx = header.index("ref_pos")
        except ValueError as exc:
            raise ValueError(f"{path} missing ref_pos column") from exc
        og_idx = {name: header.index(name) for name in outgroups if name in header}
        if not og_idx:
            raise ValueError(f"{path} has none of the outgroup columns {list(outgroups)}")
        remaining = set(wanted)
        max_pos = max(remaining)
        for line in handle:
            if not remaining:
                break
            cols = line.rstrip("\n").split("\t")
            pos = int(cols[pos_idx])
            if pos > max_pos:
                break
            if pos not in remaining:
                continue
            bases = {name: cols[idx] for name, idx in og_idx.items()}
            found[pos] = infer_ancestral_allele(bases, outgroups)
            remaining.remove(pos)
    return found


def _parse_variant_id(variant_id: str) -> tuple[int, int] | None:
    text = str(variant_id)
    if "-" not in text:
        return None
    left, right = text.split("-", 1)
    try:
        return int(left), int(right)
    except ValueError:
        return None


def collect_lof_sites(vep_path: str | Path) -> pd.DataFrame:
    """One row per unique LoF site (most severe VEP consequence, with Allele)."""

    path = Path(vep_path)
    best: dict[str, tuple[int, str, str, str]] = {}
    header: list[str] | None = None
    idx: dict[str, int] | None = None
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(VEP_META_PREFIX):
                continue
            if header is None:
                header = [part.lstrip("#") for part in line.rstrip("\n").split("\t")]
                try:
                    idx = {
                        "id": header.index("Uploaded_variation"),
                        "allele": header.index("Allele"),
                        "gene": header.index("Gene"),
                        "cons": header.index("Consequence"),
                    }
                except ValueError as exc:
                    raise ValueError(f"{path} missing VEP columns") from exc
                continue
            cols = line.rstrip("\n").split("\t")
            vid = cols[idx["id"]]
            term = most_severe_consequence(cols[idx["cons"]].replace("&", ","))
            rank = _RANK.get(term, 1000)
            prev = best.get(vid)
            if prev is None or rank < prev[0]:
                best[vid] = (rank, cols[idx["allele"]], cols[idx["gene"]], term)
    rows = [
        {"id": vid, "lof_allele": allele, "gene": gene, "consequence": term}
        for vid, (_rank, allele, gene, term) in best.items()
        if term in PLOF_HIGH_TERMS
    ]
    return pd.DataFrame(rows)


def _lookup_one_chromosome(
    chr_label: str,
    positions: list[int],
    matrix_path: str,
    subgenome: str,
) -> tuple[str, dict[int, str]]:
    return chr_label, lookup_ancestral_alleles(
        matrix_path,
        positions,
        outgroups_for_subgenome(subgenome),
    )


def run_lof_ancestry(
    vep_paths: Sequence[str | Path],
    matrix_dir: str | Path,
    output_prefix: str | Path,
    *,
    workers: int = 1,
) -> pd.DataFrame:
    """Polarize LoF sites from VEP tables against syntenic outgroup matrices."""

    tables = [collect_lof_sites(path) for path in vep_paths]
    if not tables:
        raise ValueError("No VEP tables")
    sites = pd.concat(tables, ignore_index=True)
    keys = sites["id"].map(_parse_variant_id)
    sites = sites.loc[keys.notna()].copy()
    parsed = keys.loc[sites.index]
    labels_pos = [segment_to_matrix_key(seg, pos) for seg, pos in parsed]
    sites["chr_label"] = [item[0] for item in labels_pos]
    sites["ref_pos"] = [item[1] for item in labels_pos]
    sites["subgenome"] = sites["chr_label"].str[-1]
    sites["ancestral_allele"] = ""
    matrix_root = Path(matrix_dir)
    jobs = []
    for chr_label, block in sites.groupby("chr_label", sort=False):
        matrix_path = matrix_root / f"{chr_label}.syntenic.matrix.tsv.gz"
        if not matrix_path.exists():
            continue
        jobs.append(
            (
                str(chr_label),
                block["ref_pos"].astype(int).tolist(),
                str(matrix_path),
                str(chr_label)[-1],
            )
        )
    looked_by_chr: dict[str, dict[int, str]] = {}
    if workers <= 1:
        for job in jobs:
            chr_label, mapping = _lookup_one_chromosome(*job)
            looked_by_chr[chr_label] = mapping
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(_lookup_one_chromosome, *job) for job in jobs]
            for future in as_completed(futures):
                chr_label, mapping = future.result()
                looked_by_chr[chr_label] = mapping
    for chr_label, block in sites.groupby("chr_label", sort=False):
        mapping = looked_by_chr.get(str(chr_label), {})
        sites.loc[block.index, "ancestral_allele"] = block["ref_pos"].map(mapping).fillna("")
    sites["ancestry_class"] = [
        classify_lof_ancestry(lof, anc)
        for lof, anc in zip(sites["lof_allele"], sites["ancestral_allele"])
    ]
    prefix = Path(output_prefix)
    save_df_to_tsv(sites, f"{prefix}.lof_ancestry.tsv")
    counts = sites["ancestry_class"].value_counts()
    n_derived = int(counts.get(CLASS_DERIVED, 0))
    n_ancestral = int(counts.get(CLASS_ANCESTRAL, 0))
    n_unknown = int(counts.get(CLASS_UNKNOWN, 0))
    stats = {
        "n_lof": int(len(sites)),
        "n_derived": n_derived,
        "n_ancestral": n_ancestral,
        "n_unknown": n_unknown,
        "n_remain_drop_ancestral": n_derived + n_unknown,
        "n_remain_derived_only": n_derived,
    }
    for sg in ("A", "B", "D"):
        sub = sites.loc[sites["subgenome"] == sg, "ancestry_class"].value_counts()
        stats[f"n_lof_{sg}"] = int((sites["subgenome"] == sg).sum())
        stats[f"n_derived_{sg}"] = int(sub.get(CLASS_DERIVED, 0))
        stats[f"n_ancestral_{sg}"] = int(sub.get(CLASS_ANCESTRAL, 0))
        stats[f"n_unknown_{sg}"] = int(sub.get(CLASS_UNKNOWN, 0))
    save_thresholds(stats, f"{prefix}.lof_ancestry.th.tsv")
    names = list(CLASS_ORDER)
    values = [float(counts.get(name, 0)) for name in names]
    plot_bar_chart(
        names=names,
        values=values,
        title="LoF ancestry class",
        ylabel="Sites",
        filename=f"{prefix}.lof_ancestry.class.png",
        ylim=(0.0, (max(values) * 1.15) if values else 1.05),
    )
    return sites


def _pearson_spearman(frame: pd.DataFrame, x_col: str, y_col: str) -> tuple[float, float]:
    work = frame[[x_col, y_col]].dropna()
    if len(work) < 3:
        return float("nan"), float("nan")
    return (
        float(work[x_col].corr(work[y_col], method="pearson")),
        float(work[x_col].corr(work[y_col], method="spearman")),
    )


def _gene_column(frame: pd.DataFrame) -> str:
    for name in ("gene", "gene_id", "Gene"):
        if name in frame.columns:
            return name
    raise ValueError("table needs a gene, gene_id, or Gene column")


def _flag_mask(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False)
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def _flagged_id_gene(raw: pd.DataFrame, flag_col: str) -> pd.DataFrame:
    if "ID" not in raw.columns:
        raise ValueError("flagged-site table missing ID")
    gene_col = _gene_column(raw)
    if flag_col not in raw.columns:
        raise ValueError(f"flagged-site table missing {flag_col}")
    flag = _flag_mask(raw[flag_col])
    keep = raw.loc[flag, ["ID", gene_col]].copy()
    keep = keep.rename(columns={gene_col: "gene"})
    keep["ID"] = keep["ID"].astype(str)
    keep["gene"] = keep["gene"].astype(str)
    return keep.drop_duplicates("ID")


def collect_flagged_sites(
    tables: Sequence[str | Path | pd.DataFrame],
    *,
    flag_col: str = "rare_wap",
) -> pd.DataFrame:
    """Unique variant IDs whose boolean flag is true, with the VEP gene."""

    frames: list[pd.DataFrame] = []
    for table in tables:
        if isinstance(table, pd.DataFrame):
            raw = table.copy()
        else:
            path = Path(table)
            header = pd.read_csv(path, sep="\t", nrows=0)
            wanted = [name for name in ("ID", "Gene", "gene", "gene_id", flag_col) if name in header.columns]
            parts: list[pd.DataFrame] = []
            for chunk in pd.read_csv(
                path, sep="\t", usecols=wanted, dtype=str, chunksize=500_000
            ):
                parts.append(_flagged_id_gene(chunk, flag_col))
            if not parts:
                continue
            frames.append(pd.concat(parts, ignore_index=True).drop_duplicates("ID"))
            continue
        frames.append(_flagged_id_gene(raw, flag_col))
    if not frames:
        raise ValueError("No flagged-site tables")
    out = pd.concat(frames, ignore_index=True).drop_duplicates("ID")
    return out.rename(columns={"ID": "id"})


def count_flagged_per_gene(
    tables: Sequence[str | Path | pd.DataFrame],
    *,
    flag_col: str = "rare_wap",
    count_col: str = "n_rare",
) -> pd.DataFrame:
    """One row per gene: how many unique-table rows have ``flag_col`` true."""

    counts: Counter[str] = Counter()
    empty = {"", "nan", "None", ".", "-"}
    for table in tables:
        if isinstance(table, pd.DataFrame):
            chunks: Iterable[pd.DataFrame] = (table,)
        else:
            path = Path(table)
            header = pd.read_csv(path, sep="\t", nrows=0)
            gene_names = [name for name in ("Gene", "gene", "gene_id") if name in header.columns]
            if not gene_names:
                raise ValueError("flagged-site table missing gene column")
            if flag_col not in header.columns:
                raise ValueError(f"flagged-site table missing {flag_col}")
            chunks = pd.read_csv(
                path,
                sep="\t",
                usecols=[gene_names[0], flag_col],
                chunksize=500_000,
                low_memory=False,
            )
        for chunk in chunks:
            gene_col = _gene_column(chunk)
            if flag_col not in chunk.columns:
                raise ValueError(f"flagged-site table missing {flag_col}")
            genes = chunk[gene_col].fillna("").astype(str)
            counts.update(genes[_flag_mask(chunk[flag_col])])
    for key in empty:
        counts.pop(key, None)
    if not counts:
        return pd.DataFrame(columns=["gene", count_col])
    return pd.DataFrame({"gene": list(counts), count_col: list(counts.values())})


def run_lof_vs_gene_length(
    lof_table: str | Path | pd.DataFrame,
    gff_path: str | Path,
    output_prefix: str | Path,
    *,
    count_col: str = "n_lof",
    y_noun: str = "LoF",
    file_stem: str = "lof",
) -> pd.DataFrame:
    """Count sites per gene and plot against CDS length and gene span."""

    if isinstance(lof_table, pd.DataFrame):
        lof = lof_table.copy()
    else:
        lof = pd.read_csv(lof_table, sep="\t")
    gene_col = _gene_column(lof)
    if count_col in lof.columns:
        counts = (
            pd.to_numeric(lof[count_col], errors="coerce")
            .fillna(0)
            .groupby(lof[gene_col].astype(str))
            .sum()
            .rename(count_col)
        )
    else:
        counts = lof.groupby(lof[gene_col].astype(str)).size().rename(count_col)
    genes = gene_lengths_from_gff(gff_path)
    table = genes.merge(counts, left_on="gene_id", right_index=True, how="left")
    table[count_col] = table[count_col].fillna(0).astype(int)
    table["cds_kb"] = table["cds_bp"] / 1000.0
    table["gene_span_kb"] = table["gene_span_bp"] / 1000.0
    table["n_per_kb_cds"] = [
        (n / kb) if kb > 0 else 0.0 for n, kb in zip(table[count_col], table["cds_kb"])
    ]
    coding = table.loc[table["cds_bp"] > 0].copy()
    with_hits = coding.loc[coding[count_col] > 0].copy()

    prefix = Path(output_prefix)
    save_df_to_tsv(table, f"{prefix}.{file_stem}_vs_gene_length.tsv")
    r_cds_all, rho_cds_all = _pearson_spearman(coding, "cds_bp", count_col)
    r_cds_hit, rho_cds_hit = _pearson_spearman(with_hits, "cds_bp", count_col)
    r_span_hit, rho_span_hit = _pearson_spearman(with_hits, "gene_span_bp", count_col)
    r_den, rho_den = _pearson_spearman(coding, "cds_bp", "n_per_kb_cds")
    cds_kb_sum = float(coding["cds_bp"].sum() / 1000.0)
    stats = {
        "n_genes": int(len(table)),
        "n_genes_coding": int(len(coding)),
        "n_genes_with_sites": int(len(with_hits)),
        "n_sites": int(table[count_col].sum()),
        "pearson_cds_all": r_cds_all,
        "spearman_cds_all": rho_cds_all,
        "pearson_cds_with_sites": r_cds_hit,
        "spearman_cds_with_sites": rho_cds_hit,
        "pearson_span_with_sites": r_span_hit,
        "spearman_span_with_sites": rho_span_hit,
        "pearson_density_cds": r_den,
        "spearman_density_cds": rho_den,
        "mean_n_per_kb_cds": float(table[count_col].sum() / cds_kb_sum) if cds_kb_sum else float("nan"),
    }
    save_thresholds(stats, f"{prefix}.{file_stem}_vs_gene_length.th.tsv")

    y_label = f"{y_noun} sites per gene"
    plot_correlation_with_regression(
        coding,
        "cds_kb",
        count_col,
        title=f"{y_noun} sites vs CDS length (all HC genes)",
        filename=f"{prefix}.{file_stem}_vs_cds.png",
        x_label="CDS length (kb)",
        y_label=y_label,
    )
    plot_correlation_with_regression(
        coding,
        "cds_kb",
        "n_per_kb_cds",
        title=f"{y_noun} sites per kb CDS vs CDS length",
        filename=f"{prefix}.{file_stem}_per_kb_vs_cds.png",
        x_label="CDS length (kb)",
        y_label=f"{y_noun} sites / kb CDS",
    )
    if len(with_hits) >= 2:
        plot_correlation_with_regression(
            with_hits,
            "cds_kb",
            count_col,
            title=f"{y_noun} sites vs CDS length (genes with sites)",
            filename=f"{prefix}.{file_stem}_vs_cds.nonzero.png",
            x_label="CDS length (kb)",
            y_label=y_label,
        )
        plot_correlation_with_regression(
            with_hits,
            "gene_span_kb",
            count_col,
            title=f"{y_noun} sites vs gene span (genes with sites)",
            filename=f"{prefix}.{file_stem}_vs_span.nonzero.png",
            x_label="Gene span (kb)",
            y_label=y_label,
        )
        sub = with_hits.loc[with_hits["subgenome"].isin(SUBGENOME_ORDER)].copy()
        if len(sub) >= 10:
            plot_joint_regression(
                sub,
                "cds_kb",
                count_col,
                "subgenome",
                "CDS length (kb)",
                y_label,
                f"{prefix}.{file_stem}_vs_cds.joint.png",
                title=f"{y_noun} vs CDS length",
                palette=dict(SUBGENOME_COLORS),
                hue_order=list(SUBGENOME_ORDER),
                rasterized=True,
                s=12,
                alpha=0.35,
                min_points=10,
            )
    return table

