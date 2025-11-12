#!/usr/bin/env python3
"""IO utilities for GWAS benchmark.

This module loads GWAS result tables from different software formats and harmonizes them
into a common schema:

Columns (internal representation):
    CHR: chromosome (string)
    POS: 1-based position (int)
    SNP: variant identifier (may be CHR:POS or rsid)
    PVALUE: association p-value (float)
    EFFECT: effect size or beta (float, optional)
    SE: standard error (float, optional)
    MAF: minor allele frequency (float, optional)
    INFO: imputation quality / other metric (float, optional)
    SAMPLE_SIZE: sample size (int, optional)

Detection heuristics:
    - Auto-delimiter detection (tab, comma, space) using first non-header line.
    - Column name mapping using fuzzy dictionary.
    - Fallback: if both CHR and POS missing but SNP present as 'chr:pos', split.

Supported input file types: plain text (.txt, .tsv, .csv, .gz) compressed allowed.

"""
from __future__ import annotations
import gzip
import os
import re
from typing import Iterable, Dict, List, Optional
import pandas as pd

# Fuzzy column mapping patterns
_COL_MAP = {
    "CHR": ["CHR", "CHROM", "CHROMOSOME"],
    "POS": ["POS", "BP", "POSITION"],
    "SNP": ["SNP", "MARKER", "ID", "RSID"],
    "PVALUE": ["P", "PVALUE", "PVAL", "P_VALUE", "P.VALUE", "PVALUE"],
    "EFFECT": ["BETA", "EFFECT", "B", "BETA_HAT", "BETA.HAT"],
    "SE": ["SE", "STD_ERR", "STDERR", "SE_EFFECT"],
    "MAF": ["MAF", "FREQ", "ALLELE_FREQ", "AF"],
    "INFO": ["INFO", "IMPUTE_INFO", "R2", "INFO_SCORE"],
    "SAMPLE_SIZE": ["N", "N_TOTAL", "SAMPLESIZE", "NCASES", "NCASES+NCONTROLS"],
}

def _open(path: str):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'rt')

def _detect_delimiter(first_line: str) -> Optional[str]:
    if '\t' in first_line:
        return '\t'
    if ',' in first_line:
        return ','
    # Multiple spaces -> treat as whitespace
    if re.search(r"\s{2,}", first_line):
        return None  # pandas will handle delim=None for delim_whitespace
    return '\t'  # fallback

def _map_columns(columns: Iterable[str]) -> Dict[str, str]:
    mapped = {}
    for target, alts in _COL_MAP.items():
        for c in columns:
            c_norm = re.sub(r"[^A-Za-z0-9]", "", c).upper()
            if c_norm in {re.sub(r"[^A-Za-z0-9]", "", a).upper() for a in alts}:
                mapped[target] = c
                break
    return mapped

_REQUIRED = ["PVALUE"]

def load_gwas(path: str) -> pd.DataFrame:
    """Load and harmonize a GWAS result file into standardized columns.

    Args:
        path: Path to GWAS results (can be .gz)
    Returns:
        DataFrame with standardized columns (at least PVALUE, plus CHR, POS if derivable)
    """
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    with _open(path) as fh:
        first = fh.readline()
    delim = _detect_delimiter(first)
    df = pd.read_csv(path, sep=delim, compression='infer')

    colmap = _map_columns(df.columns)
    df_std = pd.DataFrame()
    for std_name, orig in colmap.items():
        df_std[std_name] = df[orig]

    # Derive CHR, POS from SNP if necessary
    if ("CHR" not in df_std or "POS" not in df_std) and "SNP" in df_std:
        parts = df_std["SNP"].astype(str).str.extract(r"^(chr)?(?P<CHR>[0-9A-Za-z]+)[:_](?P<POS>[0-9]+)")
        if "CHR" not in df_std and "CHR" in parts:
            df_std["CHR"] = parts["CHR"]
        if "POS" not in df_std and "POS" in parts:
            df_std["POS"] = parts["POS"].astype(float).astype('Int64')

    # Ensure PVALUE exists
    missing_req = [c for c in _REQUIRED if c not in df_std]
    if missing_req:
        raise ValueError(f"Missing required columns: {missing_req} in {path}")

    # Coerce types
    if "PVALUE" in df_std:
        df_std["PVALUE"] = pd.to_numeric(df_std["PVALUE"], errors='coerce')
    if "POS" in df_std:
        df_std["POS"] = pd.to_numeric(df_std["POS"], errors='coerce')
    if "CHR" in df_std:
        df_std["CHR"] = df_std["CHR"].astype(str)

    # Create SNP synthetic if missing
    if "SNP" not in df_std and {"CHR","POS"}.issubset(df_std.columns):
        df_std["SNP"] = df_std["CHR"].astype(str) + ":" + df_std["POS"].astype(int).astype(str)

    # Drop rows with missing p-values
    df_std = df_std.dropna(subset=["PVALUE"]).reset_index(drop=True)
    return df_std

def load_truth(path: str) -> pd.DataFrame:
    """Load ground-truth causal variants or regions.

    Accepted formats:
        - columns CHR POS (single causal variants)
        - columns CHR START END (interval regions)
    Returns DataFrame standardized.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    delim = None
    with _open(path) as fh:
        head = fh.readline()
        delim = _detect_delimiter(head)
    df = pd.read_csv(path, sep=delim, compression='infer')
    cols = [c.upper() for c in df.columns]
    df.columns = cols
    if {"CHR","POS"}.issubset(cols):
        out = df[["CHR","POS"]].copy()
        out["TYPE"] = "POINT"
        return out
    elif {"CHR","START","END"}.issubset(cols):
        out = df[["CHR","START","END"]].copy()
        out["TYPE"] = "INTERVAL"
        return out
    else:
        raise ValueError("Truth file must have CHR POS or CHR START END")

def tag_hits(df: pd.DataFrame, truth: pd.DataFrame, window: int = 50000) -> pd.DataFrame:
    """Tag GWAS rows as HIT if they fall within truth regions (POINT +/- window or INTERVAL).

    Args:
        df: harmonized GWAS DataFrame (expects CHR, POS, PVALUE)
        truth: DataFrame from load_truth
        window: flank size for POINT causal variants
    Returns:
        DataFrame with an added boolean HIT column.
    """
    if not {"CHR","POS"}.issubset(df.columns):
        raise ValueError("GWAS data needs CHR and POS for hit tagging")
    df = df.copy()
    df["HIT"] = False
    # Separate point and interval truths
    point = truth[truth["TYPE"] == "POINT"] if "TYPE" in truth else truth
    interval = truth[truth["TYPE"] == "INTERVAL"] if "TYPE" in truth else pd.DataFrame(columns=["CHR","START","END"])

    # Build interval list
    intervals: List[tuple] = []
    if not interval.empty:
        intervals.extend((str(r.CHR), int(r.START), int(r.END)) for _, r in interval.iterrows())
    if not point.empty:
        intervals.extend((str(r.CHR), int(r.POS - window), int(r.POS + window)) for _, r in point.iterrows())
    # Mark hits
    by_chr: Dict[str, List[tuple]] = {}
    for chr_, s, e in intervals:
        by_chr.setdefault(chr_, []).append((s, e))
    # Simple interval check
    for chr_, group in df.groupby("CHR"):
        spans = by_chr.get(str(chr_), [])
        if not spans:
            continue
        positions = group["POS"].astype(int).values
        hit_mask = [any(s <= p <= e for s, e in spans) for p in positions]
        df.loc[group.index, "HIT"] = hit_mask
    return df
