#!/usr/bin/env python3
"""GWAS Benchmark Toolkit: Run comparison against ground truth."""
from __future__ import annotations
import argparse
import os
import re
import yaml
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Iterable
from sklearn.metrics import roc_auc_score, average_precision_score
from infra.utils import load_df_generic, plot_bar_chart, plot_gwas_qq, save_df_to_tsv

# ==================================================================================
# Metrics
# ==================================================================================

MIN_P = 1e-300

def scores_from_p(p: np.ndarray) -> np.ndarray:
    """Convert p-values to positive scores (higher=stronger signal)."""
    p = np.clip(p.astype(float), MIN_P, 1.0)
    return -np.log10(p)

def compute_roc_pr(labels: np.ndarray, scores: np.ndarray) -> Dict[str, float]:
    """Compute ROC-AUC and PR-AUC (average precision)."""
    y = labels.astype(int)
    s = scores.astype(float)
    # Guard: if only one class present, metrics undefined
    if y.max() == y.min():
        return {"roc_auc": float("nan"), "pr_auc": float("nan")}
    try:
        roc = roc_auc_score(y, s)
        ap = average_precision_score(y, s)
        return {"roc_auc": float(roc), "pr_auc": float(ap)}
    except ValueError:
        return {"roc_auc": float("nan"), "pr_auc": float("nan")}

def top_k_recall(labels: np.ndarray, scores: np.ndarray, k: int) -> float:
    """Recall at top-K (fraction of all positives captured in top K ranked by score)."""
    if k <= 0: return 0.0
    order = np.argsort(-scores)
    topk = order[:k]
    pos_total = labels.sum()
    if pos_total == 0: return float("nan")
    return float(labels[topk].sum() / pos_total)

def enrichment_at_k(labels: np.ndarray, scores: np.ndarray, k: int) -> float:
    """Fold-enrichment among top-K vs background prevalence."""
    n = len(labels)
    if k <= 0 or n == 0: return float("nan")
    prev = labels.mean()
    if prev == 0: return float("nan")
    order = np.argsort(-scores)
    topk = order[:k]
    top_rate = labels[topk].mean()
    return float(top_rate / prev)

def qq_points(p: np.ndarray, max_points: int = 5000) -> Tuple[np.ndarray, np.ndarray]:
    """Compute QQ plot expected vs observed -log10(p) values."""
    p = np.clip(p.astype(float), MIN_P, 1.0)
    p_sorted = np.sort(p)
    n = len(p_sorted)
    if n == 0:
        return np.array([]), np.array([])
    if n > max_points:
        idx = np.linspace(0, n - 1, max_points).astype(int)
        p_sorted = p_sorted[idx]
        n = len(p_sorted)
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / (n + 0.5))
    obs = -np.log10(p_sorted)
    return exp, obs

def lambda_gc(p: np.ndarray) -> float:
    """Estimate genomic inflation factor lambda_GC from p-values."""
    try:
        from scipy.stats import chi2
        p = np.clip(p.astype(float), MIN_P, 1.0)
        chisq = chi2.isf(p, 1)
        lam = np.median(chisq) / 0.454936423119572
        return float(lam)
    except (ImportError, Exception):
        return float("nan")

# ==================================================================================
# Data Loading & Processing
# ==================================================================================

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
_REQUIRED = ["PVALUE"]

def _map_columns(columns: Iterable[str]) -> Dict[str, str]:
    mapped = {}
    for target, alts in _COL_MAP.items():
        for c in columns:
            c_norm = re.sub(r"[^A-Za-z0-9]", "", c).upper()
            if c_norm in {re.sub(r"[^A-Za-z0-9]", "", a).upper() for a in alts}:
                mapped[target] = c
                break
    return mapped

def load_gwas(path: str) -> pd.DataFrame:
    """Load and harmonize a GWAS result file into standardized columns."""
    if load_df_generic is None:
        raise ImportError("infra.utils.io is missing")
        
    df = load_df_generic(path)
    if df is None:
        raise FileNotFoundError(f"Could not read {path}")

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

    # Synthetic SNP ID if missing
    if "SNP" not in df_std and {"CHR","POS"}.issubset(df_std.columns):
        df_std["SNP"] = df_std["CHR"].astype(str) + ":" + df_std["POS"].astype(int).astype(str)

    # Cleanup
    df_std = df_std.dropna(subset=["PVALUE"]).reset_index(drop=True)
    return df_std

def load_truth(path: str) -> pd.DataFrame:
    """Load ground-truth causal variants or regions."""
    if load_df_generic is None:
        raise ImportError("infra.utils.io is missing")

    df = load_df_generic(path)
    if df is None:
        raise FileNotFoundError(path)
    
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
    """Tag input GWAS rows as HIT if they overlap with truth regions."""
    if not {"CHR","POS"}.issubset(df.columns):
        raise ValueError("GWAS data needs CHR and POS for hit tagging")
    df = df.copy()
    df["HIT"] = False
    
    point = truth[truth["TYPE"] == "POINT"] if "TYPE" in truth else truth
    interval = truth[truth["TYPE"] == "INTERVAL"] if "TYPE" in truth else pd.DataFrame(columns=["CHR","START","END"])

    intervals: List[tuple] = []
    if not interval.empty:
        intervals.extend((str(r.CHR), int(r.START), int(r.END)) for _, r in interval.iterrows())
    if not point.empty:
        intervals.extend((str(r.CHR), int(r.POS - window), int(r.POS + window)) for _, r in point.iterrows())
        
    # Valid intervals map
    by_chr: Dict[str, List[tuple]] = {}
    for chr_, s, e in intervals:
        by_chr.setdefault(str(chr_), []).append((s, e))
        
    # Check
    for chr_, group in df.groupby("CHR"):
        spans = by_chr.get(str(chr_), [])
        if not spans: continue
        
        # Optimize using numpy broadcasting or simple iteration
        positions = group["POS"].astype(int).values
        # For very large GWAS and many truth regions, this can be slow. 
        # But for benchmark purposes usually fine.
        hit_mask = np.zeros(len(positions), dtype=bool)
        for s, e in spans:
            hit_mask |= ((positions >= s) & (positions <= e))
        
        df.loc[group.index, "HIT"] = hit_mask
    return df

# ==================================================================================
# Workflow / Main
# ==================================================================================

def _read_manifest(path: str) -> Dict:
    if path.lower().endswith(('.yml','.yaml')):
        with open(path, 'r') as fh:
            return yaml.safe_load(fh)
    # CSV
    df = pd.read_csv(path)
    methods = [{"name": r["name"], "file": r["file"]} for _, r in df.iterrows()]
    return {"methods": methods}

def run_benchmark(
    manifest: str, 
    truth_path: Optional[str] = None, 
    truth_window: int = 50000, 
    topk: Optional[List[int]] = None
):
    """Main execution flow for GWAS benchmarking."""
    outdir = "."
    cfg = _read_manifest(manifest)
    methods = cfg.get("methods", [])
    if not methods:
        print("[Error] No methods in manifest")
        return

    if truth_path is None:
        truth_path = cfg.get("truth")
    if not truth_path:
        print("[Error] No truth provided")
        return

    if truth_window is None:
        truth_window = int(cfg.get("truth_window", 50000))

    try:
        truth = load_truth(truth_path)
    except Exception as e:
        print(f"[Error] Failed to load truth: {e}")
        return

    summary_rows = []
    pr_metrics = {}

    for m in methods:
        name = m.get("name") or os.path.basename(m.get("file"))
        path = m.get("file")
        if not os.path.exists(path):
            print(f"[WARN] File not found for {name}: {path}")
            continue
            
        print(f"[INFO] Processing {name} ({path})")
        try:
            df = load_gwas(path)
            df = tag_hits(df, truth, window=truth_window)
            df["SCORE"] = scores_from_p(df["PVALUE"].values)
            
            # Metrics
            labels = df["HIT"].values.astype(int)
            scores = df["SCORE"].values
            base = compute_roc_pr(labels, scores)
            
            if topk is None: topk = [10, 20, 50, 100]
            topk_vals = {f"top{k}_recall": top_k_recall(labels, scores, k) for k in topk}
            enr_vals = {f"top{k}_enrichment": enrichment_at_k(labels, scores, k) for k in topk}
            lam = lambda_gc(df["PVALUE"].values)
            
            row = {"method": name, **base, **topk_vals, **enr_vals, "lambda_gc": lam, "n": len(df), "positives": int(labels.sum())}
            summary_rows.append(row)
            pr_metrics[name] = base
            
            # Save labeled
            cols_out = [c for c in df.columns if c in ("CHR","POS","SNP","PVALUE","SCORE","HIT")]
            df_out = df[cols_out]
            save_df_to_tsv(df_out, os.path.join(outdir, f"{name}.labeled.tsv"))
            
            # Plots - using infra
            if plot_gwas_qq:
                exp, obs = qq_points(df["PVALUE"].values)
                
                # Save QQ data for later plotting
                qq_df = pd.DataFrame({"Expected": exp, "Observed": obs})
                save_df_to_tsv(qq_df, os.path.join(outdir, f"{name}.qq.tsv"))

                plot_gwas_qq(
                    expected=exp, 
                    observed=obs, 
                    title=f"QQ Plot: {name}",
                    filename=os.path.join(outdir, f"qq_{name}.png")
                )
            
        except Exception as e:
            print(f"[Error] Failed to process {name}: {e}")

    # Summary
    if not summary_rows:
        print("No valid methods processed.")
        return

    summary = pd.DataFrame(summary_rows).sort_values(by=["pr_auc","roc_auc"], ascending=[False, False])
    save_df_to_tsv(summary, os.path.join(outdir, "metrics_summary.tsv"))
    
    # Combined Bar Plots
    if plot_bar_chart and pr_metrics:
        # ROC
        items = [(n, v.get("roc_auc", np.nan)) for n, v in pr_metrics.items()]
        names, vals = zip(*items)
        plot_bar_chart(names, vals, "GWAS method ROC-AUC", "ROC-AUC", os.path.join(outdir, "roc_auc_bar.png"), ylim=(0.0, 1.05))
        
        # PR
        items = [(n, v.get("pr_auc", np.nan)) for n, v in pr_metrics.items()]
        names, vals = zip(*items)
        plot_bar_chart(names, vals, "GWAS method PR-AUC", "PR-AUC", os.path.join(outdir, "pr_auc_bar.png"), color='#ff7f0e', ylim=(0.0, 1.05))

    # Best
    if not summary.empty:
        best = summary.iloc[0]
        with open(os.path.join(outdir, "best_model.txt"), 'w') as fh:
            fh.write(f"Best method by PR-AUC: {best['method']} (PR-AUC={best['pr_auc']:.3f}, ROC-AUC={best['roc_auc']:.3f})\n")
    
    print(f"[Done] Benchmark complete. Results in {outdir}")

def main():
    p = argparse.ArgumentParser(description="Benchmark GWAS methods against ground truth")
    p.add_argument("--manifest", required=True, help="YAML or CSV manifest listing methods")
    p.add_argument("--truth", help="Ground-truth file (overrides manifest)")
    p.add_argument("--truth-window", type=int, default=None, help="Window flank for point truth (bp)")
    p.add_argument("--topk", default="10,20,50,100", help="Comma list of K values")
    args = p.parse_args()

    topks = [int(x) for x in args.topk.split(',') if x]
    run_benchmark(args.manifest, truth_path=args.truth, truth_window=args.truth_window, topk=topks)

if __name__ == "__main__":
    main()
