#!/usr/bin/env python3
"""Evaluation metrics for GWAS benchmarking."""
from __future__ import annotations
from typing import Dict, Tuple
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score

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
    roc = roc_auc_score(y, s)
    ap = average_precision_score(y, s)
    return {"roc_auc": float(roc), "pr_auc": float(ap)}

def top_k_recall(labels: np.ndarray, scores: np.ndarray, k: int) -> float:
    """Recall at top-K (fraction of all positives captured in top K ranked by score)."""
    if k <= 0:
        return 0.0
    order = np.argsort(-scores)
    topk = order[:k]
    pos_total = labels.sum()
    if pos_total == 0:
        return float("nan")
    return float(labels[topk].sum() / pos_total)

def enrichment_at_k(labels: np.ndarray, scores: np.ndarray, k: int) -> float:
    """Fold-enrichment among top-K vs background prevalence."""
    n = len(labels)
    if k <= 0 or n == 0:
        return float("nan")
    prev = labels.mean()
    if prev == 0:
        return float("nan")
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
    """Estimate genomic inflation factor lambda_GC from p-values.

    Requires scipy if available; otherwise returns NaN.
    """
    try:
        from scipy.stats import chi2
    except Exception:
        return float("nan")
    p = np.clip(p.astype(float), MIN_P, 1.0)
    # Convert to 1df chi-square
    chisq = chi2.isf(p, 1)
    lam = np.median(chisq) / 0.454936423119572
    return float(lam)
