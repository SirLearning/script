#!/usr/bin/env python3
"""Plotting helpers for GWAS benchmark."""
from __future__ import annotations
import os
from typing import Dict
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_context("talk")
sns.set_style("whitegrid")

def plot_roc_pr(results: Dict[str, Dict[str, float]], outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)
    # ROC-AUC bar
    roc_items = [(name, vals.get("roc_auc", np.nan)) for name, vals in results.items()]
    pr_items = [(name, vals.get("pr_auc", np.nan)) for name, vals in results.items()]
    names, roc_vals = zip(*roc_items)
    _, pr_vals = zip(*pr_items)
    fig, ax = plt.subplots(figsize=(10,5))
    ax.bar(names, roc_vals, color="#1f77b4", alpha=0.8)
    ax.set_ylabel("ROC-AUC")
    ax.set_title("GWAS method ROC-AUC")
    ax.set_ylim(0.0, 1.05)
    for i,v in enumerate(roc_vals):
        ax.text(i, v + 0.01, f"{v:.3f}", ha='center', va='bottom', fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "roc_auc_bar.png"), dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10,5))
    ax.bar(names, pr_vals, color="#ff7f0e", alpha=0.8)
    ax.set_ylabel("PR-AUC")
    ax.set_title("GWAS method PR-AUC")
    ax.set_ylim(0.0, 1.05)
    for i,v in enumerate(pr_vals):
        ax.text(i, v + 0.01, f"{v:.3f}", ha='center', va='bottom', fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "pr_auc_bar.png"), dpi=150)
    plt.close(fig)


def plot_qq(expected: np.ndarray, observed: np.ndarray, method_name: str, outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6,6))
    ax.scatter(expected, observed, s=6, alpha=0.6)
    lim = max(expected.max(), observed.max())
    ax.plot([0, lim],[0, lim], color='red', lw=1)
    ax.set_xlabel("Expected -log10(P)")
    ax.set_ylabel("Observed -log10(P)")
    ax.set_title(f"QQ Plot: {method_name}")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"qq_{method_name}.png"), dpi=160)
    plt.close(fig)
