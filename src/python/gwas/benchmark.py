#!/usr/bin/env python3
"""CLI to benchmark GWAS methods and select best-performing models.

Example YAML manifest:

methods:
  - name: GEMMA_LMM
    file: /path/to/gemma.assoc.txt
  - name: PLINK_MLM
    file: /path/to/plink.mlm.txt
truth: /path/to/causal_truth.tsv
truth_window: 50000
phenotype: PHENO1

CSV manifest alternative (columns: name,file):
name,file
GEMMA_LMM,/path/to/gemma.assoc.txt
PLINK_MLM,/path/to/plink.mlm.txt

Usage:
  python -m gwas.benchmark --manifest bench.yaml --outdir results/bench

"""
from __future__ import annotations
import argparse
import os
import sys
from typing import Dict, List

import pandas as pd
import yaml

from .io import load_gwas, load_truth, tag_hits
from .metrics import scores_from_p, compute_roc_pr, top_k_recall, enrichment_at_k, qq_points, lambda_gc
from .plotting import plot_roc_pr, plot_qq


def _read_manifest(path: str) -> Dict:
    if path.lower().endswith(('.yml','.yaml')):
        with open(path, 'r') as fh:
            return yaml.safe_load(fh)
    # CSV
    df = pd.read_csv(path)
    methods = [{"name": r["name"], "file": r["file"]} for _, r in df.iterrows()]
    return {"methods": methods}


def run_benchmark(manifest: str, outdir: str, truth_path: str = None, truth_window: int = 50000, topk: List[int] = None) -> None:
    os.makedirs(outdir, exist_ok=True)
    cfg = _read_manifest(manifest)
    methods = cfg.get("methods", [])
    if not methods:
        raise SystemExit("No methods in manifest")

    if truth_path is None:
        truth_path = cfg.get("truth")
    if truth_path is None:
        raise SystemExit("No truth provided (use --truth or set in manifest)")

    if truth_window is None:
        truth_window = int(cfg.get("truth_window", 50000))

    truth = load_truth(truth_path)

    summary_rows = []
    pr_metrics = {}

    for m in methods:
        name = m.get("name") or os.path.basename(m.get("file"))
        path = m.get("file")
        if not os.path.exists(path):
            print(f"[WARN] File not found for {name}: {path}", file=sys.stderr)
            continue
        print(f"[INFO] Loading {name} from {path}")
        df = load_gwas(path)
        # Label hits
        df = tag_hits(df, truth, window=truth_window)
        df["SCORE"] = scores_from_p(df["PVALUE"].values)
        # Metrics
        labels = df["HIT"].values.astype(int)
        scores = df["SCORE"].values
        base = compute_roc_pr(labels, scores)
        # top-k
        if topk is None:
            topk = [10, 20, 50, 100]
        topk_vals = {f"top{k}_recall": top_k_recall(labels, scores, k) for k in topk}
        enr_vals = {f"top{k}_enrichment": enrichment_at_k(labels, scores, k) for k in topk}
        lam = lambda_gc(df["PVALUE"].values)
        row = {"method": name, **base, **topk_vals, **enr_vals, "lambda_gc": lam, "n": len(df), "positives": int(labels.sum())}
        summary_rows.append(row)
        pr_metrics[name] = base
        # Save labeled table
        df_out = df[[c for c in df.columns if c in ("CHR","POS","SNP","PVALUE","SCORE","HIT")]].copy()
        df_out.to_csv(os.path.join(outdir, f"{name}.labeled.tsv"), sep='\t', index=False)
        # QQ plot
        exp, obs = qq_points(df["PVALUE"].values)
        plot_qq(exp, obs, name, outdir)

    # Summary
    if not summary_rows:
        raise SystemExit("No valid methods processed")
    summary = pd.DataFrame(summary_rows).sort_values(by=["pr_auc","roc_auc"], ascending=[False, False])
    summary.to_csv(os.path.join(outdir, "metrics_summary.tsv"), sep='\t', index=False)
    # ROC/PR bar plots
    plot_roc_pr(pr_metrics, outdir)

    # Model selection recommendation
    best = summary.iloc[0]
    with open(os.path.join(outdir, "best_model.txt"), 'w') as fh:
        fh.write(f"Best method by PR-AUC: {best['method']} (PR-AUC={best['pr_auc']:.3f}, ROC-AUC={best['roc_auc']:.3f})\n")


def main(argv=None):
    p = argparse.ArgumentParser(description="Benchmark GWAS methods against ground truth")
    p.add_argument("--manifest", required=True, help="YAML or CSV manifest listing methods")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--truth", help="Ground-truth file (overrides manifest)")
    p.add_argument("--truth-window", type=int, default=None, help="Window flank for point truth (bp)")
    p.add_argument("--topk", default="10,20,50,100", help="Comma list of K values")
    args = p.parse_args(argv)

    topk = [int(x) for x in args.topk.split(',') if x]
    run_benchmark(args.manifest, args.outdir, truth_path=args.truth, truth_window=args.truth_window, topk=topk)

if __name__ == "__main__":
    main()
