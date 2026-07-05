"""Plot PopDep / PopDepFull / PopDepCrossChr benchmark summaries and monitor time series."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from infra.utils.graph import (
    LEGEND_FONT_SIZE,
    TICK_FONT_SIZE,
    TITLE_FONT_SIZE,
    X_LABEL_FONT_SIZE,
    Y_LABEL_FONT_SIZE,
    plot_bar_chart,
    plot_multi_line_series,
)
from infra.utils.io import load_df_generic, save_df_to_tsv

VMAP4_STATS_GENOME = Path("/data/home/tusr1/01projects/vmap4/10stats.genome")
POPDEP_BENCH_PUBLISH_DIR = Path(
    "/data1/dazheng_tusr1/vmap4.VCF.v1/benchmark/popdep_bench"
)
H2H_1M96_DATA_DIR = Path(
    "/data/home/tusr1/01projects/vmap4/00data/popdep_bench/head2head_20260619_175321"
)

H2H_1M96_SERIES = (
    ("PopDep", "monitor_popdep_rel.tsv", "#d62728"),
    ("PopDepFull", "monitor_popdepfull_rel.tsv", "#1f77b4"),
    ("PopDepCrossChr", "monitor_popdepcrosschr_rel.tsv", "#2ca02c"),
)
H2H_1M96_FILL_ALPHA = 0.32
H2H_1M96_ACTIVE_SDB_RMB_S = 50.0

MONITOR_NUMERIC = (
    "java_rss_gb",
    "samtools_n",
    "samtools_cpu_avg",
    "cpu_user",
    "cpu_iowait",
    "sdb_rMB_s",
    "sdb_wMB_s",
    "sdb_aqu_sz",
    "sdb_util",
)


@dataclass(frozen=True)
class BenchCampaign:
    campaign_id: str
    label: str
    jar: str
    bench_dir: Path
    full_summary: Path | None = None
    crosschr_summary: Path | None = None
    head2head_summary: Path | None = None
    crosschr_threads_summary: Path | None = None  # alias for 09run naming


def default_campaigns(stats_root: Path = VMAP4_STATS_GENOME) -> list[BenchCampaign]:
    root = Path(stats_root)
    return [
        BenchCampaign(
            campaign_id="n200_13run",
            label="N200 13run",
            jar="TIGER_PD_20260616",
            bench_dir=root / "13run_popdep_n200_bench",
            full_summary=root / "13run_popdep_n200_bench/run_logs/full_fork_summary.tsv",
            crosschr_summary=root / "13run_popdep_n200_bench/run_logs/crosschr_threads_summary.tsv",
        ),
        BenchCampaign(
            campaign_id="n200_14run",
            label="N200 14run",
            jar="TIGER_PD_20260619",
            bench_dir=root / "14run_popdep_n200_bench",
            full_summary=root / "14run_popdep_n200_bench/run_logs/full_fork_summary.tsv",
            crosschr_summary=root / "14run_popdep_n200_bench/run_logs/crosschr_threads_summary.tsv",
        ),
        BenchCampaign(
            campaign_id="head2head_10run",
            label="Head2head 10run (N20)",
            jar="TIGER_PD_20260616",
            bench_dir=root / "10run_popdep_head2head",
            head2head_summary=root / "10run_popdep_head2head/run_logs/head2head_summary.tsv",
        ),
        BenchCampaign(
            campaign_id="crosschr_threads_09run",
            label="CrossChr threads 09run (60 min partial)",
            jar="TIGER_PD_20260616",
            bench_dir=root / "09run_popdep_crosschr_threads_bench",
            crosschr_threads_summary=root / "09run_popdep_crosschr_threads_bench/run_logs/bench_summary.tsv",
        ),
    ]


def _jar_tag(jar: str) -> str:
    m = re.search(r"(\d{8})", jar)
    return m.group(1) if m else jar


def _parse_threads_from_run_id(run_id: str) -> int | None:
    m = re.search(r"(?:^|_)(?:e|t)(\d+)(?:_|$)", run_id)
    if m:
        return int(m.group(1))
    m = re.search(r"fork(\d+)_t(\d+)", run_id)
    if m:
        return int(m.group(2))
    return None


def _safe_float(val) -> float | None:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return None
    s = str(val).strip()
    if not s or s.upper() == "NA":
        return None
    try:
        return float(s)
    except ValueError:
        return None


def load_full_fork_summary(path: Path, campaign: BenchCampaign) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    df = load_df_generic(str(path))
    rows = []
    for _, r in df.iterrows():
        run_id = str(r.get("run_id", ""))
        threads = _safe_float(r.get("threads"))
        forks = _safe_float(r.get("max_forks"))
        rows.append(
            {
                "campaign_id": campaign.campaign_id,
                "campaign_label": campaign.label,
                "jar": campaign.jar,
                "jar_tag": _jar_tag(campaign.jar),
                "run_id": run_id,
                "app": str(r.get("app", "PopDepFull")),
                "threads": int(threads) if threads is not None else _parse_threads_from_run_id(run_id),
                "max_forks": int(forks) if forks is not None else None,
                "config_label": f"fork{int(forks)}×{int(threads)}"
                if forks is not None and threads is not None
                else run_id,
                "elapsed_min": _safe_float(r.get("elapsed_min")),
                "taxa_done": _safe_float(r.get("taxa_done")),
                "taxa_per_hour": _safe_float(r.get("taxa_per_hour")),
                "samtools_cpu_avg": None,
                "sdb_util_avg": _safe_float(r.get("sdb_util_avg")),
                "sdb_rMB_s_avg": _safe_float(r.get("sdb_rMB_s_avg")),
                "sdb_wMB_s_avg": _safe_float(r.get("sdb_wMB_s_avg")),
                "sdb_aqu_sz_avg": _safe_float(r.get("sdb_aqu_sz_avg")),
                "note": str(r.get("note", "")),
            }
        )
    return pd.DataFrame(rows)


def load_crosschr_summary(path: Path, campaign: BenchCampaign) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    df = load_df_generic(str(path))
    rows = []
    for _, r in df.iterrows():
        run_id = str(r.get("run_id", f"crosschr_e{r.get('threads', '')}"))
        threads = _safe_float(r.get("threads"))
        if threads is None:
            threads = _parse_threads_from_run_id(run_id)
        rows.append(
            {
                "campaign_id": campaign.campaign_id,
                "campaign_label": campaign.label,
                "jar": r.get("jar", campaign.jar),
                "jar_tag": _jar_tag(str(r.get("jar", campaign.jar))),
                "run_id": run_id,
                "app": "PopDepCrossChr",
                "threads": int(threads) if threads is not None else None,
                "max_forks": 1,
                "config_label": f"e{int(threads)}" if threads is not None else run_id,
                "elapsed_min": _safe_float(r.get("elapsed_min")),
                "taxa_done": _safe_float(r.get("taxa_done")),
                "taxa_per_hour": _safe_float(r.get("taxa_per_hour")),
                "samtools_cpu_avg": _safe_float(r.get("samtools_cpu_avg")),
                "sdb_util_avg": _safe_float(r.get("sdb_util_avg")),
                "sdb_rMB_s_avg": _safe_float(r.get("sdb_rMB_s_avg")),
                "sdb_wMB_s_avg": _safe_float(r.get("sdb_wMB_s_avg")),
                "sdb_aqu_sz_avg": _safe_float(r.get("sdb_aqu_sz_avg")),
                "note": str(r.get("note", "")),
            }
        )
    return pd.DataFrame(rows)


def load_head2head_summary(path: Path, campaign: BenchCampaign) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    df = load_df_generic(str(path))
    rows = []
    for _, r in df.iterrows():
        run_id = str(r.get("run_id", ""))
        app = str(r.get("app", ""))
        threads_or_forks = str(r.get("threads_or_forks", ""))
        threads, forks = None, None
        if "/" in threads_or_forks:
            parts = threads_or_forks.split("/")
            if len(parts) == 2:
                threads, forks = int(parts[0]), int(parts[1])
        rows.append(
            {
                "campaign_id": campaign.campaign_id,
                "campaign_label": campaign.label,
                "jar": campaign.jar,
                "jar_tag": _jar_tag(campaign.jar),
                "run_id": run_id,
                "app": app,
                "tier": str(r.get("tier", "")),
                "taxa": str(r.get("taxa", "")),
                "threads": threads,
                "max_forks": forks,
                "config_label": threads_or_forks,
                "elapsed_min": _safe_float(r.get("elapsed_min")),
                "taxa_done": _safe_float(r.get("taxa_done")),
                "taxa_per_hour": _safe_float(r.get("taxa_per_hour")),
                "samtools_cpu_avg": None,
                "sdb_util_avg": _safe_float(r.get("sdb_util_avg")),
                "sdb_rMB_s_avg": None,
                "sdb_wMB_s_avg": None,
                "sdb_aqu_sz_avg": None,
                "note": str(r.get("note", "")),
            }
        )
    return pd.DataFrame(rows)


def load_campaign_summaries(campaigns: list[BenchCampaign] | None = None) -> pd.DataFrame:
    campaigns = campaigns or default_campaigns()
    parts: list[pd.DataFrame] = []
    for c in campaigns:
        if c.full_summary:
            parts.append(load_full_fork_summary(c.full_summary, c))
        cross_path = c.crosschr_summary or c.crosschr_threads_summary
        if cross_path:
            parts.append(load_crosschr_summary(cross_path, c))
        if c.head2head_summary:
            parts.append(load_head2head_summary(c.head2head_summary, c))
    if not parts:
        return pd.DataFrame()
    non_empty = [p for p in parts if not p.empty]
    out = pd.concat(non_empty, ignore_index=True, sort=False)
    return out.sort_values(["campaign_id", "app", "threads"], na_position="last")


def _normalize_monitor_columns(df: pd.DataFrame) -> pd.DataFrame:
    rename = {}
    if "sdb_rMB" in df.columns and "sdb_rMB_s" not in df.columns:
        rename["sdb_rMB"] = "sdb_rMB_s"
    if "sdb_wMB" in df.columns and "sdb_wMB_s" not in df.columns:
        rename["sdb_wMB"] = "sdb_wMB_s"
    df = df.rename(columns=rename)
    if "phase" not in df.columns:
        df["phase"] = "unknown"
    for col in MONITOR_NUMERIC:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    if "timestamp" in df.columns:
        df["timestamp"] = pd.to_datetime(df["timestamp"], errors="coerce")
        t0 = df["timestamp"].min()
        df["elapsed_min"] = (df["timestamp"] - t0).dt.total_seconds() / 60.0
    return df


def load_monitor_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        df = load_df_generic(str(path))
    except Exception as exc:
        print(f"WARN: skipping malformed monitor TSV {path}: {exc}")
        return pd.DataFrame()
    if df.empty or len(df.columns) < 3:
        return pd.DataFrame()
    return _normalize_monitor_columns(df)


def discover_monitor_paths(campaign: BenchCampaign) -> list[tuple[str, str, Path]]:
    """Return (run_id, app, monitor_path) under a campaign directory."""
    found: list[tuple[str, str, Path]] = []
    bench = campaign.bench_dir
    if not bench.exists():
        return found
    # 09run summaries use synthetic run_ids; raw monitor TSVs are often malformed.
    if campaign.crosschr_threads_summary is not None:
        return found

    seen: set[Path] = set()
    for mon in sorted(bench.glob("*/run_logs/monitor.tsv")):
        if mon in seen:
            continue
        seen.add(mon)
        run_id = mon.parent.parent.name
        app = "PopDepFull" if run_id.startswith("full_") else "PopDepCrossChr"
        if run_id.startswith("t") and "_full_" in run_id:
            app = "PopDepFull"
        elif run_id.startswith("t") and "_crosschr_" in run_id:
            app = "PopDepCrossChr"
        found.append((run_id, app, mon))

    return found


def monitor_run_stats(monitor: pd.DataFrame, tiger_only: bool = True) -> dict:
    if monitor.empty:
        return {}
    df = monitor
    if tiger_only and "phase" in df.columns:
        df = df[df["phase"] == "tiger"]
        if df.empty:
            df = monitor
    stats = {}
    for col in ("samtools_cpu_avg", "samtools_n", "sdb_util", "sdb_rMB_s", "sdb_wMB_s", "cpu_iowait"):
        if col not in df.columns:
            continue
        s = df[col].dropna()
        if s.empty:
            continue
        stats[f"{col}_mean"] = float(s.mean())
        stats[f"{col}_std"] = float(s.std(ddof=0))
        stats[f"{col}_cv"] = float(s.std(ddof=0) / s.mean()) if s.mean() else None
        stats[f"{col}_max"] = float(s.max())
    return stats


def enrich_summary_with_monitor(summary: pd.DataFrame, campaigns: list[BenchCampaign] | None = None) -> pd.DataFrame:
    campaigns = campaigns or default_campaigns()
    if summary.empty:
        return summary
    out = summary.copy()
    monitor_map: dict[tuple[str, str], dict] = {}
    for c in campaigns:
        for run_id, app, mon_path in discover_monitor_paths(c):
            mon = load_monitor_tsv(mon_path)
            stats = monitor_run_stats(mon)
            if stats:
                monitor_map[(c.campaign_id, run_id)] = stats

    for idx, row in out.iterrows():
        key = (row["campaign_id"], row["run_id"])
        stats = monitor_map.get(key, {})
        if not stats:
            continue
        if pd.isna(row.get("samtools_cpu_avg")) and "samtools_cpu_avg_mean" in stats:
            out.at[idx, "samtools_cpu_avg"] = stats["samtools_cpu_avg_mean"]
        for src, dst in (
            ("sdb_util_mean", "sdb_util_avg"),
            ("sdb_rMB_s_mean", "sdb_rMB_s_avg"),
            ("sdb_wMB_s_mean", "sdb_wMB_s_avg"),
        ):
            if dst in out.columns and (pd.isna(row.get(dst)) or row.get(dst) is None) and src in stats:
                out.at[idx, dst] = stats[src]
        out.at[idx, "samtools_cpu_cv"] = stats.get("samtools_cpu_avg_cv")
        out.at[idx, "samtools_n_mean"] = stats.get("samtools_n_mean")
    return out


def _save_fig(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved plot to {path}")
    plt.close()


def plot_crosschr_throughput_by_jar(summary: pd.DataFrame, out_dir: Path) -> None:
    df = summary[
        (summary["app"] == "PopDepCrossChr")
        & summary["campaign_id"].str.startswith("n200_")
        & summary["threads"].notna()
        & summary["taxa_per_hour"].notna()
        & ~summary["run_id"].str.contains("shuffle|ckpt", regex=True)
    ].copy()
    if df.empty:
        return
    df = df.sort_values(["jar_tag", "threads"])
    plot_data = []
    for jar_tag, grp in df.groupby("jar_tag"):
        plot_data.append({"series": jar_tag, "frame": grp})
    sns.set_style("white")
    plt.figure(figsize=(10, 6))
    colors = {"20260616": "steelblue", "20260619": "darkorange"}
    for jar_tag, grp in df.groupby("jar_tag"):
        plt.plot(
            grp["threads"],
            grp["taxa_per_hour"],
            marker="o",
            linewidth=1.8,
            label=jar_tag,
            color=colors.get(jar_tag, "gray"),
        )
    plt.title("PopDepCrossChr throughput vs samtools threads (N200 truncated)", fontsize=TITLE_FONT_SIZE)
    plt.xlabel("TIGER -e (samtools threads)", fontsize=X_LABEL_FONT_SIZE)
    plt.ylabel("Taxa per hour", fontsize=Y_LABEL_FONT_SIZE)
    plt.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    plt.legend(loc="upper left", bbox_to_anchor=(0.0, -0.15), fontsize=LEGEND_FONT_SIZE, frameon=False)
    _save_fig(out_dir / "crosschr_throughput_vs_threads_by_jar.png")


def plot_crosschr_samtools_cpu(summary: pd.DataFrame, out_dir: Path) -> None:
    df = summary[
        (summary["app"] == "PopDepCrossChr")
        & summary["campaign_id"].str.startswith("n200_")
        & summary["threads"].notna()
        & summary["samtools_cpu_avg"].notna()
        & ~summary["run_id"].str.contains("shuffle|ckpt", regex=True)
    ].copy()
    if df.empty:
        return
    sns.set_style("white")
    plt.figure(figsize=(10, 6))
    colors = {"20260616": "steelblue", "20260619": "darkorange"}
    for jar_tag, grp in df.groupby("jar_tag"):
        grp = grp.sort_values("threads")
        plt.plot(
            grp["threads"],
            grp["samtools_cpu_avg"],
            marker="s",
            linewidth=1.8,
            label=jar_tag,
            color=colors.get(jar_tag, "gray"),
        )
    plt.title("Mean samtools CPU vs thread count (PopDepCrossChr, N200)", fontsize=TITLE_FONT_SIZE)
    plt.xlabel("TIGER -e (samtools threads)", fontsize=X_LABEL_FONT_SIZE)
    plt.ylabel("Mean samtools CPU (%)", fontsize=Y_LABEL_FONT_SIZE)
    plt.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    plt.legend(loc="upper left", bbox_to_anchor=(0.0, -0.15), fontsize=LEGEND_FONT_SIZE, frameon=False)
    _save_fig(out_dir / "crosschr_samtools_cpu_vs_threads_by_jar.png")


def plot_crosschr_io(summary: pd.DataFrame, out_dir: Path) -> None:
    df = summary[
        (summary["app"] == "PopDepCrossChr")
        & summary["campaign_id"].str.startswith("n200_")
        & summary["threads"].notna()
        & ~summary["run_id"].str.contains("shuffle|ckpt", regex=True)
    ].copy()
    if df.empty:
        return
    y_specs = [
        {"y_col": "sdb_util_avg", "label": "Disk %util", "color": "steelblue"},
        {"y_col": "sdb_rMB_s_avg", "label": "Read MB/s", "color": "darkgreen", "linestyle": "--"},
        {"y_col": "sdb_wMB_s_avg", "label": "Write MB/s", "color": "darkorange", "linestyle": ":"},
    ]
    for jar_tag, grp in df.groupby("jar_tag"):
        grp = grp.sort_values("threads")
        usable = [s for s in y_specs if s["y_col"] in grp.columns and grp[s["y_col"]].notna().any()]
        if not usable:
            continue
        plot_multi_line_series(
            grp,
            "threads",
            usable,
            title=f"PopDepCrossChr disk I/O vs threads ({jar_tag})",
            filename=str(out_dir / f"crosschr_io_vs_threads_{jar_tag}.png"),
            x_label="TIGER -e (samtools threads)",
            y_label="I/O metric",
        )


def plot_full_vs_crosschr_n200(summary: pd.DataFrame, out_dir: Path) -> None:
    picks = [
        ("n200_14run", "full_fork2_t32", "PopDepFull fork2×32"),
        ("n200_14run", "full_fork1_t64", "PopDepFull fork1×64"),
        ("n200_14run", "crosschr_e16", "CrossChr e16"),
        ("n200_14run", "crosschr_e64", "CrossChr e64"),
    ]
    rows = []
    for camp, run_id, label in picks:
        hit = summary[(summary["campaign_id"] == camp) & (summary["run_id"] == run_id)]
        if hit.empty:
            continue
        r = hit.iloc[0]
        rows.append({"label": label, "taxa_per_hour": r["taxa_per_hour"]})
    if not rows:
        return
    df = pd.DataFrame(rows)
    plot_bar_chart(
        df["label"].tolist(),
        df["taxa_per_hour"].tolist(),
        title="N200 throughput: PopDepFull vs PopDepCrossChr (20260619 jar)",
        ylabel="Taxa per hour",
        filename=str(out_dir / "n200_full_vs_crosschr_throughput.png"),
        ylim=(0, max(df["taxa_per_hour"]) * 1.15),
        rotate_xlabels=20,
    )


def plot_head2head_throughput(summary: pd.DataFrame, out_dir: Path) -> None:
    df = summary[summary["campaign_id"] == "head2head_10run"].copy()
    df = df[df["tier"] == "T2"]
    if df.empty:
        return
    labels = [f"{r['app']}\n({r['config_label']})" for _, r in df.iterrows()]
    values = df["taxa_per_hour"].tolist()
    plot_bar_chart(
        labels,
        values,
        title="Head2head throughput (N20, 10 Mb × 44 chr, 20260616 jar)",
        ylabel="Taxa per hour",
        filename=str(out_dir / "head2head_n20_throughput.png"),
        ylim=(0, max(values) * 1.15),
    )


def _prepare_h2h_1m96_frame(monitor: pd.DataFrame) -> pd.DataFrame:
    """Filter tiger phase with active disk read, add elapsed_sec from first sample."""
    if monitor.empty:
        return monitor
    df = monitor.copy()
    if "phase" in df.columns:
        df = df[df["phase"] == "tiger"]
    if df.empty:
        return df
    if "sdb_rMB_s" in df.columns:
        active = df[df["sdb_rMB_s"] > H2H_1M96_ACTIVE_SDB_RMB_S]
        if not active.empty:
            df = active
    if "timestamp" not in df.columns:
        return df
    df = df.sort_values("timestamp")
    t0 = df["timestamp"].min()
    df["elapsed_sec"] = (df["timestamp"] - t0).dt.total_seconds()
    return df.reset_index(drop=True)


def load_h2h_1m96_series(data_dir: Path | None = None) -> list[tuple[str, pd.DataFrame, str]]:
    """Return (label, frame, color) for each PopDep app in the 1M×96 head2head bench."""
    root = Path(data_dir or H2H_1M96_DATA_DIR)
    out: list[tuple[str, pd.DataFrame, str]] = []
    for label, fname, color in H2H_1M96_SERIES:
        mon = load_monitor_tsv(root / fname)
        frame = _prepare_h2h_1m96_frame(mon)
        if frame.empty or "elapsed_sec" not in frame.columns:
            print(f"WARN: skipping empty h2h_1M96 series {label} ({root / fname})")
            continue
        out.append((label, frame, color))
    return out


def build_h2h_1m96_summary(series: list[tuple[str, pd.DataFrame, str]]) -> pd.DataFrame:
    rows = []
    for label, df, _color in series:
        elapsed = float(df["elapsed_sec"].max())
        rows.append(
            {
                "campaign_id": "head2head_1M96",
                "campaign_label": "Head2head 1M chr × 96 taxa (Jun 2026)",
                "jar": "TIGER_dazheng.jar",
                "app": label,
                "threads": 16,
                "elapsed_sec": elapsed,
                "n_monitor_tiger": len(df),
                "samtools_cpu_avg": float(df["samtools_cpu_avg"].mean()),
                "samtools_n_mean": float(df["samtools_n"].mean()),
                "sdb_rMB_s_avg": float(df["sdb_rMB_s"].mean()) if "sdb_rMB_s" in df.columns else None,
                "sdb_wMB_s_avg": float(df["sdb_wMB_s"].mean()) if "sdb_wMB_s" in df.columns else None,
                "sdb_aqu_sz_avg": float(df["sdb_aqu_sz"].mean()) if "sdb_aqu_sz" in df.columns else None,
                "note": f"active=tiger & sdb_rMB_s>{H2H_1M96_ACTIVE_SDB_RMB_S:g}; shaded timeseries",
            }
        )
    return pd.DataFrame(rows)


def _plot_h2h_filled_timeseries(
    series: list[tuple[str, pd.DataFrame, str]],
    y_col: str,
    title: str,
    ylabel: str,
    out_path: Path,
) -> None:
    """Line + shaded area under each series (longest run drawn first)."""
    usable = [
        (label, df, color)
        for label, df, color in series
        if y_col in df.columns and df[y_col].notna().any()
    ]
    if not usable:
        return

    usable.sort(key=lambda item: item[1]["elapsed_sec"].max(), reverse=True)

    sns.set_style("white")
    fig, ax = plt.subplots(figsize=(11, 6))
    for label, df, color in usable:
        x = df["elapsed_sec"].to_numpy()
        y = df[y_col].to_numpy()
        ax.fill_between(
            x,
            y,
            0,
            color=color,
            alpha=H2H_1M96_FILL_ALPHA,
            linewidth=0,
            zorder=1,
        )
        ax.plot(
            x,
            y,
            color=color,
            linewidth=1.6,
            label=label,
            alpha=0.95,
            zorder=2,
        )

    ax.set_title(title, fontsize=TITLE_FONT_SIZE)
    ax.set_xlabel("Elapsed (s)", fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel(ylabel, fontsize=Y_LABEL_FONT_SIZE)
    ax.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.18), fontsize=LEGEND_FONT_SIZE, frameon=False)
    fig.subplots_adjust(bottom=0.22)
    _save_fig(out_path)


def plot_h2h_1m96_timeseries(
    out_dir: Path,
    info_dir: Path | None = None,
    data_dir: Path | None = None,
) -> pd.DataFrame | None:
    """PopDep vs PopDepFull vs PopDepCrossChr samtools CPU/count time series (1M×96)."""
    series = load_h2h_1m96_series(data_dir)
    if not series:
        return None

    summary = build_h2h_1m96_summary(series)
    if info_dir is not None:
        save_df_to_tsv(summary, str(info_dir / "h2h_1M96_summary.tsv"))

    prefix = "h2h_1M96"
    phase_note = "phase=tiger"
    _plot_h2h_filled_timeseries(
        series,
        "samtools_n",
        f"Concurrent samtools processes ({phase_note})",
        "samtools count",
        out_dir / f"{prefix}_samtools_n_timeseries.png",
    )
    _plot_h2h_filled_timeseries(
        series,
        "samtools_cpu_avg",
        f"Mean samtools CPU per process ({phase_note})",
        "samtools CPU avg (%)",
        out_dir / f"{prefix}_samtools_cpu_timeseries.png",
    )
    return summary


def plot_monitor_timeseries(
    campaigns: list[BenchCampaign],
    out_dir: Path,
    run_specs: list[tuple[str, str, str]],
) -> None:
    """Plot samtools CPU / count / disk util over elapsed time for selected runs."""
    series_frames: list[tuple[str, pd.DataFrame]] = []
    for campaign_id, run_id, label in run_specs:
        camp = next((c for c in campaigns if c.campaign_id == campaign_id), None)
        if camp is None:
            continue
        mon_path = camp.bench_dir / run_id / "run_logs/monitor.tsv"
        mon = load_monitor_tsv(mon_path)
        if mon.empty or "elapsed_min" not in mon.columns:
            continue
        tiger = mon[mon["phase"] == "tiger"] if "phase" in mon.columns else mon
        if tiger.empty:
            tiger = mon
        tiger = tiger.copy()
        tiger["series"] = label
        series_frames.append((label, tiger))

    if not series_frames:
        return

    metrics = [
        ("samtools_cpu_avg", "Samtools mean CPU (%)", "monitor_samtools_cpu_timeseries.png"),
        ("samtools_n", "Concurrent samtools processes", "monitor_samtools_count_timeseries.png"),
        ("sdb_util", "Disk sdb %util", "monitor_disk_util_timeseries.png"),
    ]
    sns.set_style("white")
    for col, ylabel, fname in metrics:
        plt.figure(figsize=(11, 6))
        has_data = False
        for label, df in series_frames:
            if col not in df.columns or df[col].isna().all():
                continue
            has_data = True
            plt.plot(df["elapsed_min"], df[col], linewidth=1.4, label=label, alpha=0.9)
        if not has_data:
            plt.close()
            continue
        plt.title(f"Runtime stability — {ylabel}", fontsize=TITLE_FONT_SIZE)
        plt.xlabel("Elapsed time (min)", fontsize=X_LABEL_FONT_SIZE)
        plt.ylabel(ylabel, fontsize=Y_LABEL_FONT_SIZE)
        plt.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
        plt.legend(loc="upper left", bbox_to_anchor=(0.0, -0.18), fontsize=LEGEND_FONT_SIZE, frameon=False)
        _save_fig(out_dir / fname)


def _boxplot_run_label(jar_tag: str, app: str, run_id: str) -> str:
    app_short = app.replace("PopDep", "")
    return f"{jar_tag} | {app_short} | {run_id}"


def plot_samtools_cpu_boxplot(campaigns: list[BenchCampaign], out_dir: Path) -> None:
    records = []
    for c in campaigns:
        if not c.campaign_id.startswith("n200_"):
            continue
        for run_id, app, mon_path in discover_monitor_paths(c):
            if "shuffle" in run_id or "ckpt" in run_id:
                continue
            mon = load_monitor_tsv(mon_path)
            if mon.empty:
                continue
            tiger = mon[mon["phase"] == "tiger"] if "phase" in mon.columns else mon
            if tiger.empty or "samtools_cpu_avg" not in tiger.columns:
                continue
            label = _boxplot_run_label(_jar_tag(c.jar), app, run_id)
            for val in tiger["samtools_cpu_avg"].dropna():
                records.append(
                    {
                        "label": label,
                        "samtools_cpu_avg": val,
                        "app": app,
                        "jar_tag": _jar_tag(c.jar),
                        "run_id": run_id,
                    }
                )
    if not records:
        return
    df = pd.DataFrame(records)
    order = sorted(
        df["label"].unique(),
        key=lambda x: (x.split(" | ")[0], x.split(" | ")[1], x.split(" | ")[2]),
    )
    n_cats = len(order)
    fig_w = max(22.0, n_cats * 1.55)
    sns.set_style("white")
    fig, ax = plt.subplots(figsize=(fig_w, 6.5))
    sns.boxplot(
        data=df,
        x="label",
        y="samtools_cpu_avg",
        order=order,
        color="steelblue",
        fliersize=2,
        ax=ax,
    )
    ax.set_title(
        "Samtools CPU distribution during TIGER phase (N200 bench)",
        fontsize=TITLE_FONT_SIZE,
    )
    ax.set_xlabel("Run configuration", fontsize=X_LABEL_FONT_SIZE)
    ax.set_ylabel("Samtools CPU (%) per sample", fontsize=Y_LABEL_FONT_SIZE)
    ax.tick_params(axis="x", rotation=40, labelsize=TICK_FONT_SIZE - 1)
    ax.tick_params(axis="y", labelsize=TICK_FONT_SIZE)
    for tick in ax.get_xticklabels():
        tick.set_ha("right")
    fig.subplots_adjust(bottom=0.32)
    _save_fig(out_dir / "samtools_cpu_stability_boxplot.png")


def plot_09run_io_saturation(summary: pd.DataFrame, out_dir: Path) -> None:
    """Production BAM 60-min partial sweep: samtools CPU vs disk util vs threads."""
    df = summary[
        (summary["campaign_id"] == "crosschr_threads_09run")
        & summary["note"].astype(str).str.contains("ok", na=False)
        & summary["threads"].notna()
    ].copy()
    if df.empty:
        return
    df = df.sort_values("threads")
    sns.set_style("white")
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax2 = ax1.twinx()
    ax1.plot(df["threads"], df["samtools_cpu_avg"], "o-", color="steelblue", linewidth=1.8, label="samtools CPU %")
    ax2.plot(df["threads"], df["sdb_util_avg"], "s--", color="darkorange", linewidth=1.8, label="sdb %util")
    ax1.set_xlabel("TIGER -e (samtools threads)", fontsize=X_LABEL_FONT_SIZE)
    ax1.set_ylabel("Mean samtools CPU (%)", fontsize=Y_LABEL_FONT_SIZE, color="steelblue")
    ax2.set_ylabel("Mean disk %util (sdb)", fontsize=Y_LABEL_FONT_SIZE, color="darkorange")
    ax1.tick_params(axis="y", labelcolor="steelblue", labelsize=TICK_FONT_SIZE)
    ax2.tick_params(axis="y", labelcolor="darkorange", labelsize=TICK_FONT_SIZE)
    ax1.tick_params(axis="x", labelsize=TICK_FONT_SIZE)
    ax1.set_title(
        "PopDepCrossChr I/O saturation (production BAM, 60 min partial, 20260616)",
        fontsize=TITLE_FONT_SIZE,
    )
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left", bbox_to_anchor=(0.0, -0.15), fontsize=LEGEND_FONT_SIZE, frameon=False)
    _save_fig(out_dir / "crosschr_09run_io_saturation_vs_threads.png")


def plot_n200_app_io_compare(summary: pd.DataFrame, out_dir: Path) -> None:
    """Compare disk util and samtools CPU for recommended Full vs CrossChr configs (14run)."""
    picks = [
        ("full_fork2_t32", "PopDepFull fork2×32"),
        ("crosschr_e16", "CrossChr e16"),
        ("crosschr_e64", "CrossChr e64"),
    ]
    rows = []
    for run_id, label in picks:
        hit = summary[(summary["campaign_id"] == "n200_14run") & (summary["run_id"] == run_id)]
        if hit.empty:
            continue
        r = hit.iloc[0]
        rows.append({"label": label, "sdb_util_avg": r.get("sdb_util_avg"), "samtools_cpu_avg": r.get("samtools_cpu_avg")})
    if not rows:
        return
    df = pd.DataFrame(rows)
    x = np.arange(len(df))
    width = 0.35
    sns.set_style("white")
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(x - width / 2, df["samtools_cpu_avg"], width, label="samtools CPU %", color="steelblue", alpha=0.85)
    ax.bar(x + width / 2, df["sdb_util_avg"], width, label="disk %util", color="darkorange", alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(df["label"], rotation=15, ha="right", fontsize=TICK_FONT_SIZE)
    ax.set_ylabel("Metric value", fontsize=Y_LABEL_FONT_SIZE)
    ax.set_title("N200 CPU vs disk util — PopDepFull vs CrossChr (20260619)", fontsize=TITLE_FONT_SIZE)
    ax.tick_params(axis="y", labelsize=TICK_FONT_SIZE)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.12), fontsize=LEGEND_FONT_SIZE, frameon=False)
    _save_fig(out_dir / "n200_cpu_disk_util_compare.png")


def plot_jar_evolution_elapsed(summary: pd.DataFrame, out_dir: Path) -> None:
    df = summary[
        (summary["app"] == "PopDepCrossChr")
        & summary["campaign_id"].str.startswith("n200_")
        & (summary["run_id"] == "crosschr_e64")
    ].copy()
    if len(df) < 2:
        return
    labels = [f"{r['jar_tag']}" for _, r in df.iterrows()]
    values = df["elapsed_min"].tolist()
    plot_bar_chart(
        labels,
        values,
        title="PopDepCrossChr e64 wall time — jar evolution (N200 truncated)",
        ylabel="Elapsed minutes",
        filename=str(out_dir / "crosschr_e64_elapsed_by_jar.png"),
        ylim=(0, max(values) * 1.15),
    )


def run_popdep_bench_plots(
    output_dir: str | Path,
    stats_root: str | Path = VMAP4_STATS_GENOME,
    h2h_1m96_dir: str | Path | None = None,
) -> pd.DataFrame:
    """Load bench summaries, enrich with monitor stats, write TSV + PNG figures."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    campaigns = default_campaigns(Path(stats_root))

    summary = load_campaign_summaries(campaigns)
    summary = enrich_summary_with_monitor(summary, campaigns)

    info_dir = out_dir / "info"
    plots_dir = out_dir / "plots"
    info_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    save_df_to_tsv(summary, str(info_dir / "popdep_bench_summary.tsv"))

    plot_crosschr_throughput_by_jar(summary, plots_dir)
    plot_crosschr_samtools_cpu(summary, plots_dir)
    plot_crosschr_io(summary, plots_dir)
    plot_full_vs_crosschr_n200(summary, plots_dir)
    plot_head2head_throughput(summary, plots_dir)
    plot_09run_io_saturation(summary, plots_dir)
    plot_n200_app_io_compare(summary, plots_dir)
    plot_jar_evolution_elapsed(summary, plots_dir)
    plot_samtools_cpu_boxplot(campaigns, plots_dir)

    run_specs = [
        ("n200_13run", "crosschr_e64", "13run CrossChr e64"),
        ("n200_14run", "crosschr_e64", "14run CrossChr e64"),
        ("n200_14run", "full_fork2_t32", "14run Full fork2×32"),
        ("n200_13run", "full_fork2_t32", "13run Full fork2×32"),
    ]
    plot_monitor_timeseries(campaigns, plots_dir, run_specs)
    plot_h2h_1m96_timeseries(plots_dir, info_dir, h2h_1m96_dir)

    print(f"Wrote {len(list(plots_dir.glob('*.png')))} plots under {plots_dir}")
    return summary


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Plot PopDep benchmark summaries")
    p.add_argument(
        "--output-dir",
        default=str(POPDEP_BENCH_PUBLISH_DIR),
        help="Publish root: writes plots/ and info/ (default: benchmark/popdep_bench)",
    )
    p.add_argument(
        "--stats-root",
        default=str(VMAP4_STATS_GENOME),
        help="Root of 10stats.genome bench run folders",
    )
    p.add_argument(
        "--h2h-1m96-dir",
        default=str(H2H_1M96_DATA_DIR),
        help="Monitor TSV directory for PopDep/Full/CrossChr 1M×96 head2head",
    )
    args = p.parse_args()
    run_popdep_bench_plots(args.output_dir, args.stats_root, args.h2h_1m96_dir)
