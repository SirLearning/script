import os
import numpy as np
import pandas as pd
import seaborn as sns
from infra.utils.io import load_df_from_tsv, save_df_to_tsv, save_thresholds
from infra.utils.graph import (
    plot_distribution_with_stats,
    plot_joint_regression,
    plot_binned_mean_heatmap,
    plot_loess_scatter,
    plot_loess_scatter_strata,
    plot_gam_tensor_surface,
    plot_gam_partial_curve,
    plot_contourf_panels,
    plot_line_panels,
    plot_bar_chart,
    plot_scatter_with_outliers,
    plot_scatter_continuous_color,
    plot_scatter_with_outliers_sized,
    RESIDUAL_SIGNED_CMAP,
    MISSING_RATE_CMAP,
    DEFAULT_AXIS_PADDING_FRACTION,
)
from genetics.germplasm.sample import load_df_from_tbm, anno_group
from .sample_utils import load_df_from_plink2

def coverage_dist(
    input_file="/data/home/tusr1/01projects/vmap4/00data/05taxaBamMap/all.A.taxaBamMap.txt",
    output_prefix="sample.coverage",
    tsv_file=None
):
    """
    Plots histogram of Individual Mean Coverage from .idepth file.
    Input: .idepth file
    Output: Histogram plot
    """
    if tsv_file and os.path.exists(tsv_file):
        df_coverage = load_df_from_tsv(tsv_file)
        if df_coverage is None: return
    else:
        df = load_df_from_tbm(input_file)
        if df is None: return
        # extract Sample and Coverage columns
        df_coverage = df[['Sample', 'Coverage']].copy()
        # Save processed data
        save_df_to_tsv(df_coverage, f"{output_prefix}.info.tsv")

    # get coverage thresholds
    mean_cov = df_coverage['Coverage'].mean()
    median_cov = df_coverage['Coverage'].median()
    std_cov = df_coverage['Coverage'].std()
    cov_threshold = mean_cov + 3 * std_cov
    # save thresholds
    threshold_file = f"{output_prefix}.th.tsv"
    stats_data = {
        'coverage_threshold': cov_threshold,
        'mean_coverage': mean_cov,
        'std_coverage': std_cov
    }
    save_thresholds(stats_data, threshold_file)
    print(f"[Info] Coverage Mean: {mean_cov:.4f}, Std: {std_cov:.4f}, Threshold (Mean-3SD): {cov_threshold:.4f}")
    
    threshold_to_plot = {
        'value': cov_threshold,
        'label': f'Threshold (+3SD) = {cov_threshold:.2f}',
        'color': 'purple',
        'linestyle': '-'
    }

    # Use shared plotting
    plot_distribution_with_stats(
        data=df_coverage, 
        col='Coverage',
        title="Individual Mean Coverage",
        filename=f"{output_prefix}.dist.png",
        x_label="Coverage", 
        y_label="Count",
        mean_val=mean_cov, 
        median_val=median_cov,
        std_val=std_cov,
        thresholds=[threshold_to_plot]
    )
    

def reg_missing_coverage(
    depth_file="/data/home/tusr1/01projects/vmap4/00data/05taxaBamMap/all.A.taxaBamMap.txt",
    miss_file="/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.missing.smiss",
    output_prefix="sample.coverage_vs_missing",
    group_file=None,
    tsv_file=None,
):
    if tsv_file and os.path.exists(tsv_file):
        merged = load_df_from_tsv(tsv_file)
        if merged is None: return
    else:
        # 1. Read Data
        df_coverage = load_df_from_tbm(depth_file)
        df_miss = load_df_from_plink2(miss_file)
        if df_coverage is None or df_miss is None: 
            print("[Error] Failed to load input files.")
            return
        # 2. Merge Data
        if 'Sample' in df_coverage.columns and 'Sample' in df_miss.columns:
            merged = pd.merge(df_coverage, df_miss, left_on='Sample', right_on='Sample')
        else:
            print("[Error] Columns 'Sample' (depth) or 'Sample' (miss) not found.")
            print(f"Coverage Cols: {df_coverage.columns}")
            print(f"Miss Cols: {df_miss.columns}")
            return
        # 3. Print Stats
        print(f"Coverage samples: {len(df_coverage)}")
        print(f"Missing samples: {len(df_miss)}")
        print(f"Merged samples: {len(merged)}")
        if 'Coverage' in df_coverage.columns:
            mean_coverage = df_coverage['Coverage'].mean()
            print(f"Mean Coverage (All): {mean_coverage:.4f}")

    if group_file and 'Group' not in merged.columns:
        merged = anno_group(merged, group_file, save_tsv=False)
    elif 'Group' not in merged.columns:
        merged['Group'] = 'Unknown'

    save_df_to_tsv(merged, f"{output_prefix}.miss.info.tsv")

    # 4. Regression Plot
    if 'Coverage' in merged.columns and 'F_MISS' in merged.columns:
        sns.set_theme(style="ticks")
        plot_joint_regression(
            df=merged,
            x_col='Coverage',
            y_col='F_MISS',
            group_col='Group',
            x_label='Coverage',
            y_label='Missing Rate (F_MISS)',
            y_lim=(0, 1),
            filename=f"{output_prefix}.reg_vs_miss.png",
            title="Missing Rate vs Coverage",
        )

        positive = merged[merged["Coverage"] > 0].copy()
        if len(positive) >= 10:
            positive["log10_Depth"] = np.log10(positive["Coverage"].astype(float))
            plot_joint_regression(
                df=positive,
                x_col="log10_Depth",
                y_col="F_MISS",
                group_col="Group",
                x_label="log10(Depth)",
                y_label="Missing Rate (F_MISS)",
                y_lim=(0, 1),
                filename=f"{output_prefix}.reg_vs_miss.logdepth.png",
                title="Missing Rate vs log10(Depth)",
            )

            inv = positive.copy()
            inv["inv_one_plus_miss"] = 1.0 / (1.0 + inv["F_MISS"].astype(float))
            plot_joint_regression(
                df=inv,
                x_col="log10_Depth",
                y_col="inv_one_plus_miss",
                group_col="Group",
                x_label="log10(Depth)",
                y_label="1 / (1 + Missing Rate)",
                y_lim=(0.5, 1.0),
                filename=f"{output_prefix}.reg_inv_miss.logdepth.png",
                title="1 / (1 + Missing Rate) vs log10(Depth)",
            )
        else:
            print(
                "[Warning] Too few samples with Coverage > 0 for log-depth regression plot."
            )

        # 5. Filter & Print Outliers
        high_depth = merged[merged['Coverage'] > 20]
        print(f"\n========= Samples with Coverage > 20 (Count: {len(high_depth)}) =========")
        print(high_depth[['Sample', 'F_MISS']].to_string(index=False))

        high_miss = merged[merged['F_MISS'] > 0.8]
        print(f"\n========= Samples with Missing Rate > 0.8 (Count: {len(high_miss)}) =========")
        print(high_miss[['Sample', 'F_MISS']].to_string(index=False))
    else:
        print("[Error] Creating plot: Required columns 'Coverage' or 'F_MISS' missing.")


def sample_coverage_stats_bundle(
    depth_file,
    miss_file,
    scount_file,
    output_prefix,
    group_file=None,
):
    """
    Run all sample-level coverage QC outputs for one mod (``*.sample.cov`` or ``*.sample.sg_cov``).

    Includes distribution, missing/IBS regressions, IBS×depth heatmaps, LOESS, and GAM diagnostics.
    """
    sns.set_theme(style="ticks")
    coverage_dist(depth_file, output_prefix=output_prefix)
    reg_missing_coverage(depth_file, miss_file, output_prefix, group_file=group_file)
    reg_coverage_ref_ibs(depth_file, scount_file, output_prefix, group_file=group_file)
    heatmap_ibs_depth_missing(
        depth_file,
        miss_file,
        scount_file,
        output_prefix,
        group_file=group_file,
    )


def reg_coverage_ref_ibs(
    depth_file="/data/home/tusr1/01projects/vmap4/00data/05taxaBamMap/all.A.taxaBamMap.txt",
    scount_file=None,
    output_prefix="sample.coverage",
    group_file=None,
    tsv_file=None,
    save_info=None,
):
    """
    Joint regression of mean BAM coverage vs reference IBS (from PLINK .scount).

    Emits ``{output_prefix}.ref_ibs.info.tsv`` (info) and
    ``{output_prefix}.reg_vs_ref_ibs.png`` (plot). In Nextflow, ``*.info.tsv`` is
    published to ``stats/<mod>/info/`` and ``*.png`` to ``stats/<mod>/plots/``.
    When ``tsv_file`` is supplied for plot-only reruns, info is not rewritten unless
    ``save_info=True``.
    """
    from .ref_ibs import load_and_calculate_ibs

    computed = False
    if tsv_file and os.path.exists(tsv_file):
        merged = load_df_from_tsv(tsv_file)
        if merged is None:
            return
    else:
        if scount_file is None:
            print("[Error] reg_coverage_ref_ibs: scount_file is required when tsv_file is absent.")
            return
        df_coverage = load_df_from_tbm(depth_file)
        df_ibs = load_and_calculate_ibs(scount_file)
        if df_coverage is None or df_ibs is None:
            print("[Error] Failed to load coverage or scount inputs.")
            return
        keep_cov = [c for c in ("Sample", "Coverage") if c in df_coverage.columns]
        merged = pd.merge(df_coverage[keep_cov], df_ibs, on="Sample", how="inner")
        print(f"Coverage samples: {len(df_coverage)}")
        print(f"IBS samples: {len(df_ibs)}")
        print(f"Merged samples: {len(merged)}")
        computed = True

    if group_file and "Group" not in merged.columns:
        merged = anno_group(merged, group_file, save_tsv=False)
    elif "Group" not in merged.columns:
        merged["Group"] = "Unknown"

    if save_info is None:
        save_info = computed
    if save_info:
        save_df_to_tsv(merged, f"{output_prefix}.ref_ibs.info.tsv")

    if "Coverage" in merged.columns and "IBS_Ref" in merged.columns:
        sns.set_theme(style="ticks")
        plot_joint_regression(
            df=merged,
            x_col="Coverage",
            y_col="IBS_Ref",
            group_col="Group",
            x_label="Coverage",
            y_label="IBS with Reference Genome",
            filename=f"{output_prefix}.reg_vs_ref_ibs.png",
            title="Reference IBS vs Coverage",
        )
        print(
            f"Pearson(Coverage, IBS_Ref): "
            f"{merged['Coverage'].corr(merged['IBS_Ref']):.4f}"
        )
    else:
        print("[Error] Creating plot: required columns 'Coverage' or 'IBS_Ref' missing.")


def heatmap_ibs_depth_missing(
    depth_file="/data/home/tusr1/01projects/vmap4/00data/05taxaBamMap/all.A.taxaBamMap.txt",
    miss_file=None,
    scount_file=None,
    output_prefix="sample.cov",
    group_file=None,
    tsv_file=None,
    save_info=None,
):
    """
    2D heatmap: X = reference IBS, Y = mean depth (linear or log10), color = mean missing rate.

    Emits ``{output_prefix}.ibs_depth_miss.info.tsv``,
    heatmaps, LOESS diagnostics, and GAM model comparison / surface plots.
    """
    from .ref_ibs import load_and_calculate_ibs

    computed = False
    if tsv_file and os.path.exists(tsv_file):
        merged = load_df_from_tsv(tsv_file)
        if merged is None:
            return
    else:
        if miss_file is None or scount_file is None:
            print(
                "[Error] heatmap_ibs_depth_missing: miss_file and scount_file are "
                "required when tsv_file is absent."
            )
            return
        df_coverage = load_df_from_tbm(depth_file)
        df_miss = load_df_from_plink2(miss_file)
        df_ibs = load_and_calculate_ibs(scount_file)
        if df_coverage is None or df_miss is None or df_ibs is None:
            print("[Error] Failed to load coverage, missing, or scount inputs.")
            return
        keep_cov = [c for c in ("Sample", "Coverage") if c in df_coverage.columns]
        merged = pd.merge(df_coverage[keep_cov], df_miss, on="Sample", how="inner")
        merged = pd.merge(merged, df_ibs, on="Sample", how="inner")
        print(f"Coverage samples: {len(df_coverage)}")
        print(f"Missing samples: {len(df_miss)}")
        print(f"IBS samples: {len(df_ibs)}")
        print(f"Merged samples: {len(merged)}")
        computed = True

    if group_file and "Group" not in merged.columns:
        merged = anno_group(merged, group_file, save_tsv=False)
    elif "Group" not in merged.columns:
        merged["Group"] = "Unknown"

    if save_info is None:
        save_info = computed
    if save_info:
        save_df_to_tsv(merged, f"{output_prefix}.ibs_depth_miss.info.tsv")

    required = ("Coverage", "IBS_Ref", "F_MISS")
    if not all(col in merged.columns for col in required):
        print(f"[Error] heatmap requires columns {required}; got {list(merged.columns)}")
        return

    sns.set_theme(style="ticks")
    heatmap_kw = dict(
        df=merged,
        x_col="IBS_Ref",
        y_col="Coverage",
        value_col="F_MISS",
        x_label="IBS with Reference Genome",
        cbar_label="Missing Rate",
        cmap=MISSING_RATE_CMAP,
        vmin=0.0,
        vmax=1.0,
        n_x_bins=50,
        n_y_bins=50,
    )
    plot_binned_mean_heatmap(
        **heatmap_kw,
        y_label="Depth",
        title="Missing Rate by Reference IBS and Depth",
        filename=f"{output_prefix}.heatmap_ibs_depth_miss.png",
        y_label_fmt="{:.0f}x",
    )
    plot_binned_mean_heatmap(
        **heatmap_kw,
        y_label="log10(Depth)",
        title="Missing Rate by Reference IBS and log10(Depth)",
        filename=f"{output_prefix}.heatmap_ibs_depth_miss.logdepth.png",
        y_label_fmt="{:.2f}",
        y_log10=True,
    )
    print(
        f"Pearson(Coverage, F_MISS): {merged['Coverage'].corr(merged['F_MISS']):.4f}; "
        f"Pearson(IBS_Ref, F_MISS): {merged['IBS_Ref'].corr(merged['F_MISS']):.4f}"
    )
    work = _prepare_ibs_depth_miss_work(merged)
    _loess_missing_ibs_depth_plots(work, output_prefix)
    _gam_missing_ibs_depth_analysis(work, output_prefix)


def heatmap_thin_miss_common_ibs(
    thin_info_path,
    common_info_path,
    output_prefix="sample.cov",
    save_info=True,
):
    """
    Hybrid heatmap: missing rate from ``test_thin``, IBS from ``test_common_thin``.

    Depth (``Coverage``) is sample-level and identical across both jobs; only
    ``F_MISS`` and ``IBS_Ref`` are taken from their respective info TSVs.
    When ``output_prefix`` contains ``sg_cov``, axis labels reflect subgenome
    mosdepth depth. Emits a merged info TSV and a log10(Depth) binned heatmap.
    """
    df_thin = load_df_from_tsv(thin_info_path)
    df_common = load_df_from_tsv(common_info_path)
    if df_thin is None or df_common is None:
        return

    thin_cols = ["Sample", "Coverage", "F_MISS"]
    if "Group" in df_thin.columns:
        thin_cols.append("Group")
    merged = pd.merge(
        df_thin[thin_cols],
        df_common[["Sample", "IBS_Ref"]],
        on="Sample",
        how="inner",
    )
    if "Group" not in merged.columns:
        merged["Group"] = "Unknown"

    print(
        f"Hybrid merge: thin={len(df_thin)}, common={len(df_common)}, "
        f"merged={len(merged)}"
    )
    cov_diff = (
        df_thin.set_index("Sample")["Coverage"]
        .reindex(merged["Sample"])
        .values
        - merged["Coverage"].values
    )
    if np.max(np.abs(cov_diff)) > 1e-9:
        print("[Warning] Coverage mismatch between thin and common inputs.")

    if save_info:
        save_df_to_tsv(merged, f"{output_prefix}.hybrid_thinmiss_commonibs.info.tsv")

    required = ("Coverage", "IBS_Ref", "F_MISS")
    if not all(col in merged.columns for col in required):
        print(f"[Error] hybrid heatmap requires {required}; got {list(merged.columns)}")
        return

    is_sg_cov = "sg_cov" in str(output_prefix)
    depth_y_label = "log10(Subgenome Depth)" if is_sg_cov else "log10(Depth)"
    depth_title = (
        "log10(Subgenome Depth)" if is_sg_cov else "log10(Depth)"
    )

    sns.set_theme(style="ticks")
    plot_binned_mean_heatmap(
        df=merged,
        x_col="IBS_Ref",
        y_col="Coverage",
        value_col="F_MISS",
        x_label="IBS with Reference Genome (test_common_thin)",
        y_label=depth_y_label,
        cbar_label="Missing Rate (test_thin)",
        title=f"Missing Rate (test_thin) by IBS (test_common_thin) and {depth_title}",
        filename=f"{output_prefix}.heatmap_thinmiss_commonibs.logdepth.png",
        cmap=MISSING_RATE_CMAP,
        vmin=0.0,
        vmax=1.0,
        n_x_bins=50,
        n_y_bins=50,
        y_label_fmt="{:.2f}",
        y_log10=True,
    )
    print(
        f"Pearson(Coverage, F_MISS): {merged['Coverage'].corr(merged['F_MISS']):.4f}; "
        f"Pearson(IBS_Ref, F_MISS): {merged['IBS_Ref'].corr(merged['F_MISS']):.4f}"
    )


def _prepare_ibs_depth_miss_work(merged):
    """Positive-depth subset with log10(Depth) for LOESS/GAM diagnostics."""
    work = merged.copy()
    work = work[work["Coverage"] > 0]
    work["log10_Depth"] = np.log10(work["Coverage"].astype(float))
    return work


def _loess_missing_ibs_depth_plots(work, output_prefix, loess_frac=0.2):
    """LOESS diagnostics: Missing ~ IBS, Missing ~ log10(Depth), and IBS-stratified depth curves."""
    loess_kw = dict(y_col="F_MISS", y_label="Missing Rate", y_lim=(0, 1), frac=loess_frac)

    plot_loess_scatter(
        df=work,
        x_col="IBS_Ref",
        x_label="IBS with Reference Genome",
        title="Missing Rate vs Reference IBS (LOESS)",
        filename=f"{output_prefix}.loess_miss_vs_ibs.png",
        **loess_kw,
    )
    plot_loess_scatter(
        df=work,
        x_col="log10_Depth",
        x_label="log10(Depth)",
        title="Missing Rate vs log10(Depth) (LOESS)",
        filename=f"{output_prefix}.loess_miss_vs_logdepth.png",
        **loess_kw,
    )

    ibs = work["IBS_Ref"]
    strata = [
        {"label": "IBS > 0.99", "mask": ibs > 0.99},
        {"label": "0.97 < IBS < 0.99", "mask": (ibs > 0.97) & (ibs < 0.99)},
        {"label": "IBS < 0.97", "mask": ibs < 0.97},
    ]
    plot_loess_scatter_strata(
        df=work,
        x_col="log10_Depth",
        y_col="F_MISS",
        strata=strata,
        x_label="log10(Depth)",
        y_label="Missing Rate",
        title="Missing Rate vs log10(Depth) by IBS stratum (LOESS)",
        filename=f"{output_prefix}.loess_miss_vs_logdepth.by_ibs_strata.png",
        y_lim=(0, 1),
        frac=loess_frac,
    )


def _fit_gam_te_bundle(work, n_splines=25, te_splines=20, n_grid=100):
    """
    Fit additive and tensor GAMs; return te surfaces, partial curves, and summary metrics.
    """
    from pygam import LinearGAM, s, te

    if len(work) < 50:
        return None

    x = work[["IBS_Ref", "log10_Depth"]].to_numpy(dtype=float)
    y = work["F_MISS"].to_numpy(dtype=float)
    mask = np.isfinite(x).all(axis=1) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(y) < 50:
        return None

    gam_add = LinearGAM(s(0, n_splines=n_splines) + s(1, n_splines=n_splines)).gridsearch(
        x, y, progress=False
    )
    gam_te = LinearGAM(te(0, 1, n_splines=te_splines)).gridsearch(x, y, progress=False)
    delta_aic = float(gam_te.statistics_["AIC"] - gam_add.statistics_["AIC"])
    pseudo = gam_te.statistics_.get("pseudo_r2", {})
    pseudo_r2 = pseudo.get("explained_deviance", np.nan) if isinstance(pseudo, dict) else np.nan

    x_lo, x_hi = float(x[:, 0].min()), float(x[:, 0].max())
    y_lo, y_hi = float(x[:, 1].min()), float(x[:, 1].max())
    x_pad = (x_hi - x_lo) * DEFAULT_AXIS_PADDING_FRACTION or 0.01
    y_pad = (y_hi - y_lo) * DEFAULT_AXIS_PADDING_FRACTION or 0.01
    ibs_grid = np.linspace(x_lo - x_pad, x_hi + x_pad, n_grid)
    logd_grid = np.linspace(y_lo - y_pad, y_hi + y_pad, n_grid)
    ibs_mesh, logd_mesh = np.meshgrid(ibs_grid, logd_grid)
    z_pred = np.clip(
        gam_te.predict(np.column_stack([ibs_mesh.ravel(), logd_mesh.ravel()])),
        0.0,
        1.0,
    ).reshape(ibs_mesh.shape)

    med_ibs = float(np.median(x[:, 0]))
    med_logd = float(np.median(x[:, 1]))
    ibs_line = np.linspace(x_lo - x_pad, x_hi + x_pad, 200)
    logd_line = np.linspace(y_lo - y_pad, y_hi + y_pad, 200)
    pred_ibs = np.clip(
        gam_te.predict(np.column_stack([ibs_line, np.full_like(ibs_line, med_logd)])),
        0.0,
        1.0,
    )
    pred_logd = np.clip(
        gam_te.predict(np.column_stack([np.full_like(logd_line, med_ibs), logd_line])),
        0.0,
        1.0,
    )

    ibs = work["IBS_Ref"]
    summary = {
        "n_samples": len(y),
        "ibs_min": float(ibs.min()),
        "ibs_median": float(ibs.median()),
        "ibs_max": float(ibs.max()),
        "pct_ibs_gt_0_99": float((ibs > 0.99).mean() * 100.0),
        "pearson_ibs_miss": float(work["IBS_Ref"].corr(work["F_MISS"])),
        "pearson_logdepth_miss": float(work["log10_Depth"].corr(work["F_MISS"])),
        "delta_aic_te_vs_additive": delta_aic,
        "pseudo_r2_te": float(pseudo_r2),
        "median_log10_depth": med_logd,
        "median_ibs": med_ibs,
    }
    return {
        "gam_add": gam_add,
        "gam_te": gam_te,
        "ibs_mesh": ibs_mesh,
        "logd_mesh": logd_mesh,
        "z_pred": z_pred,
        "ibs_line": ibs_line,
        "pred_ibs": pred_ibs,
        "logd_line": logd_line,
        "pred_logd": pred_logd,
        "med_ibs": med_ibs,
        "med_logd": med_logd,
        "summary": summary,
    }


def _gam_missing_ibs_depth_analysis(work, output_prefix, n_splines=25, te_splines=20):
    """
    Compare additive vs tensor-product GAMs for Missing ~ s(IBS) + s(log10 Depth) vs te(IBS, log10 Depth).

    Writes ``{output_prefix}.gam_model_compare.info.tsv`` and GAM diagnostic PNGs.
    """
    bundle = _fit_gam_te_bundle(work, n_splines=n_splines, te_splines=te_splines)
    if bundle is None:
        print(f"[Warning] Too few samples for GAM; skipping {output_prefix}.")
        return

    gam_add = bundle["gam_add"]
    gam_te = bundle["gam_te"]
    n_samples = int(bundle["summary"]["n_samples"])

    def _model_row(name, model):
        stats = model.statistics_
        pseudo = stats.get("pseudo_r2", {})
        if isinstance(pseudo, dict):
            pseudo_r2 = pseudo.get("explained_deviance", np.nan)
        else:
            pseudo_r2 = np.nan
        return {
            "model": name,
            "formula": name,
            "n_samples": n_samples,
            "AIC": stats.get("AIC", np.nan),
            "AICc": stats.get("AICc", np.nan),
            "GCV": stats.get("GCV", np.nan),
            "effective_dof": stats.get("edof", np.nan),
            "pseudo_r2_explained_deviance": pseudo_r2,
        }

    print("Fitting GAM models (additive vs tensor interaction)...")
    compare = pd.DataFrame(
        [
            _model_row("s(IBS) + s(log10_Depth)", gam_add),
            _model_row("te(IBS, log10_Depth)", gam_te),
        ]
    )
    compare["delta_AIC_vs_additive"] = compare["AIC"] - compare["AIC"].iloc[0]
    save_df_to_tsv(compare, f"{output_prefix}.gam_model_compare.info.tsv")

    best = compare.loc[compare["AIC"].idxmin(), "model"]
    print(
        f"GAM model comparison: best by AIC = {best}; "
        f"ΔAIC(te - additive) = {compare.loc[1, 'delta_AIC_vs_additive']:.2f}"
    )

    plot_gam_tensor_surface(
        x_grid=bundle["ibs_mesh"],
        y_grid=bundle["logd_mesh"],
        z_grid=bundle["z_pred"],
        x_label="IBS with Reference Genome",
        y_label="log10(Depth)",
        cbar_label="Predicted Missing Rate",
        title="GAM tensor surface: te(IBS, log10 Depth)",
        filename=f"{output_prefix}.gam_te_surface.png",
        vmin=0.0,
        vmax=1.0,
    )
    plot_gam_partial_curve(
        x_line=bundle["ibs_line"],
        y_line=bundle["pred_ibs"],
        x_label="IBS with Reference Genome",
        y_label="Predicted Missing Rate",
        title="GAM partial: Missing vs IBS (log10 Depth at median)",
        filename=f"{output_prefix}.gam_partial_ibs.png",
        y_lim=(0, 1),
        reference_note=f"log10(Depth) = {bundle['med_logd']:.2f}",
    )
    plot_gam_partial_curve(
        x_line=bundle["logd_line"],
        y_line=bundle["pred_logd"],
        x_label="log10(Depth)",
        y_label="Predicted Missing Rate",
        title="GAM partial: Missing vs log10(Depth) (IBS at median)",
        filename=f"{output_prefix}.gam_partial_logdepth.png",
        y_lim=(0, 1),
        reference_note=f"IBS = {bundle['med_ibs']:.4f}",
    )
    _gam_missing_residual_outlier_analysis(work, output_prefix, bundle)


def _positive_residual_outlier_threshold(residual, outlier_frac):
    """
    Upper-tail cutoff on signed residual: only observed > predicted (residual > 0) can be outliers.

    Uses the (1 − outlier_frac) quantile of signed residuals. If that quantile is ≤ 0, falls back
    to the smallest value among the top ``ceil(outlier_frac × n)`` strictly positive residuals.
    """
    res = residual.astype(float)
    n = len(res)
    if n == 0:
        return np.nan, pd.Series([], dtype=bool)

    threshold = float(res.quantile(1.0 - outlier_frac))
    if threshold <= 0:
        pos = res[res > 0]
        if pos.empty:
            return np.inf, pd.Series(False, index=res.index)
        n_flag = max(1, int(np.ceil(outlier_frac * n)))
        n_flag = min(n_flag, len(pos))
        threshold = float(pos.nlargest(n_flag).min())

    mask = res >= threshold
    return threshold, mask


def _gam_missing_residual_outlier_analysis(work, output_prefix, bundle, outlier_frac=0.01):
    """
    Flag samples with observed missing above the GAM prediction (positive residual upper tail).

    Residual = observed F_MISS − te(IBS, log10 Depth) prediction. Samples with residual ≤ 0
    (missing at or below prediction) are never flagged. By default the top ``outlier_frac``
    fraction of all samples by signed residual defines the cutoff (default 1%).

    Writes ``{output_prefix}.gam_residual.info.tsv``, residual distribution, group-composition bar charts, and scatter
    plots with outliers highlighted on log10(Depth) vs missing and IBS vs missing.
    """
    gam_te = bundle["gam_te"]
    is_sg_cov = "sg_cov" in str(output_prefix)
    depth_label = "log10(Subgenome Depth)" if is_sg_cov else "log10(Depth)"

    fit = work.copy()
    x = fit[["IBS_Ref", "log10_Depth"]].to_numpy(dtype=float)
    y = fit["F_MISS"].to_numpy(dtype=float)
    mask = np.isfinite(x).all(axis=1) & np.isfinite(y)
    fit = fit.loc[mask].copy()
    if len(fit) < 50:
        print(f"[Warning] Too few samples for GAM residual analysis; skipping {output_prefix}.")
        return

    pred = np.clip(
        gam_te.predict(fit[["IBS_Ref", "log10_Depth"]].to_numpy(dtype=float)),
        0.0,
        1.0,
    )
    fit["F_MISS_pred"] = pred
    fit["residual"] = fit["F_MISS"] - fit["F_MISS_pred"]
    fit["abs_residual"] = fit["residual"].abs()

    res = fit["residual"]
    threshold, fit["is_outlier"] = _positive_residual_outlier_threshold(res, outlier_frac)
    n_out = int(fit["is_outlier"].sum())
    pseudo_r2 = float(bundle["summary"].get("pseudo_r2_te", np.nan))

    out_cols = [
        "Sample",
        "Group",
        "Coverage",
        "IBS_Ref",
        "log10_Depth",
        "F_MISS",
        "F_MISS_pred",
        "residual",
        "abs_residual",
        "is_outlier",
    ]
    save_df_to_tsv(fit[out_cols], f"{output_prefix}.gam_residual.info.tsv")

    print(
        f"GAM residual outliers (high missing vs prediction): n={len(fit)}, flagged={n_out} "
        f"({100.0 * n_out / len(fit):.2f}%, residual ≥ {threshold:.5f}), pseudo_r2_te={pseudo_r2:.4f}"
    )

    if "Group" not in fit.columns:
        fit["Group"] = "Unknown"

    outlier_groups = fit.loc[fit["is_outlier"], "Group"].value_counts(normalize=True)
    group_order = outlier_groups.index.tolist()
    plot_bar_chart(
        names=group_order,
        values=[float(outlier_groups[g]) for g in group_order],
        title=(
            f"GAM high-missing outlier composition by sample group "
            f"(residual ≥ {threshold:.4f}, n={n_out})"
        ),
        ylabel="Fraction of outliers",
        filename=f"{output_prefix}.gam_residual_outlier_group_frac.bar.png",
        ylim=(0.0, max(float(outlier_groups.max()) * 1.15, 0.05)),
        rotate_xlabels=45,
    )

    group_rates = (
        fit.groupby("Group", observed=True)
        .agg(n=("Sample", "count"), n_outlier=("is_outlier", "sum"))
        .assign(outlier_rate=lambda d: d["n_outlier"] / d["n"])
        .sort_values("outlier_rate", ascending=False)
    )
    plot_bar_chart(
        names=group_rates.index.tolist(),
        values=group_rates["outlier_rate"].tolist(),
        title=(
            f"GAM high-missing outlier rate by sample group "
            f"(residual ≥ {threshold:.4f})"
        ),
        ylabel="Outlier rate",
        filename=f"{output_prefix}.gam_residual_outlier_rate_by_group.bar.png",
        ylim=(0.0, max(float(group_rates["outlier_rate"].max()) * 1.15, outlier_frac * 2)),
        rotate_xlabels=45,
    )

    _gam_residual_distribution_plot(fit, output_prefix)

    scatter_kw = dict(
        data=fit,
        y_col="F_MISS",
        outlier_col="is_outlier",
        ylabel="Missing Rate",
        alpha=0.35,
        s=18,
    )
    plot_scatter_with_outliers(
        x_col="log10_Depth",
        xlabel=depth_label,
        title=f"Missing vs {depth_label} (high missing vs GAM prediction)",
        filename=f"{output_prefix}.gam_residual_outliers_vs_logdepth.png",
        **scatter_kw,
    )
    plot_scatter_with_outliers(
        x_col="IBS_Ref",
        xlabel="IBS with Reference Genome",
        title="Missing vs IBS (high missing vs GAM prediction)",
        filename=f"{output_prefix}.gam_residual_outliers_vs_ibs.png",
        **scatter_kw,
    )
    _gam_residual_scatter_extras(fit, output_prefix)


def _gam_residual_distribution_plot(fit, output_prefix):
    """Histogram of signed GAM residuals with zero and outlier threshold markers."""
    res = fit["residual"].astype(float)
    pad = (float(res.max()) - float(res.min())) * DEFAULT_AXIS_PADDING_FRACTION
    if pad <= 0:
        pad = max(abs(float(res.max())), abs(float(res.min())), 1e-9) * DEFAULT_AXIS_PADDING_FRACTION

    thresholds = [
        {
            "value": 0.0,
            "label": "0",
            "color": "black",
            "linestyle": "-",
            "linewidth": 1.2,
        },
    ]
    if "is_outlier" in fit.columns and fit["is_outlier"].any():
        thr = float(fit.loc[fit["is_outlier"], "residual"].min())
        thresholds.append(
            {
                "value": thr,
                "label": f"outlier thr ({thr:.4f})",
                "color": "#cb181d",
                "linestyle": ":",
                "linewidth": 1.5,
            }
        )

    dist_kw = dict(
        data=fit,
        col="residual",
        mean_val=float(res.mean()),
        median_val=float(res.median()),
        std_val=float(res.std()),
        thresholds=thresholds,
        x_label="GAM residual",
        bins=80,
        xlim=(float(res.min()) - pad, float(res.max()) + pad),
    )
    plot_distribution_with_stats(
        title="GAM residual distribution (observed − predicted missing)",
        filename=f"{output_prefix}.gam_residual.dist.png",
        **dist_kw,
    )
    plot_distribution_with_stats(
        title="GAM residual distribution (observed − predicted missing, log Y)",
        filename=f"{output_prefix}.gam_residual.dist.logy.png",
        log_scale=True,
        **dist_kw,
    )


def _gam_residual_scatter_extras(fit, output_prefix):
    """Extra GAM residual scatter views: signed residual colormap and |residual|-scaled markers."""
    is_sg_cov = "sg_cov" in str(output_prefix)
    depth_label = "log10(Subgenome Depth)" if is_sg_cov else "log10(Depth)"
    residual_label = "GAM residual (observed − predicted missing)"

    scatter_y = dict(y_col="F_MISS", ylabel="Missing Rate", s=22)

    cmap_kw = dict(
        color_col="residual",
        cbar_label=residual_label,
        cmap=RESIDUAL_SIGNED_CMAP,
        center=0.0,
        center_fade_alpha=True,
        center_min_alpha=0.06,
        center_fade_power=0.5,
    )
    plot_scatter_continuous_color(
        data=fit,
        x_col="log10_Depth",
        xlabel=depth_label,
        title=f"Missing vs {depth_label} (color = GAM residual)",
        filename=f"{output_prefix}.gam_residual_outliers_vs_logdepth.residual_cmap.png",
        **cmap_kw,
        **scatter_y,
    )
    plot_scatter_continuous_color(
        data=fit,
        x_col="IBS_Ref",
        xlabel="IBS with Reference Genome",
        title="Missing vs IBS (color = GAM residual)",
        filename=f"{output_prefix}.gam_residual_outliers_vs_ibs.residual_cmap.png",
        **cmap_kw,
        **scatter_y,
    )

    size_kw = dict(
        data=fit,
        y_col="F_MISS",
        outlier_col="is_outlier",
        size_col="abs_residual",
        ylabel="Missing Rate",
        alpha=0.4,
    )
    plot_scatter_with_outliers_sized(
        x_col="log10_Depth",
        xlabel=depth_label,
        title=f"Missing vs {depth_label} (marker size = |residual|)",
        filename=f"{output_prefix}.gam_residual_outliers_vs_logdepth.outlier_size.png",
        **size_kw,
    )
    plot_scatter_with_outliers_sized(
        x_col="IBS_Ref",
        xlabel="IBS with Reference Genome",
        title="Missing vs IBS (marker size = |residual|)",
        filename=f"{output_prefix}.gam_residual_outliers_vs_ibs.outlier_size.png",
        **size_kw,
    )


def gam_residual_scatter_extras_from_info(residual_info_path, output_prefix=None):
    """Plot residual colormap / sized scatters from an existing ``*.gam_residual.info.tsv``."""
    fit = load_df_from_tsv(residual_info_path)
    if fit is None:
        return
    if output_prefix is None:
        output_prefix = residual_info_path.replace(".gam_residual.info.tsv", "")
    required = ("log10_Depth", "IBS_Ref", "F_MISS", "residual", "abs_residual", "is_outlier")
    missing = [c for c in required if c not in fit.columns]
    if missing:
        print(f"[Error] gam_residual_scatter_extras_from_info: missing columns {missing}")
        return
    _gam_residual_distribution_plot(fit, output_prefix)
    _gam_residual_scatter_extras(fit, output_prefix)


def gam_residual_outlier_diagnostics_from_info(
    info_path,
    output_prefix=None,
    outlier_frac=0.01,
    n_splines=25,
    te_splines=20,
):
    """Re-run GAM te fit + residual outlier diagnostics from an existing ibs_depth_miss info TSV."""
    merged = load_df_from_tsv(info_path)
    if merged is None:
        return
    if output_prefix is None:
        output_prefix = info_path.replace(".ibs_depth_miss.info.tsv", "")
    work = _prepare_ibs_depth_miss_work(merged)
    bundle = _fit_gam_te_bundle(work, n_splines=n_splines, te_splines=te_splines)
    if bundle is None:
        print(f"[Warning] GAM fit failed for {info_path}; skipping residual diagnostics.")
        return
    _gam_missing_residual_outlier_analysis(work, output_prefix, bundle, outlier_frac=outlier_frac)


def compare_subgenome_ibs_depth_gam(
    info_paths,
    output_prefix="ABD.sample.cov",
    mods=("A", "B", "D"),
):
    """
    Side-by-side A/B/D panels: GAM te surfaces, partial effects, and summary table.

    ``info_paths`` may be paths to ``{mod}.sample.cov.ibs_depth_miss.info.tsv`` files
    (order arbitrary; matched by subgenome prefix).
    """
    path_by_mod = {}
    for path in info_paths:
        base = os.path.basename(path)
        mod = base.split(".", 1)[0]
        if mod in mods:
            path_by_mod[mod] = path

    missing = [m for m in mods if m not in path_by_mod]
    if missing:
        print(f"[Warning] compare_subgenome_ibs_depth_gam: missing mods {missing}; skipping.")
        return

    bundles = {}
    summary_rows = []
    for mod in mods:
        merged = load_df_from_tsv(path_by_mod[mod])
        if merged is None:
            print(f"[Error] Failed to load {path_by_mod[mod]}")
            return
        work = _prepare_ibs_depth_miss_work(merged)
        bundle = _fit_gam_te_bundle(work)
        if bundle is None:
            print(f"[Warning] GAM fit failed for subgenome {mod}; skipping compare.")
            return
        bundles[mod] = bundle
        row = {"subgenome": mod, **bundle["summary"]}
        summary_rows.append(row)

    save_df_to_tsv(pd.DataFrame(summary_rows), f"{output_prefix}.gam_subgenome_summary.info.tsv")

    surface_panels = [
        {
            "title": f"{mod} subgenome",
            "x_grid": bundles[mod]["ibs_mesh"],
            "y_grid": bundles[mod]["logd_mesh"],
            "z_grid": bundles[mod]["z_pred"],
        }
        for mod in mods
    ]
    plot_contourf_panels(
        panels=surface_panels,
        x_label="IBS with Reference Genome",
        y_label="log10(Depth)",
        cbar_label="Predicted Missing Rate",
        suptitle="GAM te(IBS, log10 Depth) — A / B / D comparison",
        filename=f"{output_prefix}.gam_te_surface.panels.png",
        vmin=0.0,
        vmax=1.0,
    )

    ibs_panels = [
        {
            "title": f"{mod} (n={bundles[mod]['summary']['n_samples']})",
            "x": bundles[mod]["ibs_line"],
            "y": bundles[mod]["pred_ibs"],
            "note": f"log10(Depth)={bundles[mod]['med_logd']:.2f}",
        }
        for mod in mods
    ]
    plot_line_panels(
        panels=ibs_panels,
        x_label="IBS with Reference Genome",
        y_label="Predicted Missing Rate",
        suptitle="GAM partial: Missing vs IBS (log10 Depth at median)",
        filename=f"{output_prefix}.gam_partial_ibs.panels.png",
        y_lim=(0, 1),
    )

    logd_panels = [
        {
            "title": f"{mod} (pct IBS>0.99={bundles[mod]['summary']['pct_ibs_gt_0_99']:.1f}%)",
            "x": bundles[mod]["logd_line"],
            "y": bundles[mod]["pred_logd"],
            "note": f"IBS={bundles[mod]['med_ibs']:.4f}",
        }
        for mod in mods
    ]
    plot_line_panels(
        panels=logd_panels,
        x_label="log10(Depth)",
        y_label="Predicted Missing Rate",
        suptitle="GAM partial: Missing vs log10(Depth) (IBS at median)",
        filename=f"{output_prefix}.gam_partial_logdepth.panels.png",
        y_lim=(0, 1),
    )
    print(f"Subgenome GAM comparison saved with prefix {output_prefix}")


def select_gam_credible_deep_samples(
    gam_residual_info_paths,
    source_taxa_bam_map,
    output_prefix,
    n_top=1000,
    subgenome_mods=("A", "B", "D", "Others"),
):
    """
    Select high-WGS-depth samples whose GAM missing residuals are in the credible band.

    **Credible band:** ``is_outlier`` is False on every subgenome table (residual below
    the upper-tail cutoff from ``gam_residual_outlier_frac``, typically top 2% high-missing
    vs te(IBS, log10 subgenome depth) prediction).

    Samples are ranked by ``Coverage`` (whole-genome depth from taxaBamMap / info TSV) and
    the top ``n_top`` rows are written. A subset taxaBamMap is extracted from
    ``source_taxa_bam_map`` for PopDep via ``--popdep_taxa_bam_file``.

    ``gam_residual_info_paths`` must include one ``*.gam_residual.info.tsv`` per mod in
    ``subgenome_mods`` (basename prefix before first ``.`` = mod name).
    """
    path_by_mod = {}
    for path in gam_residual_info_paths:
        mod = os.path.basename(path).split(".", 1)[0]
        if mod in subgenome_mods:
            path_by_mod[mod] = path

    missing_mods = [m for m in subgenome_mods if m not in path_by_mod]
    if missing_mods:
        print(f"[Error] select_gam_credible_deep_samples: missing mods {missing_mods}")
        return None

    base_mod = subgenome_mods[0]
    merged = load_df_from_tsv(path_by_mod[base_mod])
    if merged is None or merged.empty:
        print(f"[Error] Failed to load {path_by_mod[base_mod]}")
        return None

    keep_cols = ["Sample", "Group", "Coverage"]
    merged = merged[[c for c in keep_cols if c in merged.columns]].copy()
    for mod in subgenome_mods:
        df = load_df_from_tsv(path_by_mod[mod])
        if df is None or df.empty:
            print(f"[Error] Failed to load {path_by_mod[mod]}")
            return None
        sub = df.set_index("Sample")[
            ["log10_Depth", "F_MISS", "F_MISS_pred", "residual", "abs_residual", "is_outlier"]
        ].add_prefix(f"{mod}_")
        merged = merged.merge(sub, left_on="Sample", right_index=True, how="inner")

    outlier_cols = [f"{m}_is_outlier" for m in subgenome_mods]
    credible = merged.loc[~merged[outlier_cols].any(axis=1)].copy()
    credible = credible.sort_values("Coverage", ascending=False, kind="mergesort")
    n_select = min(int(n_top), len(credible))
    selected = credible.head(n_select).copy()
    selected["rank_by_wgs_depth"] = range(1, n_select + 1)

    rank_cols = (
        ["rank_by_wgs_depth", "Sample", "Group", "Coverage"]
        + [f"{m}_log10_Depth" for m in subgenome_mods]
        + [f"{m}_F_MISS" for m in subgenome_mods]
        + [f"{m}_residual" for m in subgenome_mods]
    )
    save_df_to_tsv(selected[rank_cols], f"{output_prefix}.gam_credible_deep_top{n_select}.info.tsv")

    tbm = load_df_from_tbm(source_taxa_bam_map)
    if tbm is None or tbm.empty:
        print(f"[Error] Failed to load source taxaBamMap: {source_taxa_bam_map}")
        return None

    taxa_set = set(selected["Sample"])
    subset = tbm[tbm["Sample"].isin(taxa_set)].copy()
    order = {s: i for i, s in enumerate(selected["Sample"])}
    subset["_ord"] = subset["Sample"].map(order)
    subset = subset.sort_values("_ord").drop(columns="_ord")

    missing_taxa = taxa_set - set(subset["Sample"])
    if missing_taxa:
        print(
            f"[Warning] {len(missing_taxa)} selected sample(s) missing from "
            f"{source_taxa_bam_map}: {sorted(missing_taxa)[:5]}..."
        )

    map_path = f"{output_prefix}.taxaBamMap.txt"
    with open(source_taxa_bam_map, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n")
    with open(map_path, "w", encoding="utf-8") as out:
        out.write(header + "\n")
        for _, row in subset.iterrows():
            bams = row["Bam_Path"]
            if isinstance(bams, list):
                bam_field = "\t".join(str(b) for b in bams)
            else:
                bam_field = str(bams)
            out.write(f"{row['Sample']}\t{row['Coverage']}\t{bam_field}\n")

    print(
        f"GAM credible deep selection: {len(credible)} credible / {len(merged)} intersected; "
        f"wrote top {n_select} (Coverage {selected['Coverage'].min():.4f}–"
        f"{selected['Coverage'].max():.4f})"
    )
    print(f"  rank table: {output_prefix}.gam_credible_deep_top{n_select}.info.tsv")
    print(f"  taxaBamMap: {map_path}")
    return selected

