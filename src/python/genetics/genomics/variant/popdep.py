from genetics.genomics.variant import load_mq, load_vmiss
from infra.utils import load_df_from_space_sep, plot_distribution_with_stats, plot_regression_comparison, plot_mean_variance_fit, plot_qq_residuals
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats
import statsmodels.api as sm

# ==================================================================================
# Data Loading Helpers
# ==================================================================================

def load_popdep(filepath):
    """
    Loads population depth file.
    Expected columns: Position, Depth_Mean, Depth_SD (space separated)
    Calculates derived metrics: Depth_Var, Depth_CV.
    """
    print(f"[Info] Loading PopData: {filepath}")
    df = load_df_from_space_sep(filepath)
    if df is None: return None
    
    # Validation
    req_cols = ['Position', 'Depth_Mean', 'Depth_SD']
    if not all(c in df.columns for c in req_cols):
        print(f"[Error] Missing columns in popdep file. Required: {req_cols}")
        return None

    # Derived Metrics
    df['Depth_Var'] = df['Depth_SD'] ** 2
    # Avoid division by zero
    df['Depth_CV'] = np.where(
        df['Depth_Mean'] > 0, 
        df['Depth_SD'] / df['Depth_Mean'], 
        np.nan
    )
    return df

# ==================================================================================
# Analysis 1: Mean-Variance Curve Fit (Poisson vs Negative Binomial)
# ==================================================================================

def nb_variance_function(x, phi):
    """Negative Binomial Variance: Var = Mean + phi * Mean^2"""
    return x + phi * (x ** 2)

def calculate_criteria(y_true, y_pred, k, n):
    """Calculates AIC, BIC based on RSS."""
    rss = np.sum((y_true - y_pred) ** 2)
    aic = n * np.log(rss / n) + 2 * k
    bic = n * np.log(rss / n) + k * np.log(n)
    return aic, bic, rss

def analyze_curve_fit(popdep_file, output_prefix, sample_size=50000):
    """
    Fits Poisson vs Negative Binomial models to Mean-Variance relationship.
    Generates Fit Plot and QQ Plot.
    """
    df = load_popdep(popdep_file)
    if df is None: return

    # Filter valid
    df_valid = df[(df['Depth_Mean'] > 0) & (df['Depth_Var'] > 0)].dropna()
    print(f"Valid sites for curve fit: {len(df_valid)}")

    x_data = df_valid['Depth_Mean'].values
    y_data = df_valid['Depth_Var'].values
    n_samples = len(y_data)

    # 1. Poisson: Var = Mean
    y_pred_poisson = x_data

    # 2. Negative Binomial: Var = Mean + phi * Mean^2
    print("Fitting Negative Binomial Model...")
    try:
        popt, pcov = curve_fit(nb_variance_function, x_data, y_data, bounds=(0, np.inf))
        phi_est = popt[0]
        y_pred_nb = nb_variance_function(x_data, phi_est)
        print(f"Estimated Phi: {phi_est:.6f}")
    except Exception as e:
        print(f"Curve fitting failed: {e}")
        return

    # 3. Comparison
    aic_p, bic_p, rss_p = calculate_criteria(y_data, y_pred_poisson, 0, n_samples)
    aic_nb, bic_nb, rss_nb = calculate_criteria(y_data, y_pred_nb, 1, n_samples)

    print(f"{'Model':<10} | {'RSS':<15} | {'AIC':<15} | {'BIC':<15}")
    print(f"{'Poisson':<10} | {rss_p:<15.2e} | {aic_p:<15.2f} | {bic_p:<15.2f}")
    print(f"{'NegBin':<10} | {rss_nb:<15.2e} | {aic_nb:<15.2f} | {bic_nb:<15.2f}")

    # 4. Plot Mean-Variance
    x_range = np.linspace(df_valid['Depth_Mean'].min(), df_valid['Depth_Mean'].max(), 500)
    out_fit = f"{output_prefix}_mean_loss_fit.png"
    
    plot_mean_variance_fit(
        x_data=df_valid['Depth_Mean'],
        y_data=df_valid['Depth_Var'],
        x_line=x_range,
        y_line_poisson=x_range,
        y_line_nb=nb_variance_function(x_range, phi_est),
        phi_est=phi_est,
        filename=out_fit,
        sample_size=sample_size
    )

    # 5. QQ Plot (Residuals of Variance)
    resid_nb = y_data - y_pred_nb
    resid_std = (resid_nb - np.mean(resid_nb)) / np.std(resid_nb)
    
    out_qq = f"{output_prefix}_qq_residuals.png"
    plot_qq_residuals(resid_std, out_qq, title="QQ Plot of Standardized Residuals (NB Model)")

# ==================================================================================
# Analysis 2: Missing Rate vs Depth (Robust Regression)
# ==================================================================================

def analyze_missing_reg(popdep_file, vmiss_file, output_prefix, log_scale=True):
    """
    Analyzes relationship between Depth Metrics and Variant Missing Rate.
    Performs standard and log-log regression using Robust Regression (Huber).
    """
    # 1. Load & Merge
    df_dep = load_popdep(popdep_file)
    df_miss = load_vmiss(vmiss_file)
    if df_dep is None or df_miss is None: return

    merged = pd.merge(df_miss[['Position', 'F_MISS']], df_dep, on='Position', how='inner')
    if merged.empty:
        print("[Error] No overlapping positions found between PopDep and VMISS.")
        return
    print(f"Merged {len(merged)} sites.")

    # 2. Plot Distribution of F_MISS
    plot_distribution_with_stats(
        merged, 'F_MISS', 
        title='Distribution of Variant Missing Rate',
        filename=f"{output_prefix}_dist_fmiss.png",
        x_label='Missing Rate (F_MISS)',
        median_val=merged['F_MISS'].median(),
        mean_val=merged['F_MISS'].mean()
    )

    # 3. Regression Analysis
    # Define metrics to check
    # Format: (ColName, DisplayName, ApplyLog)
    metrics = [
        ('Depth_Mean', 'Mean Depth', False),
        ('Depth_SD', 'Depth SD', False),
        ('Depth_CV', 'Depth CV', False),
    ]
    
    if log_scale:
        # Create Log Columns
        for col in ['Depth_Mean', 'Depth_SD', 'Depth_CV', 'F_MISS']:
            merged[f'Log_{col}'] = np.log(merged[col].replace(0, np.nan))
        
        metrics.extend([
            ('Log_Depth_Mean', 'Log Mean Depth', True),
            ('Log_Depth_SD', 'Log Depth SD', True),
            ('Log_Depth_CV', 'Log Depth CV', True)
        ])

    for x_col, x_label, is_log in metrics:
        y_col = 'Log_F_MISS' if is_log else 'F_MISS'
        y_label = 'Log Missing Rate' if is_log else 'Missing Rate'
        
        out_file = f"{output_prefix}_reg_{y_col}_vs_{x_col}.png"
        
        plot_regression_comparison(
            merged, x_col, y_col, 
            x_label=x_label, 
            y_label=y_label, 
            filename=out_file,
            title=f"{y_label} vs {x_label}"
        )

# ==================================================================================
# Analysis 3: Mapping Quality (MQ) vs Depth
# ==================================================================================

def analyze_mq_reg(popdep_file, mq_file, output_prefix):
    """
    Analyzes MQ distribution and its relationship with Depth.
    """
    # 1. Load & Merge
    df_dep = load_popdep(popdep_file)
    df_mq = load_mq(mq_file)
    if df_dep is None or df_mq is None: return

    merged = pd.merge(df_dep, df_mq, on='Position', how='inner')
    if merged.empty:
        print("[Error] No overlapping positions found between PopDep and MQ.")
        return
    print(f"Merged {len(merged)} sites.")

    # 2. MQ Distribution
    plot_distribution_with_stats(
        merged, 'MQ',
        title='Distribution of Mapping Quality',
        filename=f"{output_prefix}_dist_mq.png",
        x_label='Mapping Quality (MQ)',
        mean_val=merged['MQ'].mean(),
        median_val=merged['MQ'].median()
    )

    # 3. Regression (Depth vs MQ)
    target_metrics = [
        ('Depth_Mean', 'Mean Depth'),
        ('Depth_SD', 'Depth SD'),
        ('Depth_CV', 'Depth CV')
    ]

    for x_col, x_label in target_metrics:
        # We model MQ as Y usually? Or MQ affecting Depth?
        # popdep_mq.py plotted Depth (X) vs MQ (Y). Let's stick to that.
        # Actually usually filtering is done on MQ to see effect on Depth consistency.
        # popdep_mq.py had plot_regression(merged, 'Depth_Mean', mq_col, ...) -> X=Depth, Y=MQ
        
        out_file = f"{output_prefix}_reg_mq_vs_{x_col}.png"
        plot_regression_comparison(
            merged, 
            x_col=x_col, 
            y_col='MQ',
            x_label=x_label, 
            y_label='Mapping Quality',
            filename=out_file
        )


# ==================================================================================
# Main CLI
# ==================================================================================

def main():
    parser = argparse.ArgumentParser(description="Population Depth Analysis Toolkit")
    subparsers = parser.add_subparsers(dest='command', help='Sub-command help')

    # Cmd 1: Curve Fit
    p_fit = subparsers.add_parser('curve_fit', help='Analyze Mean-Variance relationship')
    p_fit.add_argument('-d', '--depth', required=True, help='Population depth file')
    p_fit.add_argument('-o', '--out', default='curve_fit', help='Output prefix')

    # Cmd 2: Missing vs Depth
    p_miss = subparsers.add_parser('missing', help='Analyze Missing Rate vs Depth')
    p_miss.add_argument('-d', '--depth', required=True, help='Population depth file')
    p_miss.add_argument('-m', '--vmiss', required=True, help='Variant missingness file (.vmiss)')
    p_miss.add_argument('-o', '--out', default='missing_ana', help='Output prefix')
    p_miss.add_argument('--no_log', action='store_true', help='Disable log-scale regression')

    # Cmd 3: MQ vs Depth
    p_mq = subparsers.add_parser('mq', help='Analyze Mapping Quality vs Depth')
    p_mq.add_argument('-d', '--depth', required=True, help='Population depth file')
    p_mq.add_argument('-q', '--mq', required=True, help='Mapping Quality file (Pos MQ)')
    p_mq.add_argument('-o', '--out', default='mq_ana', help='Output prefix')

    args = parser.parse_args()

    if args.command == 'curve_fit':
        analyze_curve_fit(args.depth, args.out)
    elif args.command == 'missing':
        analyze_missing_reg(args.depth, args.vmiss, args.out, log_scale=not args.no_log)
    elif args.command == 'mq':
        analyze_mq_reg(args.depth, args.mq, args.out)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()

