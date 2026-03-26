from .variant_utils import load_df_from_plink_variant
from genetics.genomics.variant.mq import load_mq_data
from infra.utils.io import load_df_from_space_sep, save_df_to_tsv
from infra.utils.graph import plot_distribution_with_stats, plot_regression_comparison, plot_mean_variance_fit, plot_qq_residuals, plot_scatter_with_outliers
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import chi2

# ==================================================================================
# Data Loading Helpers
# ==================================================================================

def load_popdep_data(filepath):
    """
    Loads population depth file.
    Expected columns: Position, Depth_Mean, Depth_SD (space separated)
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
# Analysis 1: Mean-Variance Metrics
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

def ana_popdep_curve_fit(
    popdep_file, 
    output_prefix, 
    sample_size=50000
):
    """
    Fits Poisson vs Negative Binomial models to Mean-Variance relationship.
    """
    df = load_popdep_data(popdep_file)
    if df is None: return

    # Filter valid
    df_valid = df[(df['Depth_Mean'] > 0) & (df['Depth_Var'] > 0)].dropna()
    print(f"Valid sites for curve fit: {len(df_valid)}")

    x_data = df_valid['Depth_Mean'].values
    y_data = df_valid['Depth_Var'].values
    n_samples = len(y_data)

    # 1. Poisson
    y_pred_poisson = x_data

    # 2. Negative Binomial
    print("Fitting Negative Binomial Model...")
    phi_est = 0
    y_pred_nb = y_pred_poisson
    try:
        popt, pcov = curve_fit(nb_variance_function, x_data, y_data, bounds=(0, np.inf))
        phi_est = popt[0]
        y_pred_nb = nb_variance_function(x_data, phi_est)
        print(f"Estimated Phi: {phi_est:.6f}")
    except Exception as e:
        print(f"Curve fitting failed: {e}")
        # Continue with best effort

    # 3. Report
    aic_p, bic_p, rss_p = calculate_criteria(y_data, y_pred_poisson, 0, n_samples)
    aic_nb, bic_nb, rss_nb = calculate_criteria(y_data, y_pred_nb, 1, n_samples)
    
    print(f"Poisson | RSS: {rss_p:.2e} | AIC: {aic_p:.2f}")
    print(f"NegBin  | RSS: {rss_nb:.2e} | AIC: {aic_nb:.2f}")

    # 4. Plot Mean-Variance
    x_range = np.linspace(df_valid['Depth_Mean'].min(), df_valid['Depth_Mean'].max(), 500)
    
    plot_mean_variance_fit(
        x_data=df_valid['Depth_Mean'],
        y_data=df_valid['Depth_Var'],
        x_line=x_range,
        y_line_poisson=x_range,
        y_line_nb=nb_variance_function(x_range, phi_est),
        phi_est=phi_est,
        filename=f"{output_prefix}.mean_loss_fit.png",
        sample_size=sample_size
    )

    # 5. QQ Plot
    resid_nb = y_data - y_pred_nb
    if np.std(resid_nb) > 0:
        resid_std = (resid_nb - np.mean(resid_nb)) / np.std(resid_nb)
        plot_qq_residuals(resid_std, f"{output_prefix}.qq_residuals.png", title="QQ Plot of Residuals (NB)")

# ==================================================================================
# Analysis 2: Missing Rate vs Depth
# ==================================================================================

def ana_popdep_missing_reg(
    popdep_file, 
    vmiss_file, 
    output_prefix, 
    log_scale=True
):
    """
    Analyzes relationship between Depth Metrics and Variant Missing Rate.
    """
    # 1. Load
    df_dep = load_popdep_data(popdep_file)
    df_miss = load_df_from_plink_variant(vmiss_file) # Using renamed function
    if df_dep is None or df_miss is None: return

    # Merge
    merged = pd.merge(df_miss[['Position', 'F_MISS']], df_dep, on='Position', how='inner')
    if merged.empty:
        print("[Error] No overlapping positions found.")
        return

    # Regression Plots
    metrics = [
        ('Depth_Mean', 'Mean Depth', False),
        ('Depth_CV', 'Depth CV', False),
    ]
    
    if log_scale:
        for col in ['Depth_Mean', 'Depth_CV', 'F_MISS']:
            merged[f'Log_{col}'] = np.log(merged[col].replace(0, np.nan))
        metrics.extend([
            ('Log_Depth_Mean', 'Log Mean Depth', True),
            ('Log_Depth_CV', 'Log Depth CV', True)
        ])

    for x_col, x_label, is_log in metrics:
        y_col = 'Log_F_MISS' if is_log else 'F_MISS'
        y_label = 'Log Missing Rate' if is_log else 'Missing Rate'
        
        plot_regression_comparison(
            merged, x_col, y_col, 
            x_label=x_label, y_label=y_label, 
            filename=f"{output_prefix}.reg_{y_col}_vs_{x_col}.png",
            title=f"{y_label} vs {x_label}"
        )

# ==================================================================================
# Analysis 3: Mahalanobis Distance for Variant Depth Distribution
# ==================================================================================

def ana_popdep_mahalanobis(popdep_file, output_prefix, threshold_p_value=0.99):
    """
    Calculates Mahalanobis distance for variant depth distribution (Log_Depth_Mean vs Log_Depth_CV),
    sets an outlier threshold, saves the result, and plots it.
    """
    print(f"[Info] Running Mahalanobis distance analysis on: {popdep_file}")
    df = load_popdep_data(popdep_file)
    if df is None: return

    # Filter out missing/zero variants to perform log transform safely
    df_valid = df[(df['Depth_Mean'] > 0) & (df['Depth_CV'] > 0)].copy()
    if df_valid.empty:
        print("[Error] No valid data for Mahalanobis calculation.")
        return

    df_valid['Log_Depth_Mean'] = np.log(df_valid['Depth_Mean'])
    df_valid['Log_Depth_CV'] = np.log(df_valid['Depth_CV'])
    
    data = df_valid[['Log_Depth_Mean', 'Log_Depth_CV']].values
    
    # Fast vectorized Mahalanobis calculation
    mean_vec = np.mean(data, axis=0)
    cov_matrix = np.cov(data, rowvar=False)
    
    try:
        inv_cov_matrix = np.linalg.inv(cov_matrix)
    except np.linalg.LinAlgError:
        print("[Error] Covariance matrix is singular, cannot calculate Mahalanobis distance.")
        return

    diff = data - mean_vec
    # Vectorized computation: Mahalanobis squared distance
    mahalanobis_sq = np.sum(np.dot(diff, inv_cov_matrix) * diff, axis=1)
    
    df_valid['Mahalanobis_Sq'] = mahalanobis_sq
    df_valid['Mahalanobis'] = np.sqrt(mahalanobis_sq)

    # Threshold based on Chi-Square distribution (df=2)
    threshold = chi2.ppf(threshold_p_value, df=2)
    df_valid['Outlier'] = df_valid['Mahalanobis_Sq'] > threshold

    # Save to file
    out_file = f"{output_prefix}.info.tsv"
    out_df = df_valid[['Position', 'Depth_Mean', 'Depth_SD', 'Depth_CV', 'Mahalanobis', 'Mahalanobis_Sq', 'Outlier']]
    save_df_to_tsv(out_df, out_file)

    # Plot
    plot_file = f"{output_prefix}.mahalanobis.png"
    plot_scatter_with_outliers(
        data=df_valid,
        x_col='Log_Depth_Mean',
        y_col='Log_Depth_CV',
        outlier_col='Outlier',
        title=f'Mahalanobis Distance Outlier Detection (p={threshold_p_value})',
        filename=plot_file,
        xlabel='Log Depth Mean',
        ylabel='Log Depth CV',
        color_normal='blue',
        color_outlier='red',
        s=2
    )


