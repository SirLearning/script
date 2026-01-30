import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
import statsmodels.api as sm
from scipy.optimize import curve_fit
from scipy import stats

# Global settings
PLOT_SAMPLE_SIZE = 50000 
OUTPUT_PREFIX = "curve_fit_analysis_"

# ==========================================
# 1. Data Loading (Adapted from popdep.py)
# ==========================================
print("Processing Depth Data Analysis...")

# Define paths (Adjust these paths as necessary)
input_file = "/data/home/tusr1/01projects/vmap4/05reliable.lib/01test/chr002_popdep.txt" 
pos_file = "/data/home/tusr1/01projects/vmap4/05reliable.lib/01test/2_1_122798052.pos.txt"

# Read Depth Data
try:
    print(f"Reading {input_file}...")
    df_popdep = pd.read_csv(input_file, sep=r'\s+', usecols=['Position', 'Depth_SD', 'Depth_Mean'])
except Exception as e:
    print(f"Failed to read file: {e}")
    sys.exit(1)

# Calculate Variance
df_popdep['Depth_Var'] = df_popdep['Depth_SD'] ** 2

# Drop NaNs / Zeros
# We only care about sites with Mean > 0 for this analysis
df_valid = df_popdep[df_popdep['Depth_Mean'] > 0].dropna().copy()
print(f"Valid sites for analysis: {len(df_valid)}")

# Read Highlight Positions
df_highlight = pd.DataFrame()
try:
    df_pos_target = pd.read_csv(pos_file, sep=r'\s+', header=None, names=['Chrom', 'Position'])
    target_positions = set(df_pos_target['Position'])
    df_highlight = df_valid[df_valid['Position'].isin(target_positions)].copy()
    print(f"Found {len(df_highlight)} sites to highlight.")
except Exception as e:
    print(f"Warning: Failed to read/process position file: {e}")

# Downsample for plotting (Regression uses ALL valid data, Plotting uses sample)
if len(df_valid) > PLOT_SAMPLE_SIZE:
    df_plot = df_valid.sample(n=PLOT_SAMPLE_SIZE, random_state=42)
else:
    df_plot = df_valid

# ==========================================
# 2. Curve Fitting Logic
# ==========================================

# X: Mean, Y: Variance
x_data = df_valid['Depth_Mean'].values
y_data = df_valid['Depth_Var'].values

# --- Model 1: Poisson (Theoretical) ---
# Formula: Var = Mean
# No fitting needed, it's a fixed assumption
y_pred_poisson = x_data 

# --- Model 2: Negative Binomial (Curve Fit) ---
# Formula: Var = Mean + phi * Mean^2
def nb_variance_function(x, phi):
    return x + phi * (x**2)

print("\n--- Fitting Negative Binomial Model ---")
# Bounds: phi must be >= 0 (Overdispersion)
popt, pcov = curve_fit(nb_variance_function, x_data, y_data, bounds=(0, np.inf))
phi_est = popt[0]
y_pred_nb = nb_variance_function(x_data, phi_est)

print(f"Estimated Phi (Overdispersion Coeff): {phi_est:.6f}")

# ==========================================
# 3. Model Comparison (AIC/BIC)
# ==========================================
# Since we are doing non-linear least squares regression on summary statistics (Mean vs Var),
# we calculate AIC based on Residual Sum of Squares (RSS).
# AIC = n * ln(RSS/n) + 2k
# BIC = n * ln(RSS/n) + k * ln(n)

def calculate_criteria(y_true, y_pred, k, n):
    resid = y_true - y_pred
    rss = np.sum(resid**2)
    
    # If perfect fit or rss is 0 (unlikely), handle log domain error
    if rss <= 0: return -np.inf, -np.inf, rss

    aic = n * np.log(rss / n) + 2 * k
    bic = n * np.log(rss / n) + k * np.log(n)
    return aic, bic, rss

n_samples = len(y_data)

# Poisson: k=0 (No free parameters adjusted for the variance function itself)
aic_pois, bic_pois, rss_pois = calculate_criteria(y_data, y_pred_poisson, 0, n_samples)

# NB: k=1 (phi is estimated)
aic_nb, bic_nb, rss_nb = calculate_criteria(y_data, y_pred_nb, 1, n_samples)

print("\n--- Model Comparison (Based on RSS of Mean-Variance Fit) ---")
print(f"{'Model':<10} | {'RSS':<15} | {'AIC':<15} | {'BIC':<15}")
print("-" * 65)
print(f"{'Poisson':<10} | {rss_pois:<15.2e} | {aic_pois:<15.2f} | {bic_pois:<15.2f}")
print(f"{'NegBin':<10} | {rss_nb:<15.2e} | {aic_nb:<15.2f} | {bic_nb:<15.2f}")

diff_aic = aic_pois - aic_nb
if diff_aic > 10:
    print(f"\nConclusion: Negative Binomial fit is significantly better (Delta AIC = {diff_aic:.2f}).")
    print(f"The data exhibits Overdispersion (Phi={phi_est:.4f}).")
else:
    print("\nConclusion: Poisson model is sufficient.")

# # ==========================================
# # 4. Visualization 1: Mean-Variance Plot
# # ==========================================
# plt.figure(figsize=(10, 8))

# # 1. Plot Background data
# # Filter for log scale safety
# plt_mask = (df_plot['Depth_Mean'] > 0) & (df_plot['Depth_Var'] > 0)
# plt.scatter(df_plot.loc[plt_mask, 'Depth_Mean'], 
#             df_plot.loc[plt_mask, 'Depth_Var'], 
#             alpha=0.2, s=10, color='gray', label='Observed Data (Sampled)')

# # 2. Plot Highlight data
# if not df_highlight.empty:
#     hl_mask = (df_highlight['Depth_Mean'] > 0) & (df_highlight['Depth_Var'] > 0)
#     plt.scatter(df_highlight.loc[hl_mask, 'Depth_Mean'], 
#                 df_highlight.loc[hl_mask, 'Depth_Var'], 
#                 alpha=0.9, s=30, color='red', marker='X', edgecolors='white', linewidth=0.5, label='Highlight')

# # 3. Plot Fitted Lines
# x_range = np.linspace(df_valid['Depth_Mean'].min(), df_valid['Depth_Mean'].max(), 500)
# y_line_pois = x_range
# y_line_nb = nb_variance_function(x_range, phi_est)

# plt.plot(x_range, y_line_pois, 'b--', linewidth=2, label='Poisson (Var = Mean)')
# plt.plot(x_range, y_line_nb, 'r-', linewidth=2, label=f'NB Fit (phi={phi_est:.4f})')

# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('Mean Depth (log scale)')
# plt.ylabel('Depth Variance (log scale)')
# plt.title(f'Mean-Variance Relationship\n(Comparison of Poisson vs Negative Binomial)')
# plt.legend()
# plt.grid(True, which="both", ls="-", alpha=0.2)

# filename_fit = f"{OUTPUT_PREFIX}mean_variance_fit.png"
# plt.savefig(filename_fit, dpi=300)
# print(f"Saved fit plot to {filename_fit}")
# plt.close()

# ==========================================
# 5. Visualization 2: QQ Plot of Residuals
# ==========================================
# To validate the NB model, we check the residuals.
# "Standardized Residuals" = (Observed - Predicted) / Std_Deviation_of_Predicted
# Under the NB model assumption, Std_Dev approx sqrt(Predicted_Var) = sqrt(Mean + phi*Mean^2) ?
# Wait, curve_fit minimizes squared error of Variance.
# We look at the residuals of the VARIANCE prediction.
# Residual = Var_Obs - Var_Pred_NB

print("\n--- Generating QQ Plots ---")
resid_nb = y_data - y_pred_nb

# Standardize residuals (Z-score)
# Note: For non-linear regression on heteroscedastic data, residuals are rarely perfectly normal.
# But this helps visualize heavy tails.
resid_std = (resid_nb - np.mean(resid_nb)) / np.std(resid_nb)

plt.figure(figsize=(8, 6))
# Using statsmodels qqplot
ax = plt.gca()
sm.qqplot(resid_std, line='45', fit=True, ax=ax, alpha=0.1, markerfacecolor='teal', markeredgecolor='none')
plt.title(f'QQ Plot of Residuals (NB Model Fit)\n(Variance vs Predicted Variance)')
plt.grid(True, alpha=0.3)

# Adjust xlim based on theoretical quantiles count
# Since we are comparing to a Normal distribution (approx standard normal due to standardization)
n_points = len(resid_std)
# Calculate approximate bounds for theoretical quantiles
# Using slightly wider probability range to ensure all points are included
th_bound = abs(stats.norm.ppf(1 / (n_points + 1)))
plt.xlim(-th_bound - 0.5, th_bound + 0.5)

filename_qq = f"{OUTPUT_PREFIX}qq_plot_residuals.png"
plt.savefig(filename_qq, dpi=300)
print(f"Saved QQ plot to {filename_qq}")
plt.close()

print("Analysis Done.")
