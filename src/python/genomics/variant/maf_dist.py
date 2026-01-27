import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 1. 读取 PLINK2 的频率文件
# 使用 r'' (raw string) 来避免 SyntaxWarning: invalid escape sequence '\s'
df = pd.read_csv("/data1/dazheng_tusr1/vmap4.VCF.v1/chr002.maf.afreq", sep=r'\s+')

# 2. 计算真实的 MAF (取 Alt 频率和 Ref 频率中较小的一个)
# PLINK2 的列名可能是 'ALT_FREQS' 或 'AF'
df['MAF'] = df['ALT_FREQS'].apply(lambda x: x if x <= 0.5 else 1 - x)

# Filter data for 0-0.01
df_sub = df[df['MAF'] <= 0.01].copy()

# Calculate statistics (Global)
mean_maf_all = df['MAF'].mean()
median_maf_all = df['MAF'].median()

# Calculate count and percentage
count_sub = len(df_sub)
count_total = len(df)
percent_sub = (count_sub / count_total) * 100

print(f"Total variants: {count_total}")
print(f"Variants with MAF <= 0.01: {count_sub}")
print(f"Percentage: {percent_sub:.2f}%")
print(f"Global MAF - Mean: {mean_maf_all:.5f}, Median: {median_maf_all:.5f}")

# Plotting Function
def plot_maf_sub(log_scale, filename):
    plt.figure(figsize=(10, 6))
    
    # Aesthetic colors
    hist_color = '#4c72b0' # muted blue
    
    # Histogram
    sns.histplot(df_sub['MAF'], bins=200, kde=False, color=hist_color, edgecolor='white', linewidth=0.1, alpha=0.8)

    # Vertical Lines - Thresholds
    plt.axvline(x=0.001, color='#55a868', linestyle='--', linewidth=1.5, label='Threshold 0.001') # green
    plt.axvline(x=0.005, color='#c44e52', linestyle='--', linewidth=1.5, label='Threshold 0.005') # red

    # Vertical Lines - Stats (Global)
    plt.plot([], [], ' ', label=f'Mean (All): {mean_maf_all:.5f}')
    plt.axvline(x=median_maf_all, color='purple', linestyle='-.', linewidth=1.5, label=f'Median (All): {median_maf_all:.5f}')

    # Decorations
    scale_type = "Log" if log_scale else "Linear"
    plt.title(f'Distribution of Minor Allele Frequency (0-0.01) - {scale_type} Scale', fontsize=15)
    plt.xlabel('MAF', fontsize=12)
    plt.ylabel(f'Count of Variants', fontsize=12)
    plt.xlim(0, 0.01)
    
    if log_scale:
        plt.yscale('log')
    
    plt.legend(loc='upper right', framealpha=0.9)
    plt.grid(axis='y', linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    print(f"Plot saved to {filename}")

# 3. Plot Linear
plot_maf_sub(log_scale=False, filename='maf_distribution_0.01_linear.png')

# 4. Plot Log
plot_maf_sub(log_scale=True, filename='maf_distribution_0.01_log.png')

# 5. Plot Full Range (0-0.5) Log Scale
print("Generating Full Range MAF Plot...")
plt.figure(figsize=(10, 6))

mean_maf_all = df['MAF'].mean()
median_maf_all = df['MAF'].median()

# Histogram
sns.histplot(df['MAF'], bins=500, kde=False, color='forestgreen', edgecolor='none')

# Vertical Lines
plt.axvline(x=mean_maf_all, color='black', linestyle='-', linewidth=1.5, label=f'Mean: {mean_maf_all:.5f}')
plt.axvline(x=median_maf_all, color='purple', linestyle='-.', linewidth=1.5, label=f'Median: {median_maf_all:.5f}')

plt.title('Distribution of Minor Allele Frequency (0-0.5) - Log Scale', fontsize=15)
plt.xlabel('MAF', fontsize=12)
plt.ylabel('Count of Variants (Log Scale)', fontsize=12)
plt.xlim(0, 0.5)
plt.yscale('log')

plt.legend()
plt.grid(axis='y', alpha=0.3)
plt.tight_layout()

plt.savefig('maf_distribution_0.5_log.png', dpi=300)
print("Plot saved to maf_distribution_0.5_log.png")

# 6. Plot Range (0-0.05) Linear Scale
print("Generating Range (0-0.05) MAF Plot...")
plt.figure(figsize=(10, 6))

# Filter for plot range (optional, but good for bins)
df_sub_05 = df[df['MAF'] <= 0.05]

# Histogram (color same as 0.5 log plot: forestgreen)
sns.histplot(df_sub_05['MAF'], bins=200, kde=False, color='forestgreen', edgecolor='none')

# Vertical Lines
plt.axvline(x=0.01, color='red', linestyle='--', linewidth=1.5, label='Threshold 0.01')
plt.axvline(x=mean_maf_all, color='black', linestyle='-', linewidth=1.5, label=f'Mean (All): {mean_maf_all:.5f}')
plt.axvline(x=median_maf_all, color='purple', linestyle='-.', linewidth=1.5, label=f'Median (All): {median_maf_all:.5f}')

plt.title('Distribution of Minor Allele Frequency (0-0.05)', fontsize=15)
plt.xlabel('MAF', fontsize=12)
plt.ylabel('Count of Variants', fontsize=12)
plt.xlim(0, 0.05)

plt.legend()
plt.grid(axis='y', alpha=0.3)
plt.tight_layout()

plt.savefig('maf_distribution_0.05_linear.png', dpi=300)
print("Plot saved to maf_distribution_0.05_linear.png")
