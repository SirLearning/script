import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
# Add src/python methods
from infra.file_utils import load_df_from_space_sep

def plot_variant_missingness_distribution_plink(
    input_file = "/data1/dazheng_tusr1/vmap4.VCF.v1/check_missing.vmiss",
    output_prefix="variant_missingness"
):
    """
    Plots the distribution of variant missing rates from a .vmiss file.
    Input: .vmiss file
    Output: Histogram plot
    """
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    # 1. 读取数据 (注意：PLINK2的文件通常以 # 开头，read_csv会自动处理或需指定 sep)
    # 如果是空格分隔，使用 sep='\s+'
    df = pd.read_csv(input_file, sep='\s+')

    # 2. 检查列名，PLINK2 经常在第一列前加 '#'
    # 我们需要的是 'F_MISS' 这一列
    missing_col = 'F_MISS'

    # 3. 绘图设置
    plt.figure(figsize=(10, 6))
    sns.histplot(df[missing_col], bins=50, kde=False, color='skyblue', edgecolor='black')

    # 4. 添加标题和标签
    plt.title('Distribution of Variant Missing Rate', fontsize=15)
    plt.xlabel('Missing Rate (F_MISS)', fontsize=12)
    plt.ylabel('Count of Variants', fontsize=12)

    # 5. 可选：添加一条红线表示常用的阈值（如 0.1）
    plt.axvline(x=0.1, color='red', linestyle='--', label='Threshold = 0.1')
    plt.legend()

    plt.grid(axis='y', alpha=0.3)
    plt.savefig('missing_rate_distribution.png')
    print("Plot saved to missing_rate_distribution.png")


def plot_site_missingness_vcftools(
    input_file,
    output_prefix="site_missingness"
):
    """
    Plots histogram of Site Missingness from .lmiss file.
    Input: .lmiss file
    Output: Histogram plot
    """
    print(f"[Info] Processing Site Missingness: {input_file}")
    try:
        df = pd.read_csv(input_file, sep=r'\s+')
        if 'F_MISS' not in df.columns:
            print(f"[Error] 'F_MISS' column not found in {input_file}")
            return

        plt.figure(figsize=(10, 8))
        sns.histplot(data=df, x='F_MISS', binwidth=0.01, color="forestgreen", edgecolor="black", kde=False)
        plt.title("Site Missingness")
        plt.xlabel("Fraction of Missing Data")
        plt.ylabel("Count")
        
        outfile = f"{output_prefix}.png"
        plt.tight_layout()
        plt.savefig(outfile, dpi=300)
        plt.close()
        print(f"[Success] Plot saved to {outfile}")
    except Exception as e:
        print(f"[Error] Failed to plot site missingness: {e}")


def load_vmiss(filepath):
    """
    Loads variant missingness file (PLINK .vmiss format).
    Extracts Position from ID column if 'Position' column doesn't exist but ID does (e.g. 2-952).
    """
    print(f"[Info] Loading VMISS: {filepath}")
    df = load_df_from_space_sep(filepath)
    if df is None: return None

    # Clean header (remove #)
    df.columns = [c.replace('#', '') for c in df.columns]

    if 'Position' not in df.columns:
        if 'ID' in df.columns:
            # Try to extract position from 'ID' column (e.g., 'Chr-Pos')
            # Assuming ID format: Chr-Pos-...
            try:
                df['Position'] = df['ID'].str.split('-').str[1].astype(int)
            except Exception:
                print("[Warning] Could not extract Position from ID column. Merging might fail if 'Position' is missing.")
        else:
             print("[Error] VMISS file lacks 'Position' or 'ID' column.")

    return df

