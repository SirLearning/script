from infra.utils.io import load_df_from_tsv, load_df_from_space_sep
from infra.wheat.ref_v1 import get_chromosome, get_pos_on_chromosome

def load_df_from_plink_variant(filepath):
    """
    Loads variant missingness file (PLINK .vmiss format).
    Extracts Position from ID column if 'Position' column doesn't exist but ID does (e.g. 2-952).
    """
    print(f"[Info] Loading VMISS: {filepath}")
    df = load_df_from_tsv(filepath)
    if df is None: return None

    # Clean header (remove #)
    df.columns = [c.replace('#', '') if isinstance(c, str) else c for c in df.columns]
    
    if 'ID' in df.columns:
        # Try to extract position from 'ID' column (e.g., 'Chr-Pos')
        # Assuming ID format: Chr-Pos-...
        try:
            # Ensure ID parts are integers
            df['chr_sep_id'] = df['ID'].astype(str).str.split('-').str[0].astype(int)
            df['chr_sep_pos'] = df['ID'].astype(str).str.split('-').str[1].astype(int)

            # Get Chromosome name
            df['Chromosome'] = df['chr_sep_id'].apply(get_chromosome)
            
            # Calculate linear Position (requires passing both ID and local pos to the function)
            # get_pos_on_chromosome is not vectorized, so use apply(axis=1)
            df['Position'] = df.apply(lambda row: get_pos_on_chromosome(row['chr_sep_id'], row['chr_sep_pos']), axis=1)
        except Exception as e:
            print(f"[Warning] Could not extract Position from ID column: {e}")
    else:
            pass # Columns might just be CHROM POS etc.

    return df


def load_df_from_plink_hardy(filepath):
    print(f"[Info] Loading Hardy: {filepath}")
    df = load_df_from_space_sep(filepath)
    if df is None: return None
    df.columns = [c.replace('#', '') if isinstance(c, str) else c for c in df.columns]
    return df


def load_df_from_plink_gcount(filepath):
    print(f"[Info] Loading GCount: {filepath}")
    df = load_df_from_space_sep(filepath)
    if df is None: return None
    df.columns = [c.replace('#', '') if isinstance(c, str) else c for c in df.columns]
    return df


if __name__ == "__main__":
    # Example usage
    test_file = "/data1/dazheng_tusr1/vmap4.VCF.v1/test_plink/process/variant/A.info.vmiss"
    df = load_df_from_plink_variant(test_file)
    print(df.head())
