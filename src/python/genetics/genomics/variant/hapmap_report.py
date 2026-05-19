"""Copy precomputed hapmap tables for publish."""

from infra.utils.io import load_df_generic, save_df_to_tsv


def report_hapmap_table(hapmap_tsv_path, output_prefix):
    df = load_df_generic(hapmap_tsv_path)
    if df is None or df.empty:
        raise ValueError(f'Empty hapmap table: {hapmap_tsv_path}')
    save_df_to_tsv(df, f'{output_prefix}.hapmap.tsv')
