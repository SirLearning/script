from turtle import pd
from infra.utils import load_df_from_space_sep_no_header


def integrate_group_info(
    group_file="data/germplasm/sample_groups.txt", 
    df_merged=None
):
    df_group = load_df_from_space_sep_no_header(group_file, ['Sample', 'Group'])
    if df_group is not None:
        df_merged = pd.merge(df_merged, df_group, on='Sample', how='left')
        df_merged['Group'] = df_merged['Group'].fillna('Unknown')
    else:
        df_merged['Group'] = 'Unknown'
    return df_merged


