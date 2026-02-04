import os
import argparse
import pandas as pd
from infra.utils import load_df_generic, save_df_to_tsv, load_df_from_space_sep_no_header

# ==============================================================================
# Feature 1: Group Annotation
# ==============================================================================

def anno_group(
    input_file,
    group_file="/data1/dazheng_tusr1/vmap4.VCF.v1/sample_groups.txt"
) -> pd.DataFrame:
    print(f"[Info] Running Group Annotation...")
    print(f"  Input: {input_file}")
    print(f"  Group File: {group_file}")
    
    # Load Input
    df = load_df_generic(input_file)
    if df is None:
        print("[Error] Failed to load input file.")
        return None

    # Load Groups
    df_group = load_df_from_space_sep_no_header(group_file, ['Sample', 'Group'])
    
    if df_group is not None:
        # Check if input has 'Sample' or 'Taxa'
        join_col = 'Sample'
        if 'Sample' not in df.columns and 'Taxa' in df.columns:
             join_col = 'Taxa'
        if join_col not in df.columns: # fallback if index
             print(f"[Warning] Could not find join column (Sample/Taxa) in input. Using first column.")
             join_col = df.columns[0]

        # Rename group file column to match if needed, but it is fixed as [Sample, Group]
        # So we merge left on join_col = Sample
        
        merged = pd.merge(df, df_group, left_on=join_col, right_on='Sample', how='left')
        merged['Group'] = merged['Group'].fillna('Unknown')
        
        # Cleanup: If we joined on Taxa, we might have duplicate Sample column now
        if join_col != 'Sample':
            merged = merged.drop(columns=['Sample']) 
        
        # Generate output filename
        if input_file.lower().endswith(('.tsv', '.csv', '.txt')):
             out_path = input_file.rsplit('.', 1)[0] + '.grouped.tsv'
        else:
             out_path = input_file + '.grouped.tsv'

        save_df_to_tsv(merged, out_path)
        
        return merged
    else:
        print("[Error] Failed to load group file.")
        return None

# ==============================================================================
# Feature 2: Duplication Annotation (Consumer of sample_ana.run_dedup_analysis)
# ==============================================================================

def anno_duplication(
    input_file,
    dedup_result_file
):
    print(f"[Info] Running Duplication Annotation...")
    print(f"  Input: {input_file}")
    print(f"  Dedup Result File: {dedup_result_file}")
    
    df = load_df_generic(input_file)
    if df is None: return

    df_dedup = load_df_generic(dedup_result_file)
    if df_dedup is None: return
    
    # Dedup result usually has 'Taxa' as identifier
    # We try to join input on 'Taxa' or 'Sample' against 'Taxa' in dedup file
    
    join_col_in = 'Taxa'
    if 'Taxa' not in df.columns and 'Sample' in df.columns:
        join_col_in = 'Sample'
    
    # Fallback
    if join_col_in not in df.columns:
        # Try finding a column that looks like ID
        join_col_in = df.columns[0]
        
    join_col_dedup = 'Taxa'
    if 'Taxa' not in df_dedup.columns:
        # Maybe Sample?
        if 'Sample' in df_dedup.columns: join_col_dedup = 'Sample'
        else: join_col_dedup = df_dedup.columns[0]

    print(f"  Merging on Input[{join_col_in}] == Dedup[{join_col_dedup}]")
    
    # Columns to bring in: Dup_Group_ID, Is_Duplicate, Clean_CN, Clean_Acc
    cols_to_add = ['Dup_Group_ID', 'Is_Duplicate', 'Clean_CN', 'Clean_Acc']
    # Filter to existing
    cols_to_add = [c for c in cols_to_add if c in df_dedup.columns]
    
    # Create subset for merging
    right_df = df_dedup[[join_col_dedup] + cols_to_add].drop_duplicates()
    
    merged = pd.merge(df, right_df, left_on=join_col_in, right_on=join_col_dedup, how='left')
    
    # Cleanup duplicate key col if names differed
    if join_col_in != join_col_dedup:
        merged = merged.drop(columns=[join_col_dedup])

    # Generate output
    if input_file.lower().endswith(('.tsv', '.csv', '.txt')):
             out_path = input_file.rsplit('.', 1)[0] + '.dedup_anno.tsv'
    else:
             out_path = input_file + '.dedup_anno.tsv'
             
    save_df_to_tsv(merged, out_path)
    print(f"[Info] Annotation saved to {out_path}")


# ==============================================================================
# Feature 3: Location Annotation (Consumer of sample_ana.run_location_analysis)
# ==============================================================================

def anno_location(
    input_file,
    location_result_file
):
    print(f"[Info] Running Location Annotation...")
    print(f"  Input: {input_file}")
    print(f"  Location Result File: {location_result_file}")
    
    df = load_df_generic(input_file)
    if df is None: return

    df_loc = load_df_generic(location_result_file)
    if df_loc is None: return
    
    # Similar join logic
    join_col_in = 'Taxa'
    if 'Taxa' not in df.columns and 'Sample' in df.columns:
        join_col_in = 'Sample'
    if join_col_in not in df.columns: join_col_in = df.columns[0]
    
    join_col_loc = 'Taxa'
    # Try finding matching column in loc file
    if join_col_loc not in df_loc.columns:
        if 'Sample' in df_loc.columns: join_col_loc = 'Sample'
        else: join_col_loc = df_loc.columns[0]
        
    print(f"  Merging on Input[{join_col_in}] == Location[{join_col_loc}]")
    
    cols_to_add = ['Continent_Mapped']
    cols_to_add = [c for c in cols_to_add if c in df_loc.columns]
    
    right_df = df_loc[[join_col_loc] + cols_to_add].drop_duplicates()
    
    merged = pd.merge(df, right_df, left_on=join_col_in, right_on=join_col_loc, how='left')
    
    if join_col_in != join_col_loc:
        merged = merged.drop(columns=[join_col_loc])

    if input_file.lower().endswith(('.tsv', '.csv', '.txt')):
             out_path = input_file.rsplit('.', 1)[0] + '.loc_anno.tsv'
    else:
             out_path = input_file + '.loc_anno.tsv'
             
    save_df_to_tsv(merged, out_path)
    print(f"[Info] Annotation saved to {out_path}")


# ==============================================================================
# CLI
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(description="Germplasm Annotation Toolkit")
    subparsers = parser.add_subparsers(dest='command', help='Sub-command help')

    # 1. Group
    p_group = subparsers.add_parser('group', help='Add Group Info')
    p_group.add_argument('-i', '--input', required=True, help='Input Sample File')
    p_group.add_argument('--group-file', default='/data1/dazheng_tusr1/vmap4.VCF.v1/sample_groups.txt', help='Path to Group File')
    
    # 2. Dedup Annotation
    p_dedup = subparsers.add_parser('dedup', help='Annotate with Duplication Info')
    p_dedup.add_argument('-i', '--input', required=True, help='Input Sample File')
    p_dedup.add_argument('--dedup-file', required=True, help='Path to Dedup Analysis Result (summary.tsv)')
    
    # 3. Location Annotation
    p_loc = subparsers.add_parser('location', help='Annotate with Location Info')
    p_loc.add_argument('-i', '--input', required=True, help='Input Sample File')
    p_loc.add_argument('--loc-file', required=True, help='Path to Location Analysis Result')

    args = parser.parse_args()

    if args.command == 'group':
        anno_group(args.input, args.group_file)
        
    elif args.command == 'dedup':
        anno_duplication(args.input, args.dedup_file)
        
    elif args.command == 'location':
        anno_location(args.input, args.loc_file)
        
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
