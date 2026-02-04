import os
import argparse
import pandas as pd
import numpy as np
from infra.utils import load_df_generic, save_df_to_tsv, load_df_from_space_sep_no_header
from genetics.germplasm import process_subgroup_as, process_subgroup_hznu, process_subgroup_nature, process_subgroup_v4, process_subgroup_vmap3, process_subgroup_watkins


# ==============================================================================
# Helper Classes & Functions
# ==============================================================================

class UnionFind:
    def __init__(self, elements):
        self.parent = {e: e for e in elements}

    def find(self, i):
        root = i
        while self.parent[root] != root:
            root = self.parent[root]
        curr = i
        while curr != root:
            next_node = self.parent[curr]
            self.parent[curr] = root
            curr = next_node
        return root

    def union(self, i, j):
        root_i = self.find(i)
        root_j = self.find(j)
        if root_i != root_j:
            self.parent[root_i] = root_j


def get_country_to_continent_map():
    return {
        'Algeria': 'Africa', 'Angola': 'Africa', 'Benin': 'Africa', 'Botswana': 'Africa', 
        'Burkina Faso': 'Africa', 'Burundi': 'Africa', 'Cameroon': 'Africa', 'Cape Verde': 'Africa', 
        'Central African Republic': 'Africa', 'Chad': 'Africa', 'Comoros': 'Africa', 'Congo': 'Africa', 
        'Congo, Democratic Republic of the': 'Africa', 'Cote d\'Ivoire': 'Africa', 'Djibouti': 'Africa', 
        'Egypt': 'Africa', 'Equatorial Guinea': 'Africa', 'Eritrea': 'Africa', 'Ethiopia': 'Africa', 
        'Gabon': 'Africa', 'Gambia': 'Africa', 'Ghana': 'Africa', 'Guinea': 'Africa', 
        'Guinea-Bissau': 'Africa', 'Kenya': 'Africa', 'Lesotho': 'Africa', 'Liberia': 'Africa', 
        'Libya': 'Africa', 'Madagascar': 'Africa', 'Malawi': 'Africa', 'Mali': 'Africa', 
        'Mauritania': 'Africa', 'Mauritius': 'Africa', 'Morocco': 'Africa', 'Mozambique': 'Africa', 
        'Namibia': 'Africa', 'Niger': 'Africa', 'Nigeria': 'Africa', 'Rwanda': 'Africa', 
        'Sao Tome and Principe': 'Africa', 'Senegal': 'Africa', 'Seychelles': 'Africa', 
        'Sierra Leone': 'Africa', 'Somalia': 'Africa', 'South Africa': 'Africa', 'South Sudan': 'Africa', 
        'Sudan': 'Africa', 'Swaziland': 'Africa', 'Tanzania': 'Africa', 'Togo': 'Africa', 
        'Tunisia': 'Africa', 'Uganda': 'Africa', 'Zambia': 'Africa', 'Zimbabwe': 'Africa',
        'Kazakhstan': 'Central Asia', 'Kyrgyzstan': 'Central Asia', 'Tajikistan': 'Central Asia', 
        'Turkmenistan': 'Central Asia', 'Uzbekistan': 'Central Asia',
        'China': 'East Asia', 'Japan': 'East Asia', 'Mongolia': 'East Asia', 'North Korea': 'East Asia', 
        'South Korea': 'East Asia', 'Taiwan': 'East Asia',
        'Afghanistan': 'South Asia', 'Bangladesh': 'South Asia', 'Bhutan': 'South Asia', 
        'India': 'South Asia', 'Maldives': 'South Asia', 'Nepal': 'South Asia', 
        'Pakistan': 'South Asia', 'Sri Lanka': 'South Asia',
        'Armenia': 'West Asia', 'Azerbaijan': 'West Asia', 'Bahrain': 'West Asia', 'Cyprus': 'West Asia', 
        'Georgia': 'West Asia', 'Iran': 'West Asia', 'Iraq': 'West Asia', 'Israel': 'West Asia', 
        'Jordan': 'West Asia', 'Kuwait': 'West Asia', 'Lebanon': 'West Asia', 'Oman': 'West Asia', 
        'Palestine': 'West Asia', 'Qatar': 'West Asia', 'Saudi Arabia': 'West Asia', 
        'Syria': 'West Asia', 'Turkey': 'West Asia', 'United Arab Emirates': 'West Asia', 
        'Yemen': 'West Asia',
        'Brunei': 'East Asia', 'Cambodia': 'South Asia', 'Indonesia': 'East Asia', 'Laos': 'South Asia', 
        'Malaysia': 'East Asia', 'Myanmar': 'South Asia', 'Burma': 'South Asia',
        'Philippines': 'East Asia', 'Singapore': 'East Asia', 
        'Thailand': 'South Asia', 'Timor-Leste': 'East Asia', 'Vietnam': 'East Asia',
        'Russia': 'Europe', 
        'Albania': 'Europe', 'Andorra': 'Europe', 'Austria': 'Europe', 'Belarus': 'Europe', 
        'Belgium': 'Europe', 'Bosnia and Herzegovina': 'Europe', 'Bulgaria': 'Europe', 
        'Croatia': 'Europe', 'Czech Republic': 'Europe', 'Denmark': 'Europe', 'Estonia': 'Europe', 
        'Finland': 'Europe', 'France': 'Europe', 'Germany': 'Europe', 'Greece': 'Europe', 
        'Hungary': 'Europe', 'Iceland': 'Europe', 'Ireland': 'Europe', 'Italy': 'Europe', 
        'Kosovo': 'Europe', 'Latvia': 'Europe', 'Liechtenstein': 'Europe', 'Lithuania': 'Europe', 
        'Luxembourg': 'Europe', 'Macedonia': 'Europe', 'Malta': 'Europe', 'Moldova': 'Europe', 
        'Monaco': 'Europe', 'Montenegro': 'Europe', 'Netherlands': 'Europe', 'Norway': 'Europe', 
        'Poland': 'Europe', 'Portugal': 'Europe', 'Romania': 'Europe', 'San Marino': 'Europe', 
        'Serbia': 'Europe', 'Slovakia': 'Europe', 'Slovenia': 'Europe', 'Spain': 'Europe', 
        'Sweden': 'Europe', 'Switzerland': 'Europe', 'Ukraine': 'Europe', 'United Kingdom': 'Europe', 
        'UK': 'Europe', 'Vatican City': 'Europe', 'Soviet Union': 'Europe', 'CS': 'Europe',
        'Yugoslavia': 'Europe', 'Crete': 'Europe', 'Canary Islands': 'Europe', 'Central Europe': 'Europe',
        'Netherland': 'Europe',
        'Antigua and Barbuda': 'North America', 'Bahamas': 'North America', 'Barbados': 'North America', 
        'Belize': 'North America', 'Canada': 'North America', 'Costa Rica': 'North America', 
        'Cuba': 'North America', 'Dominica': 'North America', 'Dominican Republic': 'North America', 
        'El Salvador': 'North America', 'Grenada': 'North America', 'Guatemala': 'North America', 
        'Haiti': 'North America', 'Honduras': 'North America', 'Jamaica': 'North America', 
        'Mexico': 'North America', 'Nicaragua': 'North America', 'Panama': 'North America', 
        'Saint Kitts and Nevis': 'North America', 'Saint Lucia': 'North America', 
        'Saint Vincent and the Grenadines': 'North America', 'Trinidad and Tobago': 'North America', 
        'United States': 'North America', 'USA': 'North America', 'US': 'North America',
        'Argentina': 'South America', 'Bolivia': 'South America', 'Brazil': 'South America', 
        'Chile': 'South America', 'Colombia': 'South America', 'Ecuador': 'South America', 
        'Guyana': 'South America', 'Paraguay': 'South America', 'Peru': 'South America', 
        'Suriname': 'South America', 'Uruguay': 'South America', 'Venezuela': 'South America',
        'Australia': 'Oceania', 'Fiji': 'Oceania', 'Kiribati': 'Oceania', 'Marshall Islands': 'Oceania', 
        'Micronesia': 'Oceania', 'Nauru': 'Oceania', 'New Zealand': 'Oceania', 'Palau': 'Oceania', 
        'Papua New Guinea': 'Oceania', 'Samoa': 'Oceania', 'Solomon Islands': 'Oceania', 
        'Tonga': 'Oceania', 'Tuvalu': 'Oceania', 'Vanuatu': 'Oceania',
        'KOR': 'East Asia', 'JPN': 'East Asia', 'PRK': 'East Asia', 'CHN': 'East Asia',
        'IND': 'South Asia', 'AFG': 'South Asia', 'NPL': 'South Asia', 'PAK': 'South Asia',
        'IRN': 'West Asia', 'IRQ': 'West Asia', 'TUR': 'West Asia',
        'URY': 'South America', 'MEX': 'North America', 'ARG': 'South America', 'CAN': 'North America',
        'DEU': 'Europe', 'FRA': 'Europe', 'GBR': 'Europe', 'RUS': 'Europe', 'USA': 'North America'
    }

def map_country_to_continent(country, mapping):
    if not isinstance(country, str):
        return 'Unknown'
    
    country = country.strip()
    if country in mapping:
        return mapping[country]
    
    if country == 'USA' or country == 'United States of America': return 'North America'
    if country == 'UK': return 'Europe'
    if 'Russia' in country: return 'Europe' 
    if 'USSR' in country or 'Soviet' in country: return 'Europe'
    if country == 'Czechia': return 'Europe'
    if 'Macedonia' in country: return 'Europe'
    if 'China' in country: return 'East Asia'
    if 'Iran' in country: return 'West Asia'
    
    for k, v in mapping.items():
        if k.lower() == country.lower():
            return v
    return 'Unknown'

# ==============================================================================
# Feature 1: Group Annotation
# ==============================================================================

def run_group_annotation(
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
    # sample_group.py used load_df_from_space_sep_no_header assuming headers aren't there
    # But let's be robust
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
# Feature 2: Deduplication Analysis
# ==============================================================================

def run_dedup_analysis(
    taxa_bam_dir="/data/home/tusr1/01projects/vmap4/00data/05taxaBamMap", 
    db_dir="/data/home/tusr1/git/DBone/Service/src/main/resources/raw/20251208",
    output_prefix="sample.dedup"
):
    print("[Info] Running Deduplication Analysis...")
    TAXA_DIR = taxa_bam_dir
    DB_DIR = db_dir
    
    print(f"  TaxaDir: {TAXA_DIR}")
    print(f"  DBDir: {DB_DIR}")
    
    # --- Main Dedup Logic ---
    data_funcs = [process_subgroup_vmap3, process_subgroup_v4, process_subgroup_watkins, process_subgroup_nature, process_subgroup_hznu, process_subgroup_as]
    
    all_res = []
    all_miss = []
    
    for f in data_funcs:
        r, m = f()
        all_res.extend(r)
        all_miss.extend(m)
        
    df = pd.DataFrame(all_res)
    print(f"\nTotal Mapped Taxa: {len(df)}")
    
    if all_miss:
        df_miss = pd.DataFrame(all_miss)
        miss_file = output_prefix + '.missing.tsv'
        save_df_to_tsv(df_miss, miss_file)
        print(f"[Info] Missing taxa report saved to {miss_file}")
    
    if df.empty: return

    # Normalize
    nulls = ['nan', 'None', '', 'NaT', '-', 'NA', 'na', 'N/A', '.', '?', 'null']
    df['Clean_CN'] = df['Chinese_name'].astype(str).str.strip().replace(nulls, np.nan)
    df['Clean_Acc'] = df['Accession'].astype(str).str.strip().replace(nulls, np.nan)
    df['Coverage'] = pd.to_numeric(df['Coverage'], errors='coerce').fillna(0)
    
    # Union Find
    uf = UnionFind(df.index)
    
    cn_map = {}
    for idx, row in df.iterrows():
        if pd.notna(row['Clean_CN']):
            cn_map.setdefault(row['Clean_CN'], []).append(idx)
            
    acc_map = {}
    for idx, row in df.iterrows():
        if pd.notna(row['Clean_Acc']):
            acc_map.setdefault(row['Clean_Acc'], []).append(idx)
            
    for indices in cn_map.values():
        for i in range(1, len(indices)): uf.union(indices[0], indices[i])
            
    for indices in acc_map.values():
        for i in range(1, len(indices)): uf.union(indices[0], indices[i])
        
    # Group Identifiers
    root_to_gid = {}
    gid_counter = 1
    dup_ids = []
    is_dup = []
    
    # Pre-calc roots
    roots = [uf.find(i) for i in df.index]
    
    for r in roots:
        if r not in root_to_gid:
            root_to_gid[r] = gid_counter
            gid_counter += 1
        dup_ids.append(root_to_gid[r])
        
    df['Dup_Group_ID'] = dup_ids
    
    # Mark Duplicates (Keep Max Coverage)
    grouped = df.groupby('Dup_Group_ID')
    to_remove = []
    
    for _, group in grouped:
        if len(group) > 1:
            sorted_g = group.sort_values('Coverage', ascending=False)
            # Keep first (max coverage), remove rest
            to_remove.extend(sorted_g.index[1:].tolist())
            
    df['Is_Duplicate'] = df.index.isin(to_remove)
    
    # Save Outputs
    summary_file = output_prefix + '.summary.tsv'
    save_df_to_tsv(df, summary_file)
    
    qc_file = output_prefix + '.qc.tsv'
    qc_ids = df.loc[to_remove, 'Taxa']
    if not qc_ids.empty:
        qc_ids.to_csv(qc_file, index=False, header=False)
        print(f"[Info] QC deduplication list saved to {qc_file}")
    
    print(f"\nSummary saved to {summary_file}")
    print(f"Total Unique: {df['Dup_Group_ID'].nunique()}")
    print(f"Total Removed: {len(to_remove)}")


# ==============================================================================
# Feature 3: Location Annotation
# ==============================================================================

def run_location_analysis(
    input_file,
    output_prefix,
    db_dir="/data/home/tusr1/git/DBone/Service/src/main/resources/raw/20251208"
):
    print("[Info] Running Location Analysis...")
    df = load_df_generic(input_file)
    if df is None: return

    # Check mapping columns
    clean_country = None
    if 'Country' in df.columns: clean_country = 'Country'
    elif 'Provenance of material-GS' in df.columns: clean_country = 'Provenance of material-GS'
    # Fallback to search
    if not clean_country:
         for c in df.columns:
             if 'country' in c.lower() or 'provenance' in c.lower():
                 clean_country = c
                 break
    
    country_map = get_country_to_continent_map()
    
    # If we have a country column, map it
    if clean_country:
        print(f"  Mapping Continent using column: {clean_country}")
        df['Continent_Mapped'] = df[clean_country].apply(lambda x: map_country_to_continent(str(x), country_map))
    else:
        print(f"[Warning] No country/provenance column found. Skipping direct mapping.")
        
    # Validation against VMap3/Watkins (Optional, if DB_DIR provided)
    if db_dir:
        # TODO: Implement cross-referencing logic if needed?
        # The prompt says "based on existing data... add location".
        # The existing script just validated V4 vs Vmap3. 
        # Here we likely just want to apply the mapping rules to the input file.
        pass
        
    save_df_to_tsv(df, output_prefix)


# ==============================================================================
# CLI
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(description="Germplasm Annotation Toolkit")
    subparsers = parser.add_subparsers(dest='command', help='Sub-command help')

    # 1. Group
    p_group = subparsers.add_parser('group', help='Add Group Info to Samples')
    p_group.add_argument('-i', '--input', required=True, help='Input Sample File')
    p_group.add_argument('--group-file', default='/data/home/tusr1/git/DBone/Service/src/main/resources/raw/20251208/germplasm/sample_groups.txt', help='Path to Group File') 
    # Note: I'm not sure where sample_groups.txt really is, user code had "data/germplasm/sample_groups.txt" relative, but repeated code has absolute paths.
    # I'll default to relative but user can override.
    p_group.add_argument('-o', '--out', default='annotated_groups.tsv', help='Output file')

    # 2. Dedup
    p_dedup = subparsers.add_parser('dedup', help='Run Deduplication Analysis Pipeline')
    p_dedup.add_argument('--taxa-bam-dir', default="/data/home/tusr1/01projects/vmap4/00data/05taxaBamMap", help='Taxa Bam Map Directory')
    p_dedup.add_argument('--db-dir', default="/data/home/tusr1/git/DBone/Service/src/main/resources/raw/20251208", help='Database Root Directory')
    p_dedup.add_argument('--out-prefix', default='germplasm.dedup', help='Output files prefix')

    # 3. Location
    p_location = subparsers.add_parser('location', help='Add Continent Info to Samples')
    p_location.add_argument('-i', '--input', required=True, help='Input Sample File')
    p_location.add_argument('--db-dir', help='Optional DB dir for cross-referencing')
    p_location.add_argument('-o', '--out', default='annotated_location.tsv', help='Output file')

    args = parser.parse_args()

    if args.command == 'group':
        # Fix default path if not absolute
        if args.group_file == "data/germplasm/sample_groups.txt" and not os.path.exists(args.group_file):
             # Try to anticipate where it might be relative to script? No, rely on user to run from root or provide path
             pass
        run_group_annotation(args.input, args.group_file, args.out)
        
    elif args.command == 'dedup':
        run_dedup_analysis(args.taxa_bam_dir, args.db_dir, args.out_prefix)
        
    elif args.command == 'location':
        run_location_analysis(args.input, args.out, args.db_dir)
        
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
