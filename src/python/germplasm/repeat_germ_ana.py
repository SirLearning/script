import pandas as pd
import os
import glob
import re

# Configuration
TAXA_BAM_MAP_DIR = "/data/home/tusr1/01projects/vmap4/00data/05taxaBamMap"
DATABASE_DIR = "/data/home/tusr1/git/DBone/Service/src/main/resources/raw/20251208"

class UnionFind:
    def __init__(self, elements):
        self.parent = {e: e for e in elements}

    def find(self, i):
        # Iterative find with path compression to avoid recursion depth issues
        root = i
        while self.parent[root] != root:
            root = self.parent[root]
        
        # Path compression step
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

def read_taxa_list(filename):
    """
    Reads taxaBamMap file and returns a dictionary of {taxa: coverage}.
    """
    path = os.path.join(TAXA_BAM_MAP_DIR, filename)
    if not os.path.exists(path):
        print(f"Warning: {filename} not found.")
        return {}
    try:
        # taxaBamMap format: ID \t coverage \t path
        # Read first two columns: ID and Coverage
        df = pd.read_csv(path, sep=r'\s+', header=None, usecols=[0, 1])
        # Ensure coverage is numeric
        ids = df[0].astype(str)
        coverages = pd.to_numeric(df[1], errors='coerce').fillna(0)
        return dict(zip(ids, coverages))
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return {}

def safe_read_excel(path, sheet_name=0, header=0):
    try:
        return pd.read_excel(path, sheet_name=sheet_name, header=header)
    except Exception as e:
        print(f"Error reading {path} sheet {sheet_name}: {e}")
        return pd.DataFrame()

def process_vmap3_group():
    print("Processing VMap3 Group (D, AB, ABD)...")
    
    storage_path = os.path.join(DATABASE_DIR, "omics/LuLab/Lulab_germplasm_Storage.xlsx")
    
    # 1.BamDatabase
    print("  Loading BamDatabase...")
    df_bam_db = safe_read_excel(storage_path, sheet_name='BamDatabase')
    
    taxa_to_bampath = {}
    if not df_bam_db.empty:
        df_bam_db.columns = df_bam_db.columns.str.strip()
        if 'BamDataBaseID' in df_bam_db.columns and 'Bam-Path' in df_bam_db.columns:
            taxa_to_bampath = df_bam_db.set_index('BamDataBaseID')['Bam-Path'].to_dict()
    
    # 2. Map Sheets: AB_taxa_bam, ABD_taxa_bam, D_taxa_bam
    print("  Loading Taxa Bam Sheets...")
    map_sheets = ['AB_taxa_bam', 'ABD_taxa_bam', 'D_taxa_bam']
    bampath_to_bamtaxaname = {} # Path -> Name (Col 0)
    bampath_to_accession = {}   # Path -> Accession (New Logic)
    
    for sheet in map_sheets:
        # Read with header to identify columns
        df_sheet = safe_read_excel(storage_path, sheet_name=sheet, header=0)
        if df_sheet.empty: continue
        
        # Identify Columns
        col_acc = next((c for c in df_sheet.columns if 'Accession' in str(c)), None)
        col_name = df_sheet.columns[0] # Assuming first col is Name
        
        for idx, row in df_sheet.iterrows():
            if pd.isna(row[col_name]): continue
            name = str(row[col_name]).strip()
            
            acc_val = None
            if col_acc and pd.notna(row[col_acc]):
                acc_val = str(row[col_acc]).strip()
            
            # Paths are in cols 2, 3, 4 (indices)
            # Since header=0, iloc[2] is the 3rd column
            for col_idx in [2, 3, 4]:
                if col_idx < len(row):
                    val = row.iloc[col_idx]
                    if pd.notna(val):
                        path_val = str(val).strip()
                        bampath_to_bamtaxaname[path_val] = name
                        if acc_val:
                            bampath_to_accession[path_val] = acc_val

    # 3. Vmap3 Final
    print("  Loading Vmap3 Final Excel...")
    vmap3_path = os.path.join(DATABASE_DIR, "germplasm/LuLab/Vmap3最终版3.0.xlsx")
    df_vmap3 = safe_read_excel(vmap3_path)
    
    vmap3_lookup = {}
    vmap3_val_lookup = {} # Value -> Info
    
    if not df_vmap3.empty:
        df_vmap3.columns = df_vmap3.columns.str.strip()
        vmap3_id_col = df_vmap3.columns[0]
        col_accession_target = '合并PI(1 PI;2 Taxa;34 English_name; Sample_name)'
        col_chinese_target = 'Chinese_name'
        
        # Extended search columns
        search_cols = ['PI_accession', 'English_name', 'Sample_name', 'Chinese_name']
        valid_search_files = [c for c in search_cols if c in df_vmap3.columns]
        
        for idx, row in df_vmap3.iterrows():
            if pd.isna(row[vmap3_id_col]): continue
            rid = str(row[vmap3_id_col]).strip()
            acc = row[col_accession_target] if col_accession_target in df_vmap3.columns else None
            cn = row[col_chinese_target] if col_chinese_target in df_vmap3.columns else None
            
            info = {'Accession': acc, 'Chinese_name': cn}
            vmap3_lookup[rid] = info
            
            # Populate value lookup
            for c in valid_search_files:
                val = row[c]
                if pd.notna(val):
                    v_str = str(val).strip()
                    vmap3_val_lookup[v_str] = info

    # 4. ab_new_code.txt Loading (For NotInBamDatabase)
    print("  Loading ab_new_code.txt...")
    ab_new_path = os.path.join(DATABASE_DIR, "germplasm/LuLab/ab_new_code.txt")
    ab_new_lookup = {}
    if os.path.exists(ab_new_path):
        try:
             # Header in line 2 -> skiprows=1
             df_ab = pd.read_csv(ab_new_path, sep=',', skiprows=1)
             df_ab.columns = df_ab.columns.str.strip()
             if '测序编号' in df_ab.columns and 'Accessions' in df_ab.columns:
                 for idx, row in df_ab.iterrows():
                     k = row['测序编号']
                     v = row['Accessions'] # Fix typo assuming user meant Accessions
                     if pd.notna(k) and pd.notna(v):
                         ab_new_lookup[str(k).strip()] = str(v).strip()
        except Exception as e:
            print(f"Error reading ab_new_code.txt: {e}")

    results = []
    missing = []
    
    taxa_files_map = {
        'AB': 'AB.taxaBamMap.txt', 
        'ABD': 'ABD.taxaBamMap.txt', 
        'D': 'D.taxaBamMap.txt'
    }
    
    for subgroup, filename in taxa_files_map.items():
        taxa_map = read_taxa_list(filename)
        for taxa, coverage in taxa_map.items():
            path = taxa_to_bampath.get(taxa)
            
            # Case 1: NotInBamDatabase -> Check ab_new_code
            if not path:
                # Check taxa or taxa without 'A' suffix (if ends with A)
                found = ab_new_lookup.get(taxa)
                if not found and taxa.endswith('A'):
                    found = ab_new_lookup.get(taxa[:-1])
                
                if found:
                     results.append({
                        'Taxa': taxa,
                        'Accession': found,
                        'Chinese_name': None,
                        'Group': subgroup,
                        'Coverage': coverage
                    })
                else:
                    missing.append({'Taxa': taxa, 'Group': subgroup, 'Reason': 'NotInBamDatabase'})
                continue
                
            path_str = str(path).strip()
            bamtaxa_name = bampath_to_bamtaxaname.get(path_str)
            
            # Case 2: NotInBamSheets -> Check path->Accession loopup in Vmap3
            if not bamtaxa_name:
                mapped_acc = bampath_to_accession.get(path_str)
                if mapped_acc:
                    # Validate against Vmap3 extended columns
                    v_info = vmap3_val_lookup.get(mapped_acc)
                    if v_info:
                        results.append({
                            'Taxa': taxa,
                            'Accession': v_info['Accession'],
                            'Chinese_name': v_info['Chinese_name'],
                            'Group': subgroup,
                            'Coverage': coverage
                        })
                    else:
                        missing.append({'Taxa': taxa, 'Group': subgroup, 'Reason': 'BamSheetAccNotInVmap3'})
                else:
                    missing.append({'Taxa': taxa, 'Group': subgroup, 'Reason': 'NotInBamSheets'})
                continue
            
            # Case 3: Standard Lookup
            clean_name = bamtaxa_name
            if clean_name.startswith("V3_"):
                clean_name = clean_name[3:]
                
            info = vmap3_lookup.get(clean_name)
            if info:
                results.append({
                    'Taxa': taxa,
                    'Accession': info['Accession'],
                    'Chinese_name': info['Chinese_name'],
                    'Group': subgroup,
                    'Coverage': coverage
                })
            else:
                missing.append({'Taxa': taxa, 'Group': subgroup, 'Reason': 'NotInVmap3Final'})
            
    return results, missing

def process_v4_group():
    print("Processing V4 Group...")
    taxa_map = read_taxa_list('WAP.taxaBamMap.txt')
    
    v4_path = os.path.join(DATABASE_DIR, "germplasm/LuLab/V4_germplasm.xlsx")
    df_v4 = safe_read_excel(v4_path)
    
    v4_lookup = {}
    if not df_v4.empty:
        df_v4.columns = df_v4.columns.str.strip()
        if 'R1' in df_v4.columns:
            df_v4['R1_str'] = df_v4['R1'].astype(str).str.strip()
            if df_v4['R1_str'].duplicated().any():
                df_v4 = df_v4.drop_duplicates(subset=['R1_str'], keep='first')
            v4_lookup = df_v4.set_index('R1_str').to_dict('index')
    
    results = []
    missing = []
    col_acc = 'Accessions'
    col_cn = 'ChineseName'
    
    for t, coverage in taxa_map.items():
        t_clean = t.strip()
        row = v4_lookup.get(t_clean)
        
        if not row and t_clean.endswith('A'):
            t_mod = t_clean[:-1]
            row = v4_lookup.get(t_mod)
            
        if row:
            acc = row.get(col_acc)
            cn = row.get(col_cn)
            results.append({
                'Taxa': t,
                'Accession': acc,
                'Chinese_name': cn,
                'Group': 'WAP',
                'Coverage': coverage
            })
        else:
             missing.append({'Taxa': t, 'Group': 'WAP', 'Reason': 'NotInV4Excel'})
            
    return results, missing

def process_watkins_group():
    print("Processing Watkins Group...")
    files = ['w115.taxaBamMap.txt', 'w203.taxaBamMap.txt', 'w204.taxaBamMap.txt', 'w243.taxaBamMap.txt', 'w66.taxaBamMap.txt']
    all_taxa_map = {}
    for f in files:
        all_taxa_map.update(read_taxa_list(f))
        
    omics_path = os.path.join(DATABASE_DIR, "omics/Watkins_CRA012590.xlsx")
    
    # Run
    df_run = safe_read_excel(omics_path, sheet_name='Run')
    run_map = {}
    if not df_run.empty:
        df_run.columns = df_run.columns.str.strip()
        if 'Accession' in df_run.columns and 'Experiment accession' in df_run.columns:
             run_map = df_run.set_index('Accession')['Experiment accession'].to_dict()

    # Experiment
    df_exp = safe_read_excel(omics_path, sheet_name='Experiment')
    exp_map = {}
    if not df_exp.empty:
        df_exp.columns = df_exp.columns.str.strip()
        if 'Accession' in df_exp.columns and 'BioSample name' in df_exp.columns:
            exp_map = df_exp.set_index('Accession')['BioSample name'].to_dict()
        
    # Watkins Germplasm
    watkins_path = os.path.join(DATABASE_DIR, "germplasm/Watkins/Watkins_sm.xlsx")
    df_wat = safe_read_excel(watkins_path, header=1)
    
    wat_map = {}
    if not df_wat.empty:
        df_wat.columns = df_wat.columns.str.strip()
        col_aname = 'Accession name'
        col_gru = 'GRU StorCode'
        col_cn = 'Accession name in China (WWWG2B)'
        
        if col_aname in df_wat.columns:
            for idx, row in df_wat.iterrows():
                if pd.isna(row[col_aname]): continue
                aname = str(row[col_aname]).strip()
                key = aname
                if key.endswith('*'):
                    key = key[:-1]
                info = {
                    'Accession': row.get(col_gru),
                    'Chinese_name': row.get(col_cn)
                }
                wat_map[key] = info
                wat_map[aname] = info

    results = []
    missing = []
    
    for t, coverage in all_taxa_map.items():
        exp_acc = run_map.get(t)
        if not exp_acc: 
            missing.append({'Taxa': t, 'Group': 'Watkins', 'Reason': 'NotInRunSheet'})
            continue
        
        biosample = exp_map.get(exp_acc)
        if not biosample:
            missing.append({'Taxa': t, 'Group': 'Watkins', 'Reason': 'NotInExpSheet'})
            continue
        
        b_str = str(biosample).strip()
        info = wat_map.get(b_str)
        if info:
            results.append({
                'Taxa': t,
                'Accession': info['Accession'],
                'Chinese_name': info['Chinese_name'],
                'Group': 'Watkins',
                'Coverage': coverage
            })
        else:
             missing.append({'Taxa': t, 'Group': 'Watkins', 'Reason': 'NotInWatkinsExcel'})
            
    return results, missing

def process_ipk_group():
    print("Processing IPK Group...")
    taxa_map = read_taxa_list('Nature.taxaBamMap.txt')
    
    # ipk_sample.xlsx
    ipk_sample_path = os.path.join(DATABASE_DIR, "germplasm/ipk/ipk_sample.xlsx")
    df_sample = safe_read_excel(ipk_sample_path)
    
    wgs_to_donor = {}
    if not df_sample.empty:
        df_sample.columns = df_sample.columns.str.strip()
        col_wgs = 'WGS Biosamples ID'
        col_donor = 'Name/GBIS Donor number'
        if col_wgs in df_sample.columns and col_donor in df_sample.columns:
            wgs_to_donor = df_sample.set_index(col_wgs)[col_donor].to_dict()
    
    # ipk_germ.xlsx
    ipk_germ_path = os.path.join(DATABASE_DIR, "germplasm/ipk/ipk_germ.xlsx")
    df_germ = safe_read_excel(ipk_germ_path)
    
    valid_donors = set()
    if not df_germ.empty:
        df_germ.columns = df_germ.columns.str.strip()
        col_donor = 'Name/GBIS Donor number'
        if col_donor in df_germ.columns:
            valid_donors = set(df_germ[col_donor].dropna().astype(str).str.strip())
        
    results = []
    missing = []
    
    for t, coverage in taxa_map.items():
        donor = wgs_to_donor.get(t)
        if donor:
            d_str = str(donor).strip()
            # Verify in germ file
            if d_str in valid_donors:
                results.append({
                    'Taxa': t,
                    'Accession': donor,
                    'Chinese_name': None,
                    'Group': 'Nature',
                    'Coverage': coverage
                })
            else:
                 missing.append({'Taxa': t, 'Group': 'Nature', 'Reason': 'NotInGermExcel'})
        else:
            missing.append({'Taxa': t, 'Group': 'Nature', 'Reason': 'NotInSampleExcel'})
                
    return results, missing

def process_hznu_group():
    print("Processing HZNU Group...")
    taxa_map = read_taxa_list('HZNU.taxaBamMap.txt')
    results = []
    # No missing tracking needed here as simple mapping
    for t, coverage in taxa_map.items():
        results.append({
            'Taxa': t,
            'Accession': t, 
            'Chinese_name': None,
            'Group': 'HZNU',
            'Coverage': coverage
        })
    return results, []

def process_as_group():
    print("Processing A & S Group...")
    
    as_path = os.path.join(DATABASE_DIR, "germplasm/LuLab/A&Sgenome.xlsx")
    df_as = safe_read_excel(as_path)
    
    lookup = {}
    if not df_as.empty:
        df_as.columns = df_as.columns.str.strip()
        col_bam = next((c for c in df_as.columns if c.lower() == 'bam'), None)
        col_acc = 'Accession number'
        if col_bam and col_acc in df_as.columns:
            lookup = df_as.set_index(col_bam)[col_acc].to_dict()
        
    results = []
    missing = []
    
    taxa_files_map = {'A': 'A.taxaBamMap.txt', 'S': 'S.taxaBamMap.txt'}
    
    for subgroup, filename in taxa_files_map.items():
        taxa_map = read_taxa_list(filename)
        for t, coverage in taxa_map.items():
            # Handle A_0001 -> A0001 difference
            t_clean = t.replace('_', '')
            acc = lookup.get(t_clean)
            if not acc:
                acc = lookup.get(t)

            if acc:
                results.append({
                    'Taxa': t,
                    'Accession': acc,
                    'Chinese_name': None,
                    'Group': subgroup,
                    'Coverage': coverage
                })
            else:
                missing.append({'Taxa': t, 'Group': subgroup, 'Reason': 'NotInASGenome'})
            
    return results, missing

def main():
    all_data = []
    all_missing = []
    
    groups_funcs = [
        process_vmap3_group,
        process_v4_group,
        process_watkins_group,
        process_ipk_group,
        process_hznu_group,
        process_as_group
    ]
    
    for func in groups_funcs:
        res, mis = func()
        all_data.extend(res)
        all_missing.extend(mis)
    
    # Save Missing
    if all_missing:
        df_missing = pd.DataFrame(all_missing)
        df_missing.to_csv("germplasm_missing_taxa.csv", index=False)
        print(f"\nMissing Taxa info saved to germplasm_missing_taxa.csv. Count: {len(df_missing)}")
    else:
        print("\nNo Missing Taxa.")
        df_missing = pd.DataFrame(columns=['Taxa', 'Group', 'Reason'])

    df = pd.DataFrame(all_data)
    print(f"\nTotal Mapped Taxa: {len(df)}")
    
    if df.empty:
        print("No mapped data found.")
        return

    # Normalization
    # Define values to be treated as empty/None
    null_markers = ['nan', 'None', '', 'NaT', '-', 'NA', 'na', 'N/A', '.', '?', 'null']
    
    df['Clean_CN'] = df['Chinese_name'].astype(str).str.strip()
    df.loc[df['Clean_CN'].isin(null_markers), 'Clean_CN'] = None
    
    df['Clean_Acc'] = df['Accession'].astype(str).str.strip()
    df.loc[df['Clean_Acc'].isin(null_markers), 'Clean_Acc'] = None
    
    # Ensure coverage is numeric
    if 'Coverage' in df.columns:
        df['Coverage'] = pd.to_numeric(df['Coverage'], errors='coerce').fillna(0)
    else:
        df['Coverage'] = 0.0

    # Grouping Logic for Duplicates
    # Use index as ID
    uf = UnionFind(df.index)
    
    # Map: Value -> List of Indices
    cn_map = {}
    acc_map = {}
    
    for idx, row in df.iterrows():
        cn = row['Clean_CN']
        acc = row['Clean_Acc']
        
        if cn:
            if cn not in cn_map: cn_map[cn] = []
            cn_map[cn].append(idx)
        if acc:
            if acc not in acc_map: acc_map[acc] = []
            acc_map[acc].append(idx)
            
    # Union based on shared CN
    for cn, indices in cn_map.items():
        first = indices[0]
        for other in indices[1:]:
            uf.union(first, other)
            
    # Union based on shared Accession
    for acc, indices in acc_map.items():
        first = indices[0]
        for other in indices[1:]:
            uf.union(first, other)
            
    # Assign Group IDs
    # Find root for each, assign a group ID
    root_to_gid = {}
    next_gid = 1
    
    dup_group_ids = []
    is_duplicate = []
    
    # Calculate group sizes first
    group_sizes = {}
    for idx in df.index:
        root = uf.find(idx)
        group_sizes[root] = group_sizes.get(root, 0) + 1
        
    for idx in df.index:
        root = uf.find(idx)
        if root not in root_to_gid:
            root_to_gid[root] = f"G_{next_gid:05d}"
            next_gid += 1
        
        gid = root_to_gid[root]
        dup_group_ids.append(gid)
        
        # Duplicate if group size > 1
        is_duplicate.append('Duplicate' if group_sizes[root] > 1 else 'Unique')
        
    df['Dup_Group_ID'] = dup_group_ids
    df['Is_Duplicate'] = is_duplicate
    
    # Stats Calculation
    print("\n" + "="*80)
    print(f"{'Group':<15} | {'Unique Sample Count':<20} | {'Missing Taxa Count':<20} | {'Total Taxa (Input)':<20}")
    print("-" * 80)
    
    groups_order = ['A', 'D', 'AB', 'ABD', 'WAP', 'Watkins', 'Nature', 'HZNU', 'S']
    
    total_unique_global = df['Dup_Group_ID'].nunique()
    
    for g in groups_order:
        df_g = df[df['Group'] == g]
        df_m = df_missing[df_missing['Group'] == g]
        
        missing_count = len(df_m)
        mapped_count = len(df_g)
        input_total = mapped_count + missing_count
        
        # Unique Samples contributed by this group?
        # Simply Number of Unique Group IDs handled by this group
        unique_count_g = df_g['Dup_Group_ID'].nunique() if not df_g.empty else 0
        
        print(f"{g:<15} | {unique_count_g:<20} | {missing_count:<20} | {input_total:<20}")
        
    print("-" * 80)
    print(f"{'Total (All)':<15} | {total_unique_global:<20} | {len(df_missing):<20} | {len(df) + len(df_missing):<20}")
    print("="*80)
    
    # Output Germplasm Dedup Summary
    output_file = "germplasm_dedup_summary.csv"
    df.to_csv(output_file, index=False)
    print(f"\nDetailed mapping with Duplicate Flags saved to {output_file}")
    
    # Output Duplicate QC Removal List
    # Identity IDs to remove (Keep max coverage per group)
    to_remove = []
    
    # Iterate over existing groups
    # Group by Dup_Group_ID
    # If size > 1, sort by coverage and drop duplicates, keeping Max Coverage
    
    grouped = df.groupby('Dup_Group_ID')
    for name, group in grouped:
        if len(group) > 1:
            # Sort descending coverage
            sorted_g = group.sort_values(by='Coverage', ascending=False)
            # Keep first (highest coverage), add rest to remove list
            # iloc[1:] gets all rows after first
            to_remove.extend(sorted_g.iloc[1:]['Taxa'].tolist())
            
    # Write germ_dup_qc.id
    qc_file = "germ_dup_qc.id"
    with open(qc_file, 'w') as f:
        f.write("#IID\n")
        for t in to_remove:
            f.write(f"{t}\n")
            
    print(f"Duplicates removal list saved to {qc_file}. Count (rejected): {len(to_remove)}")

if __name__ == "__main__":
    main()

