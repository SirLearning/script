import os
import pandas as pd
import numpy as np
from infra.utils import load_df_from_excel, load_df_from_space_sep_no_header, load_df_from_tsv, load_df_from_tsv_no_header

# --- Sub-function to load taxaBamMap files to sample format ---
def load_df_from_tbm_no_header(input_file):
    """
    Loads a taxaBamMap file into a DataFrame.
    Assumes no header: Taxa, Coverage, [BAMs...]
    Returns: DataFrame with columns [Sample, Coverage, Bam_Path(list)]
    """
    print(f"[Info] Loading taxaBamMap: {input_file}")
    col_names = None  # Let it infer columns
    df = load_df_from_tsv_no_header(input_file, col_names)
    if df is None or df.empty:
        return None
    
    # Ensure at least 2 columns: Taxa and Coverage
    if df.shape[1] < 2:
        print(f"[Warning] taxaBamMap file has fewer than 2 columns.")
        return None
    
    # Rename first two columns. User requested 'Sample' for the first column.
    df = df.rename(columns={0: 'Sample', 1: 'Coverage'})
    
    # Coerce Coverage to numeric
    df['Coverage'] = pd.to_numeric(df['Coverage'], errors='coerce').fillna(0)
    
    # Collect columns 2+ into a single Bam_Path column (list)
    bam_cols = [c for c in df.columns if isinstance(c, int) and c >= 2]
    
    if bam_cols:
        df['Bam_Path'] = df[bam_cols].apply(
            lambda row: [str(x) for x in row if pd.notna(x) and str(x).strip() != ''], 
            axis=1
        )
        df = df.drop(columns=bam_cols)
    else:
        df['Bam_Path'] = [[] for _ in range(len(df))]
    
    return df


def load_df_from_tbm(input_file):
    """
    Wrapper to load taxaBamMap file, handling header if present.
    Returns: DataFrame with columns [Sample, Coverage, Bam_Path(list)]
    """
    print(f"[Info] Loading taxaBamMap with header check: {input_file}")
    # First try with header
    df = load_df_from_tsv(input_file)
    if df is not None and not df.empty and 'Taxa' in df.columns and 'Coverage' in df.columns:
        # Rename columns
        df = df.rename(columns={'Taxa': 'Sample'})
        df = df.rename(columns={'Coverage-Of-All-Bams': 'Coverage'})
        
        # Coerce Coverage to numeric
        df['Coverage'] = pd.to_numeric(df['Coverage'], errors='coerce').fillna(0)
        
        # Collect BAM columns
        bam_cols = [c for c in df.columns if c not in ['Sample', 'Coverage']]
        if bam_cols:
            df['Bam_Path'] = df[bam_cols].apply(
                lambda row: [str(x) for x in row if pd.notna(x) and str(x).strip() != ''], 
                axis=1
            )
            df = df.drop(columns=bam_cols)
        else:
            df['Bam_Path'] = [[] for _ in range(len(df))]
        
        return df
    else:
        # Fallback to no-header loading
        return load_df_from_tbm_no_header(input_file)    


def load_df_from_plink2(input_file):
    """
    Loads a PLINK2 generated file (e.g., .scount) into a DataFrame.
    """
    print(f"[Info] Loading PLINK2 file: {input_file}")
    df = load_df_from_tsv(input_file)
    if df is None or df.empty:
        return None
    # Rename columns
    df = df.rename(columns={'#IID': 'Sample'})
    return df


# --- Sub-function to read taxa list files ---
def read_taxa_list(directory, filename):
    """
    Reads taxaBamMap file and returns a dictionary of {taxa: coverage}.
    """
    path = os.path.join(directory, filename)
    # Using load_df_from_space_sep_no_header with col_names=None to infer integer columns
    df = load_df_from_space_sep_no_header(path, col_names=None)
    
    if df is None or df.empty:
        return {}
        
    try:
        if df.shape[1] < 2:
            print(f"[Warning] {filename} has fewer than 2 columns.")
            return {}
            
        ids = df.iloc[:, 0].astype(str)
        coverages = pd.to_numeric(df.iloc[:, 1], errors='coerce').fillna(0)
        return dict(zip(ids, coverages))
    except Exception as e:
        print(f"[Error] parsing {filename}: {e}")
        return {}

# --- Sub-functions adapted from repeat_germ_ana.py ---
def process_subgroup_vmap3(
    db_dir,
    taxa_dir
):
    print("  Processing VMap3 Group...")
    storage_path = os.path.join(db_dir, "omics/LuLab/Lulab_germplasm_Storage.xlsx")
    
    # BamDatabase
    df_bam_db = load_df_from_excel(storage_path, sheet_name='BamDatabase')
    taxa_to_bampath = {}
    if not df_bam_db.empty:
        df_bam_db.columns = df_bam_db.columns.str.strip()
        if 'BamDataBaseID' in df_bam_db.columns and 'Bam-Path' in df_bam_db.columns:
            taxa_to_bampath = df_bam_db.set_index('BamDataBaseID')['Bam-Path'].to_dict()

    # Sheets
    map_sheets = ['AB_taxa_bam', 'ABD_taxa_bam', 'D_taxa_bam']
    bampath_name_map = {}
    bampath_acc_map = {}
    
    for sheet in map_sheets:
        df_s = load_df_from_excel(storage_path, sheet_name=sheet)
        if df_s.empty: continue
        col_acc = next((c for c in df_s.columns if 'Accession' in str(c)), None)
        col_name = df_s.columns[0]
        for _, r in df_s.iterrows():
            if pd.isna(r[col_name]): continue
            name = str(r[col_name]).strip()
            acc = str(r[col_acc]).strip() if (col_acc and pd.notna(r[col_acc])) else None
            for cidx in [2, 3, 4]:
                if cidx < len(r) and pd.notna(r.iloc[cidx]):
                    p = str(r.iloc[cidx]).strip()
                    bampath_name_map[p] = name
                    if acc: bampath_acc_map[p] = acc
    
    # Vmap3 Final
    vmap3_path = os.path.join(db_dir, "germplasm/LuLab/Vmap3最终版3.0.xlsx")
    df_v3 = load_df_from_excel(vmap3_path)
    v3_lookup = {}
    v3_val_lookup = {}
    if not df_v3.empty:
        df_v3.columns = df_v3.columns.str.strip()
        v3_id = df_v3.columns[0]
        col_acc = '合并PI(1 PI;2 Taxa;34 English_name; Sample_name)'
        col_cn = 'Chinese_name'
        search_cols = ['PI_accession', 'English_name', 'Sample_name', 'Chinese_name']
        
        for _, r in df_v3.iterrows():
            if pd.isna(r[v3_id]): continue
            rid = str(r[v3_id]).strip()
            acc_val = r[col_acc] if col_acc in df_v3.columns else None
            cn_val = r[col_cn] if col_cn in df_v3.columns else None
            info = {'Accession': acc_val, 'Chinese_name': cn_val}
            v3_lookup[rid] = info
            for col in search_cols:
                if col in df_v3.columns and pd.notna(r[col]):
                    v3_val_lookup[str(r[col]).strip()] = info

    # AB Codes
    ab_new_path = os.path.join(db_dir, "germplasm/LuLab/ab_new_code.txt")
    ab_lookup = {}
    if os.path.exists(ab_new_path):
        try:
            df_ab = pd.read_csv(ab_new_path, sep=',', skiprows=1)
            df_ab.columns = df_ab.columns.str.strip()
            if '测序编号' in df_ab.columns and 'Accessions' in df_ab.columns:
                for _, r in df_ab.iterrows():
                    if pd.notna(r['测序编号']):
                            ab_lookup[str(r['测序编号']).strip()] = str(r['Accessions']).strip()
        except: pass

    results = []
    missing = []
    fn_map = {'AB': 'AB.taxaBamMap.txt', 'ABD': 'ABD.taxaBamMap.txt', 'D': 'D.taxaBamMap.txt'}
    
    for sub, fn in fn_map.items():
        t_map = read_taxa_list(taxa_dir, fn)
        for t, cov in t_map.items():
            path = taxa_to_bampath.get(t)
            common = {'Taxa': t, 'Group': sub, 'Coverage': cov}
            
            # Logic 1: NotInBamDB
            if not path:
                found = ab_lookup.get(t)
                if not found and t.endswith('A'): found = ab_lookup.get(t[:-1])
                if found:
                    common.update({'Accession': found, 'Chinese_name': None})
                    results.append(common)
                else:
                    missing.append({**common, 'Reason': 'NotInBamDatabase'})
                continue
            
            # Logic 2
            path_str = str(path).strip()
            name = bampath_name_map.get(path_str)
            if not name:
                mapped_acc = bampath_acc_map.get(path_str)
                if mapped_acc:
                    v_info = v3_val_lookup.get(mapped_acc)
                    if v_info:
                            results.append({**common, 'Accession': v_info['Accession'], 'Chinese_name': v_info['Chinese_name']})
                    else: missing.append({**common, 'Reason': 'BamSheetAccNotInVmap3'})
                else: missing.append({**common, 'Reason': 'NotInBamSheets'})
                continue
            
            # Logic 3
            clean = name[3:] if name.startswith("V3_") else name
            info = v3_lookup.get(clean)
            if info:
                results.append({**common, 'Accession': info['Accession'], 'Chinese_name': info['Chinese_name']})
            else: missing.append({**common, 'Reason': 'NotInVmap3Final'})
    
    return results, missing

def process_subgroup_v4(db_dir, taxa_dir):
    print("  Processing V4 Group...")
    t_map = read_taxa_list(taxa_dir, 'WAP.taxaBamMap.txt')
    df_v4 = load_df_from_excel(os.path.join(db_dir, "germplasm/LuLab/V4_germplasm.xlsx"))
    v4_lookup = {}
    if not df_v4.empty and 'R1' in df_v4.columns:
        df_v4['R1_str'] = df_v4['R1'].astype(str).str.strip()
        df_v4 = df_v4.drop_duplicates(subset=['R1_str'])
        v4_lookup = df_v4.set_index('R1_str').to_dict('index')
    
    results, missing = [], []
    for t, cov in t_map.items():
        t_cl = t.strip()
        row = v4_lookup.get(t_cl)
        if not row and t_cl.endswith('A'): row = v4_lookup.get(t_cl[:-1])
        common = {'Taxa': t, 'Group': 'WAP', 'Coverage': cov}
        if row:
            results.append({**common, 'Accession': row.get('Accessions'), 'Chinese_name': row.get('ChineseName')})
        else:
            missing.append({**common, 'Reason': 'NotInV4Excel'})
    return results, missing

def process_subgroup_watkins(db_dir, taxa_dir):
    print("  Processing Watkins Group...")
    files = ['w115.taxaBamMap.txt', 'w203.taxaBamMap.txt', 'w204.taxaBamMap.txt', 'w243.taxaBamMap.txt', 'w66.taxaBamMap.txt']
    t_map = {}
    for f in files: t_map.update(read_taxa_list(taxa_dir, f))
    
    omics_path = os.path.join(db_dir, "omics/Watkins_CRA012590.xlsx")
    df_run = load_df_from_excel(omics_path, sheet_name='Run')
    run_map = df_run.set_index('Accession')['Experiment accession'].to_dict() if not df_run.empty else {}
    
    df_exp = load_df_from_excel(omics_path, sheet_name='Experiment')
    exp_map = df_exp.set_index('Accession')['BioSample name'].to_dict() if not df_exp.empty else {}
    
    df_wat = load_df_from_excel(os.path.join(db_dir, "germplasm/Watkins/Watkins_sm.xlsx"), header=1)
    wat_map = {}
    if not df_wat.empty and 'Accession name' in df_wat.columns:
        for _, r in df_wat.iterrows():
            if pd.isna(r['Accession name']): continue
            aname = str(r['Accession name']).strip()
            k = aname[:-1] if aname.endswith('*') else aname
            info = {'Accession': r.get('GRU StorCode'), 'Chinese_name': r.get('Accession name in China (WWWG2B)')}
            wat_map[k] = info
            wat_map[aname] = info
            
    results, missing = [], []
    for t, cov in t_map.items():
        common = {'Taxa': t, 'Group': 'Watkins', 'Coverage': cov}
        exp = run_map.get(t)
        if not exp:
            missing.append({**common, 'Reason': 'NotInRunSheet'})
            continue
        bio = exp_map.get(exp)
        if not bio:
            missing.append({**common, 'Reason': 'NotInExpSheet'})
            continue
        info = wat_map.get(str(bio).strip())
        if info:
            results.append({**common, 'Accession': info['Accession'], 'Chinese_name': info['Chinese_name']})
        else:
                missing.append({**common, 'Reason': 'NotInWatkinsExcel'})
    return results, missing

def process_subgroup_nature(db_dir, taxa_dir): 
    print("  Processing Nature/IPK Group...")
    t_map = read_taxa_list(taxa_dir, 'Nature.taxaBamMap.txt')
    df_samp = load_df_from_excel(os.path.join(db_dir, "germplasm/ipk/ipk_sample.xlsx"))
    wgs_to_donor = df_samp.set_index('WGS Biosamples ID')['Name/GBIS Donor number'].to_dict() if not df_samp.empty else {}
    
    df_germ = load_df_from_excel(os.path.join(db_dir, "germplasm/ipk/ipk_germ.xlsx"))
    valid_donors = set(df_germ['Name/GBIS Donor number'].astype(str).str.strip()) if not df_germ.empty else set()
    
    results, missing = [], []
    for t, cov in t_map.items():
        common = {'Taxa': t, 'Group': 'Nature', 'Coverage': cov}
        donor = wgs_to_donor.get(t)
        if donor:
            if str(donor).strip() in valid_donors:
                    results.append({**common, 'Accession': donor, 'Chinese_name': None})
            else: missing.append({**common, 'Reason': 'NotInGermExcel'})
        else: missing.append({**common, 'Reason': 'NotInSampleExcel'})
    return results, missing

def process_subgroup_hznu(taxa_dir):
    print("  Processing HZNU Group...")
    t_map = read_taxa_list(taxa_dir, 'HZNU.taxaBamMap.txt')
    return [{'Taxa': t, 'Accession': t, 'Chinese_name': None, 'Group': 'HZNU', 'Coverage': cov} for t, cov in t_map.items()], []

def process_subgroup_as(db_dir, taxa_dir):
    print("  Processing A&S Group...")
    df_as = load_df_from_excel(os.path.join(db_dir, "germplasm/LuLab/A&Sgenome.xlsx"))
    lookup = {}
    if not df_as.empty:
        col_bam = next((c for c in df_as.columns if str(c).lower() == 'bam'), None)
        if col_bam: lookup = df_as.set_index(col_bam)['Accession number'].to_dict()
    
    results, missing = [], []
    for sub, fn in {'A': 'A.taxaBamMap.txt', 'S': 'S.taxaBamMap.txt', 'C': 'C.taxaBamMap.txt'}.items(): # Added C just in case, logic specific
        t_map = read_taxa_list(taxa_dir, fn)
        for t, cov in t_map.items():
            t_cl = t.replace('_', '')
            acc = lookup.get(t_cl, lookup.get(t))
            common = {'Taxa': t, 'Group': sub, 'Coverage': cov}
            if acc:
                    results.append({**common, 'Accession': acc, 'Chinese_name': None})
            else: missing.append({**common, 'Reason': 'NotInASGenome'})
    return results, missing
