import pandas as pd
import os

def get_country_to_continent_map():
    # A comprehensive dictionary mapping countries to continents
    # This covers many common country names found in germplasm databases
    return {
        # Africa
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
        
        # Asia - Subregions
        # Central Asia
        'Kazakhstan': 'Central Asia', 'Kyrgyzstan': 'Central Asia', 'Tajikistan': 'Central Asia', 
        'Turkmenistan': 'Central Asia', 'Uzbekistan': 'Central Asia',
        
        # East Asia
        'China': 'East Asia', 'Japan': 'East Asia', 'Mongolia': 'East Asia', 'North Korea': 'East Asia', 
        'South Korea': 'East Asia', 'Taiwan': 'East Asia',
        
        # South Asia
        'Afghanistan': 'South Asia', 'Bangladesh': 'South Asia', 'Bhutan': 'South Asia', 
        'India': 'South Asia', 'Maldives': 'South Asia', 'Nepal': 'South Asia', 
        'Pakistan': 'South Asia', 'Sri Lanka': 'South Asia',
        
        # West Asia
        'Armenia': 'West Asia', 'Azerbaijan': 'West Asia', 'Bahrain': 'West Asia', 'Cyprus': 'West Asia', 
        'Georgia': 'West Asia', 'Iran': 'West Asia', 'Iraq': 'West Asia', 'Israel': 'West Asia', 
        'Jordan': 'West Asia', 'Kuwait': 'West Asia', 'Lebanon': 'West Asia', 'Oman': 'West Asia', 
        'Palestine': 'West Asia', 'Qatar': 'West Asia', 'Saudi Arabia': 'West Asia', 
        'Syria': 'West Asia', 'Turkey': 'West Asia', 'United Arab Emirates': 'West Asia', 
        'Yemen': 'West Asia',
        
        # Southeast Asia -> Mapping to East/South Asia as per request to split only into 4 regions
        # Strategy: Myanmar, Thailand -> South Asia; Vietnam, Philippines, etc -> East Asia
        'Brunei': 'East Asia', 'Cambodia': 'South Asia', 'Indonesia': 'East Asia', 'Laos': 'South Asia', 
        'Malaysia': 'East Asia', 'Myanmar': 'South Asia', 'Burma': 'South Asia', # Burma is Myanmar
        'Philippines': 'East Asia', 'Singapore': 'East Asia', 
        'Thailand': 'South Asia', 'Timor-Leste': 'East Asia', 'Vietnam': 'East Asia',
        
        'Russia': 'Europe', 
        
        # Europe
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
        'UK': 'Europe', 'Vatican City': 'Europe', 'Soviet Union': 'Europe', 'CS': 'Europe', # Czechoslovakia
        'Yugoslavia': 'Europe', 'Crete': 'Europe', 'Canary Islands': 'Europe', 'Central Europe': 'Europe',
        'Netherland': 'Europe',
        
        # North America
        
        # North America
        'Antigua and Barbuda': 'North America', 'Bahamas': 'North America', 'Barbados': 'North America', 
        'Belize': 'North America', 'Canada': 'North America', 'Costa Rica': 'North America', 
        'Cuba': 'North America', 'Dominica': 'North America', 'Dominican Republic': 'North America', 
        'El Salvador': 'North America', 'Grenada': 'North America', 'Guatemala': 'North America', 
        'Haiti': 'North America', 'Honduras': 'North America', 'Jamaica': 'North America', 
        'Mexico': 'North America', 'Nicaragua': 'North America', 'Panama': 'North America', 
        'Saint Kitts and Nevis': 'North America', 'Saint Lucia': 'North America', 
        'Saint Vincent and the Grenadines': 'North America', 'Trinidad and Tobago': 'North America', 
        'United States': 'North America', 'USA': 'North America', 'US': 'North America',
        
        # South America
        'Argentina': 'South America', 'Bolivia': 'South America', 'Brazil': 'South America', 
        'Chile': 'South America', 'Colombia': 'South America', 'Ecuador': 'South America', 
        'Guyana': 'South America', 'Paraguay': 'South America', 'Peru': 'South America', 
        'Suriname': 'South America', 'Uruguay': 'South America', 'Venezuela': 'South America',
        
        # Oceania
        'Australia': 'Oceania', 'Fiji': 'Oceania', 'Kiribati': 'Oceania', 'Marshall Islands': 'Oceania', 
        'Micronesia': 'Oceania', 'Nauru': 'Oceania', 'New Zealand': 'Oceania', 'Palau': 'Oceania', 
        'Papua New Guinea': 'Oceania', 'Samoa': 'Oceania', 'Solomon Islands': 'Oceania', 
        'Tonga': 'Oceania', 'Tuvalu': 'Oceania', 'Vanuatu': 'Oceania',
        
        # 3-Letter Codes (ISO 3166-1 alpha-3 mostly) for IPK/Nature
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
    
    # 尝试直接匹配
    if country in mapping:
        return mapping[country]
    
    # 尝试一些常见的别名或清洗
    if country == 'USA' or country == 'United States of America': return 'North America'
    if country == 'UK': return 'Europe'
    # if 'Russia' in country: return 'Europe' # This was overriding dict.
    # The dictionary now has 'Russia': 'Asia' (which caused the issue) or we can set it to Europe in dict.
    # Re-enforcing Russia -> Europe here to fix the 10 samples
    if 'Russia' in country: return 'Europe' 
    if 'USSR' in country or 'Soviet' in country: return 'Europe'
    if country == 'Czechia': return 'Europe'
    if 'Macedonia' in country: return 'Europe'
    if 'China' in country: return 'East Asia'
    if 'Iran' in country: return 'West Asia'
    
    # 不区分大小写尝试
    for k, v in mapping.items():
        if k.lower() == country.lower():
            return v
            
    return 'Unknown'

def ck_location():
    # 1. 定义文件路径
    home_dir = os.path.expanduser("/data/home/tusr1/git/DBone/Service/src/main/resources/raw/20251208/germplasm")
    # home_dir = "/data/home/tusr1" # Explicit based on environment
    v4_path = os.path.join(home_dir, "LuLab/V4_germplasm.xlsx")
    vmap3_path = os.path.join(home_dir, "LuLab/Vmap3最终版3.0.xlsx")
    watkins_path = os.path.join(home_dir, "Watkins/Watkins_sm.xlsx")
    
    print(f"Reading V4 file: {v4_path}")
    print(f"Reading Vmap3 file: {vmap3_path}")
    print(f"Reading Watkins file: {watkins_path}")
    
    # 2. 读取 Excel 文件
    try:
        df_v4 = pd.read_excel(v4_path)
        df_vmap3 = pd.read_excel(vmap3_path)
        # Watkins 表头在第二行 (index 1)
        df_watkins = pd.read_excel(watkins_path, header=1)
        # 去除列名可能会有的前后空格
        if not df_watkins.empty:
            df_watkins.columns = df_watkins.columns.str.strip()
    except Exception as e:
        print(f"Error reading Excel files: {e}")
        # Try to install openpyxl if missing? usually installed.
        return

    # 3. 处理 V4 数据
    # 列名: 'Provenance of material-GS'
    if 'Provenance of material-GS' not in df_v4.columns:
        print("Error: 'Provenance of material-GS' column not found in V4_germplasm.xlsx")
        print("Available columns:", df_v4.columns.tolist())
        return

    # === 新增筛选步骤 ===
    print("Deduplicating V4 data...")
    print(f"Initial V4 count: {len(df_v4)}")
    
    # 查找列名: ChineseName (可能叫 Chinese_name, Chinese Name 等)
    col_chinese = None
    for c in df_v4.columns:
        if 'chinese' in c.lower() and 'name' in c.lower():
            col_chinese = c
            break
    
    # 查找列名: Accessions (可能叫 Accession, PI_accession 等)
    col_acc = None
    for c in df_v4.columns:
        if 'accession' in c.lower() or 'pi' in c.lower() and 'name' not in c.lower(): # 'Accessions' or similar
             col_acc = c
             break
    # 如果没找到模糊匹配，尝试精确匹配用户提供的
    if not col_chinese and 'ChineseName' in df_v4.columns: col_chinese = 'ChineseName'
    if not col_acc and 'Accessions' in df_v4.columns: col_acc = 'Accessions'

    # 执行去重
    if col_chinese:
        print(f"Applying deduplication on column: '{col_chinese}'")
        # 排除 NaN 值参与去重（即两个 NaN 不算重复）
        # duplicated() 返回 True 表示是重复的（保留第一个）。
        # 我们只删除那些 duplicated 为 True 且该列值不为空的行
        is_dup_name = df_v4.duplicated(subset=[col_chinese], keep='first') & df_v4[col_chinese].notna()
        df_v4 = df_v4[~is_dup_name]
        print(f"  Count after ChineseName deduplication: {len(df_v4)}")
    else:
        print("Warning: Could not find 'ChineseName' column for deduplication. Available:", df_v4.columns.tolist())

    if col_acc:
        print(f"Applying deduplication on column: '{col_acc}'")
        is_dup_acc = df_v4.duplicated(subset=[col_acc], keep='first') & df_v4[col_acc].notna()
        df_v4 = df_v4[~is_dup_acc]
        print(f"  Count after Accessions deduplication: {len(df_v4)}")
    else:
        print("Warning: Could not find 'Accessions' column for deduplication.")

    print(f"Final V4 count: {len(df_v4)}")
    # =====================

    # === 对比 V4 和 Vmap3 的重叠 ===
    print("\n" + "="*30)
    print("Checking Overlap: V4 vs Vmap3")
    print("="*30)
    
    # 1. 寻找 Vmap3 的列名
    vmap3_col_chinese = None
    vmap3_col_acc = None
    
    for c in df_vmap3.columns:
        c_str = str(c).strip()
        c_lower = c_str.lower()
        # Case insensitive matching
        if 'chinese' in c_lower and 'name' in c_lower:
            vmap3_col_chinese = c
        if 'accession' in c_lower:
            vmap3_col_acc = c
            
    # 2. 对比 Chinese Name
    if col_chinese and vmap3_col_chinese:
        print(f"Comparing columns: V4['{col_chinese}'] vs Vmap3['{vmap3_col_chinese}']")
        v4_names = set(df_v4[col_chinese].dropna().astype(str).str.strip())
        vmap3_names = set(df_vmap3[vmap3_col_chinese].dropna().astype(str).str.strip())
        overlap_names = v4_names.intersection(vmap3_names)
        print(f"  -> Overlap count (Chinese Name): {len(overlap_names)}")
    else:
        print(f"Skipping Chinese Name comparison (Cols missing). V4: {col_chinese}, Vmap3: {vmap3_col_chinese}")
        if not vmap3_col_chinese: print("  Vmap3 Columns available:", df_vmap3.columns.tolist())

    # 3. 对比 Accessions
    if col_acc and vmap3_col_acc:
        print(f"Comparing columns: V4['{col_acc}'] vs Vmap3['{vmap3_col_acc}']")
        v4_accs = set(df_v4[col_acc].dropna().astype(str).str.strip())
        vmap3_accs = set(df_vmap3[vmap3_col_acc].dropna().astype(str).str.strip())
        overlap_accs = v4_accs.intersection(vmap3_accs)
        print(f"  -> Overlap count (Accessions): {len(overlap_accs)}")
    else:
        print(f"Skipping Accession comparison (Cols missing). V4: {col_acc}, Vmap3: {vmap3_col_acc}")
        if not vmap3_col_acc and vmap3_col_chinese: # avoid printing twice if both missing
             print("  Vmap3 Columns available:", df_vmap3.columns.tolist())
    
    print("="*30 + "\n")
    # ===============================

    print("Mapping V4 countries to continents...")
    country_map = get_country_to_continent_map()
    
    # 获取 V4 对应的 State/Continent
    v4_countries = df_v4['Provenance of material-GS']
    v4_continents = v4_countries.apply(lambda x: map_country_to_continent(x, country_map))
    
    # 4. 处理 Vmap3 数据
    # 列名: 'continent' 或 'Continent(已补充)'
    target_col = 'continent'
    if target_col not in df_vmap3.columns:
        possibilities = ['Continent(已补充)', 'Ctnt(9)', 'Continent']
        found = False
        for p in possibilities:
            if p in df_vmap3.columns:
                target_col = p
                found = True
                print(f"Note: Using '{target_col}' column for Vmap3.")
                break
        
        if not found:
            print("Error: Could not find a continent column in Vmap3 file.")
            print("Available columns:", df_vmap3.columns.tolist())
            return

    vmap3_raw_continents = df_vmap3[target_col]
    
    # --- 标准化 Vmap3 大洲名称 ---
    vmap3_continents = vmap3_raw_continents.astype(str).str.strip().str.title()
    correction_map = {
        'Centralasia': 'Central Asia',
        'Eastasia': 'East Asia',
        'Northamerica': 'North America',
        'Southamerica': 'South America',
        'Southasia': 'South Asia',
        'Westasia': 'West Asia',
        '-': 'Unknown'
    }
    vmap3_continents = vmap3_continents.replace(correction_map)
    
    # 5. 处理 Watkins 数据
    print("Mapping Watkins countries to continents...")
    if 'Country of origin' not in df_watkins.columns:
        print("Error: 'Country of origin' column not found in Watkins_sm.xlsx")
        return

    watkins_countries = df_watkins['Country of origin']
    watkins_continents = watkins_countries.apply(lambda x: map_country_to_continent(x, country_map))

    # ===============================
    # 5b. Processing IPK/Nature Data
    # ===============================
    print("\n" + "="*30)
    print("Processing Nature/IPK Data")
    print("="*30)
    
    nature_path = "/data/home/tusr1/01projects/vmap4/00data/05taxaBamMap/Nature.taxaBamMap.txt"
    ipk_sample_path = os.path.join(home_dir, "ipk/ipk_sample.xlsx")
    ipk_germ_path = os.path.join(home_dir, "ipk/ipk_germ.xlsx")
    ipk_ena_path = os.path.join(home_dir, "ipk/ipk_ena.xlsx")
    
    nature_continents_list = []
    
    # Check files exist
    if os.path.exists(nature_path) and os.path.exists(ipk_sample_path) and os.path.exists(ipk_germ_path):
        try:
            # 1. Read Nature list (Tab separated, no header)
            # Column 0 is the ID (e.g., SAMEA...)
            df_nature = pd.read_csv(nature_path, sep=r'\s+', header=None, usecols=[0])
            nature_ids = df_nature[0].astype(str).str.strip().unique()
            print(f"Loaded {len(nature_ids)} samples from Nature.taxaBamMap.txt")
            
            # --- Debug Lists ---
            debug_asia_origins = set()
            debug_americas_origins = set()
            missing_in_ipk_sample = [] # List for finding missing samples
            # -------------------

            # 2. Read IPK Sample Map
            df_ipk_sample = pd.read_excel(ipk_sample_path)
            # Clean columns
            df_ipk_sample.columns = df_ipk_sample.columns.str.strip()
            
            # Find key columns
            col_wgs_id = next((c for c in df_ipk_sample.columns if 'WGS Biosamples ID' in c), None)
            col_donor_num_sample = next((c for c in df_ipk_sample.columns if 'Name/GBIS Donor number' in c), None)
            
            if not col_wgs_id or not col_donor_num_sample:
                print(f"Error: Columns missing in ipk_sample.xlsx. Found: {df_ipk_sample.columns.tolist()}")
            else:
                # Create Dictionary: ID -> Donor Number
                print("Checking ipk_sample.xlsx for duplicate mappings...")
                wgs_counts = df_ipk_sample[col_wgs_id].value_counts()
                dups = wgs_counts[wgs_counts > 1]
                if len(dups) > 0:
                     print(f"Warning: {len(dups)} WGS IDs appear multiple times in ipk_sample.xlsx.")
                     
                sample_map = df_ipk_sample.set_index(col_wgs_id)[col_donor_num_sample].to_dict()
                
                # 3. Read IPK Germ Map
                df_ipk_germ = pd.read_excel(ipk_germ_path)
                df_ipk_germ.columns = df_ipk_germ.columns.str.strip()
                
                col_donor_num_germ = next((c for c in df_ipk_germ.columns if 'Name/GBIS Donor number' in c), None)
                # Look for columns flexibly
                col_continent = next((c for c in df_ipk_germ.columns if 'Continent' in c), None)
                col_macro = next((c for c in df_ipk_germ.columns if 'Macroregin' in c or 'Macroregion' in c), None)
                col_origin = next((c for c in df_ipk_germ.columns if 'Origin' in c), None)
                
                # 3.5 Read IPK ENA (Fallback)
                ena_map = {}
                if os.path.exists(ipk_ena_path):
                    print(f"Reading IPK ENA file: {ipk_ena_path}")
                    try:
                        df_ena = pd.read_excel(ipk_ena_path)
                        df_ena.columns = df_ena.columns.str.strip()
                        col_acc_ena = next((c for c in df_ena.columns if 'sample_accession' in c.lower()), None)
                        col_alias_ena = next((c for c in df_ena.columns if 'sample_alias' in c.lower()), None)
                        if col_acc_ena and col_alias_ena:
                            ena_map = df_ena.set_index(col_acc_ena)[col_alias_ena].to_dict()
                            print(f"Loaded {len(ena_map)} ENA mappings.")
                        else:
                            print(f"Warning: Columns missing in ipk_ena.xlsx (Need sample_accession, sample_alias). Found: {df_ena.columns.tolist()}")
                    except Exception as e:
                        print(f"Error reading ipk_ena.xlsx: {e}")
                
                 # 3.6 Read IPK ENA 2 (Fallback 2)
                ipk_ena2_path = os.path.join(home_dir, "ipk_ena_2.xlsx")
                if os.path.exists(ipk_ena2_path):
                    print(f"Reading IPK ENA 2 file: {ipk_ena2_path}")
                    try:
                        df_ena2 = pd.read_excel(ipk_ena2_path)
                        df_ena2.columns = df_ena2.columns.str.strip()
                        # Reuse same column names/patterns if possible, or adapt if file structure differs
                        col_acc_ena2 = next((c for c in df_ena2.columns if 'sample_accession' in c.lower()), None)
                        col_alias_ena2 = next((c for c in df_ena2.columns if 'sample_alias' in c.lower()), None)
                        
                        if col_acc_ena2 and col_alias_ena2:
                            # Update dictionary. If duplicates exist between file 1 and 2, file 2 overwrites (or use separate details)
                            # We can merge into ena_map directly
                            count_new = 0
                            for idx, row in df_ena2.iterrows():
                                acc = row[col_acc_ena2]
                                alias = row[col_alias_ena2]
                                if acc not in ena_map:
                                    ena_map[acc] = alias
                                    count_new += 1
                            print(f"Added {count_new} new mappings from ENA 2.")
                        else:
                             print(f"Warning: Columns missing in ipk_ena_2.xlsx. Found: {df_ena2.columns.tolist()}")
                    except Exception as e:
                         print(f"Error reading ipk_ena_2.xlsx: {e}")

                if not col_donor_num_germ:
                    print(f"Error: Columns missing in ipk_germ.xlsx. Found: {df_ipk_germ.columns.tolist()}")
                else:
                    # Create lookup dataframe indexed by Donor Number
                    df_germ_lookup = df_ipk_germ.set_index(col_donor_num_germ)
                    
                    # --- DEDUPLICATION & MERGE LOGIC ---
                    seen_donors = {} 
                    merged_samples_log = []
                    unique_nature_ids_to_process = []
                    
                    # Pre-pass to handle deduplication logic
                    # We need to know the donor for every ID first (direct or via ENA)
                    
                    # Helper to get donor num from ID
                    def get_donor_num(nid):
                        # 1. Try ipk_sample
                        d_num = sample_map.get(nid)
                        if d_num: return d_num, 'Direct'
                        
                        # 2. Try ipk_ena -> parse alias
                        alias = ena_map.get(nid)
                        if alias and isinstance(alias, str):
                            parts = alias.strip().split('_')
                            if len(parts) > 0:
                                if parts[0] == 'TRI':
                                    if len(parts) > 1:
                                        return f"TRI {parts[1]}", 'ENA'
                                elif parts[0] == 'B': # Special case for B_1, B_100, etc.
                                    if len(parts) > 1:
                                        return f"B {parts[1]}", 'ENA'
                                else:
                                    # Fallback for like PF096_leer -> PF096
                                    return parts[0], 'ENA'
                        return None, None
                    
                    # Store ENA Matches for reporting
                    ena_matches_log = []

                    for nid in nature_ids:
                        donor_num, source = get_donor_num(nid)
                        
                        if donor_num:
                            if source == 'ENA':
                                ena_matches_log.append(f"Sample: {nid} -> Found via ENA -> Donor: {donor_num}")

                            if donor_num in seen_donors:
                                existing_nid = seen_donors[donor_num]
                                merged_samples_log.append(f"Donor: {donor_num} | Skipping {nid} (Exists as {existing_nid})")
                                continue # Skip processing this sample
                            seen_donors[donor_num] = nid
                        else:
                            # No donor found, but we still process it to count as Unknown or log it
                            pass
                        
                        unique_nature_ids_to_process.append(nid)
                        
                    if merged_samples_log:
                        print(f"\n[INFO] Merged {len(merged_samples_log)} duplicate samples sharing the same Donor:")
                        # for log in merged_samples_log[:5]:
                        #    print(f"  - {log}")
                    
                    # Report ENA matches
                    if ena_matches_log:
                        print(f"\n[INFO] Found {len(ena_matches_log)} samples via ENA lookup:")
                        for msg in ena_matches_log:
                            print(f"  {msg}")

                    # 4. Process each Nature ID
                    mapped_count = 0
                    for nid in unique_nature_ids_to_process:
                        final_cont = 'Unknown'
                        
                        donor_num, source = get_donor_num(nid)
                        
                        if not donor_num:
                            missing_in_ipk_sample.append(nid)
                        
                        if donor_num:
                            # Step B: Donor Number -> Geography
                            if donor_num in df_germ_lookup.index:
                                row = df_germ_lookup.loc[donor_num]
                                # Handle duplicates in index by taking first
                                if isinstance(row, pd.DataFrame):
                                    row = row.iloc[0]
                                
                                # Priority 1: Continent column
                                if col_continent and pd.notna(row[col_continent]):
                                    c_val = str(row[col_continent]).strip()
                                    if c_val.lower() != 'nan':
                                        final_cont = c_val
                                
                                # Priority 2: Macroregion -> Clean manually if needed
                                if final_cont == 'Unknown' and col_macro and pd.notna(row[col_macro]):
                                     m_val = str(row[col_macro]).strip()
                                     if 'Asia' in m_val: final_cont = 'Asia'
                                     elif 'Europe' in m_val: final_cont = 'Europe'
                                     elif 'America' in m_val: final_cont = 'Americas'
                                     elif 'Africa' in m_val: final_cont = 'Africa'
                                
                                # Priority 3: Origin -> use map_country_to_continent
                                if final_cont == 'Unknown' and col_origin and pd.notna(row[col_origin]):
                                    origin_val = str(row[col_origin]).strip()
                                    final_cont = map_country_to_continent(origin_val, country_map)
                                
                                # === Refinement Step for Generic Continents ===
                                if final_cont in ['Asia', 'Americas', 'America', 'Unknown']:
                                    potential_refined = 'Unknown'
                                    
                                    # Try Origin first
                                    if col_origin and pd.notna(row[col_origin]):
                                        origin_val = str(row[col_origin]).strip()
                                        mapped_origin = map_country_to_continent(origin_val, country_map)
                                        if mapped_origin != 'Unknown':
                                            potential_refined = mapped_origin
                                    
                                    # If Origin didn't give specific (still 'Unknown' or same generic), try detailed Macroregion
                                    if potential_refined == 'Unknown' or potential_refined == final_cont:
                                        if col_macro and pd.notna(row[col_macro]):
                                            m_val = str(row[col_macro]).strip()
                                            if 'East Asia' in m_val: potential_refined = 'East Asia'
                                            elif 'West Asia' in m_val: potential_refined = 'West Asia'
                                            elif 'South Asia' in m_val: potential_refined = 'South Asia'
                                            elif 'Central Asia' in m_val: potential_refined = 'Central Asia'
                                            elif 'North America' in m_val: potential_refined = 'North America'
                                            elif 'South America' in m_val: potential_refined = 'South America'
                                            elif 'Europe' in m_val: potential_refined = 'Europe' # Standardize if mixed

                                    # Update if better
                                    if potential_refined != 'Unknown':
                                        # Only update if it makes it more specific, or if original was Unknown
                                        if final_cont == 'Unknown':
                                            final_cont = potential_refined
                                        elif final_cont == 'Asia' and 'Asia' in potential_refined:
                                            final_cont = potential_refined
                                        elif final_cont in ['Americas', 'America'] and 'America' in potential_refined:
                                            final_cont = potential_refined
                                
                                # Debugging Collection
                                if final_cont == 'Asia' and col_origin and pd.notna(row[col_origin]):
                                     debug_asia_origins.add(str(row[col_origin]).strip())
                                if final_cont in ['Americas', 'America'] and col_origin and pd.notna(row[col_origin]):
                                     debug_americas_origins.add(str(row[col_origin]).strip())
                        
                        nature_continents_list.append(final_cont)
                        if final_cont != 'Unknown':
                            mapped_count += 1
                            
                    print(f"Mapped {mapped_count} out of {len(nature_ids)} Nature samples to a continent.")                    
                    
                    if missing_in_ipk_sample:
                        print(f"\n[INFO] Samples NOT found (even with ENA lookup): {len(missing_in_ipk_sample)}")
                        print(f"List of Missing Samples: {missing_in_ipk_sample}")
                        
                    if debug_asia_origins:
                        print("\n[DEBUG] Origins remaining as generic 'Asia' (Could not map to subregion):")
                        print(debug_asia_origins)
                    if debug_americas_origins:
                        print("\n[DEBUG] Origins remaining as generic 'Americas' (Could not map to subregion):")
                        print(debug_americas_origins)

        except Exception as e:
            print(f"Error processing IPK files: {e}")
            import traceback
            traceback.print_exc()
    else:
        print("Warning: One or more IPK/Nature files not found. Skipping this step.")

    nature_continents = pd.Series(nature_continents_list)
    # Normalize names
    nature_continents = nature_continents.astype(str).str.strip().str.title()
    nature_continents = nature_continents.replace(correction_map)

    # ===============================
    # 6. Overlap Analysis & Deduplication (Merger)
    # Target: Base on Vmap3, add unique V4, add Watkins
    # ===============================
    print("\n" + "="*30)
    print("Analyzing Overlap (V4 vs Vmap3)")
    print("="*30)

    # A. Identify Linking Columns in Vmap3
    vmap3_col_chinese = None
    vmap3_col_acc = None
    
    for c in df_vmap3.columns:
        c_str = str(c).strip()
        c_lower = c_str.lower()
        if 'chinese' in c_lower and 'name' in c_lower:
            vmap3_col_chinese = c
        if 'accession' in c_lower:
            vmap3_col_acc = c

    # B. Build Lookup Dictionaries for Vmap3
    # We map Key -> Index to update Vmap3 if needed
    vmap3_cn_idx_map = {}
    vmap3_acc_idx_map = {}

    # Need align indices
    # df_vmap3 and vmap3_continents share index
    # Iterate and populate
    if vmap3_col_chinese:
        print(f"Using Vmap3 column '{vmap3_col_chinese}' for Name matching.")
        # Zip inputs: names, indices
        for idx, name in df_vmap3[vmap3_col_chinese].items():
            if pd.isna(name): continue
            name_str = str(name).strip()
            vmap3_cn_idx_map[name_str] = idx
    
    if vmap3_col_acc:
        print(f"Using Vmap3 column '{vmap3_col_acc}' for Accession matching.")
        for idx, acc in df_vmap3[vmap3_col_acc].items():
            if pd.isna(acc): continue
            acc_str = str(acc).strip()
            vmap3_acc_idx_map[acc_str] = idx

    # C. Check V4 against Vmap3
    v4_unique_continents = []
    
    overlap_count = 0
    consistent_count = 0
    conflict_count = 0
    
    # Iterate V4
    # df_v4 and v4_continents share index
    # Use reset_index to iterate safely or zip
    
    print("Checking V4 samples against Vmap3 database...")
    
    # Helper to check if item is in Vmap3
    # Returns (is_overlap, vmap3_idx)
    def check_overlap_idx(row):
        # Check Name
        if col_chinese and pd.notna(row[col_chinese]):
            name = str(row[col_chinese]).strip()
            if name in vmap3_cn_idx_map:
                return True, vmap3_cn_idx_map[name]
                
        # Check Accession
        if col_acc and pd.notna(row[col_acc]):
            acc = str(row[col_acc]).strip()
            if acc in vmap3_acc_idx_map:
                return True, vmap3_acc_idx_map[acc]
        
        return False, None

    for idx, row in df_v4.iterrows():
        v4_cont_val = v4_continents.loc[idx]
        is_overlap, vmap3_idx = check_overlap_idx(row)
        
        if is_overlap:
            overlap_count += 1
            # Check consistency
            vmap3_cont_val = vmap3_continents.loc[vmap3_idx]
            
            # simple string comparison
            if str(v4_cont_val).lower() == str(vmap3_cont_val).lower():
                consistent_count += 1
            else:
                # Resolve Conflict
                vmap3_curr = str(vmap3_cont_val).lower()
                v4_curr = str(v4_cont_val).lower()
                invalids = ['unknown', '-', 'nan']
                
                # Case 1: Vmap3 is Unknown, V4 is known -> Update Vmap3
                # (Assuming 'Unknown' or '-' are invalid)
                if vmap3_curr in invalids and v4_curr not in invalids:
                    print(f"  [Update] Vmap3 {vmap3_idx} updated from {vmap3_cont_val} to {v4_cont_val} based on V4 {row.get(col_chinese)}")
                    vmap3_continents.at[vmap3_idx] = v4_cont_val
                    consistent_count += 1 
                # Case 2: V4 is Unknown, Vmap3 is known -> Trust Vmap3 (Do nothing to data, but count as consistent)
                elif v4_curr in invalids and vmap3_curr not in invalids:
                     # Just log locally if needed, but treat as resolved.
                     # Since we use Vmap3 list for counts, the final result will be the Known value from Vmap3.
                     # print(f"  [Resolved] V4 is Unknown, keeping Vmap3 {vmap3_cont_val}")
                     consistent_count += 1
                else:
                    conflict_count += 1
                    if conflict_count <= 5: # Print first 5 conflicts
                       print(f"  [Conflict] V4 Name/Acc: {row.get(col_chinese)}/{row.get(col_acc)} | V4 Cont: {v4_cont_val} != Vmap3 Cont: {vmap3_cont_val}")
        else:
            # Unique to V4, keep it
            v4_unique_continents.append(v4_cont_val)
            
    print("-" * 20)
    print(f"Total V4 Samples: {len(df_v4)}")
    print(f"Overlapping with Vmap3: {overlap_count}")
    print(f"  - Consistent Continent: {consistent_count}")
    print(f"  - Conflicting Continent: {conflict_count}")
    print(f"Unique V4 Samples (Added to total): {len(v4_unique_continents)}")
    print("="*30 + "\n")

    # 7. 合并数据 & 统计
    print("Combining data (Vmap3 + Unique V4 + Watkins + Nature)...")
    
    final_list = list(vmap3_continents)
    final_list.extend(v4_unique_continents)
    final_list.extend(watkins_continents)
    final_list.extend(nature_continents)
    
    all_continents = pd.Series(final_list)
    # Ensure Title Case
    all_continents = all_continents.astype(str).str.title()
    
    # 统计
    counts = all_continents.value_counts()
    
    # 8. 输出结果
    print("\n" + "="*40)
    print("Final Sample Counts (Deduplicated + Nature)")
    print("="*40)
    
    # Custom sort order
    def sort_key(name):
        name = name.lower()
        if 'asia' in name:
            return f"0_{name}" 
        elif 'america' in name:
            return f"1_{name}" 
        elif 'europe' in name:
            return f"2_{name}"
        elif 'africa' in name:
            return f"3_{name}"
        elif 'oceania' in name:
            return f"4_{name}"
        else:
            return f"5_{name}" 

    sorted_cats = sorted(list(counts.index), key=sort_key)
    
    print(f"{'Continent':<20} | {'Count':<10}")
    print("-" * 33)
    total = 0
    for cont in sorted_cats:
        if cont.lower() == 'nan': continue
        c = counts[cont]
        total += c
        print(f"{cont:<20} | {c:<10}")
    print("-" * 33)
    print(f"{'Total':<20} | {total:<10}")
    
    # Breakdown is tricky now because we merged datasets.
    # We can reconstruct roughly:
    # Vmap3 Total, V4 Unique Total, Watkins Total.
    print("\n" + "="*40)
    print("Dataset Contribution")
    print("="*40)
    print(f"Vmap3 (Total): {len(vmap3_continents)}")
    print(f"V4 (Unique added): {len(v4_unique_continents)} (Overlap: {overlap_count})")
    print(f"Watkins (Total): {len(watkins_continents)}")
    print(f"Nature (Total): {len(nature_continents)}")
    print(f"Sum: {len(vmap3_continents) + len(v4_unique_continents) + len(watkins_continents) + len(nature_continents)}")


