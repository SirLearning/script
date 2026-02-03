import pandas as pd

def save_thresholds(stats_dict, output_file):
    """
    Saves thresholds/stats to a TSV file.
    Format:
    Header1\tHeader2...
    Value1\tValue2...
    """
    try:
        # Create DataFrame for easy TSV writing
        df = pd.DataFrame([stats_dict])
        df.to_csv(output_file, sep='\t', index=False)
        print(f"Thresholds saved to {output_file}")
    except Exception as e:
        print(f"[Error] Failed to save thresholds: {e}")
        
def load_thresholds(input_file):
    """
    Loads thresholds/stats from a TSV file into a dictionary.
    """
    try:
        df = pd.read_csv(input_file, sep='\t')
        if df.empty:
            print(f"[Warning] No data found in {input_file}")
            return {}
        # Convert first row to dictionary
        stats_dict = df.iloc[0].to_dict()
        return stats_dict
    except Exception as e:
        print(f"[Error] Failed to load thresholds: {e}")
        return {}

