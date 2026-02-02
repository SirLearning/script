from turtle import pd


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
