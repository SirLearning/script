import argparse
import pandas as pd
from infra.utils.io import load_df_generic, save_df_to_tsv


def run_tagsnp_selection(input_file, output_prefix, max_tags=1000, ld_threshold=0.8):
    df = load_df_generic(input_file)
    if df is None or df.empty:
        raise ValueError(f'No genotype matrix loaded from: {input_file}')
    if 'Sample' not in df.columns:
        raise ValueError('Input requires Sample column')

    x = df.drop(columns=['Sample']).apply(pd.to_numeric, errors='coerce').fillna(0.0)
    marker_maf = x.mean(axis=0).apply(lambda v: min(v, 1.0 - v))
    ordered_markers = marker_maf.sort_values(ascending=False).index.tolist()

    selected = []
    corr = x.corr().abs()
    for mk in ordered_markers:
        if len(selected) >= int(max_tags):
            break
        if not selected:
            selected.append(mk)
            continue

        max_ld = corr.loc[mk, selected].max()
        if pd.isna(max_ld) or max_ld < float(ld_threshold):
            selected.append(mk)

    out = pd.DataFrame({'TagSNP': selected})
    save_df_to_tsv(out, f'{output_prefix}.tagsnp.tsv')


def main():
    p = argparse.ArgumentParser(description='WWWG2B-style tagSNP selection')
    p.add_argument('--input', required=True)
    p.add_argument('--output-prefix', required=True)
    p.add_argument('--max-tags', type=int, default=1000)
    p.add_argument('--ld-threshold', type=float, default=0.8)
    args = p.parse_args()
    run_tagsnp_selection(args.input, args.output_prefix, max_tags=args.max_tags, ld_threshold=args.ld_threshold)


if __name__ == '__main__':
    main()
