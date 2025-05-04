import sys
import pandas as pd

def count_line(file):
    germ = pd.read_csv(file, sep=',', header=0)
    germ['count'] = germ['ChineseName'].map(germ['ChineseName'].value_counts())
    return germ


def align(germ_file, ref_file):
    germ = pd.read_csv(germ_file, sep=',', header=0)
    ref = pd.read_csv(ref_file, sep=',', header=0)
    ref = ref.merge(germ, on='ChineseName').drop_duplicates()
    return ref


def main():
    germ_file = sys.argv[1]
    ref_file = sys.argv[2]
    out_file = sys.argv[3]
    # out_file = sys.argv[2]
    # germ = count_line(germ_file)
    ref = align(germ_file, ref_file)
    ref.to_csv(out_file, sep=',', index=False)

if __name__ == '__main__':
    main()