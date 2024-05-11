import sys
import pandas as pd


def fai_to_bed(fai_name):
    fai = pd.read_table(fai_name, header=None)
    fai.columns = ['chr', 'size', 'start', 'line', 'width']
    fai['start'] = 1
    fai['end'] = fai['size']
    bed = fai[['chr', 'start', 'end']]
    return bed


def main():
    bed = fai_to_bed(sys.argv[1])
    bed.to_csv(sys.argv[2], sep='\t', header=False, index=False)

if __name__ == '__main__':
    main()