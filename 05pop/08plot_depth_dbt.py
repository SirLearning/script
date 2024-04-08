import collections
from operator import itemgetter
import itertools as it
import matplotlib.pyplot as plt
import pandas as pd


def main():
    args = "data/test.txt"
    summary_name = "data/test"
    prb = (x.rstrip().split("\t") for x in open(args))
    dtb = collections.defaultdict(list)
    summary = pd.read_table(summary_name, sep="\t", header=None)
    summary.columns = ['chrom', 'length', 'bases', 'mean', 'min', 'max']
    for chrom, data in it.groupby(prb, itemgetter(0)):
        data_list = list(data)
        depth, dens = [], []
        step = float(data_list[0][1]) / 20
        tab = 0
        xb = float(data_list[0][1])
        yb = 0
        found = False
        v50 = 0
        depth.append(xb)
        dens.append(yb)
        count = 0
        for _, x, y in data_list:
            count += 1
            tab = xb - float(x)
            if tab > step:
                depth.append((xb + float(x)) / 2)
                xb = float(x)
                dens.append(float(y) - yb)
                yb = float(y)
                tab = 0
            if not found and float(y) > 0.5:
                v50 = x
                found = True
        depth.append(0)
        dens.append(0)

        exp = sum([d * p for d, p in zip(dens, depth)])
        print(exp, count, summary[summary['chrom'] == chrom]['length'].values[0])

        dtb[chrom].append({
            'depth': [round(x, 3) for x in depth],
            'dens': [round(y, 3) for y in dens],
            'mode': 'lines',
            'name': "1B" + (" (%.1f)" % float(v50))
        })

        plt.plot(depth, dens, label=chrom)
        plt.title('Probability Density Plot for ' + chrom)
        plt.xlabel('Depth')
        plt.ylabel('Probability Density')
        plt.legend([v50])
        plt.show()
    print([8.01, 31.31, 517.14, 3445.41, 3437.80, 12.20, 17.69])

if __name__ == '__main__':
    main()
