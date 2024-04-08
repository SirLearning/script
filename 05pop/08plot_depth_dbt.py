import collections
from operator import itemgetter
import itertools as it

import matplotlib.pyplot as plt


def main():
    args = "data/test.txt"
    prb = (x.rstrip().split("\t") for x in open(args))
    dtb = collections.defaultdict(list)
    for chrom, data in it.groupby(prb, itemgetter(0)):
        data_list = list(data)
        depth, dens = [], []
        step = float(data_list[0][1]) / 20
        tab = 0
        xb = float(data_list[0][1])
        yb = 0
        found = False
        v50 = 0

        for _, x, y in data_list:
            tab = xb - float(x)
            print(step, tab)
            if tab > step:
                depth.append(xb)
                xb = float(x)
                dens.append(float(y) - yb)
                yb = float(y)
                tab = 0
            if not found and float(y) > 0.5:
                v50 = x
                found = True
        depth.append(0)
        dens.append(0)

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
        plt.show()

if __name__ == '__main__':
    main()
