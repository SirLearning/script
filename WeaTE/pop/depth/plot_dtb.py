import collections
from operator import itemgetter
import itertools as it
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.signal import find_peaks_cwt
from scipy.signal import find_peaks

def main():
    summary = pd.concat([chr_main('1A'), chr_main('1B'), chr_main('1D')])
    summary.reset_index(drop=True, inplace=True)
    summary = summary.sort_values(by='type')

    fig, ax = plt.subplots()
    # plt.rcParams['font.size'] = 24
    ax.figure.set_size_inches(18, 9)

    sns.boxplot(x='type', y='cv', hue='chr', data=summary, showfliers=False,
                palette='Set3', linewidth=2, dodge=True, ax=ax, width=0.8)
    plt.axhline(1, color='r', linestyle='--')
    plt.title('Variation of Depth Distribution')
    plt.ylabel('Coefficient of Variation')
    plt.legend(title='Name')
    plt.show()

    with open('data/vu_reads_depth/threshold.txt', 'w') as f:
        for i in range(len(summary)):
            if summary.loc[i, 'threshold']:
                f.write(summary.loc[i, 'chrom'] + '\t' + summary.loc[i, 'chr'] + '\n')


def peak(data):
    peak_idx = []
    flag = False
    bias0 = 0
    for i in range(len(data) - 1):
        bias1 = data[i] - data[i + 1]
        if i < len(data) - 2:
            bias2 = data[i + 1] - data[i + 2]
        else:
            bias2 = 0
        if flag:
            flag = False
            continue
        if bias1 > min(0.1, max(data)/3) and bias0 <= 0:
            peak_idx.append(i)
        elif bias1 < -min(0.1, max(data)/3) and bias2 >= 0:
            flag = True
            peak_idx.append(i + 1)
        bias0 = bias1
    return peak_idx


def chr_main(chr):
    summary_name = "data/vu_reads_depth/" + chr + ".mosdepth.summary.txt"
    summary = pd.read_table(summary_name, sep="\t", header=None)
    summary.columns = ['chrom', 'length', 'bases', 'mean', 'min', 'max']

    args = "data/depth_dtb/" + chr + ".mosdepth.global.dist.txt"
    prb = (x.rstrip().split("\t") for x in open(args))
    summary = dtb_form(prb, summary)

    # TE type
    summary['type'] = summary['chrom'].str.split('#').str[1]
    tecode = pd.read_table("data/TEcode", sep=",", header=None)
    tecode.columns = ['cls', 'new_cls']
    for i in range(0, len(tecode)):
        summary.loc[summary['type'] == tecode['cls'][i], 'type'] = tecode['new_cls'][i]
    summary['ref_depth'] = summary['mean']
    summary['chr'] = chr
    return summary


def dtb_form(prb, summary):
    times = 0
    summary['threshold'] = False
    for chrom, data in it.groupby(prb, itemgetter(0)):
        # if "LTR" not in chrom: continue
        times += 1
        data_list = list(data)
        depth, dens = [], []
        sep = 60
        step = float(data_list[0][1]) / min(sep, len(data_list))
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
            if tab > step or float(x) == 0:
                depth.append(float(x))
                xb = float(x)
                dens.append(float(y) - yb)
                yb = float(y)
                tab = 0
            if not found and float(y) > 0.5:
                v50 = x
                found = True

        std_dev = 0
        mean_value = 0
        for i in range(1, len(dens)):
            mean_value = float(summary.loc[summary['chrom'] == chrom, 'mean'].iloc[0])
            std_dev += dens[i] * (depth[i] - mean_value) ** 2
        std_dev = np.sqrt(std_dev)
        summary.loc[summary['chrom'] == chrom, 'cv'] = std_dev / mean_value
        if std_dev / mean_value > 1:
            summary.loc[summary['chrom'] == chrom, 'threshold'] = True
        summary.loc[summary['chrom'] == chrom, 'v50'] = v50

        # # to output
        # dtb[chrom].append({
        #     'ref_depth': [round(x, 3) for x in depth],
        #     'dens': [round(y, 3) for y in dens],
        #     'name': chr + (" (%.1f)" % float(v50)),
        # })

        # # test mean depth
        # exp = sum([d * p for d, p in zip(dens, depth)])
        # print(exp, count, summary[summary['chrom'] == chrom]['length'].values[0])

        dens_array = np.array(dens)
        peak_index = peak(dens_array)
        peak_index_new = find_peaks_cwt(dens_array, 5)

        # if times > 1000:
        #     plt.figure(figsize=(14, 7))
        #     plt.plot(depth, dens, label=chrom)
        #     plt.scatter(depth, dens, alpha=0.5, s=30)
        #     for i in peak_index:
        #         plt.scatter(depth[i], dens[i], color='red', s=50)
        #     for i in peak_index_new:
        #         plt.scatter(depth[i], dens[i], color='orange', alpha=0.5, s=50)
        #     plt.title('Probability Density Plot for ' + chrom)
        #     plt.xlabel('Depth')
        #     plt.ylabel('Probability')
        #     plt.legend([v50])
        #     plt.show()
        # if times == 1010:
        #     break
    return summary


# def peak(data):
#     p_data = np.zeros_like(data, dtype=np.int32)
#     count = data.shape[0]
#     arr_rowsum = []
#     for k in range(1, count // 2 + 1):
#         row_sum = 0
#         for i in range(k, count - k):
#             if data[i] > data[i - k] and data[i] > data[i + k]:
#                 row_sum -= 1
#         arr_rowsum.append(row_sum)
#     min_index = np.argmin(arr_rowsum)
#     max_window_length = min_index
#     for k in range(1, max_window_length + 1):
#         for i in range(k, count - k):
#             if data[i] > data[i - k] and data[i] > data[i + k]:
#                 p_data[i] += 1
#     return np.where(p_data == max_window_length)[0]


if __name__ == '__main__':
    main()
