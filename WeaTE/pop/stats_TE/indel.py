import matplotlib.pyplot as plt
import pandas as pd


def mod(name):
    ins_ana = pd.read_csv("WeaTE/data/indel/insertion_reads_kg." + str(name) + ".txt", header=None, sep='\t')
    del_ana = pd.read_csv("WeaTE/data/indel/deletion_reads_kg." + str(name) + ".txt", header=None, sep='\t')
    ins_ana.columns = ['ID', 'reads']
    del_ana.columns = ['ID', 'reads']
    ins_ana['test'] = str(name)
    del_ana['test'] = str(name)
    return ins_ana, del_ana


def num_plot(ins_ana, del_ana):
    fig, ax = plt.subplots()
    ins_ana.groupby('test')['reads'].sum().plot(kind='bar', ax=ax, color='b', alpha=0.7, label='insertion')
    del_ana.groupby('test')['reads'].sum().plot(kind='bar', ax=ax, color='r', alpha=0.7, label='deletion')
    plt.show()


def main():
    for i in range(1, 5):
        if i == 1:
            ins_ana, del_ana = mod(i)
        else:
            ins_, del_ = mod(i)
            ins_ana = pd.concat([ins_ana, ins_])
            del_ana = pd.concat([del_ana, del_])
    num_plot(ins_ana, del_ana)


if __name__ == "__main__":
    main()