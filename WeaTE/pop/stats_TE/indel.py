import matplotlib.pyplot as plt
import pandas as pd


def mod(name):
    ins_ana = pd.read_csv("WeaTE/transposon/indel/insertion_reads_kg." + str(name) + ".txt", header=None, sep='\t')
    del_ana = pd.read_csv("WeaTE/transposon/indel/deletion_reads_kg." + str(name) + ".txt", header=None, sep='\t')
    ins_ana.columns = ['insertion', 'reads']
    del_ana.columns = ['deletion', 'reads']
    ins_ana['test'] = str(name)
    del_ana['test'] = str(name)
    ins_ana.drop(columns='reads', inplace=True)
    del_ana.drop(columns='reads', inplace=True)
    return ins_ana, del_ana


def num_plot(ins_ana, del_ana):
    # plt.style.use('seaborn')
    fig, ax = plt.subplots()
    ins_ana.groupby('test').count().plot(kind='bar', ax=ax, color='#f38121', alpha=0.8, label='insertion')
    del_ana.groupby('test').count().plot(kind='bar', ax=ax, color='#1b78b3', alpha=0.6, label='deletion')
    plt.ylabel('Number')
    plt.title('Number of insertion and deletion detected')
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