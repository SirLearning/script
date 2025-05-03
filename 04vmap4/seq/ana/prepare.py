import pandas as pd
import numpy as np


def size_determine(size_file, size):
    size_df = pd.read_csv(size_file, sep='\t', header=None)
    size_df.columns = ['size', 'file']
    size_df['Gb'] = size_df['size'].str.split('G').str[0].astype(float)
    sum_size = 0
    file_name = []
    for i in range(size_df.shape[0]):
        sum_size += size_df['Gb'][i]
        if sum_size >= size:
            break
        file_name.append(size_df['file'][i].split('/')[1])
    return file_name


def print_out_size(number):
    file_name = size_determine('seq/data/fq.size', number)
    for i in range(len(file_name)):
        print(file_name[i])
    for i in range(len(file_name)):
        if i == 0:
            print('(' + file_name[i], end=' ')
        elif i == len(file_name) - 1:
            print(file_name[i] + ')')
        else:
            print(file_name[i], end=' ')


def check_bam(name1, name2):
    ref_md5 = pd.read_csv('seq/data/' + name1 + '.md5', sep=r'\s+', header=None)
    ref_md5.columns = ['md5', 'file']
    files = pd.read_csv('seq/data/' + name2 + '.md5', sep=r'\s+', header=None)
    files.columns = ['md5', 'file']
    for i in range(files.shape[0]):
        if files['md5'][i] != ref_md5.loc[ref_md5['file'] == files['file'][i], 'md5'].values:
            print(files['file'][i] + ' is incorrect')


def check_md5(name):
    # define ref_md5
    ## check .fq
    ref_md5 = pd.read_csv('seq/data/CRA012590.csv', sep=',', header=0)
    ref_md5['file1'] = ref_md5['Read filename 1'].str.split().str[0]
    ref_md5['file2'] = ref_md5['Read filename 2'].str.split().str[0]
    # define files
    # files = pd.read_csv('seq/data/' + name + '.md5', sep=r'\s+', header=None)
    files = pd.read_csv('seq/data/' + name + '.md5', sep=r'\s+', header=None)
    files.columns = ['md5', 'file']
    files['file'] = files['file'].str.split('/').str[-1]
    files['to_ref'] = np.zeros(files.shape[0])
    files['url'] = ''
    pass_md5 = 0
    # print(ref_md5.loc[ref_md5['file1'] == 'CRR876744_f1.fastq.gz', 'Read file1 MD5'].values)
    line_md5 = ''
    for i in range(files.shape[0]):
        if files['file'][i].endswith('f1.fastq.gz') or files['file'][i].endswith('f1.fq.gz'):
            line_md5 = ref_md5.loc[ref_md5['file1'] == files['file'][i], 'Read file1 MD5'].values
            line_sample = ref_md5.loc[ref_md5['file1'] == files['file'][i], 'Run title'].values
            # files.loc[i, 'url'] = ref_md5.loc[ref_md5['file1'] == files['file'][i], 'DownLoad Read file1'].values
            # print(ref_md5.loc[ref_md5['file1'] == files['file'][i], 'Read file1 MD5'].values)
            # print(files['file'][i])
            # print(ref_md5['file1'][i])
            # line_md5 = match_site['Read file1 MD5']
        else:
            line_md5 = ref_md5.loc[ref_md5['file2'] == files['file'][i], 'Read file2 MD5'].values
            line_sample = ref_md5.loc[ref_md5['file2'] == files['file'][i], 'Run title'].values
            # files.loc[i, 'url'] = ref_md5.loc[ref_md5['file2'] == files['file'][i], 'DownLoad Read file2'].values
            # line_md5 = match_site['Read file2 MD5']
        if line_md5.size > 0 and line_md5 == files['md5'][i]:
            # print(files['file'][i] + ' pass')
            files.loc[i, 'to_ref'] = 1
        else:
            print(files['file'][i] + ' is incorrect')
            # with open('seq/data/wrong_'+name+'.txt', 'a') as f:
            #     f.write(files['url'][i] + '\t' + files['file'][i].split('_')[0] + '\n')
            #     f.write(files['file'][i].split('_')[0] + '\n')
                # f.write(files['file'][i].split('_')[0] + '\t' + line_sample[0] + '\n')
        if i%2 == 1:
            if files['to_ref'][i] + files['to_ref'][i-1] == 2:
                with open('seq/data/pass_'+name+'.txt', 'a') as f:
                # with open('seq/data/pass_md5_usb.txt', 'a') as f:
                    f.write(files['file'][i].split('_')[0] + '\t' + line_sample[0] + '\n')
                pass_md5 += 1
            elif files['to_ref'][i] + files['to_ref'][i-1] == 0:
                print('Both files are incorrect')
            #     with open('seq/data/wrong_'+name+'.txt', 'a') as f:
            #         f.write(files['file'][i].split('_')[0] + '\n')
            # elif files['to_ref'][i] + files['to_ref'][i-1] == 1 and files['to_ref'][i] == 0:
            #     with open('seq/data/wrong_'+name+'.txt', 'a') as f:
            #         f.write(files['file'][i].split('_')[0] + '/' + files['file'][i] + '\n')
            # elif files['to_ref'][i] + files['to_ref'][i-1] == 1 and files['to_ref'][i-1] == 0:
            #     with open('seq/data/wrong_'+name+'.txt', 'a') as f:
            #         f.write(files['file'][i-1].split('_')[0] + '/' + files['file'][i-1] + '\n')
        line_md5 = ''
        line_sample = ''
    print('Rate = ' + str(pass_md5/files.shape[0]*200) + '%')
    # Rate_data1 = 38.429752066115704%
    # Rate_data2 = 27.300613496932513%
    # Rate_115 = 100.0%


def choose_v1():
    v1_size = pd.read_csv('seq/data/V1_size.csv', sep=',', header=0)
    v1_size['class'] = v1_size['Sample'].str.extract(r'([A-Za-z]+)', expand=False)
    v1_size['size'] = v1_size['Bytes'].str.split().str[0].astype(float)
    # print(v1_size.columns)  # Debugging line to check columns
    # print(v1_size.tail())  # Debugging line to check the first few rows
    v1 = v1_size.groupby('class')['size'].sum()
    print(v1.head())
    class_sum = pd.DataFrame({'A': [0], 'B': [0], 'D': [0], 'TW': [0]})
    old_class = ''
    # for i in range(v1_size.shape[0]):
    #     class_sum.loc[0, v1_size['class'][i]] += v1_size['size'][i]
    #     if class_sum.loc[0, v1_size['class'][i]] >= 400 and v1_size['class'][i] != old_class:
    #         print(v1_size['Sample'][i])
    #         old_class = v1_size['class'][i]
    list_name = ["A001", "A002", "A003", "A004", "A005", "A006", "A007", "A008", "A017", "A019", "A031", "A043", "A057", "A058", "A059", "A060", "A061", "A063", "A088", "A089", "A090", "A091", "B001", "B002", "B003", "B004", "B005", "B006", "B007", "B008", "B009", "B010", "B011", "B012", "B013", "B014", "B015", "B016", "B017", "B018", "B019", "B020", "B021", "B022", "B023", "B024", "B025", "TW001", "TW002", "TW003", "TW004", "TW005", "TW006", "TW007", "TW008", "TW009", "TW010", "TW011", "TW012", "TW013", "TW014", "TW015", "TW016", "TW017", "TW018", "TW019", "TW020", "TW021", "TW022", "D001", "D002", "D003", "D004", "D005", "D006", "D007", "D008", "D009", "D010", "D011", "D012", "D013", "D014", "D015", "D016", "D017", "D018", "D019", "D020", "D021", "D022", "D023", "D024", "D025", "D026"]
    sum_size = 0.0
    for i in list_name:
        # print(v1_size.loc[v1_size['Sample'] == i, 'size'])
        sum_size += v1_size.loc[v1_size['Sample'] == i, 'size'].values[0]
    print("final size = " + str(sum_size))


def count_bam():
    ref_md5 = pd.read_csv('seq/data/CRA012590.csv', sep=',', header=0)
    for name in ['34', '35', '36', '37', '38', '39', '40', '41']:
        with open('seq/data/pass_' + name + '.txt', 'r') as f:
            lines = [line.strip() for line in f.readlines() if not line.startswith('#')]
            bams = pd.DataFrame(lines)
            bams = bams[0].str.split('\t', expand=True).dropna()
            bams.columns = ['fq', 'sample']
            bams['finished'] = bams['sample'].str.split(r'\s+').str[-1]
            for index, row in bams.iterrows():
                if row['finished'] == "#":
                    ref_md5.loc[ref_md5['Accession'] == row['fq'], 'finished'] = "done"
                    ref_md5.loc[ref_md5['Accession'] == row['fq'], 'hardware'] = name
    with open('seq/data/count_bam.txt', 'w', newline='') as f:
        f.write(ref_md5[['Accession', 'finished', 'hardware']].to_csv(sep='\t', index=False, header=False))
    print('Finished number = ' + str(ref_md5['finished'].notna().sum()))
    with open('seq/data/ignored_bam.txt', 'w', newline='') as f:
        ignored_bam = ref_md5.loc[ref_md5['finished'].isna(), 'Accession'].sort_values()
        # print(ignored_bam.iloc[1])
        f.write(ignored_bam.to_csv(sep='\t', index=False, header=False))

def main():
    # print_out_size(3300)
    # check_md5('40.2')
    # check_bam('usb', 'old')
    # choose_v1()
    count_bam()


if __name__ == '__main__':
    main()

