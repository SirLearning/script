import os

import pandas as pd
from cookiecutter.prompt import read_user_dict


def org_China(China_data, file_name):
    file_data = pd.read_csv('phenotype/data/Watseq/China_phenotype/' + file_name, sep='\t', header=0)
    if China_data.empty:
        China_data = file_data
    else:
        China_data = pd.merge(China_data, file_data, on='IID')
    return China_data


def check_trait(China_data):
    ref_trait = pd.read_csv('phenotype/data/Watseq/trait_information.csv', sep=',', header=0)
    # print(ref_trait.columns)
    measured_num = 0
    for column in China_data.columns:
        if column.split('_')[0] not in ref_trait['Abbreviation used in plots'].values:
            print(column + ' is not in the reference trait list')
    for trait in ref_trait['Abbreviation used in plots']:
        if trait not in [col.split('_')[0] for col in China_data.columns]:
            print(trait + ' is not in the China data')
            continue
        measured_num += 1
    print('The number of measured traits is ' + str(measured_num))


def main():
    China_data = pd.DataFrame()
    # print(len(os.listdir('phenotype/data/Watseq/China_phenotype')))
    for file_name in os.listdir('phenotype/data/Watseq/China_phenotype'):
        China_data = org_China(China_data, file_name)
    # China_data.to_csv('phenotype/data/Watseq/China_data_output.csv', index=False)
    check_trait(China_data)


if __name__ == '__main__':
    main()

