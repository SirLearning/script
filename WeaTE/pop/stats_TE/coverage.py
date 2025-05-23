import pandas as pd


def main():
    cover_name = 'transposon/6A.lib.fa.out'
    fai_1A_name = 'transposon/1A.lib.fa.fai'
    fai_6A_name = 'transposon/6A.lib.fa.fai'
    output_name = 'transposon/same_high.txt'

    fai_1A = pd.read_table(fai_1A_name, sep='\t', header=None)
    fai_1A.columns = ['name', 'length', 'offset', 'linebases', 'linewidth']
    fai_6A = pd.read_table(fai_6A_name, sep='\t', header=None)
    fai_6A.columns = ['name', 'length', 'offset', 'linebases', 'linewidth']

    cover = pd.read_table(cover_name, sep='\s+', header=None)
    cover.columns = ['SW score', 'perc div.', 'perc del.', 'perc ins.', 'Query', 'begin', 'end', '(left)',
                     'C+', 'Match', 'match_class', 'begin', 'end', 'left', 'ID']

    cover['query_class'] = cover['Query'].str.split('#').str[1]
    cover['Match'] = cover['Match'] + '#' + cover['match_class']
    cover['class'] = 'diff'
    cover['cover'] = 'low'
    threshold = 0.8

    for i in cover.index:
        print(fai_6A.loc[fai_1A['name'] == cover.loc[i, 'Query'], 'length'])
        if cover.loc[i, 'query_class'] == cover.loc[i, 'match_class']:
            cover.loc[i, 'class'] = 'same'
        cover.loc[i, 'query_cover'] = (cover.loc[i, 'end'] - cover.loc[i, 'begin'] + 1) / \
                                      fai_6A.loc[fai_6A['name'] == cover.loc[i, 'Query'], 'length']
        cover.loc[i, 'match_cover'] = (cover.loc[i, 'end'] - cover.loc[i, 'begin'] + 1) / \
                                      fai_1A.loc[fai_1A['name'] == cover.loc[i, 'Match'], 'length']
        if cover.loc[i, 'query_cover'] > threshold and cover.loc[i, 'match_cover'] > threshold:
            cover.loc[i, 'cover'] = 'high'

    with open(output_name, 'w') as f:
        f.write('Query\tMatch\tclass\tcover\n')
        for i in cover.index:
            if cover.loc[i, 'class'] == 'same' and cover.loc[i, 'cover'] == 'high':
                f.write(
                    cover.loc[i, 'Query'] + '\t' + cover.loc[i, 'Match'] + '\t' + cover.loc[i, 'class'] + '\t' + cover.loc[
                        i, 'cover'] + '\n')

    grouped = cover.groupby(['Query', 'Match', 'class', 'cover']).size()
    print(grouped)


if __name__ == '__main__':
    main()
