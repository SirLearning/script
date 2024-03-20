import matplotlib.pyplot as plt
import pandas as pd

plt.style.use('bmh')

# 假设我们有一个包含染色体长度的数据框
seq_df = pd.DataFrame({
    'sup_fam': ['RLG', 'RLX'],
    'length': [3120054, 3120054]
})

# 假设我们有一个包含基因位置的数据框
TE_dtb = pd.read_table('data/cut.gff3', sep='\s+', header=None)
TE_dtb.columns = ['chrom', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes', 'sup_fam']

# 过滤出sup_fam列为'RLG'或'RLX'的行
TE_dtb = TE_dtb[TE_dtb['sup_fam'].isin(['RLG', 'RLX'])]

fig, ax = plt.subplots()
ax.figure.set_size_inches(20, 10)

# 绘制染色体
ax.barh(seq_df['sup_fam'], seq_df['length'], color='gray', alpha=0.5)

# 绘制TE
# for te in seq_df['sup_fam'].unique():
#     subset = TE_dtb[TE_dtb['sup_fam'] == te]
#     for i in range(0, len(subset)-1):
#         # 获取染色体的长度
#         te_length = 3120054
#         # 计算TE的相对位置
#         start_pos = (int(subset['start'].iloc[i]) - 92963774)
#         end_pos = (int(subset['end'].iloc[i]) - 92963774)
#         ax.vlines(x=[start_pos, end_pos], ymin=te, ymax=te, color='red')

# 绘制TE
for idx, te in enumerate(seq_df['sup_fam'].unique()):
    subset = TE_dtb[TE_dtb['sup_fam'] == te]
    for i in range(len(subset)):
        # 获取染色体的长度
        te_length = seq_df[seq_df['sup_fam'] == te]['length'].values[0]
        # 计算TE的相对位置
        start_pos = subset['start'].iloc[i] / te_length
        end_pos = subset['end'].iloc[i] / te_length
        ax.vlines(x=[start_pos, end_pos], ymin=idx, ymax=idx, color='red')

plt.xlabel('Relative Position')
plt.show()