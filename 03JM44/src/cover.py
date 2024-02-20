import matplotlib.pyplot as plt
import numpy as np

# 提供数据
categories = ['DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'DTX', 'DHH', 'DXX', 'RLC', 'RLG', 'RLX', 'RIX', 'RSX', 'XXX', 'no_match', 'Ontology']
countsCS = [53, 29750, 1861, 3457, 6945, 4331, 27, 482, 22463, 60531, 5918, 4437, 271, 7655, 1832, 52847]
countsEDTA = [1854, 4400, 4777, 5392, 11134, 0, 16490, 0, 4947, 11215, 7257, 1306, 0, 0, 0, 0]

index = np.arange(len(categories))
width = 0.8

# 计算百分比
total_counts1 = sum(countsCS)
percentages1 = [count / total_counts1 * 100 for count in countsCS]

total_counts2 = sum(countsEDTA)
percentages2 = [count / total_counts2 * 100 for count in countsEDTA]

# 创建水平条形图
fig, ax = plt.subplots()

# 绘制条形图
bars1 = ax.barh(index, countsCS, width, label='countsCS', color='#6780B6')
bars2 = ax.barh(index, countsEDTA, width, label='countsEDTA', color='#A3A1FB', alpha=0.5)

# 添加标签，显示百分比
for bar, count, percent in zip(bars1, countsCS, percentages1):
    ax.text(bar.get_width(), bar.get_y() + bar.get_height()/2, f'{percent:.1f}%', ha='left', va='center', fontsize=10)

# 添加标题和标签
ax.set_title('Category Percentages for CS and EDTA on chr1A')
ax.set_xlabel('Counts')
ax.set_ylabel('Categories')

# 设置y轴的刻度标签
ax.set_yticks(index)
ax.set_yticklabels(categories)

ax.legend()

# 显示图表
plt.show()