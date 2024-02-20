import matplotlib.pyplot as plt

# 提供数据
categories = ['DHH', 'DTA', 'DTC', 'DTH', 'DTM', 'DTT', 'DTX', 'DXX', 'RIX', 'RLC', 'RLG', 'RLX', 'XXX']
counts = [735, 1194, 1095976, 48828, 95785, 184898, 108199, 12509, 109873, 596077, 1586131, 147763, 285813]

# 计算比例
total_counts = sum(counts)
percentages = [count / total_counts * 100 for count in counts]

# 创建饼图
plt.pie(percentages, labels=categories, autopct='%1.1f%%', startangle=140)

# 添加标题
plt.title('Category Percentages')

# 显示图表
plt.show()
