import matplotlib.pyplot as plt

# 计算占比
total_TE = 150013
overlap_TE = 21364
non_overlap_TE = total_TE - overlap_TE

# 创建一个列表，包含两部分的值
sizes = [overlap_TE, non_overlap_TE]

# 创建饼图
fig1, ax1 = plt.subplots()
ax1.pie(sizes, labels=['Overlap TE', 'Non-overlap TE'], autopct='%1.1f%%',
        shadow=True, startangle=90)

# Equal aspect ratio ensures that pie is drawn as a circle.
ax1.axis('equal')

plt.show()
