import matplotlib.pyplot as plt
import numpy as np

# 计算半径
total_TE = 150013
overlap_TE = 21364
overlap_TE_same_type = 10752

radius_total = np.sqrt(total_TE / np.pi)
radius_overlap = np.sqrt(overlap_TE / np.pi)
radius_overlap_same_type = np.sqrt(overlap_TE_same_type / np.pi)

# 创建图形和轴
fig, ax = plt.subplots()

# 创建三个圆
circle_total = plt.Circle((0, 0), radius_total, color='b', fill=False)
circle_overlap = plt.Circle((0, 0), radius_overlap, color='r', fill=False)
circle_overlap_same_type = plt.Circle((0, 0), radius_overlap_same_type, color='g', fill=False)

# 将圆添加到图形中
ax.add_patch(circle_total)
ax.add_patch(circle_overlap)
ax.add_patch(circle_overlap_same_type)

# 设置轴的比例和限制
ax.set_aspect('equal', adjustable='datalim')
ax.plot()  # Causes an autoscale update.
plt.show()