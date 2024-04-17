import numpy as np
import matplotlib.pyplot as plt

# 设置指数分布的参数，即平均值或者说是事件发生的间隔
scale = 1.0
n = 4
score = 2
ptl = []
for i in range(1, 8):
    mean = 0.85**i*(2+0.5*i)
    print(mean)


# # 生成10000个符合指数分布的随机数
# data = np.random.exponential(scale, 10000)
# print(data)
#
# # 绘制直方图，bins参数表示条形的数量
# plt.hist(data, bins=100, density=True)
# plt.show()