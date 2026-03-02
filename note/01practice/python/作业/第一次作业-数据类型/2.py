import random
import math
import time
DARTS = int(2**20)
hits = 0

############补充代码#############
print("\n{:^50}".format("开始计算"))
#################################

start_time = time.perf_counter()
for i in range(1,DARTS):
    x, y = random.random(), random.random()
    dist = math.sqrt(x**2 + y**2)
    if dist <= 1.0:
        hits = hits + 1
    #进度百分比
    percentage = round(i/DARTS*100)
    star = "*"*round(i/DARTS*50)
    dot = "."*(50-round(i/DARTS*50))
    during_time = time.perf_counter() - start_time
    ############补充代码#############
    print("\r{:3}%[{}{}-> {:.2f}s]".format(percentage, star, dot, during_time), end="")
    #################################

############补充代码#############
print("\n{:^50}".format("计算结束"))
#################################

pi = 4 * (hits/DARTS)
print("Pi的值是", pi)
during_time = time.perf_counter() - start_time
print("程序运行时间是", during_time)

