import matplotlib.pyplot as plt
import numpy as np

plt.style.use('seaborn-v0_8-deep')

fthrsd = np.arange(1, 0, -0.05)
nav = [54542, 37598, 32698, 29581, 26905, 24715, 22680, 20930, 19507, 17839,
       16150, 14810, 13661, 12625, 11502, 10648, 9688, 8736, 7687, 6784]
anv = [53448, 36507, 31609, 28498, 25841, 23659, 21718, 20087, 18687, 17124,
       15396, 14204, 13021, 11880, 10962, 10113, 9102, 8219, 7466, 6596]

sum_num = (320051+318956)/2

nav = [i/sum_num for i in nav]
anv = [i/sum_num for i in anv]

fig, ax = plt.subplots()
ax.figure.set_size_inches(12, 8)
ax.plot(fthrsd, nav, label='Only in N-sliced chr1A', alpha=0.8)
ax.plot(fthrsd, anv, label='Only in intact chr1A', alpha=0.8)
ax.axvline(x=0.8, color='r', linestyle='--')

ax.annotate(f'{26905}', fontsize=18, color='darkred', xy=(0.8, 26905/sum_num), xytext=(0.68, 29005/sum_num), textcoords='data', xycoords='data', arrowprops=dict(facecolor='black', shrink=0.05, color='darkred'))
ax.annotate(f'{25841}', fontsize=18, color='darkred', xy=(0.8, 25841/sum_num), xytext=(0.84, 23041/sum_num), textcoords='data', xycoords='data', arrowprops=dict(facecolor='black', shrink=0.05, color='darkred'))

ax.set_xlabel('threshold', fontsize=20)
ax.set_ylabel('TE_differ proportion', fontsize=20)
ax.set_title('Threshod & TE difference', fontsize=28)
ax.legend(fontsize=18, framealpha=0.5)
plt.show()