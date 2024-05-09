import matplotlib.pyplot as plt
import numpy as np

plt.style.use('seaborn-v0_8-deep')

fthrsd = [0.04, 0.2, 0.4, 1, 2, 3]
mapped = [58.31, 58.31, 58.31, 58.31, 58.30, 58.31]
nodp_mapped = [57.16, 55.79, 54.76, 52.61, 50.13, 48.27]
nodp_primary = [56.15, 54.72, 53.64, 51.37, 48.77, 46.80]
nodp_properly = [46.80, 46.06, 45.16, 43.03, 40.50, 38.59]

fig, ax = plt.subplots()
ax.figure.set_size_inches(12, 8)
ax.plot(fthrsd, mapped, label='mapped', alpha=0.8, linewidth=4)
ax.plot(fthrsd, nodp_mapped, label='mapped no duplication', alpha=0.8, linewidth=4)
ax.plot(fthrsd, nodp_primary, label='primary mapped no duplication', alpha=0.8, linewidth=4)
ax.plot(fthrsd, nodp_properly, label='properly mapped no duplication', alpha=0.8, linewidth=4)

ax.set_xlabel('depth', fontsize=20)
ax.set_ylabel('proportion of mapped reads (%)', fontsize=20)
ax.set_title('Depth & reads mapping', fontsize=28)
ax.legend(fontsize=18, framealpha=0.5)
plt.show()