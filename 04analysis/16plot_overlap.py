import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# overlap in mean 80%
total_TE = 150013

overlap = pd.read_table('data/overlap.txt', sep='\s+', header=None)
overlap.columns = ['method', 'overlap', 'same type']

overlap['overlap'] = overlap['overlap'] / total_TE
overlap['same type'] = overlap['same type'] / total_TE

print(overlap['same type'])

# plot
plt.figure(figsize=(12, 8))
sns.barplot(data=overlap, x='method', y='overlap', label='overlap', width=0.2)
sns.barplot(data=overlap, x='method', y='same type', label='same type', width=0.2, alpha=0.6)
plt.title('Number of TEs in 80% overlap')
plt.legend()
plt.show()
