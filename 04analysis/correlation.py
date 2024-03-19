import matplotlib.pyplot as plt
import pandas as pd

# 假设你的DataFrame是这样的：
df = pd.DataFrame({
    'column1': [1, 2, 3, 4, 5],
    'column2': [2, 3, 4, 5, 6]
})

plt.scatter(df['column1'], df['column2'])
plt.xlabel('column1')
plt.ylabel('column2')
plt.title('Scatter plot of column1 vs column2')
plt.show()