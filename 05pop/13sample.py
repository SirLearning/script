import numpy as np

for i in np.arange(0.2, 3.2, 0.2):
    cal = i/15.79*10000
    print(f'{i:.1f} {cal:.0f}')