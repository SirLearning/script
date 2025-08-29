# process the initial transposon
# TE annotation analysis
import os
import pandas as pd

data = pd.read_csv("new.txt", sep = '\t')
length_counts = data.groupby("length").size().reset_pwdindex(name='count')
print(length_counts)
length_counts.to_csv("length_counts.txt", sep = '\t', index = False)