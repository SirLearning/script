from collections import defaultdict
from pandas import DataFrame, Series
import pandas as pd
import numpy as np

path = '../resources/transposon/JM44.repeat.masked.gff'
attribute = pd.read_table('../resources/transposon/JM44.repeat.masked.gff', sep ='\t', header = None, usecols = [8])
Classification = pd.read_table('../resources/transposon/JM44.repeat.masked.gff', sep =';', header = None, usecols = [2])
print(Classification)