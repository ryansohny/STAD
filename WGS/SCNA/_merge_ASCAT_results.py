from glob import glob
import pandas as pd

files = sorted(glob('*_T/Results_table_*csv'))

dfList = list()

for file in files:
        dfList.append(pd.read_csv(file, index_col=0))

df = pd.concat(dfList)

df.to_csv("ASCAT_Results_N84.csv")
