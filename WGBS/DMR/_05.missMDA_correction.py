dfh = open("DMR_mat_imputed.csv", 'r')
rfh = open("DMR_mat_imputed_corrected.csv", 'w')

rfh.write(dfh.readline())
for i in dfh:
        line = i.strip('\n').split(',')
        met = list(map(lambda x: float(x), line[1:-1]))
        corrected_met = list(map(lambda x: '0' if x < 0 else ('100' if x > 100 else str(x)), met))
        rfh.write(line[0] + ',' + ','.join(corrected_met) + ',' + line[-1] + '\n')
        rfh.flush()

import pandas as pd
dmr = pd.read_csv("DMR_mat_imputed_corrected.csv", index_col=0)
dmr.iloc[84:, :-1].to_csv("DMR_mat_imputed_corrected_tumor.csv")
dmr.iloc[:84, :-1].to_csv("DMR_mat_imputed_corrected_normal.csv")
