import scanpy as sc
import pandas as pd
dmr = pd.read_csv("DMR_mat_imputed.csv", index_col=0).iloc[:,:-1]
dmr = sc.AnnData(dmr)
dmr.raw = dmr
dmr.layers['Percent_met'] = dmr.X
dmr.obs['TN'] = pd.DataFrame(['Normal']*84 + ['Tumor']*84, index=dmr.obs.index)

sc.pp.scale(dmr)
sc.tl.pca(dmr, n_comps=100, zero_center=True)
sc.pl.pca(dmr, color='TN')
