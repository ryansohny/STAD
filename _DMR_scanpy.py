import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors
sc.settings.verbosity = 3
plt.rcParams['figure.figsize'] = (5,5)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])
%matplotlib

dmr = pd.read_csv("DMR_mat_imputed.csv", index_col=0).iloc[:,:-1]
dmr = sc.AnnData(dmr)
dmr.raw = dmr
dmr.layers['Percent_met'] = dmr.X
dmr.obs['TN'] = pd.DataFrame(['Normal']*84 + ['Tumor']*84, index=dmr.obs.index)
##### 다른 clinical data 집어넣은다음에 umap leiden에 넣어서 봐볼것!

sc.pp.scale(dmr)
sc.tl.pca(dmr, n_comps=100, zero_center=True)
sc.pl.pca(dmr, color='TN')
#sc.pl.pca_variance_ratio(dmr, log=True)

### Tests for varying degrees of n_neighbors needed ###
sc.pp.neighbors(dmr, n_neighbors=3, n_pcs=6)
sc.tl.umap(dmr, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')

sc.pp.neighbors(dmr, n_neighbors=15, n_pcs=6)
sc.tl.umap(dmr, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')


sc.pl.umap(dmr, color='TN')


sc.pp.neighbors(dmr, n_neighbors=10, n_pcs=10) ################################################################ 이게 좋은 것 같음
sc.tl.leiden(dmr, resolution=1.0, key_added='leiden_r1') ## 오우 이거 괜찮


# https://github.com/aertslab/pySCENIC/issues/357 이거랑 (Vascular Aging)
# https://doi.org/10.1016/j.celrep.2018.10.045 이거 archiving 해놓을것!! (Vascular Aging)
# https://github.com/aertslab/pySCENIC/issues/136 이것도 archiving ==> pySCENIC on bulk RNA-seq DATA!!!!!!!!!!!!!
# https://github.com/aertslab/pySCENIC/issues/169 이것도 archiving ==> multiple pySCENIC RUN?????? (Vascular Aging)
# https://github.com/aertslab/pySCENIC/find/master 여기에 .ipynb 들 
