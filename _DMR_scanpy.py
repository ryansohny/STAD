import scanpy as sc
import pandas as pd
import numpy as np
import seaborn sns
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
clinic_info = pd.read_csv('/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/2021_WC300_clinical_information_Xadded_NormalCombined_new.csv', index_col='ID')
dmr.obs = clinic_info
##### 다른 clinical data 집어넣은다음에 umap leiden에 넣어서 봐볼것!

sc.pp.scale(dmr)
sc.tl.pca(dmr, n_comps=100, zero_center=True)
sc.pl.pca(dmr, color='TN')
sc.pl.pca(dmr, color='TN', add_outline=True, size=100, palette={'Normal':'Blue', 'Tumor':'Red'})
#sc.pl.pca_variance_ratio(dmr, log=True)
pca_variance = pd.DataFrame(dmr.uns['pca']['variance_ratio'], index=list(map(lambda x: 'PC' + str(x), list(range(1,101)))), columns=['Variance_ratio'])
np.sum(pca_variance.values.flatten()[:11])
# 0.7016306 ==> PC1~11까지가 variance 70% 설명

sc.pp.neighbors(dmr, n_neighbors=10, n_pcs=11)
sc.tl.leiden(dmr, resolution=1.0, key_added='leiden_r1') 
sc.tl.umap(dmr, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
sc.pl.umap(dmr, color='TN', add_outline=True, size=100, palette={'Normal':'Blue', 'Tumor':'Red'})


dmr_t = sc.AnnData(pd.read_csv("DMR_mat_tumor_imputed.csv", index_col=0))
dmr_t.raw = dmr
dmr_t.layers['Percent_met'] = dmr_t.X
clinic_info = pd.read_csv('/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/2021_WC300_clinical_information_Xadded.csv', index_col='ID')
clinic_info[['T.stage', 'N.stage', 'M.stage', 'TNM.stage', 'HER2IHC']] = clinic_info[['T.stage', 'N.stage', 'M.stage', 'TNM.stage', 'HER2IHC']].astype(str)
clinic_info = clinic_info.iloc[:,np.r_[0:17, -1]]
dmr_t.obs = clinic_info
sc.pp.scale(dmr_t)
sc.tl.pca(dmr_t, n_comps=83, zero_center=True)
sc.pl.pca(dmr_t, color='Lauren')
# sc.pl.pca_variance_ratio(dmr_t, log=True)

#pca_variance_t = pd.DataFrame(dmr_t.uns['pca']['variance_ratio'], index=list(map(lambda x: 'PC' + str(x), list(range(1,84)))), columns=['Variance_ratio'])
# np.sum(pca_variance_t.values.flatten()[:12])

sc.pp.neighbors(dmr_t, n_neighbors=5, n_pcs=12)
sc.tl.leiden(dmr_t, resolution=1.0, key_added='leiden_r1')
sc.tl.umap(dmr_t, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
sc.pl.umap(dmr_t, color='Lauren', add_outline=False)


df = dmr_t.obs[['WHO', 'leiden_r1']]
pd.crosstab(df['WHO'], df['leiden_r1'], normalize=1).T.plot.bar(stacked=True)

# https://github.com/aertslab/pySCENIC/issues/357 이거랑 (Vascular Aging)
# https://doi.org/10.1016/j.celrep.2018.10.045 이거 archiving 해놓을것!! (Vascular Aging)
# https://github.com/aertslab/pySCENIC/issues/136 이것도 archiving ==> pySCENIC on bulk RNA-seq DATA!!!!!!!!!!!!!
# https://github.com/aertslab/pySCENIC/issues/169 이것도 archiving ==> multiple pySCENIC RUN?????? (Vascular Aging)
# https://github.com/aertslab/pySCENIC/find/master 여기에 .ipynb 들 
