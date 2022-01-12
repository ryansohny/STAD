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

# DMR table for ALL (tumor and normal)
dmr = pd.read_csv("DMR_mat_imputed.csv", index_col=0).iloc[:,:-1]
dmr = sc.AnnData(dmr)
dmr.raw = dmr
dmr.layers['Percent_met'] = dmr.X
# np.ndarray.min(dmr.raw.X) ==> 전체 table에서 minimum value 값 (maximum은 min==> max)
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



# DMR table for tumor only
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

sc.pp.neighbors(dmr_t, n_neighbors=15, n_pcs=12)
sc.tl.leiden(dmr_t, resolution=1.0, key_added='leiden_r1')
sc.tl.umap(dmr_t, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
sc.pl.umap(dmr_t, color='Lauren', add_outline=False)


df = dmr_t.obs[['WHO', 'leiden_r1']]
pd.crosstab(df['WHO'], df['leiden_r1'], normalize=1).T.plot.bar(stacked=True)


dmr_met = pd.DataFrame(dmr_t.raw.X.T, index=dmr_t.var.index, columns=dmr_t.obs.index)

col_colors = list(dict(zip(list(dmr_t.obs['leiden_r1'].value_counts().index), dmr_t.uns['leiden_r1_colors']))[x] for x in dmr_t.obs['leiden_r1'])
#sns.clustermap(dmr_met, method='ward', metric='euclidean', z_score=None, standard_scale=0, cmap=cmap, xticklabels=True)
g = sns.clustermap(dmr_met, method='ward', metric='euclidean', z_score=None, standard_scale=0, cmap=cmap, xticklabels=True, col_colors=col_colors)
g = sns.clustermap(dmr_met, method='ward', metric='euclidean', z_score=None, standard_scale=None, cmap=cmap, xticklabels=True, col_colors=col_colors)
dmr_t_colorder = g.dendrogram_col.reordered_ind
dmr_t_roworder = g.dendrogram_row.reordered_ind

df = sc.get.obs_df(dmr_t, keys=['EpiBurden', 'leiden_r1'], use_raw=True)
ax = sns.boxplot(x='leiden_r1', y='EpiBurden', data=df)
ax = sns.stripplot(x='leiden_r1', y='EpiBurden', data=df, color=".3")

# DMR table for normal only (for clustermap)
dmr_n = sc.AnnData(pd.read_csv("DMR_mat_normal_imputed.csv", index_col=0))
dmr_n.raw = dmr_n
dmr_n.layers['Percent_met'] = dmr_n.X
dmr_n_met = pd.DataFrame(dmr_n.raw.X.T, index=dmr_n.var.index, columns=dmr_n.obs.index)
dmr_n_met.iloc[dmr_t_roworder, dmr_t_colorder]
sns.heatmap(dmr_n_met.iloc[dmr_t_roworder, dmr_t_colorder], cmap=cmap, xticklabels=True)


# RNA-seq
rna = sc.AnnData(pd.read_csv("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/STAD_SNUH_Tumor_leiden_vst.csv", index_col=0).T)
rna.raw = rna
rna.layers['vst'] = rna.X
#rna.obs = dmr_t.obs # 이렇게 하면, rna.obs에 뭔가 update할 때마다 dmr_t.obs에도 그게 똑같이 들어감
rna.obs = dmr_t.obs.copy()
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps=83, zero_center=True)
sc.pl.pca(rna, color='leiden_r1')
# pca_variance_rna = pd.DataFrame(rna.uns['pca']['variance_ratio'], index=list(map(lambda x: 'PC' + str(x), list(range(1,84)))), columns=['Variance_ratio'])
# np.sum(pca_variance_rna.values.flatten()[:25]) 

sc.pp.neighbors(rna, n_neighbors=20, n_pcs=15)
sc.tl.leiden(rna, resolution=0.5, key_added='rna_leiden_r05')
sc.tl.leiden(rna, resolution=0.75, key_added='rna_leiden_r075')
sc.tl.umap(rna, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
sc.pl.umap(rna, color=['leiden_r1', 'rna_leiden_r05', 'rna_leiden_r075'], add_outline=False)


# https://github.com/aertslab/pySCENIC/issues/357 이거랑 (Vascular Aging)
# https://doi.org/10.1016/j.celrep.2018.10.045 이거 archiving 해놓을것!! (Vascular Aging)
# https://github.com/aertslab/pySCENIC/issues/136 이것도 archiving ==> pySCENIC on bulk RNA-seq DATA!!!!!!!!!!!!!
# https://github.com/aertslab/pySCENIC/issues/169 이것도 archiving ==> multiple pySCENIC RUN?????? (Vascular Aging)
# https://github.com/aertslab/pySCENIC/find/master 여기에 .ipynb 들 
# UMAP 에 대한 영상인데, UMAP 만든 사람이 좋다고 함 https://www.youtube.com/watch?v=6BPl81wGGP8
