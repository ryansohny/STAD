import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors
sc.settings.verbosity = 3
plt.rcParams['figure.figsize'] = (5,5)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])
%matplotlib

###### DMR table for ALL (tumor and normal)
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



###### DMR table for tumor only
dmr_t = sc.AnnData(pd.read_csv("DMR_mat_tumor_imputed.csv", index_col=0))
dmr_t.raw = dmr_t
dmr_t.layers['Percent_met'] = dmr_t.X
clinic_info = pd.read_csv('/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/2021_WC300_clinical_information_Xadded_new.csv', index_col='ID')
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
sc.tl.leiden(dmr_t, resolution=0.5, key_added='leiden_r05')
sc.tl.umap(dmr_t, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
sc.pl.umap(dmr_t, color='Lauren', add_outline=False)

leiden_name_change_dict1 = {'4': 'leiden_A',
                           '1': 'leiden_B',
                           '0': 'leiden_C',
                           '3': 'leiden_D',
                           '2': 'leiden_E'}
leiden_name_change_dict2 = {'1': 'leiden_A',
                           '0': 'leiden_B',
                           '2': 'leiden_C'}

dmr_t.obs['DMR_leiden'] = dmr_t.obs['leiden_r1'].map(lambda x: leiden_name_change_dict1[x]).astype('category')
dmr_t.obs['DMR_leiden2'] = dmr_t.obs['leiden_r05'].map(lambda x: leiden_name_change_dict2[x]).astype('category')
sc.pl.umap(dmr_t, color='DMR_leiden', add_outline=False, legend_loc='on data')

lin = tuple(sorted(list(dmr_t.obs['DMR_leiden'].values.unique())))
dmr_t.obs['DMR_leiden'] = dmr_t.obs['DMR_leiden'].cat.reorder_categories(list(lin), ordered=True)
color_dict = {
    "leiden_A": "#fdbf6f",
    "leiden_B": "#ff7f0e",
    "leiden_C": "#1f77b4",
    "leiden_D": "#d62728",
    "leiden_E": '#2ca02c',
}

lin = tuple(sorted(list(dmr_t.obs['DMR_leiden2'].values.unique())))
dmr_t.obs['DMR_leiden2'] = dmr_t.obs['DMR_leiden2'].cat.reorder_categories(list(lin), ordered=True)
color_dict = {
    "leiden_A": "#9467bd",
    "leiden_B": "#a6cee3",
    "leiden_C": "#b15928",
}

sns.boxplot(data=dmr_t.obs, x='DMR_leiden', y='EpiBurden', palette=color_dict)
sns.stripplot(data=dmr_t.obs, x='DMR_leiden', y='EpiBurden', jitter=True, marker='o', color='black', alpha=0.8)

df = dmr_t.obs[['WHO', 'Lauren', 'DMR_leiden', 'Cellularity']]

# DMR_leiden clusters and WHO classification
ax = pd.crosstab(df['WHO'], df['DMR_leiden'], normalize=1).T.plot.bar(stacked=True, color=dict(zip(dmr_t.obs['WHO'].unique(), sns.color_palette("Set2", len(dmr_t.obs['WHO'].unique())))))
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))
plt.tight_layout()

# DMR_leiden clusters and Lauren's classification
ax = pd.crosstab(df['Lauren'], df['DMR_leiden'], normalize=1).T.plot.bar(stacked=True, color=dict(zip(dmr_t.obs['Lauren'].unique(), sns.color_palette("Accent", len(dmr_t.obs['Lauren'].unique())))))
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))
plt.tight_layout()

# DMR_leiden clusters2 and WHO classification
ax = pd.crosstab(dmr_t.obs['WHO'], dmr_t.obs['DMR_leiden2'], normalize=1).T.plot.bar(stacked=True, color=dict(zip(dmr_t.obs['WHO'].unique(), sns.color_palette("Set2", len(dmr_t.obs['WHO'].unique())))))
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))
plt.tight_layout()

# DMR_leiden clusters2 and Lauren's classification
ax = pd.crosstab(dmr_t.obs['Lauren'], dmr_t.obs['DMR_leiden2'], normalize=1).T.plot.bar(stacked=True, color=dict(zip(dmr_t.obs['Lauren'].unique(), sns.color_palette("Accent", len(dmr_t.obs['Lauren'].unique())))))
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))
plt.tight_layout()



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

###### DMR table for normal only (for clustermap)
dmr_n = sc.AnnData(pd.read_csv("DMR_mat_normal_imputed.csv", index_col=0))
dmr_n.raw = dmr_n
dmr_n.layers['Percent_met'] = dmr_n.X
dmr_n_met = pd.DataFrame(dmr_n.raw.X.T, index=dmr_n.var.index, columns=dmr_n.obs.index)
dmr_n_met.iloc[dmr_t_roworder, dmr_t_colorder]
sns.heatmap(dmr_n_met.iloc[dmr_t_roworder, dmr_t_colorder], cmap=cmap, xticklabels=True)


###### RNA-seq
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

sc.pp.neighbors(rna, n_neighbors=15, n_pcs=18)
sc.tl.leiden(rna, resolution=0.5, key_added='rna_leiden_r05')
sc.tl.leiden(rna, resolution=0.6, key_added='rna_leiden_r06')
sc.tl.leiden(rna, resolution=0.75, key_added='rna_leiden_r075')
sc.tl.umap(rna, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
sc.pl.umap(rna, color=['leiden_r1', 'rna_leiden_r05', 'rna_leiden_r06', 'rna_leiden_r075'], add_outline=False)




###### DMR table for tumor-normal
dmr_tn = sc.AnnData(pd.read_csv("DMR_mat_tumor-normal_imputed.csv", index_col=0))
dmr_tn.raw = dmr_tn
dmr_tn.layers['met_diff'] = dmr_tn.X
clinic_info = pd.read_csv('/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/2021_WC300_clinical_information_Xadded_new.csv', index_col='ID')
clinic_info[['T.stage', 'N.stage', 'M.stage', 'TNM.stage', 'HER2IHC']] = clinic_info[['T.stage', 'N.stage', 'M.stage', 'TNM.stage', 'HER2IHC']].astype(str)
clinic_info = clinic_info.iloc[:,np.r_[0:17, -1]]
dmr_tn.obs = clinic_info
sc.pp.scale(dmr_tn)
sc.tl.pca(dmr_tn, n_comps=83, zero_center=True)
sc.pl.pca(dmr_tn, color='Lauren')
# sc.pl.pca_variance_ratio(dmr_tn, log=True)

#pca_variance_tn = pd.DataFrame(dmr_tn.uns['pca']['variance_ratio'], index=list(map(lambda x: 'PC' + str(x), list(range(1,84)))), columns=['Variance_ratio'])
# np.sum(pca_variance_tn.values.flatten()[:12])

sc.pp.neighbors(dmr_tn, n_neighbors=15, n_pcs=17)
sc.tl.leiden(dmr_tn, resolution=1.0, key_added='tn_leiden_r1')
sc.tl.umap(dmr_tn, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
sc.pl.umap(dmr_tn, color='tn_leiden_r1', add_outline=False)


# https://github.com/aertslab/pySCENIC/issues/357 이거랑 (Vascular Aging)
# https://doi.org/10.1016/j.celrep.2018.10.045 이거 archiving 해놓을것!! (Vascular Aging)
# https://github.com/aertslab/pySCENIC/issues/136 이것도 archiving ==> pySCENIC on bulk RNA-seq DATA!!!!!!!!!!!!!
# https://github.com/aertslab/pySCENIC/issues/169 이것도 archiving ==> multiple pySCENIC RUN?????? (Vascular Aging)
# https://github.com/aertslab/pySCENIC/find/master 여기에 .ipynb 들 
# UMAP 에 대한 영상인데, UMAP 만든 사람이 좋다고 함 https://www.youtube.com/watch?v=6BPl81wGGP8
# sc.pp.highly_variable_genes 의 batch_key 항목 꼭 check!!!
