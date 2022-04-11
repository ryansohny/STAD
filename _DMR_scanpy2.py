import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors
from statsmodels.stats.multitest import multipletests
sc.settings.verbosity = 3
plt.rcParams['figure.figsize'] = (5,5)
sns.set(font="Arial", font_scale=1, style='ticks')
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])
cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#0057b7", "#000000", "#ffd700"])
cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#0057B8", "#000000", "#ffd700"])
%matplotlib
%autoindent

clinic_info = pd.read_csv('/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/2022_WC300_clinical_information_Xadded.csv', index_col='Sample')

###### DMR table for ALL (tumor and normal)
dmr = pd.read_csv("DMR_abs15_ALL_imputed_corrected.csv", index_col=0).iloc[:,:-1] # Type column removal
dmr = sc.AnnData(dmr)
dmr.raw = dmr
dmr.layers['Percent_met'] = dmr.X
# np.ndarray.min(dmr.raw.X) ==> 전체 table에서 minimum value 값 (maximum은 min==> max)

dmr.obs = clinic_info.copy() # copy 반드시 집어넣어야함

##### 다른 clinical data 집어넣은다음에 umap leiden에 넣어서 봐볼것!


sc.pp.scale(dmr)
sc.tl.pca(dmr, n_comps=100, zero_center=True)
#sc.pl.pca(dmr, color='TN', add_outline=True, size=100, palette={'Normal':'Blue', 'Tumor':'Red'})
#sc.pl.pca_variance_ratio(dmr, log=True)
pca_variance = pd.DataFrame(dmr.uns['pca']['variance_ratio'], index=list(map(lambda x: 'PC' + str(x), list(range(1,101)))), columns=['Variance_ratio'])
np.sum(pca_variance.values.flatten()[:11])
# 0.7058488


sc.pp.neighbors(dmr, n_neighbors=13, n_pcs=11) # 13 ==> Good
sc.tl.leiden(dmr, resolution=0.5, key_added='leiden_r05')
sc.tl.leiden(dmr, resolution=1.0, key_added='leiden_r1')
sc.tl.umap(dmr, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
sc.pl.umap(dmr, color=['TN'], palette={'Normal':'Blue', 'Tumor':'Red'})
sc.pl.umap(dmr, color=['TN', 'leiden_r05', 'leiden_r1', 'PercentMet_COV5'], add_outline=False, size=400, color_map=cmap)

dmr_met = pd.DataFrame(dmr.raw.X.T, index=dmr.var.index, columns=dmr.obs.index)

col_colors1= list(dict(zip(list(dmr.obs['TN'].value_counts().index), ['#C0C0C0', '#000000']))[x] for x in dmr.obs['TN'])
col_colors2 = list(dict(zip(list(dmr.obs['leiden_r1'].value_counts().index), dmr.uns['leiden_r1_colors']))[x] for x in dmr.obs['leiden_r1'])
#sns.clustermap(dmr_met, method='ward', metric='euclidean', z_score=None, standard_scale=0, cmap=cmap, xticklabels=True)
sns.clustermap(dmr_met,
               method='ward',
               metric='euclidean',
               z_score=None,
               standard_scale=0,
               cmap=cmap,
               xticklabels=True,
               yticklabels=False,
               col_colors=[col_colors1, col_colors2])

#### Tumor - Normal Difference
dmr_met = pd.DataFrame(dmr.raw.X.T, index=dmr.var.index, columns=dmr.obs.index)

dmr_tn_met = dmr_met.iloc[:, 84:] - dmr_met.iloc[:,:84].values
col_colors3 = list(dict(zip(['0', '1', '2', '3', '4', '5'], ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']))[x] for x in dmr[84:].obs['leiden_r1'])

sns.clustermap(dmr_tn_met,
               method='ward',
               metric='euclidean',
               z_score=None,
               standard_scale=0,
               cmap=cmap,
               xticklabels=True,
               yticklabels=False,
               col_colors=[col_colors3])


############################################################################################################
############################################################################################################
############################################################################################################
#################### 2022-03-27 ####################

dmr_met = pd.read_csv("DMR_abs15_ALL_imputed_corrected.csv", index_col=0).iloc[:,:-1].T
dmr_met.index = list(map(lambda x: x.split('.')[0] + ':' + x.split('.')[1] + '-' + x.split('.')[2], dmr_met.index))

dmr_info = pd.read_table("DMR_abs15_Hyper-Hypo_annotation.txt", index_col=0)
dmr_info['Norm_CpGdensity'] = (np.max(dmr_info['CpGdensity']) - dmr_info['CpGdensity']) / np.max(dmr_info['CpGdensity'])

row_colors_dmr1 = list(dict(zip(['Hypo', 'Hyper'], ['#6C8EAD', '#A23E48']))[x] for x in dmr_info['Type'])


g = sns.clustermap(dmr_met,
                   method='ward',
                   metric='euclidean',
                   z_score=None,
                   standard_scale=0,
                   cmap=cmap,
                    robust=True,
                    col_colors=[col_colors1],
               row_colors=[row_colors_dmr1],
                   xticklabels=False,
                    yticklabels=False,
                     cbar_kws={'label': 'DNA methylation'})

g.cax.set_visible(False)


p = sns.boxplot(data=dmr_info, x='Type', y='Norm_CpGdensity', palette={'Hyper': '#A23E48', 'Hypo': '#6C8EAD'}, showfliers = False)
p = sns.stripplot(data=dmr_info, x='Type', y='Norm_CpGdensity', jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
p.set_ylabel("Normalized CpG density of DMR")
p.set_xlabel("Types of DMR")
p.set_xticklabels(['Hypomethylation', 'Hypermethylation'])
plt.tight_layout()

###### DMR table for tumor only
dmr_t = sc.AnnData(pd.read_csv("DMR_abs15_ALL_imputed_corrected.csv", index_col=0).iloc[84:,:-1])
dmr_t.var_names = list(map(lambda x: x.split('.')[0] + ':' + x.split('.')[1] + '-' + x.split('.')[2], dmr_t.var_names))
dmr_t.raw = dmr_t
dmr_t.layers['Percent_met'] = dmr_t.X
#clinic_info = pd.read_csv('/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/2022_WC300_clinical_information_Xadded.csv', index_col='Sample')
#clinic_info[['T.stage', 'N.stage', 'M.stage', 'TNM.stage', 'HER2IHC']] = clinic_info[['T.stage', 'N.stage', 'M.stage', 'TNM.stage', 'HER2IHC']].astype(str)
#clinic_info = clinic_info.iloc[:,np.r_[0:17, -1]] ?????????????????????????????????????????????????????????????#??? 이거 왜 했지
dmr_t.obs = clinic_info.iloc[84:].copy()
sc.pp.scale(dmr_t)
sc.tl.pca(dmr_t, n_comps=83, zero_center=True)
dmr_t.obs.rename(columns={"EpiBurden": "Epimutation Burden"}, inplace=True)
sc.pl.pca(dmr_t, color='Epimutation Burden', add_outline=False, size=200, color_map=cmap, annotate_var_explained=True)
# sc.pl.pca_variance_ratio(dmr_t, log=True)

#pca_variance_t = pd.DataFrame(dmr_t.uns['pca']['variance_ratio'], index=list(map(lambda x: 'PC' + str(x), list(range(1,84)))), columns=['Variance_ratio'])
# np.sum(pca_variance_t.values.flatten()[:12])

#sc.pp.neighbors(dmr_t, n_neighbors=15, n_pcs=12)
#sc.pp.neighbors(dmr_t, n_neighbors=11, n_pcs=12)
sc.pp.neighbors(dmr_t, n_neighbors=20, n_pcs=12)
sc.tl.leiden(dmr_t, resolution=1.0, key_added='leiden_r1')
sc.tl.leiden(dmr_t, resolution=0.75, key_added='leiden_r075')
sc.tl.leiden(dmr_t, resolution=0.5, key_added='leiden_r05')
sc.tl.leiden(dmr_t, resolution=0.3, key_added='leiden_r03')
sc.tl.umap(dmr_t, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
sc.pl.umap(dmr_t, color=['leiden_r1', 'leiden_r075', 'leiden_r03', 'EpiBurden'], add_outline=False, color_map=cmap)

leiden_name_change_dict1 = {'3': 'leiden_A',
                           '1': 'leiden_B',
                           '0': 'leiden_C',
                           '2': 'leiden_D'} # leiden_r075

leiden_name_change_dict2 = {'1': 'leiden_A',
                           '0': 'leiden_B',
                           '2': 'leiden_C'}

dmr_t.obs['DMR_leiden'] = dmr_t.obs['leiden_r075'].map(lambda x: leiden_name_change_dict1[x]).astype('category')
#dmr_t.obs['DMR_leiden'] = dmr_t.obs['leiden_r05'].map(lambda x: leiden_name_change_dict2[x]).astype('category')
sc.pl.umap(dmr_t, color='DMR_leiden', add_outline=False, legend_loc='on data', palette='tab10')
dmr_t.obs.rename(columns={"DMR_leiden": "DMR Clusters"}, inplace=True)

sc.pl.umap(dmr_t, color='DMR Clusters', add_outline=False, legend_loc='on data', palette='Dark2')
dmr_t.obs.rename(columns={"DMR Clusters": "DMR_leiden"}, inplace=True)

lin2 = tuple(sorted(list(dmr_t.obs['DMR_leiden'].values.unique())))
dmr_t.obs['DMR_leiden'] = dmr_t.obs['DMR_leiden'].cat.reorder_categories(list(lin2), ordered=True)
color_dict2 = {
    "leiden_A": "#1b9e77",
    "leiden_B": "#7570b3",
    "leiden_C": "#e6ab02",
    "leiden_D": "#666666"} # DMR_leiden (leiden_r075)

p = sns.boxplot(data=dmr_t.obs, x='DMR_leiden', y='EpiBurden', palette=color_dict2)
p = sns.stripplot(data=dmr_t.obs, x='DMR_leiden', y='EpiBurden', jitter=True, marker='o', color='black', alpha=0.8)
p.set_ylabel('DMR Clusters')
p.set_xlabel('DMR Clusters')
p.set_ylabel('Epimutation Burden')
plt.tight_layout()
sns.despine()

import statsmodels.api as sm
from statsmodels.formula.api import ols
df = dmr_t.obs[['EpiBurden', 'DMR_leiden']]
df_lm = ols('EpiBurden ~ C(DMR_leiden)', data=df).fit() # C() ==> categorical data (not necessary here because batch is already categorical)
print(sm.stats.anova_lm(df_lm, typ=2))

from statsmodels.stats.multicomp import pairwise_tukeyhsd
print(pairwise_tukeyhsd(df['EpiBurden'], df['DMR_leiden'], alpha=0.05))

# Diffusion pseudotime
sc.tl.diffmap(dmr_t)
sc.pl.diffmap(dmr_t, color=['DMR_leiden'], add_outline=False, color_map=cmap)

start_cell = np.isin(dmr_t.obs['DMR_leiden'], 'leiden_A') # leiden_r1의 3번 cluster
max_start_id = np.argmax(dmr_t.obsm['X_diffmap'][start_cell, 1])
root_id = np.arange(len(start_cell))[start_cell][max_start_id]
dmr_t.uns['iroot'] = root_id
sc.tl.dpt(dmr_t, n_branchings=0, n_dcs=15)

p = sns.lmplot(data=dmr_t.obs, x='dpt_pseudotime', y='EpiBurden')
p.set_xlabels('Diffusion Pseudotime')
p.set_ylabels('Epimutation Burden')

stats.spearmanr(dmr_t.obs['dpt_pseudotime'], dmr_t.obs['EpiBurden'])

lin = tuple(sorted(list(dmr_t.obs['DMR_leiden'].values.unique())))
dmr_t.obs['DMR_leiden'] = dmr_t.obs['DMR_leiden'].cat.reorder_categories(list(lin), ordered=True)

# DMR_leiden clusters and WHO classification
ax = pd.crosstab(dmr_t.obs['WHO'], dmr_t.obs['DMR_leiden'], normalize=1).T.plot.bar(stacked=True, color=dict(zip(dmr_t.obs['WHO'].unique(), sns.color_palette("Set2", len(dmr_t.obs['WHO'].unique())))))
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))
plt.tight_layout()

# DMR_leiden clusters and Lauren's classification
ax = pd.crosstab(dmr_t.obs['Lauren'], dmr_t.obs['DMR_leiden'], normalize=1).T.plot.bar(stacked=True, color=dict(zip(dmr_t.obs['Lauren'].unique(), sns.color_palette("Accent", len(dmr_t.obs['Lauren'].unique())))))
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))
plt.tight_layout()

# DMR_leiden clusters and T.stage
ax = pd.crosstab(dmr_t.obs['T.stage'], dmr_t.obs['DMR_leiden'], normalize=1).T.plot.bar(stacked=True, color=dict(zip(dmr_t.obs['T.stage'].unique(), sns.color_palette("tab10", len(dmr_t.obs['T.stage'].unique())))))
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))
plt.tight_layout()

# DMR_leiden clusters and N.stage
ax = pd.crosstab(dmr_t.obs['N.stage'], dmr_t.obs['DMR_leiden'], normalize=1).T.plot.bar(stacked=True, color=dict(zip(dmr_t.obs['N.stage'].unique(), sns.color_palette("tab10", len(dmr_t.obs['N.stage'].unique())))))
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))
plt.tight_layout()

# DMR_leiden clusters and M.stage
ax = pd.crosstab(dmr_t.obs['M.stage'], dmr_t.obs['DMR_leiden'], normalize=1).T.plot.bar(stacked=True, color=dict(zip(dmr_t.obs['M.stage'].unique(), sns.color_palette("tab10", len(dmr_t.obs['M.stage'].unique())))))
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))
plt.tight_layout()


# dmr_t_met clustering
dmr_t_met = pd.DataFrame(dmr_t.raw.X.T, index=dmr_t.var.index, columns=dmr_t.obs.index)

col_colors_DMR_leiden = list(dict(zip(sorted(list(dmr_t.obs['DMR_leiden'].value_counts().index)), \
                                      dmr_t.uns['DMR_leiden_colors']))[x] for x in dmr_t.obs['DMR_leiden'])
#col_colors2 = list(dict(zip(list(dmr.obs['TN'].value_counts().index), ['#C0C0C0', '#000000']))[x] for x in dmr.obs['TN']
#sns.clustermap(dmr_met, method='ward', metric='euclidean', z_score=None, standard_scale=0, cmap=cmap, xticklabels=True)

g = sns.clustermap(dmr_t_met,
                   method='ward',
                   metric='euclidean',
                   z_score=None,
                   standard_scale=0,
                   cmap=cmap,
                   xticklabels=False,
                   yticklabels=False,
                   col_colors=[col_colors_DMR_leiden],
                   row_colors=[row_colors_dmr1],
                   cbar_kws={'label': 'DNA methylation'})

g.ax_row_dendrogram.set_visible(False)
g.cax.set_visible(False)

dmr_t_met_colorder = g.dendrogram_col.reordered_ind
dmr_t_met_roworder = g.dendrogram_row.reordered_ind

dmr_n_met = dmr_met.iloc[:,:84]
sns.clustermap(dmr_n_met.iloc[dmr_t_met_roworder, dmr_t_met_colorder],
               col_cluster=False,
               row_cluster=False,
                   z_score=None,
                   standard_scale=None,
                   cmap=cmap,
                   xticklabels=False,
                   yticklabels=False)

##### CIMP-H and CIMP-L
cpgi_met = pd.read_table("CpGi_ALL.txt", index_col=0)
cpgi_met.columns = list(map(lambda x: 'X'+x, cpgi_met.columns))
cpgi_met.index = list(map(lambda x: '/'.join(x.split('/')[:2]), cpgi_met.index))
# Normal DNA methylation  < 40%
cpgi_met = cpgi_met[cpgi_met.isna().sum(axis=1) == 0][cpgi_met[cpgi_met.isna().sum(axis=1) == 0].iloc[:, :84].mean(axis=1) < 40]

# Promoter CpGi
cpgi_pls = pd.read_table("PLS_CpGi.txt", index_col=0)
cpgi_pls = cpgi_pls[cpgi_pls.index.isin(list(map(lambda x: '/'.join(x.split('/')[:2]), cpgi_met.index)))].index
cpgi_met = cpgi_met[cpgi_met.index.isin(cpgi_pls)]

# mean(Tumor - Normal) >= 10%
cpgi_tn_met = cpgi_met.iloc[:, 84:] - cpgi_met.iloc[:, :84].values
cpgi_tn_met = cpgi_tn_met[cpgi_tn_met.mean(axis=1) >= 10]

g = sns.clustermap(
    cpgi_tn_met,
    method='ward',
    metric='euclidean',
    z_score=None,
    standard_scale=None,
    cmap=cmap,
    xticklabels=False,
    yticklabels=False,
    col_colors=[col_colors_DMR_leiden])

dmr_t.obs['CIMP'] = list(map(lambda x: 'CIMP(+)' if x == True else 'CIMP(-)', dmr_t.obs.index.isin(cpgi_tn_met.iloc[:,g.dendrogram_col.reordered_ind[:37]].columns)))
dmr_t.obs[['DMR Clusters', 'CIMP']].value_counts().sort_index()
#DMR Clusters  CIMP
#leiden_A      CIMP(-)    10
#leiden_B      CIMP(+)     2
#              CIMP(-)    20
#leiden_C      CIMP(+)    22
#              CIMP(-)    13
#leiden_D      CIMP(+)    13
#              CIMP(-)     4
#dtype: int64
sc.pl.umap(dmr_t, color=['CIMP'], add_outline=False, legend_loc='right margin', palette={'CIMP(+)': 'darkred', 'CIMP(-)': 'navy'})
sc.pl.umap(dmr_t, color=['DMR Clusters', 'CIMP'], add_outline=False, legend_loc='right margin')
col_colors_CIMP = list(dict(zip(sorted(['CIMP(-)', 'CIMP(+)']), \
                                      ['#8b0000ff', '#000080ff']))[x] for x in dmr_t.obs['CIMP'])

g = sns.clustermap(
    cpgi_tn_met,
    method='ward',
    metric='euclidean',
    z_score=None,
    standard_scale=None,
    cmap=cmap,
    xticklabels=False,
    yticklabels=False,
    col_colors=[col_colors_DMR_leiden, col_colors_CIMP])


ax = (pd.crosstab(dmr_t.obs['DMR Clusters'], dmr_t.obs['CIMP'], normalize=0)*100).plot.bar(stacked=True, color=['#8b0000ff', '#000080ff'], rot=0)
plt.ylabel("Proportion (%)")
ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1.0))
plt.tight_layout()

g = sns.clustermap(dmr_t_met,
                   method='ward',
                   metric='euclidean',
                   z_score=None,
                   standard_scale=0,
                   cmap=cmap,
                   xticklabels=False,
                   yticklabels=False,
                   col_colors=[col_colors_DMR_leiden, col_colors_CIMP],
                   row_colors=[row_colors_dmr1],
                   cbar_kws={'label': 'DNA methylation'})
g.ax_row_dendrogram.set_visible(False)
g.cax.set_visible(False)


###### RNA-seq
rna = pd.read_table("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/STAD_SNUH_vst_new.txt", index_col=0, sep=' ')
rna.columns = list(map(lambda x: 'X'+x, rna.columns))
rna = rna[rna.index.isin(list(map(lambda x: x.split('/')[2], pls.index)))]

rna_tonly = pd.read_table("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/STAD_SNUH_vst_new_Tumor.txt", index_col=0, sep=' ')
rna_tonly.columns = list(map(lambda x: 'X'+x, rna_tonly.columns))
#rna_tonly = rna_tonly[rna_tonly.index.isin(list(map(lambda x: x.split('/')[2], pls.index)))]

deg_tn = pd.read_csv("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/Tumor.Normal.compare.csv", index_col=0)
deg_tn = deg_tn[deg_tn.index.isin(list(map(lambda x: x.split('/')[2], pls.index)))]

deg_tn_uplist = deg_tn[(deg_tn['baseMean'] >=10) & (deg_tn['padj'] < 0.05) & (deg_tn['log2FoldChange'] > 1)].index
deg_tn_downlist = deg_tn[(deg_tn['baseMean'] >=10) & (deg_tn['padj'] < 0.05) & (deg_tn['log2FoldChange'] < -1)].index

# Gene division by expression
np.log10(deg_tn['baseMean']).describe()
p = sns.histplot(data=pd.DataFrame(np.log10(deg_tn['baseMean'])), x='baseMean', kde=True, stat='count', bins=100)
p.set_xticks([0,1,2,3,4,5,6])
p.set_xlabel("Log10 mean normalized read counts across samples")
sns.despine()




deg_da = pd.read_csv("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/Leiden_D.A.compare.csv", index_col=0)
deg_da_uplist = deg_da[(deg_da['baseMean'] >=10) & (deg_da['padj'] < 0.05) & (deg_da['log2FoldChange'] > 1)].index
deg_da_downlist = deg_da[(deg_da['baseMean'] >=10) & (deg_da['padj'] < 0.05) & (deg_da['log2FoldChange'] < -1)].index

# Tumor vs Normal DEG (log2FC >=1 & baseMean >=10, padj < 0.05) ==> complete linkage (and correlation)
sns.clustermap(
    rna[rna.index.isin(list(deg_tn_uplist) + list(deg_tn_downlist))],
    method='complete',
    metric='correlation',
    z_score=0,
    standard_scale=None,
    cmap=cmap3,
    col_colors=[col_colors1],
    xticklabels=False,
    yticklabels=False, vmin=-2.5, vmax=2.5)

sns.clustermap(
    rna[rna.index.isin(list(deg_tn_uplist) + list(deg_tn_downlist))].iloc[:, 84:],
    method='complete',
    metric='correlation',
    z_score=0,
    standard_scale=None,
    cmap=cmap3,
    col_colors=[col_colors_DMR_leiden, col_colors_CIMP],
    xticklabels=False,
    yticklabels=False, vmin=-2.5, vmax=2.5)

sns.clustermap(
    rna_tonly[rna_tonly.index.isin(list(deg_da_uplist) + list(deg_da_downlist))],
    method='ward',
    metric='euclidean',
    z_score=0,
    standard_scale=None,
    cmap=cmap3,
    col_colors=[col_colors_DMR_leiden, col_colors_CIMP],
    xticklabels=False,
    yticklabels=False, vmin=-2.5, vmax=2.5)

###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### Promoter Methylation (PLS)
pls = pd.read_table("PLS_ALL.txt", index_col=0)
pls_tn = pls.iloc[:,84:] - pls.iloc[:,:84].values
pls_tn[pls_tn.isna().sum(axis=1) <= 21] # 21 / 84 == 25%
pls[pls.index.isin(pls_tn[pls_tn.isna().sum(axis=1) <= 21].index)].to_csv("PLS_ALL_isna25percent.txt", sep="\t")


pls = pd.read_table("PLS_ALL.txt", index_col=0)
pls_index = pls[pls.isna().sum(axis=1) <= 10].index
pls = pd.read_csv("PLS_ALL_isna10_imputed_corrected.csv", index_col=0).iloc[:,:-1].T
pls.loc['chrX:155997348-155997895/EH38D6141605/IL9R/ENSG00000124334/ENST00000369423/protein_coding/6'] = pls.loc['chrX.155997348.155997895.EH38D6141605.IL9R.ENSG00000124334.ENST00000244174.protein_coding.6'] # imputation할 때 빠뜨렸었음
pls.index = pls_index

pls_info = pd.read_table("PLS_annotated_table_full.txt", index_col=0)

pls = pls.loc[pls_info[pls_info.index.isin(pls.index)].index]
pls_info = pls_info[pls_info.index.isin(pls.index)]

pls_info = pls_info[pls_info['GeneID'].isin(rna.index)]
pls = pls[pls.index.isin(pls_info.index)]
pls = pls[pls.index.isin(pls_info.index)]

pls_deg_da_up[pls_deg_da_up.isna().sum(axis=1) == 0].iloc[:, 84:]
sns.clustermap(pls_deg_da_up[pls_deg_da_up.isna().sum(axis=1) == 0].iloc[:,84:],
                   method='ward',
                   metric='euclidean',
                   z_score=0,
                   standard_scale=None,
                   cmap=cmap,
                    robust=True,
                    col_colors=[col_colors_DMR_leiden],
               row_colors=None,
                   xticklabels=False,
                    yticklabels=False, vmin=-1, vmax=1)

# No expression (log10 baseMean < 1)
np.log10(deg_tn['baseMean']).describe(percentiles=[.20, .25, .50, .75, .90])
#count    20762.000000
#mean         1.616227
#std          1.220618
#min         -1.526428
#20%          0.447281
#25%          0.688768
#50%          1.906520 ==> low (below), intermediate (upper)
#75%          2.608915
#90%          2.978760 ==> high expression
#max          6.937551

percentile_50 = np.log10(deg_tn['baseMean']).quantile([.5]).values[0]
percentile_90 = np.log10(deg_tn['baseMean']).quantile([.9]).values[0]

a = pd.DataFrame(pls_info[pls_info['GeneID'].isin(deg_tn[np.log10(deg_tn['baseMean']) < 1].index)].iloc[:,-8:]['K36me3'].value_counts())

b = pd.DataFrame(pls_info[pls_info['GeneID'].isin(deg_tn[(np.log10(deg_tn['baseMean']) >= 1) & (np.log10(deg_tn['baseMean']) < percentile_50)].index)].iloc[:-8:]['K36me3'].value_counts())

c = pd.DataFrame(pls_info[pls_info['GeneID'].isin(deg_tn[(np.log10(deg_tn['baseMean']) >= percentile_50) & (np.log10(deg_tn['baseMean']) < percentile_90)].index)].iloc[:-8:]['K36me3'].value_counts())

d = pd.DataFrame(pls_info[pls_info['GeneID'].isin(deg_tn[np.log10(deg_tn['baseMean']) >= percentile_90].index)].iloc[:,-8:]['K36me3'].value_counts())

abcd = pd.concat([a,b,c,d], axis=1)
abcd.columns = ['No expression', 'Low expression', 'Intermediate expression', 'High expression']

ax = abcd.plot.pie(subplots=True, colors = ['lightgrey', 'aqua'], figsize=(15,8), legend=False, counterclock=False, fontsize=13, autopct='%0.01f%%')

######### Total Promoter #########
# Total Promoter DNA methylation
a = pd.DataFrame(pls[pls.isna().sum(axis=1) == 0].iloc[:,:84].mean(axis=1), columns=['Normal'])
b = pd.DataFrame(pls[pls.isna().sum(axis=1) == 0].iloc[:,84:].mean(axis=1), columns=['Tumor'])
ab = pd.concat([a,b], axis=1)
#p = sns.violinplot(data=ab, palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale='width')
p = sns.boxplot(data=ab, palette={'Normal':'navy', 'Tumor':'darkred'}, showfliers = False)
p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
p.set_ylabel("Promoter DNA methylation (%)")
p.set_title("All Promoters")
plt.tight_layout()
stats.ttest_rel(ab['Tumor'], ab['Normal'])

######### CpGi 유무 ##########
# CpGi None (Promoter DNA methylation)
a = pd.DataFrame(pls[pls.index.isin(pls_info[pls_info['CpGi'] == 'na'].index)].iloc[:, :84].mean(axis=1), columns=['Normal'])
b = pd.DataFrame(pls[pls.index.isin(pls_info[pls_info['CpGi'] == 'na'].index)].iloc[:, 84:].mean(axis=1), columns=['Tumor'])
ab = pd.concat([a,b], axis=1)
#p = sns.violinplot(data=ab, palette={'Normal':'navy', 'Tumor':'darkred'}, scale='width')
p = sns.boxplot(data=ab, palette={'Normal':'navy', 'Tumor':'darkred'}, showfliers = False)
p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
p.set_ylabel("Promoter DNA methylation (%)")
p.set_title("Non-CpGi Promoter (N=17,949)")
plt.tight_layout()
stats.ttest_rel(ab['Tumor'], ab['Normal'])

# CpGi None (RNA expression)
a = pd.DataFrame(rna[rna.index.isin(list(set(list(map(lambda x: x.split('/')[2], pls_info[pls_info['CpGi'] == 'na'].index)))))].iloc[:, :84].mean(axis=1), columns=['Normal'])
b = pd.DataFrame(rna[rna.index.isin(list(set(list(map(lambda x: x.split('/')[2], pls_info[pls_info['CpGi'] == 'na'].index)))))].iloc[:, 84:].mean(axis=1), columns=['Tumor'])
ab = pd.concat([a,b], axis=1)
#p = sns.violinplot(data=ab, palette={'Normal':'navy', 'Tumor':'darkred'}, scale='width')
p = sns.boxplot(data=ab, palette={'Normal':'navy', 'Tumor':'darkred'}, showfliers = False)
p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
p.set_ylabel("RNA expression")
p.set_title("Non-CpGi Promoter (N=10,940)")
plt.tight_layout()
stats.ttest_rel(ab['Tumor'], ab['Normal'])

# CpGi Yes (Promoter DNA methylation)
a = pd.DataFrame(pls[pls.index.isin(pls_info[pls_info['CpGi'] == 'Yes'].index)].iloc[:, :84].mean(axis=1), columns=['Normal'])
b = pd.DataFrame(pls[pls.index.isin(pls_info[pls_info['CpGi'] == 'Yes'].index)].iloc[:, 84:].mean(axis=1), columns=['Tumor'])
ab = pd.concat([a,b], axis=1)
#p = sns.violinplot(data=ab, palette={'Normal':'navy', 'Tumor':'darkred'}, scale='width')
p = sns.boxplot(data=ab, palette={'Normal':'navy', 'Tumor':'darkred'}, showfliers = False)
p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
p.set_ylabel("Promoter DNA methylation (%)")
p.set_title("CpGi Promoter (N=28,633)")
plt.tight_layout()
stats.ttest_rel(ab['Tumor'], ab['Normal'])

# CpGi Yes (RNA expression)
a = pd.DataFrame(rna[rna.index.isin(list(set(list(map(lambda x: x.split('/')[2], pls_info[pls_info['CpGi'] == 'Yes'].index)))))].iloc[:, :84].mean(axis=1), columns=['Normal'])
b = pd.DataFrame(rna[rna.index.isin(list(set(list(map(lambda x: x.split('/')[2], pls_info[pls_info['CpGi'] == 'Yes'].index)))))].iloc[:, 84:].mean(axis=1), columns=['Tumor'])
ab = pd.concat([a,b], axis=1)
#p = sns.violinplot(data=ab, palette={'Normal':'navy', 'Tumor':'darkred'}, scale='width')
p = sns.boxplot(data=ab, palette={'Normal':'navy', 'Tumor':'darkred'}, showfliers = False)
p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
p.set_ylabel("RNA expression")
p.set_title("CpGi Promoter (N=16,019)")
plt.tight_layout()
stats.ttest_rel(ab['Tumor'], ab['Normal'])

# Promoter methylation boxplot
color_dict3 = {
    "Normal": "#FFFAFA",
    "leiden_A": "#1b9e77",
    "leiden_B": "#7570b3",
    "leiden_C": "#e6ab02",
    "leiden_D": "#666666"}

gene = 'TWIST1'
df = pd.concat([ pls.loc[pls_info[pls_info['GeneID'] == gene].index].iloc[0, 84:], pd.DataFrame(dmr_t.obs['DMR Clusters']) ], axis=1)
sns.boxplot(data=df, x='DMR Clusters', y=df.columns[0], palette=color_dict2)
sns.swarmplot(data=df, x='DMR Clusters', y=df.columns[0], color=".2").set(ylabel=df.columns[0].split('/')[2] + ' promoter methylation')

# WIth normal Promoter methylation
df = pd.concat([pls.loc[pls_info[pls_info['GeneID'] == gene].index].iloc[0, 84:], pd.DataFrame(dmr_t.obs['DMR Clusters'])], axis=1)
df2 = pd.concat([pls.loc[pls_info[pls_info['GeneID'] == gene].index].iloc[0, :84], pd.DataFrame(["Normal"]*84, columns=["DMR Clusters"], index=pls.iloc[:,:84].columns)], axis=1)
df = pd.concat([df, df2])
sns.boxplot(data=df, x='DMR Clusters', y=df.columns[0], palette=color_dict3, order=["Normal", "leiden_A", "leiden_B", "leiden_C", "leiden_D"], showfliers=False)
sns.stripplot(data=df, x='DMR Clusters', y=df.columns[0], order=["Normal", "leiden_A", "leiden_B", "leiden_C", "leiden_D"], color=".2", size=2, alpha=0.4).set(ylabel=df.columns[0].split('/')[2] + ' promoter methylation')
stats.mannwhitneyu(df[df['DMR Clusters'] == 'Normal'].iloc[:,0], df[df['DMR Clusters'] == 'leiden_D'].iloc[:,0])

# RNA expression boxplot
df = pd.concat([pd.DataFrame(rna.iloc[:,84:].loc[gene]), pd.DataFrame(dmr_t.obs['DMR Clusters'])], axis=1)
p = sns.boxplot(data=pd.concat([pd.DataFrame(rna.iloc[:,84:].loc[gene]), pd.DataFrame(dmr_t.obs['DMR Clusters'])], axis=1), x='DMR Clusters', y=gene, palette=color_dict2, showfliers=False)
p = sns.stripplot(data=pd.concat([pd.DataFrame(rna.iloc[:,84:].loc[gene]), pd.DataFrame(dmr_t.obs['DMR Clusters'])], axis=1), x='DMR Clusters', y=gene, color=".2").set(ylabel=gene + ' expression')

df = pd.concat([pd.DataFrame(rna_tonly.loc[gene]), pd.DataFrame(dmr_t.obs['DMR Clusters'])], axis=1)
p = sns.boxplot(data=df, x='DMR Clusters', y=gene, palette=color_dict2, showfliers=False)
p = sns.stripplot(data=df, x='DMR Clusters', y=df.columns[0], color=".2", size=2, alpha=0.4).set(ylabel=df.columns[0] + ' expression')

# WIth normal RNA expression boxplot
df = pd.concat([pd.DataFrame(rna.iloc[:,84:].loc[gene]), pd.DataFrame(dmr_t.obs['DMR Clusters'])], axis=1)
df2 = pd.concat([pd.DataFrame(rna.iloc[:,:84].loc[gene]), pd.DataFrame(["Normal"]*84, columns=["DMR Clusters"], index=rna.iloc[:,:84].columns)], axis=1)
df = pd.concat([df, df2])
p = sns.boxplot(data=df, x='DMR Clusters', y=gene, palette=color_dict3, order=["Normal", "leiden_A", "leiden_B", "leiden_C", "leiden_D"], width=0.6, showfliers=False)
p = sns.stripplot(data=df, x='DMR Clusters', y=gene, order=["Normal", "leiden_A", "leiden_B", "leiden_C", "leiden_D"], color=".2", size=2, alpha=0.4).set(ylabel=gene + ' expression')
stats.mannwhitneyu(df[df['DMR Clusters'] == 'Normal'].iloc[:,0], df[df['DMR Clusters'] == 'leiden_D'].iloc[:,0])

stats.pearsonr(rna.loc['ESRRG'].iloc[84:], dmr_t.obs['dpt_pseudotime'])
stats.pearsonr(pls.loc[pls_info[pls_info['GeneID'] == 'ESRRG'].index].iloc[8, 84:], dmr_t.obs['dpt_pseudotime'])



for i in deg_up:
    try:
        pearson = stats.pearsonr(np.delete(pro_met.loc[pro_met_info[pro_met_info['GeneID'] == i].index].iloc[:,84:].values[0], np.argwhere(np.isnan(pro_met.loc[pro_met_info[pro_met_info['GeneID'] == i].index].iloc[:,84:].values[0]))), np.delete(deg_tn_protein.loc[i].values, np.argwhere(np.isnan(pro_met.loc[pro_met_info[pro_met_info['GeneID'] == i].index].iloc[:,84:].values[0]))))
        if pearson[1] < 0.05: # p value만 가지고 해야 함 나중에는
            if pearson[0] < 0:
                deg_up_sig_cor.append(i)
            print(i + '\t' + str(pearson[0]) + '\t' + str(pearson[1]))
    except (IndexError, ValueError):
        pass

for i in deg_down:
    try:
        pearson = stats.pearsonr(np.delete(pro_met.loc[pro_met_info[pro_met_info['GeneID'] == i].index].iloc[:,84:].values[0], np.argwhere(np.isnan(pro_met.loc[pro_met_info[pro_met_info['GeneID'] == i].index].iloc[:,84:].values[0]))), np.delete(deg_tn_protein.loc[i].values, np.argwhere(np.isnan(pro_met.loc[pro_met_info[pro_met_info['GeneID'] == i].index].iloc[:,84:].values[0]))))
        if pearson[1] < 0.05:
            if pearson[0] < 0: # p value만 가지고 해야 함 나중에는
                deg_down_sig_cor.append(i)
            print(i + '\t' + str(pearson[0]) + '\t' + str(pearson[1]))
    except (IndexError, ValueError):
        pass



#### Genebody Methylation ####
gb = pd.read_csv("GeneBody_ALL_isna10_imputed_corrected.csv", index_col=0).iloc[:,:-1].T
def R_FormattedIndex_correction(string):
    return string[:string.index('ENSG')].split('.')[0]  + ':' + string[:string.index('ENSG')].split('.')[1] + '-' + string[:string.index('ENSG')].split('.')[2] + '/' + '.'.join(string[:string.index('ENSG')].split('.')[3:-1]) + '/' + '/'.join(string[string.index('ENSG'):].split('.'))

gb.index = list(map(lambda x: R_FormattedIndex_correction(x), gb.index))


pearsonr_gb_pseudotime = pd.DataFrame(gb.iloc[:, 84:].corrwith(dmr_t.obs['dpt_pseudotime'], axis=1, method='pearson'), columns=['PearsonR'])

# WIth normal Genebody methylation
gene = 'HOXA9'
df = pd.concat([gb[gb.index.isin(gb_info[gb_info['GeneID'] == 'HOXA9'].index)].iloc[1, 84:], pd.DataFrame(dmr_t.obs['DMR Clusters'])], axis=1)
df2 = pd.concat([gb[gb.index.isin(gb_info[gb_info['GeneID'] == 'HOXA9'].index)].iloc[1, :84], pd.DataFrame(["Normal"]*84, columns=["DMR Clusters"], index=pls.iloc[:,:84].columns)], axis=1)
df = pd.concat([df, df2])
p = sns.boxplot(data=df, x='DMR Clusters', y=df.columns[0], palette=color_dict3, order=["Normal", "leiden_A", "leiden_B", "leiden_C", "leiden_D"], showfliers=False)
p = sns.stripplot(data=df, x='DMR Clusters', y=df.columns[0], order=["Normal", "leiden_A", "leiden_B", "leiden_C", "leiden_D"], color=".2", size=2, alpha=0.4).set(ylabel=df.columns[0].split('/')[1] + ' genebody methylation')
plt.tight_layout()

stats.mannwhitneyu(df[df['DMR Clusters'] == 'Normal'].iloc[:,0], df[df['DMR Clusters'] == 'leiden_D'].iloc[:,0])

# WIth normal RNA expression boxplot
df = pd.concat([pd.DataFrame(rna.iloc[:,84:].loc[gene]), pd.DataFrame(dmr_t.obs['DMR Clusters'])], axis=1)
df2 = pd.concat([pd.DataFrame(rna.iloc[:,:84].loc[gene]), pd.DataFrame(["Normal"]*84, columns=["DMR Clusters"], index=rna.iloc[:,:84].columns)], axis=1)
df = pd.concat([df, df2])
sns.boxplot(data=df, x='DMR Clusters', y=gene, palette=color_dict3, order=["Normal", "leiden_A", "leiden_B", "leiden_C", "leiden_D"], width=0.6, showfliers=False)
sns.stripplot(data=df, x='DMR Clusters', y=gene, order=["Normal", "leiden_A", "leiden_B", "leiden_C", "leiden_D"], color=".2", size=2, alpha=0.4).set(ylabel=gene + ' expression')
stats.mannwhitneyu(df[df['DMR Clusters'] == 'Normal'].iloc[:,0], df[df['DMR Clusters'] == 'leiden_D'].iloc[:,0])





























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
