# source activate complexheatmap
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import scanpy as sc
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
sc.settings.verbosity = 3
plt.rcParams['figure.figsize'] = (5,5)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])
cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#0057b7", "#000000", "#ffd700"])
cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#0057B8", "#000000", "#ffd700"])

%matplotlib
%autoindent


# Step 1. Create Tumor Normal imputed matrix
import pandas as pd
dmr = pd.read_csv("DMR_mat_imputed_corrected.csv", index_col=0)
dmr.iloc[84:, :-1].to_csv("DMR_mat_imputed_corrected_tumor.csv")
dmr.iloc[:84, :-1].to_csv("DMR_mat_imputed_corrected_normal.csv")

# Step 2-1.
dmr_t = sc.AnnData(pd.read_csv("DMR_mat_imputed_corrected_tumor.csv", index_col=0))
dmr_t.raw = dmr_t
dmr_t.layers['Percent_met'] = dmr_t.X
clinic_info = pd.read_csv('/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/2022_WC300_clinical_information_Xadded.csv', index_col=0)
clinic_info[['T.stage', 'N.stage', 'M.stage', 'TNM.stage', 'HER2IHC']] = clinic_info[['T.stage', 'N.stage', 'M.stage', 'TNM.stage', 'HER2IHC']].astype(str)
#clinic_info = clinic_info.iloc[:,np.r_[0:17, -1]]
dmr_t.obs = clinic_info
sc.pp.scale(dmr_t)
sc.tl.pca(dmr_t, n_comps=83, zero_center=True)
#sc.pl.pca(dmr_t, color='Lauren')
# sc.pl.pca_variance_ratio(dmr_t, log=True)

pca_variance_t = pd.DataFrame(dmr_t.uns['pca']['variance_ratio'], index=list(map(lambda x: 'PC' + str(x), list(range(1,84)))), columns=['Variance_ratio'])
np.sum(pca_variance_t.values.flatten()[:13])
# 0.7076739
#sc.pp.neighbors(dmr_t, n_neighbors=11, n_pcs=13) ###
#sc.pp.neighbors(dmr_t, n_neighbors=20, n_pcs=13)
sc.pp.neighbors(dmr_t, n_neighbors=11, n_pcs=13)
sc.tl.leiden(dmr_t, resolution=1.0, key_added='leiden_r1')
sc.tl.leiden(dmr_t, resolution=0.5, key_added='leiden_r05')
sc.tl.umap(dmr_t, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
sc.pl.umap(dmr_t, color=['leiden_r1', 'leiden_r05', 'EpiBurden'], add_outline=False, color_map=cmap)

leiden_name_change_dict1 = {'3': 'leiden_A',
                           '1': 'leiden_B',
                           '0': 'leiden_C',
                           '2': 'leiden_D'}
dmr_t.obs['DMR_leiden'] = dmr_t.obs['leiden_r05'].map(lambda x: leiden_name_change_dict1[x]).astype('category')

sc.tl.diffmap(dmr_t)
sc.pl.diffmap(dmr_t, color=['DMR_leiden'], add_outline=False, color_map=cmap)

start_cell = np.isin(dmr_t.obs['DMR_leiden'], 'leiden_A') # leiden_r1의 3번 cluster
max_start_id = np.argmax(dmr_t.obsm['X_diffmap'][start_cell, 1])
root_id = np.arange(len(start_cell))[start_cell][max_start_id]
dmr_t.uns['iroot'] = root_id
sc.tl.dpt(dmr_t, n_branchings=0, n_dcs=15)

sns.lmplot(data=dmr_t.obs, x='dpt_pseudotime', y='EpiBurden')

lin = tuple(sorted(list(dmr_t.obs['DMR_leiden'].values.unique())))
dmr_t.obs['DMR_leiden'] = dmr_t.obs['DMR_leiden'].cat.reorder_categories(list(lin), ordered=True)

color_dict = {
    "leiden_A": "#d62728",
    "leiden_B": "#ff7f0e",
    "leiden_C": "#1f77b4",
    "leiden_D": "#2ca02c"
} # equivalent to dict(zip(list(dmr_t.obs['DMR_leiden'].value_counts().index), dmr_t.uns['DMR_leiden_colors']))

sns.boxplot(data=dmr_t.obs, x='DMR_leiden', y='EpiBurden', palette=color_dict)
sns.swarmplot(data=dmr_t.obs, x='DMR_leiden', y='EpiBurden', color=".2")

rna = pd.read_table("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/STAD_SNUH_vst.txt", index_col=0, sep=' ')
rna = pd.read_table("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/STAD_SNUH_vst.txt", index_col=0, sep=' ')
rna.columns = list(map(lambda x: "X" + x, rna.columns))

deg_tn_protein = pd.read_table("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/STAD_SNUH_Tumor_leiden_vst_DEG_Leiden_A_D_protein.txt", index_col="ID")
deg_tn_protein.columns = list(map(lambda x: "X" + x, deg_tn_protein.columns))



pro_met = pd.read_table("Promoter_up500down500_ALL.txt", index_col="ID")
pro_met = pd.read_table("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/03.WGBS/Promoter_cCRE_ALL.txt", index_col="ID")
pro_met = pd.read_table("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/03.WGBS/NEW/Promoter_cCRE_ALL.txt", index_col="ID")
pro_met.columns = list(map(lambda x: "X" + x, pro_met.columns))
pro_met_info = pd.DataFrame(list(zip(list(map(lambda x: x.split('/')[0], pro_met.index)), list(map(lambda x: x.split('/')[1], pro_met.index)), list(map(lambda x: x.split('/')[2], pro_met.index)), list(map(lambda x: x.split('/')[-1], pro_met.index)))), columns=['Loc', 'GeneID', 'EnsemblID', 'CpG'], index=pro_met.index)
pro_met.loc[pro_met_info[pro_met_info['GeneID'] == 'OR4F5'].index]

deg_up = list(map(lambda x: x.strip('\n'), open("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/deg_up_protein.txt", 'r').readlines()))
deg_down = list(map(lambda x: x.strip('\n'), open("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/deg_down_protein.txt", 'r').readlines()))

tfs = list(map(lambda x: x.strip('\n'), open("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/hs_hgnc_curated_tfs.txt", 'r').readlines()))
tfs_in_deg =
deg_rowcolors1 = list(map(lambda x: '#800000' if x in deg_up else '#000080', list(deg_tn_protein.loc[(deg_tn_protein.index.isin(deg_up)) | (deg_tn_protein.index.isin(deg_down))].index)))
deg_rowcolors2 = list(map(lambda x: '#000000' if x in tfs else '#FFFFFF', list(deg_tn_protein.loc[(deg_tn_protein.index.isin(deg_up)) | (deg_tn_protein.index.isin(deg_down))].index)))

sns.clustermap(
    rna.loc[(rna.index.isin(deg_up)) | (rna.index.isin(deg_down))].iloc[:,84:],
    method="ward",
    metric="euclidean",
    z_score=None,
    standard_scale=0,
    row_cluster=True,
    col_cluster=True,
    xticklabels=True,
    yticklabels=False,
    col_colors=[col_colors1],
    row_colors=deg_rowcolors1,
    cmap=cmap2
)

deg_tn_protein.loc[:, (deg_tn_protein.columns.isin(dmr_t.obs['DMR_leiden'][dmr_t.obs['DMR_leiden'] == 'leiden_A'].index)) | (deg_tn_protein.columns.isin(dmr_t.obs['DMR_leiden'][dmr_t.obs['DMR_leiden'] == 'leiden_D'].index)) ].loc[deg_tn_protein.index.isin(deg_up)]

#rfh = open("DMR_tumor_pseudotime_and_DEG_TN_protein_Correlation.txt", 'w')
#rfh.write("ID\tPearsonR\tPvalue\n")

r_list, p_list = list(), list()
for i in range(len(deg_tn_protein.index)):
    r, p = stats.pearsonr(dmr_t.obs['dpt_pseudotime'], deg_tn_protein.iloc[i])
    r_list.append(r)
    p_list.append(p)

multipletests(pvals=p_list, alpha=0.05, method='fdr_bh')[1]
pt_degtn = pd.DataFrame([r_list, multipletests(pvals=p_list, alpha=0.05, method='fdr_bh')[1]], columns=deg_tn_protein.index, index=['PearsonR', 'Padj']).T

leiden_name_change_dict1 = {'3': 'leiden_A',
                           '1': 'leiden_B',
                           '0': 'leiden_C',
                           '2': 'leiden_D'}
dmr_t.obs['DMR_leiden'] = dmr_t.obs['leiden_r05'].map(lambda x: leiden_name_change_dict1[x]).astype('category')
col_colors1 = list(dict(zip(list(dmr_t.obs['DMR_leiden'].value_counts().index), dmr_t.uns['DMR_leiden_colors']))[x] for x in dmr_t.obs['DMR_leiden'])
sample_cluster_color = pd.concat([dmr_t.obs[['DMR_leiden']], pd.DataFrame(col_colors1, index=dmr_t.obs.index, columns=['Colors'])], axis=1)

#col_colors1 = list(dict(zip(list(dmr_t.obs['DMR_leiden'].value_counts().index), dmr_t.uns['DMR_leiden_colors']))[x] for x in dmr_t.obs['DMR_leiden']) * 2

dmr_t_met = pd.DataFrame(dmr_t.raw.X.T, index=dmr_t.var.index, columns=dmr_t.obs.index)
g = sns.clustermap(dmr_t_met, method='ward', metric='euclidean', z_score=None, standard_scale=None, cmap=cmap, xticklabels=True, yticklabels=False, col_colors=[col_colors1])

deg_up_sig_cor = list()
for i in deg_up:
    try:
        pearson = stats.pearsonr(np.delete(pro_met.loc[pro_met_info[pro_met_info['GeneID'] == i].index].iloc[:,84:].values[0], np.argwhere(np.isnan(pro_met.loc[pro_met_info[pro_met_info['GeneID'] == i].index].iloc[:,84:].values[0]))), np.delete(deg_tn_protein.loc[i].values, np.argwhere(np.isnan(pro_met.loc[pro_met_info[pro_met_info['GeneID'] == i].index].iloc[:,84:].values[0]))))
        if pearson[1] < 0.05: # p value만 가지고 해야 함 나중에는
            if pearson[0] < 0:
                deg_up_sig_cor.append(i)
            print(i + '\t' + str(pearson[0]) + '\t' + str(pearson[1]))
    except (IndexError, ValueError):
        pass

deg_up_sig_cor = list()
for i in deg_up:
    try:
        spearman = stats.spearmanr(np.delete(pro_met.loc[pro_met_info[pro_met_info['GeneID'] == i].index].iloc[:,84:].values[0], np.argwhere(np.isnan(pro_met.loc[pro_met_info[pro_met_info['GeneID'] == i].index].iloc[:,84:].values[0]))), np.delete(deg_tn_protein.loc[i].values, np.argwhere(np.isnan(pro_met.loc[pro_met_info[pro_met_info['GeneID'] == i].index].iloc[:,84:].values[0]))))
        if spearman[1] < 0.05:
            if spearman[0] < 0: # p value만 가지고 해야 함 나중에는
                deg_up_sig_cor.append(i)
            print(i + '\t' + str(spearman[0]) + '\t' + str(spearman[1]))
    except (IndexError, ValueError):
        pass

deg_down_sig_cor = list()
for i in deg_down:
    try:
        pearson = stats.pearsonr(np.delete(pro_met.loc[pro_met_info[pro_met_info['GeneID'] == i].index].iloc[:,84:].values[0], np.argwhere(np.isnan(pro_met.loc[pro_met_info[pro_met_info['GeneID'] == i].index].iloc[:,84:].values[0]))), np.delete(deg_tn_protein.loc[i].values, np.argwhere(np.isnan(pro_met.loc[pro_met_info[pro_met_info['GeneID'] == i].index].iloc[:,84:].values[0]))))
        if pearson[1] < 0.05:
            if pearson[0] < 0: # p value만 가지고 해야 함 나중에는
                deg_down_sig_cor.append(i)
            print(i + '\t' + str(pearson[0]) + '\t' + str(pearson[1]))
    except (IndexError, ValueError):
        pass

deg_down_sig_cor = list()
for i in deg_down:
    try:
        spearman = stats.spearmanr(np.delete(pro_met.loc[pro_met_info[pro_met_info['GeneID'] == i].index].iloc[:,84:].values[0], np.argwhere(np.isnan(pro_met.loc[pro_met_info[pro_met_info['GeneID'] == i].index].iloc[:,84:].values[0]))), np.delete(deg_tn_protein.loc[i].values, np.argwhere(np.isnan(pro_met.loc[pro_met_info[pro_met_info['GeneID'] == i].index].iloc[:,84:].values[0]))))
        if spearman[1] < 0.05:
            if spearman[0] < 0: # p value만 가지고 해야 함 나중에는
                deg_down_sig_cor.append(i)
            print(i + '\t' + str(spearman[0]) + '\t' + str(spearman[1]))
    except (IndexError, ValueError):
        pass

# Promoter methylation boxplot
gene = 'SOX2'
df = pd.concat([pro_met.loc[pro_met_info[pro_met_info['GeneID'] == gene].index].iloc[:,84:].T, pd.DataFrame(dmr_t.obs['DMR_leiden2'])], axis=1)
sns.boxplot(data=df, x='DMR_leiden2', y=df.columns[0], palette=color_dict2)
sns.swarmplot(data=df, x='DMR_leiden2', y=df.columns[0], color=".2").set(ylabel=df.columns[0].split('/')[1] + ' promoter methylation')

# Gene boxplot
# (i) 이렇게 하거나
gene = 'TWIST1'
sns.boxplot(data=pd.concat([pd.DataFrame(deg_tn_protein.loc[gene]), pd.DataFrame(dmr_t.obs['DMR_leiden'])], axis=1), x='DMR_leiden', y=gene, palette=color_dict)
sns.swarmplot(data=pd.concat([pd.DataFrame(deg_tn_protein.loc[gene]), pd.DataFrame(dmr_t.obs['DMR_leiden'])], axis=1), x='DMR_leiden', y=gene, color=".2").set(ylabel=gene + ' expression')

# (ii) or 이렇게 해
sns.boxplot(data=pd.concat([pd.DataFrame(rna.iloc[:,84:].loc[gene]), pd.DataFrame(dmr_t.obs['DMR_leiden2'])], axis=1), x='DMR_leiden2', y=gene, palette=color_dict2)
sns.swarmplot(data=pd.concat([pd.DataFrame(rna.iloc[:,84:].loc[gene]), pd.DataFrame(dmr_t.obs['DMR_leiden2'])], axis=1), x='DMR_leiden2', y=gene, color=".2").set(ylabel=gene + ' expression')

# 굳이 둘 간의 mann-whitney U test가 하고 싶다면
df = pd.concat([pd.DataFrame(rna.iloc[:,84:].loc['ESRRG']), pd.DataFrame(dmr_t.obs['DMR_leiden'])], axis=1)
print(stats.mannwhitneyu(df[df['DMR_leiden'] == 'leiden_D']['ESRRG'], df[df['DMR_leiden'] == 'leiden_A']['ESRRG'])[1])


# Leiden D vs A 에서 up 되고 down 된 gene 을 row로 두고 leiden_A, D 에 속하는 샘플들만 hierarchical clustering
fuck1 = deg_tn_protein.loc[:, (deg_tn_protein.columns.isin(dmr_t.obs['DMR_leiden'][dmr_t.obs['DMR_leiden'] == 'leiden_A'].index)) | (deg_tn_protein.columns.isin(dmr_t.obs['DMR_leiden'][dmr_t.obs['DMR_leiden'] == 'leiden_D'].index)) ].loc[deg_tn_protein.index.isin(deg_up)]

fuck2 = deg_tn_protein.loc[:, (deg_tn_protein.columns.isin(dmr_t.obs['DMR_leiden'][dmr_t.obs['DMR_leiden'] == 'leiden_A'].index)) | (deg_tn_protein.columns.isin(dmr_t.obs['DMR_leiden'][dmr_t.obs['DMR_leiden'] == 'leiden_D'].index)) ].loc[deg_tn_protein.index.isin(deg_down)]

sns.clustermap(fuck1,
    method="ward",
    metric="euclidean",
    z_score=None,
    standard_scale=0,
    row_cluster=True,
    col_cluster=True,
    xticklabels=True,
    yticklabels=False,
    col_colors=list(sample_cluster_color.loc[fuck2.columns]['Colors'].values),
    cmap=cmap
)
sns.clustermap(fuck2,
    method="ward",
    metric="euclidean",
    z_score=None,
    standard_scale=0,
    row_cluster=True,
    col_cluster=True,
    xticklabels=True,
    yticklabels=False,
    col_colors=list(sample_cluster_color.loc[fuck2.columns]['Colors'].values),
    cmap=cmap
)


# Leiden D vs A 에서 up 되고 down 된 gene의 Promoter 을 row로 두고 leiden_A, D 에 속하는 샘플들만 hierarchical clustering
df = pro_met.loc[pro_met_info[pro_met_info['GeneID'].isin(deg_up)].index].iloc[:,84:]
df = df.loc[:, (df.columns.isin(dmr_t.obs['DMR_leiden'][dmr_t.obs['DMR_leiden'] == 'leiden_A'].index)) | (df.columns.isin(dmr_t.obs['DMR_leiden'][dmr_t.obs['DMR_leiden'] == 'leiden_D'].index))]
df.index = list(map(lambda x: x.split('/')[1], list(df.index)))

""""""
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
col_dism = 1 - df.corr() # samples
row_dism = 1 - df.T.corr() # promoters
row_linkage = hc.linkage(sp.distance.squareform(row_dism), method='complete')
col_linkage = hc.linkage(sp.distance.squareform(col_dism), method='complete')


sns.clustermap(df,
    row_linkage=row_linkage,
    col_linkage=col_linkage,
    row_cluster=True,
    col_cluster=True,
    z_score=None,
    standard_scale=None,
    xticklabels=True,
    yticklabels=True,
    col_colors=list(sample_cluster_color.loc[df.columns]['Colors'].values),
    cmap=cmap
)
""""""


sns.clustermap(df.fillna(0.000000001),
    metric='euclidean',
    method='ward',
    row_cluster=True,
    col_cluster=True,
    z_score=None,
    standard_scale=None,
    xticklabels=True,
    yticklabels=False,
    col_colors=list(sample_cluster_color.loc[df.columns]['Colors'].values),
    cmap=cmap
)

df = pro_met.loc[pro_met_info[pro_met_info['GeneID'].isin(deg_down)].index].iloc[:,84:]
df = df.loc[:, (df.columns.isin(dmr_t.obs['DMR_leiden'][dmr_t.obs['DMR_leiden'] == 'leiden_A'].index)) | (df.columns.isin(dmr_t.obs['DMR_leiden'][dmr_t.obs['DMR_leiden'] == 'leiden_D'].index))]
df.index = list(map(lambda x: x.split('/')[1], list(df.index)))
sns.clustermap(df.fillna(0.000000001),
    metric='euclidean',
    method='ward',
    row_cluster=True,
    col_cluster=True,
    z_score=None,
    standard_scale=None,
    xticklabels=True,
    yticklabels=False,
    col_colors=list(sample_cluster_color.loc[df.columns]['Colors'].values),
    cmap=cmap
)

df = pro_met.loc[pro_met_info[pro_met_info['GeneID'].isin(deg_up_sig_cor)].index].iloc[:,84:]
df = df.loc[:, (df.columns.isin(dmr_t.obs['DMR_leiden'][dmr_t.obs['DMR_leiden'] == 'leiden_A'].index)) | (df.columns.isin(dmr_t.obs['DMR_leiden'][dmr_t.obs['DMR_leiden'] == 'leiden_D'].index))]
df.index = list(map(lambda x: x.split('/')[1], list(df.index)))
sns.clustermap(df.fillna(0.000000001),
    metric='euclidean',
    method='ward',
    row_cluster=True,
    col_cluster=True,
    z_score=None,
    standard_scale=None,
    xticklabels=True,
    yticklabels=True,
    col_colors=list(sample_cluster_color.loc[df.columns]['Colors'].values),
    cmap=cmap
)

df = pro_met.loc[pro_met_info[pro_met_info['GeneID'].isin(deg_down_sig_cor)].index].iloc[:,84:]
df = df.loc[:, (df.columns.isin(dmr_t.obs['DMR_leiden'][dmr_t.obs['DMR_leiden'] == 'leiden_A'].index)) | (df.columns.isin(dmr_t.obs['DMR_leiden'][dmr_t.obs['DMR_leiden'] == 'leiden_D'].index))]
df.index = list(map(lambda x: x.split('/')[1], list(df.index)))
sns.clustermap(df.fillna(0.000000001),
    metric='euclidean',
    method='ward',
    row_cluster=True,
    col_cluster=True,
    z_score=None,
    standard_scale=None,
    xticklabels=True,
    yticklabels=True,
    col_colors=list(sample_cluster_color.loc[df.columns]['Colors'].values),
    cmap=cmap
)

# DMR_leiden clusters and WHO classification
ax = pd.crosstab(dmr_t.obs['WHO'], dmr_t.obs['leiden_r05'], normalize=1).T.plot.bar(stacked=True, color=dict(zip(dmr_t.obs['WHO'].unique(), sns.color_palette("Set2", len(dmr_t.obs['WHO'].unique())))))
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))
plt.tight_layout()

# DMR_leiden clusters and Lauren's classification
ax = pd.crosstab(df['Lauren'], df['leiden_r05'], normalize=1).T.plot.bar(stacked=True, color=dict(zip(dmr_t.obs['Lauren'].unique(), sns.color_palette("Accent", len(dmr_t.obs['Lauren'].unique())))))
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))
plt.tight_layout()
