import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors
import scipy.stats as stats
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd

sc.settings.verbosity = 3
plt.rcParams['figure.figsize'] = (5,5)
sns.set(font="Arial", font_scale=1, style='ticks')
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])
cmap4 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#191970", "#ffdab9", "#8B0000"])
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#000000", "#8b0a50"])
cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#0057b7", "#000000", "#ffd700"])
cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#0057B8", "#000000", "#ffd700"])
%matplotlib
%autoindent

clinic_info = pd.read_csv('/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/2022_WC300_clinical_information_Xadded.csv', index_col='Sample')

# Average methylation for Tumor vs Normal (Figure 1)
p = sns.violinplot(data=clinic_info.iloc[84:][['PercentMet_COV5_Normal', 'PercentMet_COV5_Tumor']], palette={'PercentMet_COV5_Normal':'midnightblue', 'PercentMet_COV5_Tumor':'darkred'}, cut=0, scale="count")
p = sns.stripplot(data=clinic_info.iloc[84:][['PercentMet_COV5_Normal', 'PercentMet_COV5_Tumor']], color="black")
p.set_xticklabels(['Normal (N=84)', 'Tumor (N=84)'])
p.set_ylabel("Average methylation (%)")
sns.despine()

# PMD methylation
sns.clustermap(pmd_met,
                   col_cluster=False,
                   method='ward',
                   metric='euclidean',
                   z_score=None,
                   standard_scale=None,
                   cmap=cmap4,
                   xticklabels=False,
                   yticklabels=False,
                   col_colors=[col_colors1],
                   row_colors=None,
                   cbar_kws={'label': 'DNA methylation'}, vmin=30, vmax=100)

# PMD fraction
pmdf = pd.read_table("PMD_fraction.txt", index_col=0, sep='\t')
p = sns.violinplot(data=pmdf, x='TN', y='PMD_Fraction', palette={'Normal':'midnightblue', 'Tumor':'darkred'}, cut=0, scale="count")
p = sns.stripplot(data=pmdf, x='TN', y='PMD_Fraction', color=".3", alpha=0.5)
p.set_ylabel("Fraction of PMDs in the genome")
p.set_xlabel("")
sns.despine()
plt.tight_layout()

# TCGA-STAD 450K comparison
tcga = pd.read_table("./TCGA/TCGA-STAD_analysis_ready_woY.tsv", index_col=0)
snu = pd.read_table("./TCGA/DNA_methylation_NT84_modified.bg", index_col=0)
tcga = tcga[~tcga.index.duplicated(keep='first')]

comb_tcga_snu = pd.concat([snu, tcga], axis=1)
corr_tcga_snu = comb_tcga_snu.corr()

colors = ['#E6E6FA']*84 + ['#800000']*84 + ['#00008B']*2 + ['#F08080']*395

sns.clustermap(corr_tcga_snu,
method='ward',
annot=False,
cmap=cmap,
row_colors=colors,
col_colors=colors,
xticklabels=False,
yticklabels=False,
vmin=0.7,
vmax=1)



# Call Imputed DMR for all samples
dmr = pd.read_csv("DMR_abs15_ALL_imputed_corrected.csv", index_col=0).iloc[:,:-1] # Type column removal
dmr = sc.AnnData(dmr)
dmr.raw = dmr
dmr.layers['Percent_met'] = dmr.X
# np.ndarray.min(dmr.raw.X) ==> 전체 table에서 minimum value 값 (maximum은 min==> max)

dmr.obs = clinic_info.copy() # copy 반드시 집어넣어야함

sc.pp.scale(dmr)
sc.tl.pca(dmr, n_comps=100, zero_center=True)
sc.pl.pca(dmr, color='TN', palette={'Normal':'midnightblue', 'Tumor':'darkred'}, annotate_var_explained=True, size=100)
sns.despine()
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


#### Tumor Analysis Start ####
# Imputed DMR
dmr_met = pd.read_csv("DMR_abs15_ALL_imputed_corrected.csv", index_col=0).iloc[:,:-1].T

# Row name change
dmr_met.index = list(map(lambda x: x.split('.')[0] + ':' + x.split('.')[1] + '-' + x.split('.')[2], dmr_met.index))

# DMR annotation table
dmr_info = pd.read_table("DMR_abs15_Hyper-Hypo_annotation.txt", index_col=0)

# Normalized CpG density to DMR annotation table
dmr_info['Norm_CpGdensity'] = (np.max(dmr_info['CpGdensity']) - dmr_info['CpGdensity']) / np.max(dmr_info['CpGdensity'])

# CpG density boxplot between Hyper-DMR and Hypo-DMR
p = sns.boxplot(data=dmr_info, x='Type', y='Norm_CpGdensity', order=['Hyper', 'Hypo'], palette={'Hyper': '#A23E48', 'Hypo': '#6C8EAD'}, showfliers = False)
p = sns.stripplot(data=dmr_info, x='Type', y='Norm_CpGdensity', order=['Hyper', 'Hypo'], jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
p.set_ylabel("Normalized CpG density of DMR")
p.set_xlabel("Types of DMR")
p.set_xticklabels(['Hypermethylation', 'Hypomethylation'])
plt.tight_layout()

# DMR annotation colormap for sns.clustermap
row_colors_dmr1 = list(dict(zip(['Hypo', 'Hyper'], ['#6C8EAD', '#A23E48']))[x] for x in dmr_info['Type'])

# Clustering
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
g.cax.set_visible(False) # Legend removal


# DNA methylation of DMR for Tumor-only
dmr_t = sc.AnnData(pd.read_csv("DMR_abs15_ALL_imputed_corrected.csv", index_col=0).iloc[84:,:-1])

# Row name change
dmr_t.var_names = list(map(lambda x: x.split('.')[0] + ':' + x.split('.')[1] + '-' + x.split('.')[2], dmr_t.var_names))

# Percentage methylation stored on raw and layer attribute
dmr_t.raw = dmr_t
dmr_t.layers['Percent_met'] = dmr_t.X

# clinic info attachment
dmr_t.obs = clinic_info.iloc[84:].copy()

# RNA batch info
rna_batch = pd.read_table("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/RNA_batch.txt", index_col=0)
dmr_t.obs['RNA_batch'] = rna_batch.iloc[84:]

# Scaling and PCA
sc.pp.scale(dmr_t)
sc.tl.pca(dmr_t, n_comps=83, zero_center=True)

# Check for PCA numbers
# pca_variance_t = pd.DataFrame(dmr_t.uns['pca']['variance_ratio'], index=list(map(lambda x: 'PC' + str(x), list(range(1,84)))), columns=['Variance_ratio'])
# np.sum(pca_variance_t.values.flatten()[:12])

# Neighbor Graph construction and leiden community detection
sc.pp.neighbors(dmr_t, n_neighbors=20, n_pcs=12)
sc.tl.leiden(dmr_t, resolution=0.75, key_added='leiden_r075')
sc.tl.leiden(dmr_t, resolution=0.5, key_added='leiden_r05')
sc.tl.umap(dmr_t, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
sc.pl.umap(dmr_t, color=['leiden_r075', 'leiden_r05', 'EpiBurden'], add_outline=False, legend_loc='right margin', color_map=cmap)

# Leiden cluster name change
leiden_name_change_dict1 = {'3': 'leiden_A',
                           '1': 'leiden_B',
                           '0': 'leiden_C',
                           '2': 'leiden_D'} # leiden_r075

leiden_name_change_dict2 = {'1': 'leiden_A',
                           '0': 'leiden_B',
                           '2': 'leiden_C'} # leiden_r05

# Leiden cluster name change: leiden_r075 ==> DMR_leiden
dmr_t.obs['DMR Clusters'] = dmr_t.obs['leiden_r075'].map(lambda x: leiden_name_change_dict1[x]).astype('category')
dmr_t.obs['DMR Clusters2'] = dmr_t.obs['leiden_r05'].map(lambda x: leiden_name_change_dict2[x]).astype('category')

# UMAP projection plot for DMR Clusters and Epimutation Burden
dmr_t.obs['Epimutation Burden'] = dmr_t.obs['EpiBurden']
sc.pl.umap(dmr_t, color=['DMR Clusters', 'Epimutation Burden'], add_outline=False, legend_loc='on data', palette='Dark2')
sns.despine()

# Ordering of DMR Clusters leiden_A to D
lin = tuple(sorted(list(dmr_t.obs['DMR Clusters'].values.unique())))
dmr_t.obs['DMR Clusters'] = dmr_t.obs['DMR Clusters'].cat.reorder_categories(list(lin), ordered=True)
color_dict2 = {
    "leiden_A": "#1b9e77",
    "leiden_B": "#7570b3",
    "leiden_C": "#e6ab02",
    "leiden_D": "#666666"} # DMR Clusters (leiden_r075)
lin = tuple(sorted(list(dmr_t.obs['DMR Clusters2'].values.unique())))
dmr_t.obs['DMR Clusters2'] = dmr_t.obs['DMR Clusters2'].cat.reorder_categories(list(lin), ordered=True)

color_dict = {
    "leiden_A": "#a6cee3",
    "leiden_B": "#fdbf6f",
    "leiden_C": "#b15928"} # DMR Clusters (leiden_r05)


# Distribution of Epimutation burden along DMR Clusters
p = sns.boxplot(data=dmr_t.obs, x='DMR Clusters2', y='EpiBurden', palette=color_dict)
p = sns.stripplot(data=dmr_t.obs, x='DMR Clusters2', y='EpiBurden', jitter=True, marker='o', color='black', alpha=0.8)
p.set_xlabel('DMR Clusters')
p.set_ylabel('Epimutation Burden')
plt.tight_layout()
sns.despine()

# ANOVA for Epimutation burden along DMR Clusters
df = dmr_t.obs[['EpiBurden', 'leiden_r05']]
df_lm = ols('EpiBurden ~ C(leiden_r05)', data=df).fit() # C() ==> categorical data (not necessary here because batch is already categorical)
print(sm.stats.anova_lm(df_lm, typ=2))
print(pairwise_tukeyhsd(df['EpiBurden'], df['leiden_r05'], alpha=0.05))

# Diffusion maps
sc.tl.diffmap(dmr_t)
sc.pl.diffmap(dmr_t, color=['DMR Clusters'], add_outline=False, color_map=cmap)

# Diffusion pseudotime calculation
start_cell = np.isin(dmr_t.obs['DMR Clusters'], 'leiden_A')
max_start_id = np.argmax(dmr_t.obsm['X_diffmap'][start_cell, 1])
root_id = np.arange(len(start_cell))[start_cell][max_start_id]
dmr_t.uns['iroot'] = root_id
sc.tl.dpt(dmr_t, n_branchings=0, n_dcs=15)

# Plot scatterplot and regression model fits for Diffusion pseudotime and Epimutation burden

p = sns.lmplot(data=dmr_t.obs, x='dpt_pseudotime', y='Epimutation Burden')
p.set_xlabels('Diffusion Pseudotime')
p.set_ylabels('Epimutation Burden')
stats.spearmanr(dmr_t.obs['dpt_pseudotime'], dmr_t.obs['EpiBurden'])



# DMR Clusters and WHO classification (LI, VI, PI, Lauren's...etc)
ax = pd.crosstab(dmr_t.obs['WHO'], dmr_t.obs['DMR_leiden'], normalize=1).T.plot.bar(stacked=True, color=dict(zip(dmr_t.obs['WHO'].unique(), sns.color_palette("Set2", len(dmr_t.obs['WHO'].unique())))))
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))
plt.tight_layout()

# Hierarchical clustering of DNA methylation level of DMR for Tumors
dmr_t_met = pd.DataFrame(dmr_t.raw.X.T, index=dmr_t.var.index, columns=dmr_t.obs.index)

col_colors_DMR_leiden = list(dict(zip(sorted(list(dmr_t.obs['DMR Clusters'].value_counts().index)), \
                                      dmr_t.uns['DMR Clusters_colors']))[x] for x in dmr_t.obs['DMR Clusters'])
col_colors_DMR_leiden2 = list(dict(zip(sorted(list(dmr_t.obs['DMR Clusters2'].value_counts().index)), \
                                      dmr_t.uns['DMR Clusters2_colors']))[x] for x in dmr_t.obs['DMR Clusters2'])
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

# Heatmap of Normal DNA methylation in DMR
dmr_t_met_colorder = g.dendrogram_col.reordered_ind
dmr_t_met_roworder = g.dendrogram_row.reordered_ind
dmr_n_met = dmr_met.iloc[:,:84]
g = sns.clustermap(dmr_n_met.iloc[dmr_t_met_roworder, dmr_t_met_colorder],
                   col_cluster=False,
                   row_cluster=False,
                   cmap=cmap,
                   xticklabels=False,
                   yticklabels=False)


## Determining CIMP Tumors

# Call CpGi methylation
#cpgi_met = pd.read_table("CpGi_ALL.txt", index_col=0)
cpgi_met = pd.read_table("CpGi_smooth.txt", index_col=0)
cpgi_met = cpgi_met * 100
# Column name change
cpgi_met.columns = list(map(lambda x: 'X'+x, cpgi_met.columns))
#cpgi_met.index = list(map(lambda x: '/'.join(x.split('/')[:2]), cpgi_met.index))

# Discard missing CpGi DNA methylation rows & pick CpGi sites where Normal DNA methylation < 40
cpgi_met = cpgi_met[cpgi_met.isna().sum(axis=1) == 0][cpgi_met[cpgi_met.isna().sum(axis=1) == 0].iloc[:, :84].mean(axis=1) < 40]

# Call Promoter CpGi
cpgi_pls = list(map(lambda x: x.strip('\n').split('/')[0], open("PLS_CpGi.txt", 'r').readlines()))

# Select Promoter CpGi
cpgi_met = cpgi_met[cpgi_met.index.isin(cpgi_pls)]

# mean(Tumor - Normal) >= 10%
cpgi_tn_met = cpgi_met.iloc[:, 84:] - cpgi_met.iloc[:, :84].values
cpgi_tn_met = cpgi_tn_met[cpgi_tn_met.mean(axis=1) >= 10]

# Hierarchical clustering
g = sns.clustermap(cpgi_tn_met,
                   method='ward',
                   metric='euclidean',
                   z_score=None,
                   standard_scale=None,
                   cmap=cmap,
                   xticklabels=False,
                   yticklabels=False,
                   col_colors=[col_colors_DMR_leiden])

# Attach CIMP annotation for Tumor DNA methylation
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


# Deposit CIMP information inside Tumor DNA methylation anndata
dmr_t.obs['CIMP'] = list(map(lambda x: 'CIMP(+)' if x == True else 'CIMP(-)', dmr_t.obs.index.isin(cpgi_tn_met.iloc[:,g.dendrogram_col.reordered_ind[:35]].columns)))
dmr_t.obs[['DMR Clusters2', 'CIMP']].value_counts().sort_index()
#DMR Clusters2  CIMP
#leiden_A       CIMP(+)     2
#               CIMP(-)    30
#leiden_B       CIMP(+)    22
#               CIMP(-)    13
#leiden_C       CIMP(+)    11
#               CIMP(-)     6
#dtype: int64


# CIMP proportional plot
ax = (pd.crosstab(dmr_t.obs['DMR Clusters'], dmr_t.obs['CIMP'], normalize=0)*100).plot.bar(stacked=True, color=['#8b0000ff', '#000080ff'], rot=0)
plt.ylabel("Proportion (%)")
ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1.0))
plt.tight_layout()
sns.despine()

sc.pl.umap(dmr_t, color=['CIMP'], add_outline=False, legend_loc='right margin', palette={'CIMP(-)':'navy', 'CIMP(+)':'darkred'})
sns.despine()

## Promoter DNA methylation

# Call PLS DNA methylation files
pls = pd.read_table("PLS_ALL.txt", index_col=0)
pls_tn = pls.iloc[:,84:] - pls.iloc[:,:84].values
pls_tn[pls_tn.isna().sum(axis=1) <= 21] # 21 / 84 == 25%
pls_index = pls_tn[pls_tn.isna().sum(axis=1) <= 21].index

pls = pd.read_csv("PLS_ALL_isna25percent_imputed.csv", index_col=0).iloc[:,:-1].T
pls.index = pls_index
del pls_index

# Call PLS information files (contain histone modification and etc)
pls_info = pd.read_table("PLS_annotated_table_full.txt", index_col=0)

# Pick PLS only if it contains PLS information on pls_info

# pls[~pls.index.isin(pls_info[pls_info.index.isin(pls.index)].index)]
# ==> THRA1/BTR 유전자 transcript 두 가지가 나온다. 이들은 제외!
pls = pls.drop(labels=['chr17:48293862-48294411/EH38D2938225/THRA1/BTR/ENSG00000235300/ENST00000621191/antisense/11', 'chr17:48293862-48294411/EH38D2938225/THRA1/BTR/ENSG00000235300/ENST00000604191/antisense/11'], axis=0)
pls = pls.loc[pls_info[pls_info.index.isin(pls.index)].index]

# pick PLS info only if it contains PLS.index
pls_info = pls_info[pls_info.index.isin(pls.index)]

# Promoter methylation boxplot
color_dict3 = {
    "Normal": "#FFFAFA",
    "leiden_A": "#1b9e77",
    "leiden_B": "#7570b3",
    "leiden_C": "#e6ab02",
    "leiden_D": "#666666"}

# Further filtering PLS using DEG genes
pls_info = pls_info[pls_info['GeneID'].isin(deg_tn.index)]
pls = pls.loc[pls_info[pls_info.index.isin(pls.index)].index]

pls_dmr_abs10_info = pd.read_table("PLS_annotated_table_full_DMR_abs10.txt", index_col=0)
pls_dmr_abs10_info = pls_dmr_abs10_info[pls_dmr_abs10_info.index.isin(pls.index)]

pls_dmr_info = pd.read_table("PLS_annotated_table_full_DMR.txt", index_col=0)
pls_dmr_info = pls_dmr_info[pls_dmr_info.index.isin(pls.index)]

### Regulatory element enrichment
fe = pd.read_table("FoldEnrichment_Regulatory_Element.txt", index_col=0)
fe = fe.loc[['Promoter', 'CpGi', 'H3K4me3', 'H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K9me3', 'Super Enhancer', 'TFBS' ]]
ax = fe.plot.bar(rot=30, color = ['lightgrey', 'darkred', 'navy'], figsize=(15,8), fontsize=13, legend=True)
ax.legend(loc='upper left', bbox_to_anchor=(0.88, 0.9), fontsize=20)
plt.xlabel("")
plt.ylabel("Fold Enrichment (x)", fontsize=13)
plt.tight_layout()
sns.despine()

### PLS division

# CpGi effect
a = pd.DataFrame(rna_norm_log2[rna_norm_log2.index.isin(list(pls_info[pls_info['CpGi'] == 'Yes']['GeneID'].values))].iloc[:,:84].mean(axis=1), columns=['CpGi Promoter'])
b = pd.DataFrame(rna_norm_log2[rna_norm_log2.index.isin(list(pls_info[pls_info['CpGi'] == 'na']['GeneID'].values))].iloc[:,:84].mean(axis=1), columns=['Non-CpGi Promoter'])
ab = pd.concat([a,b], axis=1)
p = sns.boxplot(data=ab, palette={'CpGi Promoter':'#00203FFF', 'Non-CpGi Promoter':'#ADEFD1FF'}, width=0.5, showfliers = False)
p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
p.set_ylabel("Gene expression")
p.set_title("CpG islands present in Promoter or not")
plt.tight_layout()
sns.despine()
stats.ttest_ind(a, b)

# Open chromatin effect
a = pd.DataFrame(rna_norm_log2[rna_norm_log2.index.isin(list(pls_info[pls_info['ATAC'] == 'Yes']['GeneID'].values))].iloc[:,:84].mean(axis=1), columns=['Open Promoter'])
b = pd.DataFrame(rna_norm_log2[rna_norm_log2.index.isin(list(pls_info[pls_info['ATAC'] == 'na']['GeneID'].values))].iloc[:,:84].mean(axis=1), columns=['Closed Promoter'])
ab = pd.concat([a,b], axis=1)
p = sns.boxplot(data=ab, palette={'Open Promoter':'#B1624EFF', 'Closed Promoter':'#5CC8D7FF'}, width=0.5, showfliers = False)
p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
p.set_ylabel("Gene expression")
p.set_title("Promoter Open or closed")
plt.tight_layout()
sns.despine()
stats.ttest_ind(a, b)

# K4me3 effect
a = pd.DataFrame(rna_norm_log2[rna_norm_log2.index.isin(list(pls_info[pls_info['K4me3'] == 'Yes']['GeneID'].values))].iloc[:,:84].mean(axis=1), columns=['H3K4me3 Promoter'])
b = pd.DataFrame(rna_norm_log2[rna_norm_log2.index.isin(list(pls_info[pls_info['K4me3'] == 'na']['GeneID'].values))].iloc[:,:84].mean(axis=1), columns=['Non-H3K4me3 Promoter'])
ab = pd.concat([a,b], axis=1)
p = sns.boxplot(data=ab, palette={'H3K4me3 Promoter':'#89ABE3FF', 'Non-H3K4me3 Promoter':'#FCF6F5FF'}, width=0.5, showfliers = False)
p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
p.set_ylabel("Gene expression")
p.set_title("H3K4me3 present in Promoter or not")
plt.tight_layout()
sns.despine()
stats.ttest_ind(a, b)

# K27ac effect
a = pd.DataFrame(rna_norm_log2[rna_norm_log2.index.isin(list(pls_info[pls_info['K27ac'] == 'Yes']['GeneID'].values))].iloc[:,:84].mean(axis=1), columns=['H3K27ac Promoter'])
b = pd.DataFrame(rna_norm_log2[rna_norm_log2.index.isin(list(pls_info[pls_info['K27ac'] == 'na']['GeneID'].values))].iloc[:,:84].mean(axis=1), columns=['Non-H3K27ac Promoter'])
ab = pd.concat([a,b], axis=1)
p = sns.boxplot(data=ab, palette={'H3K27ac Promoter':'#F2AA4CFF', 'Non-H3K27ac Promoter':'#FCF6F5FF'}, width=0.5, showfliers = False)
p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
p.set_ylabel("Gene expression")
p.set_title("H3K27ac present in Promoter or not")
plt.tight_layout()
sns.despine()
stats.ttest_ind(a, b)

# K4me1 effect
a = pd.DataFrame(rna_norm_log2[rna_norm_log2.index.isin(list(pls_info[pls_info['K4me1'] == 'Yes']['GeneID'].values))].iloc[:,:84].mean(axis=1), columns=['H3K4me1 Promoter'])
b = pd.DataFrame(rna_norm_log2[rna_norm_log2.index.isin(list(pls_info[pls_info['K4me1'] == 'na']['GeneID'].values))].iloc[:,:84].mean(axis=1), columns=['Non-H3K4me1 Promoter'])
ab = pd.concat([a,b], axis=1)
p = sns.boxplot(data=ab, palette={'H3K4me1 Promoter':'#D7C49EFF', 'Non-H3K4me1 Promoter':'#343148FF'}, width=0.5, showfliers = False)
p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
p.set_ylabel("Gene expression")
p.set_title("H3K4me1 present in Promoter or not")
plt.tight_layout()
sns.despine()
stats.ttest_ind(a, b)

# K27me3 effect
a = pd.DataFrame(rna_norm_log2[rna_norm_log2.index.isin(list(pls_info[pls_info['K27me3'] == 'Yes']['GeneID'].values))].iloc[:,:84].mean(axis=1), columns=['H3K27me3 Promoter'])
b = pd.DataFrame(rna_norm_log2[rna_norm_log2.index.isin(list(pls_info[pls_info['K27me3'] == 'na']['GeneID'].values))].iloc[:,:84].mean(axis=1), columns=['Non-H3K27me3 Promoter'])
ab = pd.concat([a,b], axis=1)
p = sns.boxplot(data=ab, palette={'H3K27me3 Promoter':'#E94B3CFF', 'Non-H3K27me3 Promoter':'#FCF6F5FF'}, width=0.5, showfliers = False)
p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
p.set_ylabel("Gene expression")
p.set_title("H3K27me3 present in Promoter or not")
plt.tight_layout()
sns.despine()
stats.ttest_ind(a, b)

# K36me3 effect (Normal)
a = pd.DataFrame(rna_norm_log2[rna_norm_log2.index.isin(list(pls_info[pls_info['K36me3'] == 'Yes']['GeneID'].values))].iloc[:,:84].mean(axis=1), columns=['H3K36me3 Genes'])
b = pd.DataFrame(rna_norm_log2[rna_norm_log2.index.isin(list(pls_info[pls_info['K36me3'] == 'na']['GeneID'].values))].iloc[:,:84].mean(axis=1), columns=['Non-H3K36me3 Genes'])
ab = pd.concat([a,b], axis=1)
p = sns.boxplot(data=ab, palette={'H3K36me3 Genes':'#990011FF', 'Non-H3K36me3 Genes':'#FCF6F5FF'}, width=0.5, showfliers = False)
p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
p.set_ylabel("Gene expression")
p.set_title("H3K36me3 present in Gene Body or not")
plt.tight_layout()
sns.despine()
stats.ttest_ind(a, b)

# K36me3 effect (Tumor)
a = pd.DataFrame(rna_norm_log2[rna_norm_log2.index.isin(list(pls_info[pls_info['K36me3'] == 'Yes']['GeneID'].values))].iloc[:,84:].mean(axis=1), columns=['H3K36me3 Genes'])
b = pd.DataFrame(rna_norm_log2[rna_norm_log2.index.isin(list(pls_info[pls_info['K36me3'] == 'na']['GeneID'].values))].iloc[:,84:].mean(axis=1), columns=['Non-H3K36me3 Genes'])
ab = pd.concat([a,b], axis=1)
p = sns.boxplot(data=ab, palette={'H3K36me3 Genes':'#990011FF', 'Non-H3K36me3 Genes':'#FCF6F5FF'}, width=0.5, showfliers = False)
p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
p.set_ylabel("Gene expression")
p.set_title("H3K36me3 present in Gene Body or not")
plt.tight_layout()
sns.despine()
stats.ttest_ind(a, b)


pls_dmr_abs10_info.loc[:, ['K4me3','K27ac', 'K27me3']].value_counts()

#K4me3  K27ac  K27me3
#Yes    Yes    na        1746
#na     na     na        1011
#Yes    na     Yes        550
#       Yes    Yes        380
#na     Yes    na         377
#       na     Yes        174
#Yes    na     na         118
#na     Yes    Yes         15

pls_dmr_info.loc[:, ['K4me3','K27ac', 'K27me3']].value_counts()



# Tumor > Normal DEG에서 DMR이 겹치는 Promoter 중에
pls_dmr_abs10_info[pls_dmr_abs10_info['GeneID'].isin(list(deg_tn_uplist))].iloc[:, [-5, -4, -2]].value_counts()

pls_dmr_abs10_info[pls_dmr_abs10_info['GeneID'].isin(list(deg_tn_uplist))].iloc[:, -7:].value_counts()

pls_dmr_deg_tn_up = pls_dmr_abs10_info[pls_dmr_abs10_info['GeneID'].isin(list(deg_tn_uplist))]

pls[pls.index.isin(pls_dmr_abs10_info[(pls_dmr_abs10_info['K4me3'] == 'na') & (pls_dmr_abs10_info['K27ac'] == 'Yes') & (pls_dmr_abs10_info['K27me3'] == 'Yes')].index)].iloc[:, :84]



########################################################################################################################################################
# K4me3 Yes / K27ac Yes / K27me3 Yes
total = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'Yes') & (pls_dmr_info['K27ac'] == 'Yes') & (pls_dmr_info['K27me3'] == 'Yes')].index)]
N = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'Yes') & (pls_dmr_info['K27ac'] == 'Yes') & (pls_dmr_info['K27me3'] == 'Yes')].index)].iloc[:, :84]
T = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'Yes') & (pls_dmr_info['K27ac'] == 'Yes') & (pls_dmr_info['K27me3'] == 'Yes')].index)].iloc[:, 84:]

a = pd.DataFrame(N.mean(axis=1), columns=['Normal'])
b = pd.DataFrame(T.mean(axis=1), columns=['Tumor'])
ab = pd.concat([a,b], axis=1)
#p = sns.boxplot(data=ab, palette={'Normal':'midnightblue', 'Tumor':'darkred'}, width=0.5, showfliers = False)
#p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
#p.set_ylabel("Promoter Methylation")
#p.set_title("CpG islands present in Promoter or not")
p = sns.kdeplot(data=ab, x="Normal", y="Tumor", cmap="afmhot", shade=True, thresh=0)
p.set_xlabel("Normal promoter methylation (%)")
p.set_ylabel("Tumor promoter methylation (%)")
plt.tight_layout()
sns.despine()
plt.xlim((0,85))
plt.ylim((0,85))

stats.ttest_rel(a, b)

p_list = list()
t_list = list()
for i in list(range(len(N.index))):
    ttest = stats.ttest_rel(T.iloc[i], N.iloc[i])
    t_list.append(ttest[0])
    p_list.append(ttest[1])

from statsmodels.stats.multitest import multipletests
t_fdrp = pd.DataFrame([t_list, multipletests(pvals=p_list, alpha=0.05, method='fdr_bh')[1]],
                      columns=N.index, index=['Tstatistics', 'Padj']).T
p = len(t_fdrp[(t_fdrp['Tstatistics'] > 0) & (t_fdrp['Padj'] < 0.05)])
m = len(t_fdrp[(t_fdrp['Tstatistics'] < 0) & (t_fdrp['Padj'] < 0.05)])


print(str(len(N.index)) + '\t' + str(p) + '\t' + str(m))
# 380	260	42

sns.histplot(np.log10(deg_tn[deg_tn.index.isin(list(set(pls_dmr_info[pls_dmr_info.index.isin(total.index)]['GeneID'].values)))]['baseMean']))

p = sns.histplot(t_fdrp[t_fdrp['Padj'] < 0.05]['Tstatistics'])
p.set_xlabel('Paired t-test t-statistics (Tumor-Normal)')
plt.tight_layout()

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


pls_tumorup = list(set(pls_info[pls_info.index.isin(t_fdrp[(t_fdrp['Padj'] < 0.05) & (t_fdrp['Tstatistics'] > 0)].index)]['GeneID'].values)) # 570
pls_tumordown = list(set(pls_info[pls_info.index.isin(t_fdrp[(t_fdrp['Padj'] < 0.05) & (t_fdrp['Tstatistics'] < 0)].index)]['GeneID'].values)) # 35
a = intersection(pls_tumorup, deg_tn_uplist) # 9 ==> ETS2 (Proto oncogene)가 있네
b = intersection(pls_tumorup, deg_tn_downlist) # 78
c = intersection(pls_tumordown, deg_tn_uplist) # 30
d = intersection(pls_tumordown , deg_tn_downlist) # 17

print(len(a))
print(len(b))
print(len(c))
print(len(d))

e = pls_dmr_info[pls_dmr_info['GeneID'].isin(a)]
pls_e = pls[pls.index.isin(e[(e['K4me3'] == 'Yes') & (e['K27ac'] == 'Yes') & (e['K27me3'] == 'Yes')].index)]
rna_e = rna_norm_trans_log2[rna_norm_trans_log2.index.isin(list(e[(e['K4me3'] == 'Yes') & (e['K27ac'] == 'Yes') & (e['K27me3'] == 'Yes')]['ENSTID'].values))]


########################################################################################################################################################
# K4me3 Yes / K27ac Yes / K27me3 No
total = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'Yes') & (pls_dmr_info['K27ac'] == 'Yes') & (pls_dmr_info['K27me3'] == 'na')].index)]
N = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'Yes') & (pls_dmr_info['K27ac'] == 'Yes') & (pls_dmr_info['K27me3'] == 'na')].index)].iloc[:, :84]
T = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'Yes') & (pls_dmr_info['K27ac'] == 'Yes') & (pls_dmr_info['K27me3'] == 'na')].index)].iloc[:, 84:]

a = pd.DataFrame(N.mean(axis=1), columns=['Normal'])
b = pd.DataFrame(T.mean(axis=1), columns=['Tumor'])
ab = pd.concat([a,b], axis=1)
#p = sns.boxplot(data=ab, palette={'Normal':'midnightblue', 'Tumor':'darkred'}, width=0.5, showfliers = False)
#p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
#p.set_ylabel("Promoter Methylation")
#p.set_title("CpG islands present in Promoter or not")
p = sns.kdeplot(data=ab, x="Normal", y="Tumor", cmap="afmhot", shade=True, thresh=0)
p.set_xlabel("Normal promoter methylation (%)")
p.set_ylabel("Tumor promoter methylation (%)")
plt.tight_layout()
sns.despine()
plt.xlim((0,100))
plt.ylim((0,100))
stats.ttest_rel(a, b)

p_list = list()
t_list = list()
for i in list(range(len(N.index))):
    ttest = stats.ttest_rel(T.iloc[i], N.iloc[i])
    t_list.append(ttest[0])
    p_list.append(ttest[1])

from statsmodels.stats.multitest import multipletests
t_fdrp = pd.DataFrame([t_list, multipletests(pvals=p_list, alpha=0.05, method='fdr_bh')[1]],
                      columns=N.index, index=['Tstatistics', 'Padj']).T
p = len(t_fdrp[(t_fdrp['Tstatistics'] > 0) & (t_fdrp['Padj'] < 0.05)])
m = len(t_fdrp[(t_fdrp['Tstatistics'] < 0) & (t_fdrp['Padj'] < 0.05)])


print(str(len(N.index)) + '\t' + str(p) + '\t' + str(m))
# 1746  331 245

sns.histplot(np.log10(deg_tn[deg_tn.index.isin(list(set(pls_dmr_info[pls_dmr_info.index.isin(total.index)]['GeneID'].values)))]['baseMean']))

p = sns.histplot(t_fdrp[t_fdrp['Padj'] < 0.05]['Tstatistics'])
p.set_xlabel('Paired t-test t-statistics (Tumor-Normal)')
plt.tight_layout()

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


pls_tumorup = list(set(pls_info[pls_info.index.isin(t_fdrp[(t_fdrp['Padj'] < 0.05) & (t_fdrp['Tstatistics'] > 0)].index)]['GeneID'].values)) # 570
pls_tumordown = list(set(pls_info[pls_info.index.isin(t_fdrp[(t_fdrp['Padj'] < 0.05) & (t_fdrp['Tstatistics'] < 0)].index)]['GeneID'].values)) # 35
a = intersection(pls_tumorup, deg_tn_uplist) # 35
b = intersection(pls_tumorup, deg_tn_downlist) # 223
c = intersection(pls_tumordown, deg_tn_uplist) # 68
d = intersection(pls_tumordown , deg_tn_downlist) # 65

print(len(a))
print(len(b))
print(len(c))
print(len(d))

e = pls_dmr_info[pls_dmr_info['GeneID'].isin(a)]
pls_e = pls[pls.index.isin(e[(e['K4me3'] == 'Yes') & (e['K27ac'] == 'Yes') & (e['K27me3'] == 'na')].index)]
rna_e = rna_norm_trans_log2[rna_norm_trans_log2.index.isin(list(e[(e['K4me3'] == 'Yes') & (e['K27ac'] == 'Yes') & (e['K27me3'] == 'na')]['ENSTID'].values))]


sns.clustermap(rna_e,
                   col_cluster=False,
                   row_cluster=False,
                   method='complete',
                   metric='correlation',
                   z_score=0,
                   standard_scale=None,
                   cmap='CMRmap',
                   col_colors=[col_colors1],
                   xticklabels=False,
                   yticklabels=False, vmin=-2.5, vmax=3.5)

sns.clustermap(pls_e,
                   col_cluster=False,
                   row_cluster=False,
                   method='complete',
                   metric='correlation',
                   z_score=None,
                   standard_scale=0,
                   cmap=cmap,
                   col_colors=[col_colors1],
                   xticklabels=False,
                   yticklabels=False)





########################################################################################################################################################
# K4me3 Yes / K27ac No / K27me3 Yes
total = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'Yes') & (pls_dmr_info['K27ac'] == 'na') & (pls_dmr_info['K27me3'] == 'Yes')].index)]
N = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'Yes') & (pls_dmr_info['K27ac'] == 'na') & (pls_dmr_info['K27me3'] == 'Yes')].index)].iloc[:, :84]
T = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'Yes') & (pls_dmr_info['K27ac'] == 'na') & (pls_dmr_info['K27me3'] == 'Yes')].index)].iloc[:, 84:]

a = pd.DataFrame(N.mean(axis=1), columns=['Normal'])
b = pd.DataFrame(T.mean(axis=1), columns=['Tumor'])
ab = pd.concat([a,b], axis=1)
#p = sns.boxplot(data=ab, palette={'Normal':'midnightblue', 'Tumor':'darkred'}, width=0.5, showfliers = False)
#p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
#p.set_ylabel("Promoter Methylation")
#p.set_title("CpG islands present in Promoter or not")
p = sns.kdeplot(data=ab, x="Normal", y="Tumor", cmap="afmhot", shade=True, thresh=0)
p.set_xlabel("Normal promoter methylation (%)")
p.set_ylabel("Tumor promoter methylation (%)")
plt.tight_layout()
sns.despine()
plt.xlim((0,85))
plt.ylim((0,85))
stats.ttest_rel(a, b)

p_list = list()
t_list = list()
for i in list(range(len(N.index))):
    ttest = stats.ttest_rel(T.iloc[i], N.iloc[i])
    t_list.append(ttest[0])
    p_list.append(ttest[1])

from statsmodels.stats.multitest import multipletests
t_fdrp = pd.DataFrame([t_list, multipletests(pvals=p_list, alpha=0.05, method='fdr_bh')[1]],
                      columns=N.index, index=['Tstatistics', 'Padj']).T
p = len(t_fdrp[(t_fdrp['Tstatistics'] > 0) & (t_fdrp['Padj'] < 0.05)])
m = len(t_fdrp[(t_fdrp['Tstatistics'] < 0) & (t_fdrp['Padj'] < 0.05)])


print(str(len(N.index)) + '\t' + str(p) + '\t' + str(m))

sns.histplot(np.log10(deg_tn[deg_tn.index.isin(list(set(pls_dmr_info[pls_dmr_info.index.isin(total.index)]['GeneID'].values)))]['baseMean']))

p = sns.histplot(t_fdrp[t_fdrp['Padj'] < 0.05]['Tstatistics'])
p.set_xlabel('Paired t-test t-statistics (Tumor-Normal)')
plt.tight_layout()

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


pls_tumorup = list(set(pls_info[pls_info.index.isin(t_fdrp[(t_fdrp['Padj'] < 0.05) & (t_fdrp['Tstatistics'] > 0)].index)]['GeneID'].values)) # 570
pls_tumordown = list(set(pls_info[pls_info.index.isin(t_fdrp[(t_fdrp['Padj'] < 0.05) & (t_fdrp['Tstatistics'] < 0)].index)]['GeneID'].values)) # 35
a = intersection(pls_tumorup, deg_tn_uplist) # 9 ==> ETS2 (Proto oncogene)가 있네
b = intersection(pls_tumorup, deg_tn_downlist) # 78
c = intersection(pls_tumordown, deg_tn_uplist) # 30
d = intersection(pls_tumordown , deg_tn_downlist) # 17

print(len(a))
print(len(b))
print(len(c))
print(len(d))

e = pls_dmr_info[pls_dmr_info['GeneID'].isin(a)]
pls_e = pls[pls.index.isin(e[(e['K4me3'] == 'Yes') & (e['K27ac'] == 'na') & (e['K27me3'] == 'Yes')].index)]
rna_e = rna_norm_trans_log2[rna_norm_trans_log2.index.isin(list(e[(e['K4me3'] == 'Yes') & (e['K27ac'] == 'na') & (e['K27me3'] == 'Yes')]['ENSTID'].values))]

########################################################################################################################################################
# K4me3 No / K27ac No / K27me3 Yes
total = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'na') & (pls_dmr_info['K27ac'] == 'na') & (pls_dmr_info['K27me3'] == 'Yes')].index)]
N = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'na') & (pls_dmr_info['K27ac'] == 'na') & (pls_dmr_info['K27me3'] == 'Yes')].index)].iloc[:, :84]
T = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'na') & (pls_dmr_info['K27ac'] == 'na') & (pls_dmr_info['K27me3'] == 'Yes')].index)].iloc[:, 84:]

a = pd.DataFrame(N.mean(axis=1), columns=['Normal'])
b = pd.DataFrame(T.mean(axis=1), columns=['Tumor'])
ab = pd.concat([a,b], axis=1)
#p = sns.boxplot(data=ab, palette={'Normal':'midnightblue', 'Tumor':'darkred'}, width=0.5, showfliers = False)
#p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
#p.set_ylabel("Promoter Methylation")
#p.set_title("CpG islands present in Promoter or not")
p = sns.kdeplot(data=ab, x="Normal", y="Tumor", cmap="afmhot", shade=True, thresh=0)
p.set_xlabel("Normal promoter methylation (%)")
p.set_ylabel("Tumor promoter methylation (%)")
plt.tight_layout()
sns.despine()
plt.xlim((0,100))
plt.ylim((0,100))
stats.ttest_rel(a, b)

p_list = list()
t_list = list()
for i in list(range(len(N.index))):
    ttest = stats.ttest_rel(T.iloc[i], N.iloc[i])
    t_list.append(ttest[0])
    p_list.append(ttest[1])

from statsmodels.stats.multitest import multipletests
t_fdrp = pd.DataFrame([t_list, multipletests(pvals=p_list, alpha=0.05, method='fdr_bh')[1]],
                      columns=N.index, index=['Tstatistics', 'Padj']).T
p = len(t_fdrp[(t_fdrp['Tstatistics'] > 0) & (t_fdrp['Padj'] < 0.05)])
m = len(t_fdrp[(t_fdrp['Tstatistics'] < 0) & (t_fdrp['Padj'] < 0.05)])


print(str(len(N.index)) + '\t' + str(p) + '\t' + str(m))

sns.histplot(np.log10(deg_tn[deg_tn.index.isin(list(set(pls_dmr_info[pls_dmr_info.index.isin(total.index)]['GeneID'].values)))]['baseMean']))

p = sns.histplot(t_fdrp[t_fdrp['Padj'] < 0.05]['Tstatistics'])
p.set_xlabel('Paired t-test t-statistics (Tumor-Normal)')
plt.tight_layout()


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


pls_tumorup = list(set(pls_info[pls_info.index.isin(t_fdrp[(t_fdrp['Padj'] < 0.05) & (t_fdrp['Tstatistics'] > 0)].index)]['GeneID'].values)) # 570
pls_tumordown = list(set(pls_info[pls_info.index.isin(t_fdrp[(t_fdrp['Padj'] < 0.05) & (t_fdrp['Tstatistics'] < 0)].index)]['GeneID'].values)) # 35
a = intersection(pls_tumorup, deg_tn_uplist) # 9 ==> ETS2 (Proto oncogene)가 있네
b = intersection(pls_tumorup, deg_tn_downlist) # 78
c = intersection(pls_tumordown, deg_tn_uplist) # 30
d = intersection(pls_tumordown , deg_tn_downlist) # 17

print(len(a))
print(len(b))
print(len(c))
print(len(d))

e = pls_dmr_info[pls_dmr_info['GeneID'].isin(a)]
pls_e = pls[pls.index.isin(e[(e['K4me3'] == 'na') & (e['K27ac'] == 'na') & (e['K27me3'] == 'Yes')].index)]
rna_e = rna_norm_trans_log2[rna_norm_trans_log2.index.isin(list(e[(e['K4me3'] == 'na') & (e['K27ac'] == 'na') & (e['K27me3'] == 'Yes')]['ENSTID'].values))]

########################################################################################################################################################
# K4me3 Yes / K27ac No / K27me3 No
total = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'Yes') & (pls_dmr_info['K27ac'] == 'na') & (pls_dmr_info['K27me3'] == 'na')].index)]
N = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'Yes') & (pls_dmr_info['K27ac'] == 'na') & (pls_dmr_info['K27me3'] == 'na')].index)].iloc[:, :84]
T = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'Yes') & (pls_dmr_info['K27ac'] == 'na') & (pls_dmr_info['K27me3'] == 'na')].index)].iloc[:, 84:]

a = pd.DataFrame(N.mean(axis=1), columns=['Normal'])
b = pd.DataFrame(T.mean(axis=1), columns=['Tumor'])
ab = pd.concat([a,b], axis=1)
#p = sns.boxplot(data=ab, palette={'Normal':'midnightblue', 'Tumor':'darkred'}, width=0.5, showfliers = False)
#p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
#p.set_ylabel("Promoter Methylation")
#p.set_title("CpG islands present in Promoter or not")
p = sns.kdeplot(data=ab, x="Normal", y="Tumor", cmap="afmhot", shade=True, thresh=0)
p.set_xlabel("Normal promoter methylation (%)")
p.set_ylabel("Tumor promoter methylation (%)")
plt.tight_layout()
sns.despine()
plt.xlim((0,100))
plt.ylim((0,100))
stats.ttest_rel(a, b)

p_list = list()
t_list = list()
for i in list(range(len(N.index))):
    ttest = stats.ttest_rel(T.iloc[i], N.iloc[i])
    t_list.append(ttest[0])
    p_list.append(ttest[1])

from statsmodels.stats.multitest import multipletests
t_fdrp = pd.DataFrame([t_list, multipletests(pvals=p_list, alpha=0.05, method='fdr_bh')[1]],
                      columns=N.index, index=['Tstatistics', 'Padj']).T
p = len(t_fdrp[(t_fdrp['Tstatistics'] > 0) & (t_fdrp['Padj'] < 0.05)])
m = len(t_fdrp[(t_fdrp['Tstatistics'] < 0) & (t_fdrp['Padj'] < 0.05)])

print(str(len(N.index)) + '\t' + str(p) + '\t' + str(m))

sns.histplot(np.log10(deg_tn[deg_tn.index.isin(list(set(pls_dmr_info[pls_dmr_info.index.isin(total.index)]['GeneID'].values)))]['baseMean']))

p = sns.histplot(t_fdrp[t_fdrp['Padj'] < 0.05]['Tstatistics'])
p.set_xlabel('Paired t-test t-statistics (Tumor-Normal)')
plt.tight_layout()


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


pls_tumorup = list(set(pls_info[pls_info.index.isin(t_fdrp[(t_fdrp['Padj'] < 0.05) & (t_fdrp['Tstatistics'] > 0)].index)]['GeneID'].values)) # 570
pls_tumordown = list(set(pls_info[pls_info.index.isin(t_fdrp[(t_fdrp['Padj'] < 0.05) & (t_fdrp['Tstatistics'] < 0)].index)]['GeneID'].values)) # 35
a = intersection(pls_tumorup, deg_tn_uplist) # 9 ==> ETS2 (Proto oncogene)가 있네
b = intersection(pls_tumorup, deg_tn_downlist) # 78
c = intersection(pls_tumordown, deg_tn_uplist) # 30
d = intersection(pls_tumordown , deg_tn_downlist) # 17

print(len(a))
print(len(b))
print(len(c))
print(len(d))

e = pls_dmr_info[pls_dmr_info['GeneID'].isin(a)]
pls_e = pls[pls.index.isin(e[(e['K4me3'] == 'Yes') & (e['K27ac'] == 'na') & (e['K27me3'] == 'na')].index)]
rna_e = rna_norm_trans_log2[rna_norm_trans_log2.index.isin(list(e[(e['K4me3'] == 'Yes') & (e['K27ac'] == 'na') & (e['K27me3'] == 'na')]['ENSTID'].values))]

########################################################################################################################################################
# K4me3 No / K27ac Yes / K27me3 No
total = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'na') & (pls_dmr_info['K27ac'] == 'Yes') & (pls_dmr_info['K27me3'] == 'na')].index)]
N = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'na') & (pls_dmr_info['K27ac'] == 'Yes') & (pls_dmr_info['K27me3'] == 'na')].index)].iloc[:, :84]
T = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'na') & (pls_dmr_info['K27ac'] == 'Yes') & (pls_dmr_info['K27me3'] == 'na')].index)].iloc[:, 84:]

a = pd.DataFrame(N.mean(axis=1), columns=['Normal'])
b = pd.DataFrame(T.mean(axis=1), columns=['Tumor'])
ab = pd.concat([a,b], axis=1)
#p = sns.boxplot(data=ab, palette={'Normal':'midnightblue', 'Tumor':'darkred'}, width=0.5, showfliers = False)
#p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
#p.set_ylabel("Promoter Methylation")
#p.set_title("CpG islands present in Promoter or not")
p = sns.kdeplot(data=ab, x="Normal", y="Tumor", cmap="afmhot", shade=True, thresh=0, kde_kws={'clip': (0.0, 100)})
p.set_xlabel("Normal promoter methylation (%)")
p.set_ylabel("Tumor promoter methylation (%)")
plt.tight_layout()
sns.despine()
plt.xlim((0,100))
plt.ylim((0,100))
stats.ttest_rel(a, b)

p_list = list()
t_list = list()
for i in list(range(len(N.index))):
    ttest = stats.ttest_rel(T.iloc[i], N.iloc[i])
    t_list.append(ttest[0])
    p_list.append(ttest[1])

from statsmodels.stats.multitest import multipletests
t_fdrp = pd.DataFrame([t_list, multipletests(pvals=p_list, alpha=0.05, method='fdr_bh')[1]],
                      columns=N.index, index=['Tstatistics', 'Padj']).T
p = len(t_fdrp[(t_fdrp['Tstatistics'] > 0) & (t_fdrp['Padj'] < 0.05)])
m = len(t_fdrp[(t_fdrp['Tstatistics'] < 0) & (t_fdrp['Padj'] < 0.05)])


print(str(len(N.index)) + '\t' + str(p) + '\t' + str(m))

sns.histplot(np.log10(deg_tn[deg_tn.index.isin(list(set(pls_dmr_info[pls_dmr_info.index.isin(total.index)]['GeneID'].values)))]['baseMean']))

p = sns.histplot(t_fdrp[t_fdrp['Padj'] < 0.05]['Tstatistics'])
p.set_xlabel('Paired t-test t-statistics (Tumor-Normal)')
plt.tight_layout()


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


pls_tumorup = list(set(pls_info[pls_info.index.isin(t_fdrp[(t_fdrp['Padj'] < 0.05) & (t_fdrp['Tstatistics'] > 0)].index)]['GeneID'].values)) # 570
pls_tumordown = list(set(pls_info[pls_info.index.isin(t_fdrp[(t_fdrp['Padj'] < 0.05) & (t_fdrp['Tstatistics'] < 0)].index)]['GeneID'].values)) # 35
a = intersection(pls_tumorup, deg_tn_uplist) # 9 ==> ETS2 (Proto oncogene)가 있네
b = intersection(pls_tumorup, deg_tn_downlist) # 78
c = intersection(pls_tumordown, deg_tn_uplist) # 30
d = intersection(pls_tumordown , deg_tn_downlist) # 17

print(len(a))
print(len(b))
print(len(c))
print(len(d))

e = pls_dmr_info[pls_dmr_info['GeneID'].isin(a)]
pls_e = pls[pls.index.isin(e[(e['K4me3'] == 'na') & (e['K27ac'] == 'Yes') & (e['K27me3'] == 'na')].index)]
rna_e = rna_norm_trans_log2[rna_norm_trans_log2.index.isin(list(e[(e['K4me3'] == 'na') & (e['K27ac'] == 'Yes') & (e['K27me3'] == 'na')]['ENSTID'].values))]

########################################################################################################################################################
# K4me3 No / K27ac No / K27me3 No
total = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'na') & (pls_dmr_info['K27ac'] == 'na') & (pls_dmr_info['K27me3'] == 'na')].index)]
N = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'na') & (pls_dmr_info['K27ac'] == 'na') & (pls_dmr_info['K27me3'] == 'na')].index)].iloc[:, :84]
T = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'na') & (pls_dmr_info['K27ac'] == 'na') & (pls_dmr_info['K27me3'] == 'na')].index)].iloc[:, 84:]

a = pd.DataFrame(N.mean(axis=1), columns=['Normal'])
b = pd.DataFrame(T.mean(axis=1), columns=['Tumor'])
ab = pd.concat([a,b], axis=1)
#p = sns.boxplot(data=ab, palette={'Normal':'midnightblue', 'Tumor':'darkred'}, width=0.5, showfliers = False)
#p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
#p.set_ylabel("Promoter Methylation")
#p.set_title("CpG islands present in Promoter or not")
p = sns.kdeplot(data=ab, x="Normal", y="Tumor", cmap="afmhot", shade=True, thresh=0)
p.set_xlabel("Normal promoter methylation (%)")
p.set_ylabel("Tumor promoter methylation (%)")
plt.tight_layout()
sns.despine()
plt.xlim((0,100))
plt.ylim((0,100))
stats.ttest_rel(a, b)

p_list = list()
t_list = list()
for i in list(range(len(N.index))):
    ttest = stats.ttest_rel(T.iloc[i], N.iloc[i])
    t_list.append(ttest[0])
    p_list.append(ttest[1])

from statsmodels.stats.multitest import multipletests
t_fdrp = pd.DataFrame([t_list, multipletests(pvals=p_list, alpha=0.05, method='fdr_bh')[1]],
                      columns=N.index, index=['Tstatistics', 'Padj']).T
p = len(t_fdrp[(t_fdrp['Tstatistics'] > 0) & (t_fdrp['Padj'] < 0.05)])
m = len(t_fdrp[(t_fdrp['Tstatistics'] < 0) & (t_fdrp['Padj'] < 0.05)])


print(str(len(N.index)) + '\t' + str(p) + '\t' + str(m))


sns.histplot(np.log10(deg_tn[deg_tn.index.isin(list(set(pls_dmr_info[pls_dmr_info.index.isin(total.index)]['GeneID'].values)))]['baseMean']))

p = sns.histplot(t_fdrp[t_fdrp['Padj'] < 0.05]['Tstatistics'])
p.set_xlabel('Paired t-test t-statistics (Tumor-Normal)')
plt.tight_layout()


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


pls_tumorup = list(set(pls_info[pls_info.index.isin(t_fdrp[(t_fdrp['Padj'] < 0.05) & (t_fdrp['Tstatistics'] > 0)].index)]['GeneID'].values)) # 570
pls_tumordown = list(set(pls_info[pls_info.index.isin(t_fdrp[(t_fdrp['Padj'] < 0.05) & (t_fdrp['Tstatistics'] < 0)].index)]['GeneID'].values)) # 35
a = intersection(pls_tumorup, deg_tn_uplist) # 6 ==> ETS2 (Proto oncogene)가 있네
b = intersection(pls_tumorup, deg_tn_downlist) # 10
c = intersection(pls_tumordown, deg_tn_uplist) # 79
d = intersection(pls_tumordown , deg_tn_downlist) # 102

print(len(a))
print(len(b))
print(len(c))
print(len(d))

e = pls_dmr_info[pls_dmr_info['GeneID'].isin(d)]
pls_e = pls[pls.index.isin(e[(e['K4me3'] == 'na') & (e['K27ac'] == 'na') & (e['K27me3'] == 'na')].index)]
rna_e = rna_combat_tx_log2[rna_combat_tx_log2.index.isin(list(e[(e['K4me3'] == 'na') & (e['K27ac'] == 'na') & (e['K27me3'] == 'na')]['ENSTID'].values))]


#########################################################################################################################################################
total1 = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'Yes') & (pls_dmr_info['K27ac'] == 'Yes') & (pls_dmr_info['K27me3'] == 'Yes')].index)]
total2 = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'Yes') & (pls_dmr_info['K27ac'] == 'Yes') & (pls_dmr_info['K27me3'] == 'na')].index)]
total3 = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'Yes') & (pls_dmr_info['K27ac'] == 'na') & (pls_dmr_info['K27me3'] == 'Yes')].index)]
total4 = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'Yes') & (pls_dmr_info['K27ac'] == 'na') & (pls_dmr_info['K27me3'] == 'na')].index)]
total5 = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'na') & (pls_dmr_info['K27ac'] == 'Yes') & (pls_dmr_info['K27me3'] == 'na')].index)]
total6 = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'na') & (pls_dmr_info['K27ac'] == 'na') & (pls_dmr_info['K27me3'] == 'Yes')].index)]
total7 = pls[pls.index.isin(pls_dmr_info[(pls_dmr_info['K4me3'] == 'na') & (pls_dmr_info['K27ac'] == 'na') & (pls_dmr_info['K27me3'] == 'na')].index)]

p = sns.histplot(np.log10(deg_tn[deg_tn.index.isin(list(set(pls_dmr_info[pls_dmr_info.index.isin(total1.index)]['GeneID'].values)))]['baseMean']), kde=True, stat='count', color=sns.color_palette("Accent", 7)[0])
p = sns.histplot(np.log10(deg_tn[deg_tn.index.isin(list(set(pls_dmr_info[pls_dmr_info.index.isin(total2.index)]['GeneID'].values)))]['baseMean']), kde=True, stat='count', color=sns.color_palette("Accent", 7)[1])
p = sns.histplot(np.log10(deg_tn[deg_tn.index.isin(list(set(pls_dmr_info[pls_dmr_info.index.isin(total3.index)]['GeneID'].values)))]['baseMean']), kde=True, stat='count', color=sns.color_palette("Accent", 7)[2])
p = sns.histplot(np.log10(deg_tn[deg_tn.index.isin(list(set(pls_dmr_info[pls_dmr_info.index.isin(total4.index)]['GeneID'].values)))]['baseMean']), kde=True, stat='count', color=sns.color_palette("Accent", 7)[3])
p = sns.histplot(np.log10(deg_tn[deg_tn.index.isin(list(set(pls_dmr_info[pls_dmr_info.index.isin(total5.index)]['GeneID'].values)))]['baseMean']), kde=True, stat='count', color=sns.color_palette("Accent", 7)[4])
p = sns.histplot(np.log10(deg_tn[deg_tn.index.isin(list(set(pls_dmr_info[pls_dmr_info.index.isin(total6.index)]['GeneID'].values)))]['baseMean']), kde=True, stat='count', color=sns.color_palette("Accent", 7)[5])
p = sns.histplot(np.log10(deg_tn[deg_tn.index.isin(list(set(pls_dmr_info[pls_dmr_info.index.isin(total7.index)]['GeneID'].values)))]['baseMean']), kde=True, stat='count', color=sns.color_palette("Accent", 7)[6])
p.set_xlabel("Log10 mean normalized RNA read counts across samples")
sns.despine()


### TFBS enrichment
hypertf = pd.read_table("Enrichment_TFBS_HyperDMR.txt", index_col=0)
hypotf = pd.read_table("Enrichment_TFBS_HypoDMR.txt", index_col=0)
hypertf.sort_values(by="Enrichment", ascending=False)
hypotf.sort_values(by="Enrichment", ascending=False)
combined_tf = pd.concat([hypertf, hypotf], axis=1)
combined_tf.columns = ['Hyper DMR', 'Hypo DMR']

ax = combined_tf.sort_values(ascending=False, by='Hyper DMR').iloc[:10].plot.bar(rot=30, color = ['darkred', 'navy'], figsize=(15,8), fontsize=13, legend=True)
ax.legend(loc='upper left', bbox_to_anchor=(0.7, 0.9), fontsize=20)
ax.set_xticklabels(labels=list(combined_tf.sort_values(ascending=False, by='Hyper DMR').iloc[:10].index), fontstyle='italic')
plt.xlabel("")
plt.ylabel("Fold Enrichment (x)", fontsize=13)
plt.tight_layout()
sns.despine()

ax = combined_tf.sort_values(ascending=False, by='Hypo DMR').iloc[:10].plot.bar(rot=30, color = ['darkred', 'navy'], figsize=(15,8), fontsize=13, legend=True)
ax.legend(loc='upper left', bbox_to_anchor=(0.7, 0.9), fontsize=20)
ax.set_xticklabels(labels=list(combined_tf.sort_values(ascending=False, by='Hypo DMR').iloc[:10].index), fontstyle='italic')
plt.xlabel("")
plt.ylabel("Fold Enrichment (x)", fontsize=13)
plt.tight_layout()
sns.despine()

combined_tf_diff = pd.DataFrame(combined_tf['Hyper DMR'] - combined_tf['Hypo DMR'].values)
combined_tf_diff.columns = ['Difference']







### cCREs enrichment
hyperccres = pd.read_table("Enrichment_cCREs_HyperDMR.txt", index_col=0)
hypoccres = pd.read_table("Enrichment_cCREs_HypoDMR.txt", index_col=0)
hyperccres.sort_values(by="Enrichment", ascending=False)
hypoccres.sort_values(by="Enrichment", ascending=False)
combined_ccres = pd.concat([hyperccres, hypoccres ], axis=1)
combined_ccres.columns = ['Hyper DMR', 'Hypo DMR']

ax = combined_ccres.sort_index().plot.bar(rot=30, color = ['darkred', 'navy'], figsize=(15,8), fontsize=13, legend=True)
ax.legend(loc='upper left', bbox_to_anchor=(0.8, 1.0), fontsize=20)
ax.set_xticklabels(labels=list(combined_ccres.sort_index().index), fontstyle=None)
plt.xlabel("")
plt.ylabel("Fold Enrichment (x)", fontsize=13)
plt.tight_layout()
sns.despine()







## HOX Clusters DNA methylation

# Call HOX Cluster DNA methylation files

hox = pd.read_table("HOX_Clusters_ALL.txt", index_col=0)
hox_fillna0 = hox.fillna(0)
hoxa = hox_fillna0.iloc[100:210].copy()
hoxb = hox_fillna0.iloc[330:].copy()
hoxc = hox_fillna0.iloc[210:330].copy()
hoxd = hox_fillna0.iloc[:100].copy()
hox_fillna0_new = pd.concat([hoxa, hoxb, hoxc, hoxd])
sns.clustermap(hoxa,
                   method='ward',
                   metric='euclidean',
                   col_cluster=False,
                   row_cluster=False,
                   z_score=None,
                   standard_scale=0,
                   cmap='gnuplot2',
                   xticklabels=False,
                   yticklabels=False,
                   col_colors=[col_colors1],
                   row_colors=None,
                   cbar_kws={'label': 'DNA methylation'})


#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################
#############################################################################################################################################################
###### RNA-seq
rna = pd.read_table("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/STAD_SNUH_vst_new.txt", index_col=0, sep=' ')
rna.columns = list(map(lambda x: 'X'+x, rna.columns))
rna_norm = pd.read_table("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/STAD_SNUH_norm_counts.txt", index_col=0, sep=' ')
rna_norm.columns = list(map(lambda x: 'X'+x, rna_norm.columns))
rna_norm_log2 = np.log2(rna_norm + 1)
rna_norm_log2.columns = list(map(lambda x: 'X'+x, rna_norm_log2.columns))

rna = rna[rna.index.isin(list(map(lambda x: x.split('/')[2], pls.index)))]

rna_norm_trans = pd.read_table("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/STAD_SNUH_transcripts_norm_counts.txt", index_col=0, sep=' ')
rna_norm_trans_log2 = np.log2(rna_norm_trans + 1)
rna_norm_trans_log2.columns = list(map(lambda x: 'X'+x, rna_norm_trans_log2.columns))
rna_norm_trans_log2 = rna_norm_trans_log2[rna_norm_trans_log2.index.isin(list(map(lambda x: x.split('/')[4], pls.index)))]

rna_combat_tx = pd.read_table("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/STAD_SNUH_tx_combat_counts.txt", index_col=0, sep=' ')
rna_combat_tx_log2 = np.log2(rna_combat_tx + 1)
rna_combat_tx_log2.columns = list(map(lambda x: 'X'+x, rna_combat_tx_log2.columns))
rna_combat_tx_log2 = rna_combat_tx_log2[rna_combat_tx_log2.index.isin(list(map(lambda x: x.split('/')[4], pls.index)))]

# RNA batch info
rna_batch = pd.read_table("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/RNA_batch.txt", index_col=0)
col_colors_rna_batch = list(dict(zip(sorted(list(set(rna_batch.values.flatten()))), sns.color_palette("Dark2", 8)))[x] for x in rna_batch['RNA_batch'])

from matplotlib.patches import Patch
lut = dict(zip(sorted(list(set(rna_batch.values.flatten()))), sns.color_palette("Dark2", 8)))
handles = [Patch(facecolor=lut[name]) for name in lut]
# seaborn clustermap
plt.legend(handles, lut, title='RNA batch', bbox_to_anchor=(1,1), bbox_transform=plt.gcf().transFigure, loc='upper right')



txdb = pd.read_table("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/txdb_geneSymbol.txt", index_col=0)

deg_tn = pd.read_csv("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/Tumor.Normal.compare.csv", index_col=0)
deg_tn = deg_tn[deg_tn.index.isin(list(map(lambda x: x.split('/')[2], pls.index)))]

# Gene separation by RNA expression
p = sns.histplot(data=pd.DataFrame(np.log10(deg_tn['baseMean'])), x='baseMean', kde=True, stat='count', bins=100)
p.set_xticks([-1,0,1,2,3,4,5,6])
p.set_xlabel("Log10 mean normalized read counts across samples")
sns.despine()
plt.tight_layout()

# No expression (log10 baseMean < 1)
np.log10(deg_tn['baseMean']).describe(percentiles=[.20, .25, .50, .75, .90])
#count    23268.000000
#mean         1.631751
#std          1.210938
#min         -1.564640
#20%          0.476809
#25%          0.713531
#50%          1.918940
#75%          2.612718
#90%          2.983458
#max          6.937551
#Name: baseMean, dtype: float64

percentile_50 = np.log10(deg_tn['baseMean']).quantile([.5]).values[0]
percentile_90 = np.log10(deg_tn['baseMean']).quantile([.9]).values[0]

a = pd.DataFrame(pls_info[pls_info['GeneID'].isin(deg_tn[np.log10(deg_tn['baseMean']) < 1].index)].iloc[:,-8:]['K36me3'].value_counts())

b = pd.DataFrame(pls_info[pls_info['GeneID'].isin(deg_tn[(np.log10(deg_tn['baseMean']) >= 1) & (np.log10(deg_tn['baseMean']) < percentile_50)].index)].iloc[:-8:]['K36me3'].value_counts())

c = pd.DataFrame(pls_info[pls_info['GeneID'].isin(deg_tn[(np.log10(deg_tn['baseMean']) >= percentile_50) & (np.log10(deg_tn['baseMean']) < percentile_90)].index)].iloc[:-8:]['K36me3'].value_counts())

d = pd.DataFrame(pls_info[pls_info['GeneID'].isin(deg_tn[np.log10(deg_tn['baseMean']) >= percentile_90].index)].iloc[:,-8:]['K36me3'].value_counts())

abcd = pd.concat([a,b,c,d], axis=1)
abcd.columns = ['No expression', 'Low expression', 'Intermediate expression', 'High expression']

ax = abcd.plot.pie(subplots=True, colors = ['lightgrey', 'aqua'], figsize=(15,8), legend=False, counterclock=False, fontsize=13, autopct='%0.01f%%')

# For clustering purpose only
deg_tn_uplist = deg_tn[(deg_tn['baseMean'] >=10) & (deg_tn['padj'] < 0.01) & (deg_tn['log2FoldChange'] > 1)].index
deg_tn_downlist = deg_tn[(deg_tn['baseMean'] >=10) & (deg_tn['padj'] < 0.01) & (deg_tn['log2FoldChange'] < -1)].index

g = sns.clustermap(rna_norm_log2[rna_norm_log2.index.isin(list(deg_tn_uplist) + list(deg_tn_downlist))],
                   method='complete',
                   metric='correlation',
                   z_score=0,
                   standard_scale=None,
                   cmap=cmap3,
                    col_colors=[col_colors1],
                   xticklabels=False,
                    yticklabels=False, vmin=-3.5, vmax=3)

deg_tn_uplist = deg_tn[(deg_tn['baseMean'] >=10) & (deg_tn['padj'] < 0.005) & (deg_tn['log2FoldChange'] > 0)].index
deg_tn_downlist = deg_tn[(deg_tn['baseMean'] >=10) & (deg_tn['padj'] < 0.005) & (deg_tn['log2FoldChange'] < 0)].index


deg_da = pd.read_csv("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/Leiden_D.A.compare.csv", index_col=0)
deg_da = deg_da[deg_da.index.isin(list(map(lambda x: x.split('/')[2], pls.index)))]
deg_da_uplist = deg_da[(deg_da['baseMean'] >=10) & (deg_da['padj'] < 0.005) & (deg_da['log2FoldChange'] > 0)].index
deg_da_downlist = deg_da[(deg_da['baseMean'] >=10) & (deg_da['padj'] < 0.005) & (deg_da['log2FoldChange'] < 0)].index

### Gene-set enrichment analysis on Tumor - Normal

hallmark_up = gp.enrichr(gene_list=list(deg_tn_uplist),
			 gene_sets=['MSigDB_Hallmark_2020'],
			 organism='Human',
			 description='Hallmark_deg_tn_up',
			 outdir='./Hallmark_deg_tn_up',
			 cutoff=0.05,
			 no_plot=True)
gp.dotplot(hallmark_up.res2d, title='HALLMARK DEG UP', column='Adjusted P-value', cutoff=0.05, top_term=10, cmap='viridis_r', ofname='./Hallmark_deg_tn_up/Dotplot_HALLMARK_DEG_TN_UP.pdf')

hallmark_down = gp.enrichr(gene_list=list(deg_tn_downlist),
			 gene_sets=['MSigDB_Hallmark_2020'],
			 organism='Human',
			 description='Hallmark_deg_tn_down',
			 outdir='./Hallmark_deg_tn_down',
			 cutoff=0.05, no_plot=True)
gp.dotplot(hallmark_down.res2d, title='HALLMARK DEG DOWN', column='Adjusted P-value', cutoff=0.05, top_term=10, cmap='viridis_r', ofname='./Hallmark_deg_tn_down/Dotplot_HALLMARK_DEG_TN_DOWN.pdf')

gobp_up = gp.enrichr(gene_list=list(deg_tn_uplist),
		     gene_sets=['GO_Biological_Process_2021'],
		     organism='Human',
		     description='GOBP_deg_tn_up',
		     outdir='./GOBP_deg_tn_up',
		     cutoff=0.05, no_plot=True)
gp.dotplot(gobp_up.res2d, title='GOBP DEG UP', column='Adjusted P-value', cutoff=0.05, top_term=10, cmap='viridis_r', ofname='./GOBP_deg_tn_up/Dotplot_GOBP_DEG_TN_UP.pdf')

gobp_down = gp.enrichr(gene_list=list(deg_tn_downlist),
		       gene_sets=['GO_Biological_Process_2021'],
		       organism='Human',
		       description='GOBP_deg_tn_down',
		       outdir='./GOBP_deg_tn_down',
		       cutoff=0.05, no_plot=True)
gp.dotplot(gobp_down.res2d, title='GOBP DEG DOWN', column='Adjusted P-value', cutoff=0.05, top_term=10, cmap='viridis_r', ofname='./GOBP_deg_tn_down/Dotplot_GOBP_DEG_TN_DOWN.pdf')

kegg_up = gp.enrichr(gene_list=list(deg_tn_uplist),
		     gene_sets=['KEGG_2021_Human'],
		     organism='Human',
		     description='KEGG_deg_tn_up',
		     outdir='./KEGG_deg_tn_up',
		     cutoff=0.05, no_plot=True)
gp.dotplot(kegg_up.res2d, title='KEGG DEG UP', column='Adjusted P-value', cutoff=0.05, top_term=10, cmap='viridis_r', ofname='./KEGG_deg_tn_up/Dotplot_KEGG_DEG_TN_UP.pdf')

kegg_down = gp.enrichr(gene_list=list(deg_tn_downlist),
		       gene_sets=['KEGG_2021_Human'],
		       organism='Human',
		       description='KEGG_deg_tn_down',
		       outdir='./KEGG_deg_tn_down',
		       cutoff=0.05, no_plot=True)
gp.dotplot(kegg_down.res2d, title='KEGG DEG DOWN', column='Adjusted P-value', cutoff=0.05, top_term=10, cmap='viridis_r', ofname='./KEGG_deg_tn_down/Dotplot_KEGG_DEG_TN_DOWN.pdf')



### Gene-set enrichment analysis on Tumor - Normal

hallmark_up = gp.enrichr(gene_list=list(deg_da_uplist),
			 gene_sets=['MSigDB_Hallmark_2020'],
			 organism='Human',
			 description='Hallmark_deg_da_up',
			 outdir='./Hallmark_deg_da_up',
			 cutoff=0.05,
			 no_plot=True)
gp.dotplot(hallmark_up.res2d, title='HALLMARK DEG UP', column='Adjusted P-value', cutoff=0.05, top_term=10, cmap='viridis_r', ofname='./Hallmark_deg_da_up/Dotplot_HALLMARK_DEG_DA_UP.pdf')

hallmark_down = gp.enrichr(gene_list=list(deg_da_downlist),
			 gene_sets=['MSigDB_Hallmark_2020'],
			 organism='Human',
			 description='Hallmark_deg_da_down',
			 outdir='./Hallmark_deg_da_down',
			 cutoff=0.05, no_plot=True)
gp.dotplot(hallmark_down.res2d, title='HALLMARK DEG DOWN', column='Adjusted P-value', cutoff=0.05, top_term=10, cmap='viridis_r', ofname='./Hallmark_deg_da_down/Dotplot_HALLMARK_DEG_DA_DOWN.pdf')

gobp_up = gp.enrichr(gene_list=list(deg_da_uplist),
		     gene_sets=['GO_Biological_Process_2021'],
		     organism='Human',
		     description='GOBP_deg_da_up',
		     outdir='./GOBP_deg_da_up',
		     cutoff=0.05, no_plot=True)
gp.dotplot(gobp_up.res2d, title='GOBP DEG UP', column='Adjusted P-value', cutoff=0.05, top_term=10, cmap='viridis_r', ofname='./GOBP_deg_da_up/Dotplot_GOBP_DEG_DA_UP.pdf')

gobp_down = gp.enrichr(gene_list=list(deg_da_downlist),
		       gene_sets=['GO_Biological_Process_2021'],
		       organism='Human',
		       description='GOBP_deg_da_down',
		       outdir='./GOBP_deg_da_down',
		       cutoff=0.05, no_plot=True)
gp.dotplot(gobp_down.res2d, title='GOBP DEG DOWN', column='Adjusted P-value', cutoff=0.05, top_term=10, cmap='viridis_r', ofname='./GOBP_deg_da_down/Dotplot_GOBP_DEG_DA_DOWN.pdf')

kegg_up = gp.enrichr(gene_list=list(deg_da_uplist),
		     gene_sets=['KEGG_2021_Human'],
		     organism='Human',
		     description='KEGG_deg_da_up',
		     outdir='./KEGG_deg_da_up',
		     cutoff=0.05, no_plot=True)
gp.dotplot(kegg_up.res2d, title='KEGG DEG UP', column='Adjusted P-value', cutoff=0.05, top_term=10, cmap='viridis_r', ofname='./KEGG_deg_da_up/Dotplot_KEGG_DEG_DA_UP.pdf')

kegg_down = gp.enrichr(gene_list=list(deg_da_downlist),
		       gene_sets=['KEGG_2021_Human'],
		       organism='Human',
		       description='KEGG_deg_da_down',
		       outdir='./KEGG_deg_da_down',
		       cutoff=0.05, no_plot=True)
gp.dotplot(kegg_down.res2d, title='KEGG DEG DOWN', column='Adjusted P-value', cutoff=0.05, top_term=10, cmap='viridis_r', ofname='./KEGG_deg_da_down/Dotplot_KEGG_DEG_DA_DOWN.pdf')



### EMT index Tumor - Normal Difference

# (i) VST
emtindex = pd.read_table("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/EMT_index_vst_new.txt", index_col=0)
emtindex['EMT_index_diff'] = emtindex['Tumor'] - emtindex['Normal'].values
emtindex.columns = ['EMT_index_Normal', 'EMT_index_Tumor', 'EMT_index_Diff']
dmr_t.obs[['EMT_index_Normal', 'EMT_index_Tumor', 'EMT_index_Diff']] = emtindex

# (ii) Norm counts
emtindex = pd.read_table("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/EMT_index_norm_counts.txt", index_col=0)
emtindex['EMT_index_diff'] = emtindex['Tumor'] - emtindex['Normal'].values
emtindex.columns = ['EMT_index_Normal', 'EMT_index_Tumor', 'EMT_index_Diff']
dmr_t.obs[['EMT_index_Normal', 'EMT_index_Tumor', 'EMT_index_Diff']] = emtindex

# (iii) TPM
emtindex = pd.read_table("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/EMT_index_TPM.txt", index_col=0)
emtindex['EMT_index_diff'] = emtindex['Tumor'] - emtindex['Normal'].values
emtindex.columns = ['EMT_index_Normal', 'EMT_index_Tumor', 'EMT_index_Diff']
dmr_t.obs[['EMT_index_Normal', 'EMT_index_Tumor', 'EMT_index_Diff']] = emtindex

# Distribution of EMT index difference along DMR Clusters
# (i) Mean = 0, SD = 1 (zero mean, unit variance scaling)
from sklearn.preprocessing import StandardScaler
df = pd.DataFrame(dmr_t.obs['EMT_index_Diff'].copy())
df_scaled = StandardScaler().fit_transform(df.to_numpy())
df_scaled = pd.DataFrame(df_scaled, columns=['EMT_index_Diff_scaled'], index=dmr_t.obs.index)
dmr_t.obs[['EMT_index_Diff_scaled']] = df_scaled

p = sns.boxplot(data=dmr_t.obs, x='DMR Clusters2', y='EMT_index_Diff_scaled', palette=color_dict2, showfliers = False)
p = sns.stripplot(data=dmr_t.obs, x='DMR Clusters2', y='EMT_index_Diff_scaled', jitter=True, marker='o', color='black', size=1.5, alpha=0.5)
p.set_xlabel('DMR Clusters')
p.set_ylabel('EMT index difference (Tumor-Normal)')
plt.tight_layout()
sns.despine()

# (ii) Min-Max Standardization
df = pd.DataFrame(dmr_t.obs['EMT_index_Diff'].copy())
df_norm = (df - df.min()) / (df.max() - df.min())
dmr_t.obs[['EMT_index_Diff_norm']] = df_norm

p = sns.boxplot(data=dmr_t.obs, x='DMR Clusters2', y='EMT_index_Diff_norm', palette=color_dict, showfliers = False)
p = sns.stripplot(data=dmr_t.obs, x='DMR Clusters2', y='EMT_index_Diff_norm', jitter=True, marker='o', color='black', size=1.5, alpha=0.5)
p.set_xlabel('DMR Clusters')
p.set_ylabel('Normalized EMT index difference (Tumor-Normal)')
plt.tight_layout()
sns.despine()

# Kruskal-Wallis test for EMT index difference among DMR Clusters
df = dmr_t.obs[['DMR Clusters2', 'EMT_index_Diff_norm']]
stats.kruskal(df[df['DMR Clusters2'] == 'leiden_A']['EMT_index_Diff_norm'], df[df['DMR Clusters2'] == 'leiden_B']['EMT_index_Diff_norm'], df[df['DMR Clusters2'] == 'leiden_C']['EMT_index_Diff_norm'])

# ANOVA analysis for EMT index difference among DMR Clusters
df = dmr_t.obs[['leiden_r05', 'EMT_index_Diff_norm']]
df_lm = ols('EMT_index_Diff_norm ~ C(leiden_r05)', data=df).fit()
print(sm.stats.anova_lm(df_lm, typ=2))
print(pairwise_tukeyhsd(df['EMT_index_Diff_norm'], df['leiden_r05'], alpha=0.05))


# Association between Diffusion pseudotime and EMT index difference
stats.pearsonr(dmr_t.obs['EMT_index_Diff_norm'], dmr_t.obs['dpt_pseudotime'])
p = sns.lmplot(data=dmr_t.obs[['dpt_pseudotime', 'EMT_index_Diff_norm']], x='dpt_pseudotime', y='EMT_index_Diff_norm')
p.set_xlabels('Diffusion Pseudotime')
p.set_ylabels('Normalized EMT index difference (Tumor-Normal)')





sns.boxplot(data=dmr_t.obs[['EMT_index_Normal', 'EMT_index_Tumor']], palette={'EMT_index_Normal': 'navy', 'EMT_index_Tumor':'darkred'})






### Super Enhancer Overlapped Genes
se_genes = pd.read_table("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/SE_DMR_abs15_hyper_overlapped_genes.txt")
se_genes = list(se_genes.values.flatten()) # 100개

deg_tn[deg_tn.index.isin(se_genes)] # ==> 70개

# deg_tn_uplist[deg_tn_uplist.isin(se_genes)] ==> 10개
# deg_tn_downlist[deg_tn_downlist.isin(se_genes)] ==> 20개
se_deg = list( deg_tn_uplist[deg_tn_uplist.isin(se_genes)] ) + list( deg_tn_downlist[deg_tn_downlist.isin(se_genes)] )

g = sns.clustermap(rna_norm_log2[rna_norm_log2.index.isin( se_deg )],
                   col_cluster=False,
                    method='ward',
                   metric='euclidean',
                   z_score=0,
                   standard_scale=None,
                   cmap=cmap3,
                    col_colors=[col_colors1],
                   xticklabels=False,
                    yticklabels=True, vmin=-2.5, vmax=2.5)
g.ax_heatmap.set_yticklabels(labels=g.ax_heatmap.get_yticklabels(), fontstyle='italic')

se_deg_rna_roworder = g.dendrogram_row.reordered_ind

#  Array of codes for making SE_DEG_ALL.txt (refer to /mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/DMR/DMR_min55_new/Enrichment/Stomach_SE)
reorder = list(rna_norm_log2[rna_norm_log2.index.isin( se_deg )].iloc[se_deg_rna_roworder].index)

se_deg_met = pd.read_table("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/03.WGBS/NEW/SE_DEG_ALL.txt", index_col=0)
se_deg_met.columns = list(map(lambda x: 'X'+x, se_deg_met.columns))

a = pd.DataFrame(list(map(lambda x: x.split('/')[0], se_deg_met.index)), columns=['Region'], index=se_deg_met.index)
b = pd.DataFrame(list(map(lambda x: x.split('/')[1], se_deg_met.index)), columns=['GeneID'], index=se_deg_met.index)
c = pd.DataFrame(list(map(lambda x: x.split('/')[2], se_deg_met.index)), columns=['CpG'], index=se_deg_met.index)
se_deg_met_info = pd.concat([a,b,c], axis=1)

reorder_iloc = list()

for i in reorder:
    reorder_iloc.append(list(se_deg_met_info['GeneID'].values).index(i))

# DNA methylation of Super enhancers overlapped with DEG
g = sns.clustermap(se_deg_met.iloc[reorder_iloc],
                   col_cluster=False,
                   row_cluster=False,
                   cmap='Spectral_r',
                   z_score=None,
                   standard_scale=0,
                   col_colors=[col_colors1],
                   xticklabels=False,
                   yticklabels=False)



normal_cpgi = cpgi_met.iloc[:, :84].mean()
tumor_cpgi = cpgi_met.iloc[:, 84:].mean()

pmd_met = pd.read_table("PMD_ALL.txt", index_col=0)
pmd_met.columns = list(map(lambda x: 'X'+x, pmd_met.columns))
normal_pmd = pmd_met.iloc[:, :84].mean()
tumor_pmd = pmd_met.iloc[:, 84:].mean()

a = pd.concat([normal_cpgi, normal_pmd, pd.Series(['Normal']*84, index=pmd_met.columns[:84])], axis=1)
a.columns = ['CpGi', 'PMD', 'Type']
b = pd.concat([tumor_cpgi, tumor_pmd, dmr_t.obs['CIMP']], axis=1)
b.columns = ['CpGi', 'PMD', 'Type']

c = pd.concat([a,b])

ax = sns.scatterplot(data=c, x='CpGi', y='PMD', hue='Type', linewidth=0, palette={'Normal': 'darkblue', 'CIMP(-)': 'salmon', 'CIMP(+)': 'maroon'}, s=50)
ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1.0))
sns.despine()
plt.tight_layout()
ax.set_xlabel('CIMP-CGI DNA methylation (%)')
ax.set_ylabel('PMD DNA methylation (%)')











# Imputed DMR
dmr_met = pd.read_csv("DMR_abs15_ALL_imputed_corrected.csv", index_col=0).iloc[:,:-1].T

# Row name change
dmr_met.index = list(map(lambda x: x.split('.')[0] + ':' + x.split('.')[1] + '-' + x.split('.')[2], dmr_met.index))
dmr_tn = dmr_met.iloc[:,84:] - dmr_met.iloc[:,:84].values

# DNA methylation of Tumor-Normal
dmr_tn = sc.AnnData(dmr_tn.T)

# Percentage methylation stored on raw and layer attribute
dmr_tn.raw = dmr_tn
dmr_tn.layers['TN_met'] = dmr_tn.X

# clinic info attachment
dmr_tn.obs = clinic_info.iloc[84:].copy()

# Scaling and PCA
sc.pp.scale(dmr_tn)
sc.tl.pca(dmr_tn, n_comps=83, zero_center=True)

# Check for PCA numbers
# pca_variance_tn = pd.DataFrame(dmr_tn.uns['pca']['variance_ratio'], index=list(map(lambda x: 'PC' + str(x), list(range(1,84)))), columns=['Variance_ratio'])
# print( np.sum(pca_variance_tn.values.flatten()[:16]) )
# PC1 ~ 16 ==> 70.98%

# Neighbor Graph construction and leiden community detection
sc.pp.neighbors(dmr_tn, n_neighbors=10, n_pcs=16)
sc.tl.leiden(dmr_tn, resolution=0.75, key_added='leiden_r075')
sc.tl.leiden(dmr_tn, resolution=0.5, key_added='leiden_r05')
sc.tl.umap(dmr_tn, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
sc.pl.umap(dmr_tn, color=['leiden_r075', 'leiden_r05', 'EpiBurden'], add_outline=False, legend_loc='right margin', color_map=cmap)

# Leiden cluster name change
leiden_name_change_dict1 = {'3': 'leiden_A',
                           '1': 'leiden_B',
                           '0': 'leiden_C',
                           '2': 'leiden_D'} # leiden_r075

leiden_name_change_dict2 = {'1': 'leiden_A',
                           '0': 'leiden_B',
                           '2': 'leiden_C'} # leiden_r05

# Leiden cluster name change: leiden_r075 ==> DMR_leiden
dmr_t.obs['DMR Clusters'] = dmr_t.obs['leiden_r075'].map(lambda x: leiden_name_change_dict1[x]).astype('category')

# UMAP projection plot for DMR Clusters and Epimutation Burden
dmr_t.obs['Epimutation Burden'] = dmr_t.obs['EpiBurden']
sc.pl.umap(dmr_t, color=['DMR Clusters', 'Epimutation Burden'], add_outline=False, legend_loc='on data', palette='Dark2')
sns.despine()





