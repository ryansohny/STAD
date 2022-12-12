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
%matplotlib
%autoindent

clinic_info = pd.read_csv('/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/2022_WC300_clinical_information_Xadded_ver2.0.csv', index_col='Sample')

# Average methylation for Tumor vs Normal (Figure 1)
p = sns.violinplot(data=clinic_info.iloc[84:][['PercentMet_COV5_Normal', 'PercentMet_COV5_Tumor']], palette={'PercentMet_COV5_Normal':'midnightblue', 'PercentMet_COV5_Tumor':'darkred'}, cut=0, scale="count")
p = sns.stripplot(data=clinic_info.iloc[84:][['PercentMet_COV5_Normal', 'PercentMet_COV5_Tumor']], color="black")
p.set_xticklabels(['Normal (N=84)', 'Tumor (N=84)'])
p.set_ylabel("Average methylation (%)")
sns.despine()

# Partially Methylated Domains (PMDs)
pmd_met = pd.read_table("PMD_ALL.txt", index_col=0)
pmd_met.columns = list(map(lambda x: 'X'+x, pmd_met.columns))
normal_pmd = pmd_met.iloc[:, :84].mean()
tumor_pmd = pmd_met.iloc[:, 84:].mean()


## Determining CIMP Tumors

# Call CpGi methylation
cpgi_met = pd.read_table("CpGi_smooth.txt", index_col=0)
cpgi_met = cpgi_met * 100

# Column name change
cpgi_met.columns = list(map(lambda x: 'X'+x, cpgi_met.columns))

# Discard missing CpGi DNA methylation rows & pick CpGi sites where Normal DNA methylation < 40
cpgi_met = cpgi_met[cpgi_met.iloc[:, :84].mean(axis=1) < 40] # (i)
#cpgi_met = cpgi_met[cpgi_met.iloc[:, :84].mean(axis=1) < 20] # (ii)
# Call Promoter CpGi
cpgi_pls = list(map(lambda x: x.strip('\n').split('/')[0], open("PLS_CpGi.txt", 'r').readlines()))

# Select Promoter CpGi
cpgi_pls_met = cpgi_met[cpgi_met.index.isin(cpgi_pls)]

# mean(Tumor - Normal) >= 10%
cpgi_pls_tn_met = cpgi_pls_met.iloc[:, 84:] - cpgi_pls_met.iloc[:, :84].values
cpgi_pls_tn_met = cpgi_pls_tn_met[cpgi_pls_tn_met.mean(axis=1) >= 10]

# Hierarchical clustering

g = sns.clustermap(cpgi_pls_tn_met,
                   method='ward',
                   metric='euclidean',
                   z_score=None,
                   standard_scale=None,
                   cmap='RdYlBu_r',
                   xticklabels=False,
                   yticklabels=False,
                   col_colors=None)


cimp_positive_samples = list(cpgi_pls_tn_met.iloc[:, g.dendrogram_col.reordered_ind[:33]].columns) # (i)
#cimp_positive_samples = list(cpgi_pls_tn_met.iloc[:, g.dendrogram_col.reordered_ind[:32]].columns) # (i)

total_sample_cimp_info = list(map(lambda x: 'CIMP(+) tumor (N=33)' if x in cimp_positive_samples else ('Normal (N = 84)' if x[-1] == 'N' else 'CIMP(-) tumor (N = 51)'), cpgi_pls_met.columns))
total_sample_cimp_info = pd.Series(dict(zip(list(cpgi_pls_met.columns), total_sample_cimp_info)))

# Association between PMD methylation and CIMP-CGI DNA methylation

cimp_pmd = pd.concat([cpgi_pls_met[cpgi_pls_met.index.isin(cpgi_pls_tn_met.index)].mean(), pmd_met.mean(), total_sample_cimp_info], axis=1)
cimp_pmd.columns = ['CpGimet', 'PMDmet', 'CIMPtype']
rho = round(stats.spearmanr(cimp_pmd['CpGimet'], cimp_pmd['PMDmet'])[0], 3)

fig, ax = plt.subplots(1,1, figsize=(7,7))
sns.scatterplot(data=cimp_pmd, x='CpGimet', y='PMDmet', hue='CIMPtype', linewidth=0, palette={'Normal (N = 84)': 'darkblue', 'CIMP(-) tumor (N = 51)': 'salmon', 'CIMP(+) tumor (N=33)': 'maroon'}, s=50, ax=ax)
handles, labels = ax.get_legend_handles_labels()
order = [0, 2, 1]
ax.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='upper left', bbox_to_anchor=(0.01, 0.18), frameon=True, edgecolor='black', fancybox=False)
ax.set_xlabel('CIMP-CGI DNA methylation (%)')
ax.set_ylabel('PMD DNA methylation (%)')
ax.text(0.47, 0.1, f"Spearman's Rho: {rho}", size=11, weight='bold', transform=ax.transAxes)
ax.set_xlim((0, 75))
ax.set_ylim((48, 83))
sns.despine(ax=ax)

g = sns.lmplot(data=cimp_pmd, x='CpGimet', y='PMDmet', hue='CIMPtype', palette={'Normal (N = 84)': 'darkblue', 'CIMP(-) tumor (N = 51)': 'salmon', 'CIMP(+) tumor (N=33)': 'maroon'})
g.ax.set_xlabel('CIMP-CGI DNA methylation (%)')
g.ax.set_ylabel('PMD DNA methylation (%)')
g.ax.set_xlim((0, 75))
g.ax.set_ylim((48, 83))
g.ax.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='upper left', bbox_to_anchor=(0.01, 0.212), frameon=True, edgecolor='black', fancybox=False)

# CIMP proportional plot
#ax = (pd.crosstab(dmr_t.obs['DMR Clusters'], dmr_t.obs['CIMP'], normalize=0)*100).plot.bar(stacked=True, color=['#8b0000ff', '#000080ff'], rot=0)
#plt.ylabel("Proportion (%)")
#ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1.0))
#plt.tight_layout()
#sns.despine()


# DMR
dmr_met = pd.read_table("DMR_abs10_smooth.txt", index_col=0)
dmr_met.columns = list(map(lambda x: 'X'+x, dmr_met.columns))
dmr_met = dmr_met*100

dmr_info = pd.read_table("DMR_abs10_Hyper-Hypo_annotation.txt", index_col=0)

# Remove PMD-overlapped Hypo-DMR
with open("DMR_abs10_hypo_wPMD.index", 'r') as dfh:
    hypodmr_wPMD = list(map(lambda x: x.strip(), dfh.readlines()))

dmr_met = dmr_met[~dmr_met.index.isin(hypodmr_wPMD)]
dmr_info = dmr_info[~dmr_info.index.isin(hypodmr_wPMD)]

del hypodmr_wPMD

col_colors_tn = ['#C0C0C0']*84 + ['#000000']*84
col_colors_cimp = list(dict(zip(['Normal (N = 84)', 'CIMP(-) tumor (N = 51)', 'CIMP(+) tumor (N=33)'], ['darkblue', 'salmon', 'maroon']))[x] for x in cimp_pmd['CIMPtype'])
row_colors_dmr = list(dict(zip(['Hypo', 'Hyper'], ['#6C8EAD', '#A23E48']))[x] for x in dmr_info['Type'])


g = sns.clustermap(dmr_met,
                   method='ward',
                   metric='euclidean',
                   z_score=None,
                   standard_scale=None,
                   cmap='RdYlBu_r',
                   xticklabels=False,
                   yticklabels=False,
                   vmin=0,
                   vmax=100,
                   col_colors=[col_colors_tn, col_colors_cimp],
                   row_colors=row_colors_dmr) # standard scale: 0 (rows) or 1 columns (subtract min and divide by max)
g.ax_heatmap.set_ylabel('')

# Only Tumor    
g = sns.clustermap(dmr_met,
                   method='ward',
                   metric='euclidean',
                   z_score=None,
                   standard_scale=None,
                   cmap='RdYlBu_r',
                   xticklabels=False,
                   yticklabels=False,
                   col_colors=['#C0C0C0']*84 + ['#000000']*84,
                   row_colors=row_colors_dmr)
g.ax_heatmap.set_ylabel('')

# Normalized CpG density to DMR annotation table
dmr_info['Norm_CpGdensity'] = (dmr_info['CpGdensity'] - np.min(dmr_info['CpGdensity'])) / (np.max(dmr_info['CpGdensity']) - np.min(dmr_info['CpGdensity']))
print(stats.ttest_ind(dmr_info[dmr_info['Type'] == 'Hyper']['Norm_CpGdensity'], dmr_info[dmr_info['Type'] == 'Hypo']['Norm_CpGdensity'], equal_var=False))
#Ttest_indResult(statistic=34.445773291849996, pvalue=1.1938040649902871e-228)

# CpG density kdeplot between Hyper-DMR and Hypo-DMR
fig, ax = plt.subplots(1,1, figsize=(7,7))
sns.kdeplot(dmr_info[dmr_info['Type'] == 'Hypo']['Norm_CpGdensity'], color='#6C8EAD', fill=True, ax=ax)
sns.kdeplot(dmr_info[dmr_info['Type'] == 'Hyper']['Norm_CpGdensity'], color='#A23E48', fill=True, ax=ax)
ax.set_xlim((0, 1))
ax.set_ylim((0, 10))
sns.despine(ax=ax)


# DMR annotation colormap for sns.clustermap
row_colors_dmr1 = list(dict(zip(['Hypo', 'Hyper'], ['#6C8EAD', '#A23E48']))[x] for x in dmr_info['Type'])
row_colors_dmr2 = list(dict(zip(['Hypo', 'Hyper'], ['#6C8EAD', '#A23E48']))[x] for x in dmr_info2['Type'])
col_colors1 = ['#C0C0C0']*84 + ['#000000']*84
# Clustering 1
g = sns.clustermap(dmr_met.loc[dmr_info2.index],
                   method='ward',
                   metric='euclidean',
                   z_score=None,
                   standard_scale=None,
                   cmap=cmap,
                   robust=True,
                   col_colors=[col_colors1],
                   row_colors=[row_colors_dmr2],
                   xticklabels=False,
                   yticklabels=False,
                   cbar_kws={'label': 'DNA methylation'})
g.cax.set_visible(False) # Legend removal



# LOLA Functional Enrichment 
input_dir="/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/03.WGBS/NEW/Functional_Enrichment/LOLA/"
lola_dmr = pd.read_table(f"{input_dir}LOLA_2_allEnrichments.tsv")
lola_hyper_dmr = lola_dmr[lola_dmr['userSet'] == 1]
lola_hypo_dmr = lola_dmr[lola_dmr['userSet'] == 2]







# DMR UMAP Projection

dmr_met_adata = sc.AnnData(dmr_met.T)

# Percentage methylation stored on raw and layer attribute
dmr_met_adata.raw = dmr_met_adata
dmr_met_adata.layers['Percent_met'] = dmr_met_adata.X.copy()

# clinic info attachment
dmr_met_adata.obs = clinic_info.copy()

# RNA batch info
rna_batch = pd.read_table("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/RNA_batch.txt", index_col=0)
dmr_met_adata.obs['RNA_batch'] = rna_batch.copy()
dmr_met_adata.obs['CIMP'] = cimp_pmd['CIMPtype']

# Scaling (optional) and PCA
sc.pp.scale(dmr_met_adata) # z-score scaling, which is (X-mean)/std. So the mean of the each variable becomes 0 (almost zero) and gets to have a unit variance.
sc.tl.pca(dmr_met_adata, n_comps=83, zero_center=True) # zero_center=True => compute standard PCA from covariance matrix
sc.pl.pca(dmr_met_adata, color=['TN', 'EpiBurden', 'CIMP'], add_outline=False, legend_loc='right margin', color_map=cmap, use_raw=True, annotate_var_explained=True, size=100, components=['1,2'])

# Check for PCA numbers
pca_variance = pd.DataFrame(dmr_met_adata.uns['pca']['variance_ratio'], index=list(map(lambda x: 'PC' + str(x), list(range(1,84)))), columns=['Variance_ratio'])
np.sum(pca_variance.values.flatten()[:12])

# Neighbor Graph construction and leiden community detection
sc.pp.neighbors(dmr_met_adata, n_neighbors=10, n_pcs=9)
sc.tl.umap(dmr_met_adata, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
sc.pl.umap(dmr_met_adata, color=['TN', 'EpiBurden', 'CIMP'], add_outline=False, legend_loc='right margin', color_map=cmap)


sc.tl.leiden(dmr_met_adata, resolution=0.75, key_added='leiden_r075')
sc.tl.leiden(dmr_met_adata, resolution=0.5, key_added='leiden_r05')
sc.tl.umap(dmr_met_adata, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
sc.pl.umap(dmr_met_adata, color=['leiden_r075', 'leiden_r05', 'EpiBurden'], add_outline=False, legend_loc='right margin', color_map=cmap)