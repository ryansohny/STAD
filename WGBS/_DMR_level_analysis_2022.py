# source activate wc300
# ipython --profile=wc300

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
sns.set(font="Arial", font_scale=1.2, style='ticks')
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])
cmap4 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#191970", "#ffdab9", "#8B0000"])
%matplotlib
%autoindent

# Clinical information
clinic_info = pd.read_csv('/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/2022_WC300_clinical_information_Xadded.csv', index_col='Sample')
#clinic_info = pd.read_csv('/home/mhryan/Workspace/02.Projects/02.WC300/2022_WC300_clinical_information_Xadded.csv', index_col='Sample')

# Imputed DMR
dmr_met = pd.read_table("DMR_abs10_smooth.txt", index_col=0)
dmr_met.columns = list(map(lambda x: 'X'+x, dmr_met.columns))

# DMR annotation table
dmr_info2 = pd.read_table("DMR_abs15_Hyper-Hypo_annotation.txt", index_col=0)
dmr_info = pd.read_table("DMR_abs10_Hyper-Hypo_annotation.txt", index_col=0)

# Normalized CpG density to DMR annotation table
dmr_info['Norm_CpGdensity'] = (dmr_info['CpGdensity'] - np.min(dmr_info['CpGdensity'])) / (np.max(dmr_info['CpGdensity']) - np.min(dmr_info['CpGdensity']))
print(stats.ttest_ind(dmr_info[dmr_info['Type'] == 'Hyper']['Norm_CpGdensity'], dmr_info[dmr_info['Type'] == 'Hypo']['Norm_CpGdensity'], equal_var=False))
#Ttest_indResult(statistic=34.445773291849996, pvalue=1.1938040649902871e-228)

# CpG density boxplot between Hyper-DMR and Hypo-DMR 
p = sns.boxplot(data=dmr_info, x='Type', y='Norm_CpGdensity', order=['Hyper', 'Hypo'], palette={'Hyper': '#A23E48', 'Hypo': '#6C8EAD'}, showfliers = False)
#p = sns.violinplot(data=dmr_info, x='Type', y='Norm_CpGdensity', order=['Hyper', 'Hypo'], palette={'Hyper': '#A23E48', 'Hypo': '#6C8EAD'}, cut=0, scale='area')
p = sns.stripplot(data=dmr_info, x='Type', y='Norm_CpGdensity', order=['Hyper', 'Hypo'], jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
p.set_ylabel("Normalized CpG density of DMR")
p.set_xlabel("Types of DMR")
p.set_xticklabels(['Hypermethylation', 'Hypomethylation'])
plt.tight_layout()
sns.despine()

# CpG density kdeplot between Hyper-DMR and Hypo-DMR
ax = plt.subplot(1,1,1)
sns.kdeplot(dmr_info[dmr_info['Type'] == 'Hypo']['Norm_CpGdensity'], color='#6C8EAD', fill=True, ax=ax)
sns.kdeplot(dmr_info[dmr_info['Type'] == 'Hyper']['Norm_CpGdensity'], color='#A23E48', fill=True, ax=ax)
ax.set_xlabel("Min-Max Normalized CpG Density of DMR")
plt.xlim((0, 1))
plt.ylim((0, 20))
sns.despine()
plt.tight_layout()

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

# Clustering 2
g = sns.clustermap(dmr_met.loc[dmr_info2[dmr_info2['Type'] == 'Hyper'].index],
                   method='ward',
                   metric='euclidean',
                   z_score=None,
                   standard_scale=None,
                   cmap=cmap,
                   robust=True,
                   col_colors=[col_colors1],
                   row_colors=None,
                   xticklabels=False,
                   yticklabels=False,
                   cbar_kws={'label': 'DNA methylation'})
g.cax.set_visible(False) # Legend removal

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



# Call Smooth DMR for all samples
dmr = sc.AnnData(dmr_met.T)
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

