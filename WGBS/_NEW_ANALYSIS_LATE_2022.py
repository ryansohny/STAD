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
#cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#104e8b", "#ffdab9", "#8b0a50"])
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["darkblue", "#ffdab9", "maroon"])
%matplotlib
%autoindent

clinic_info = pd.read_csv('/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/2022_WC300_clinical_information_Xadded_ver2.0.csv', index_col='Sample')

pmd_frac = pd.read_table('PMD_fraction.txt', index_col=0)
pmd_frac.index = list(map(lambda x: 'X'+x, pmd_frac.index))
clinic_info['PMD_Fraction'] = pmd_frac['PMD_Fraction'].copy()

# Average methylation for Tumor vs Normal (Figure 1A)
p = sns.violinplot(data=clinic_info.iloc[84:][['PercentMet_COV5_Normal', 'PercentMet_COV5_Tumor']], palette={'PercentMet_COV5_Normal':'midnightblue', 'PercentMet_COV5_Tumor':'darkred'}, cut=0, scale="area")
p = sns.stripplot(data=clinic_info.iloc[84:][['PercentMet_COV5_Normal', 'PercentMet_COV5_Tumor']], color="black")
p.set_xticklabels(['Normal (N=84)', 'Tumor (N=84)'])
p.set_ylabel("Average methylation (%)")
sns.despine()

# PMD Fraction for Tumor vs Normal (Figure 1B)
p = sns.violinplot(data=clinic_info, x='TN', y='PMD_Fraction', palette={'Normal':'midnightblue', 'Tumor':'darkred'}, cut=0, scale="width")
p = sns.stripplot(data=clinic_info, x='TN', y='PMD_Fraction', color="black")
p.set_xticklabels(['Normal (N=84)', 'Tumor (N=84)'])
p.set_ylabel("Fraction of PMDs in the Genome")
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
g.ax_heatmap.set_ylabel('')

cimp_positive_samples = list(cpgi_pls_tn_met.iloc[:, g.dendrogram_col.reordered_ind[:33]].columns) # (i)
#cimp_positive_samples = list(cpgi_pls_tn_met.iloc[:, g.dendrogram_col.reordered_ind[:32]].columns) # (i)

total_sample_cimp_info = list(map(lambda x: 'CIMP(+) tumor (N=33)' if x in cimp_positive_samples else ('Normal (N=84)' if x[-1] == 'N' else 'CIMP(-) tumor (N=51)'), cpgi_pls_met.columns))
total_sample_cimp_info = pd.Series(dict(zip(list(cpgi_pls_met.columns), total_sample_cimp_info)))
clinic_info['CIMP'] = total_sample_cimp_info.values

# Association between PMD methylation and CIMP-CGI DNA methylation

cimp_pmd = pd.concat([cpgi_pls_met[cpgi_pls_met.index.isin(cpgi_pls_tn_met.index)].mean(), pmd_met.mean(), total_sample_cimp_info], axis=1)
cimp_pmd.columns = ['CpGimet', 'PMDmet', 'CIMPtype']
rho = round(stats.spearmanr(cimp_pmd['CpGimet'], cimp_pmd['PMDmet'])[0], 3)

fig, ax = plt.subplots(1,1, figsize=(7,7))
sns.scatterplot(data=cimp_pmd, x='CpGimet', y='PMDmet', hue='CIMPtype', linewidth=0, palette={'Normal (N=84)': 'darkblue', 'CIMP(-) tumor (N=51)': 'salmon', 'CIMP(+) tumor (N=33)': 'maroon'}, s=50, ax=ax)
handles, labels = ax.get_legend_handles_labels()
order = [0, 2, 1]
ax.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='upper left', bbox_to_anchor=(0.01, 0.18), frameon=True, edgecolor='black', fancybox=False)
ax.set_xlabel('CIMP-CGI DNA methylation (%)')
ax.set_ylabel('PMD DNA methylation (%)')
ax.text(0.47, 0.1, f"Spearman's Rho: {rho}", size=11, weight='bold', transform=ax.transAxes)
ax.set_xlim((0, 75))
ax.set_ylim((48, 83))
sns.despine(ax=ax)

g = sns.lmplot(data=cimp_pmd, x='CpGimet', y='PMDmet', hue='CIMPtype', palette={'Normal (N=84)': 'darkblue', 'CIMP(-) tumor (N=51)': 'salmon', 'CIMP(+) tumor (N=33)': 'maroon'})
g.ax.set_xlabel('CIMP-CGI DNA methylation (%)')
g.ax.set_ylabel('PMD DNA methylation (%)')
g.ax.set_xlim((0, 75))
g.ax.set_ylim((48, 83))
g.ax.legend([handles[idx] for idx in order], [labels[idx] for idx in order], loc='upper left', bbox_to_anchor=(0.01, 0.212), frameon=True, edgecolor='black', fancybox=False)

# 
g = sns.lmplot(data=clinic_info.iloc[84:, :], x="Age", y="EpiBurden", hue="CIMP", palette={'CIMP(+) tumor (N=33)': 'maroon', 'CIMP(-) tumor (N=51)': 'salmon'})
g.ax.set_xlabel('Age')
g.ax.set_ylabel('Epimutation Burden')



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
'''
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
'''

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
ax.set_xlabel("Min-Max Normalized CpG Density of DMR")
sns.despine(ax=ax)

# CpG density histplot between Hyper-DMR and Hypo-DMR - ver2
fig, ax = plt.subplots(1, 1, figsize=(7, 7))
sns.histplot(dmr_info[dmr_info['Type'] == 'Hypo']['Norm_CpGdensity'], kde=True, stat='count', color='#6C8EAD', ax=ax)
sns.histplot(dmr_info[dmr_info['Type'] == 'Hyper']['Norm_CpGdensity'], kde=True, stat='count', color='#A23E48', ax=ax)
ax.set_xlim((0, 1))


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
lola_hyper_dmr = lola_dmr[(lola_dmr['userSet'] == 1) & (lola_dmr['collection'] == 'encode_tfbs')]
lola_hypo_dmr = lola_dmr[(lola_dmr['userSet'] == 2) & (lola_dmr['collection'] == 'encode_tfbs')]

lola_hyper_dmr['Antibody from ENCODE cell type'] = lola_hyper_dmr[['antibody', 'cellType']].apply(lambda x: ' from '.join(x), axis=1)
lola_hypo_dmr['Antibody from ENCODE cell type'] = lola_hypo_dmr[['antibody', 'cellType']].apply(lambda x: ' from '.join(x), axis=1)


lola_hypo_dmr.loc[:, ('antibody', 'cellType', 'treatment')].apply(lambda x: ', '.join(x), axis=1)
test = pd.concat([lola_hyper_dmr.iloc[:10, :], lola_hypo_dmr.iloc[:10, :]])
fig, ax = plt.subplots(1, 1, figsize=(10,6), constrained_layout=True)
sns.scatterplot(data=test.sort_values(by='pValueLog', ascending=False), x="pValueLog", y="Antibody from ENCODE cell type", hue='pValueLog', size="oddsRatio", style='userSet', sizes=(20,200), ax=ax)





# DMR UMAP Projection

dmr_met_adata = sc.AnnData(dmr_met.iloc[:, 84:].T)

# Percentage methylation stored on raw and layer attribute
dmr_met_adata.raw = dmr_met_adata
dmr_met_adata.layers['Percent_met'] = dmr_met_adata.X.copy()

# clinic info attachment
dmr_met_adata.obs = clinic_info.iloc[84:, :].copy()

# RNA batch info
rna_batch = pd.read_table("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/RNA_batch.txt", index_col=0)
dmr_met_adata.obs['RNA_batch'] = rna_batch.iloc[84:].copy()
dmr_met_adata.obs[['CpGimet', 'PMDmet', 'CIMPtype']] = cimp_pmd.copy()

# Scaling (optional) and PCA
sc.pp.scale(dmr_met_adata) # z-score scaling, which is (X-mean)/std. So the mean of the each variable becomes 0 (almost zero) and gets to have a unit variance.
sc.tl.pca(dmr_met_adata, n_comps=83, zero_center=True) # zero_center=True => compute standard PCA from covariance matrix
sc.pl.pca(dmr_met_adata, color=['EpiBurden', 'CIMP'], add_outline=False, legend_loc='right margin', color_map=cmap, use_raw=True, annotate_var_explained=True, size=200, components=['1,2'])

# Check for PCA numbers
pca_variance = pd.DataFrame(dmr_met_adata.uns['pca']['variance_ratio'], index=list(map(lambda x: 'PC' + str(x), list(range(1,84)))), columns=['Variance_ratio'])
np.sum(pca_variance.values.flatten()[:14]) # 0.7040434

# Neighbor Graph construction and leiden community detection
sc.pp.neighbors(dmr_met_adata, n_neighbors=10, n_pcs=14)
sc.tl.umap(dmr_met_adata, min_dist=0.5, spread=1.0, n_components=2, alpha=1.0, gamma=1.0, init_pos='spectral', method='umap')
sc.pl.umap(dmr_met_adata, color=['EpiBurden', 'CIMP'], add_outline=False, legend_loc='right margin', color_map=cmap, size=1000)
sns.despine()


sc.tl.leiden(dmr_met_adata, resolution=0.75, key_added='leiden_r075')
sc.tl.leiden(dmr_met_adata, resolution=0.5, key_added='leiden_r05')
sc.pl.umap(dmr_met_adata, color=['leiden_r075', 'leiden_r05', 'CIMP', 'EpiBurden'], add_outline=False, legend_loc='right margin', color_map=cmap)
sc.pl.umap(dmr_met_adata, color=['leiden_r075', 'Age', 'CIMP', 'EpiBurden'], add_outline=False, legend_loc='right margin', color_map=cmap)

from statsmodels.stats.anova import anova_lm
from statsmodels.formula.api import ols
model = ols('EpiBurden ~ C(leiden_r075)', dmr_met_adata.obs[['leiden_r075', 'EpiBurden']]).fit()
print(sm.stats.anova_lm(model))

from statsmodels.stats.multicomp import pairwise_tukeyhsd
print(pairwise_tukeyhsd(dmr_met_adata.obs[['leiden_r075', 'EpiBurden']]['EpiBurden'], dmr_met_adata.obs[['leiden_r075', 'EpiBurden']]['leiden_r075'], alpha=0.05)) # 1과 나머지는 모두 통과

# CIMP proportional plot
order = ['3', '2', '0', '1']
ax = (pd.crosstab(dmr_met_adata.obs['leiden_r075'], dmr_met_adata.obs['CIMP'], normalize=0).loc[order]*100).plot.bar(stacked=True, color=['maroon', 'salmon'], rot=0)
ax.set_ylabel("Proportion (%)")
ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1.0))
plt.tight_layout()
sns.despine()


# Diffusion pseudotime
sc.tl.diffmap(dmr_met_adata)
sc.pl.diffmap(dmr_met_adata, color=['leiden_r075'], add_outline=False, color_map=cmap)

start_cell = np.isin(dmr_met_adata.obs['leiden_r075'], '3') # leiden_r075의 3번 cluster
max_start_id = np.argmin(dmr_met_adata.obsm['X_diffmap'][start_cell, 2]) # 2 ==> DC1
root_id = np.arange(len(start_cell))[start_cell][max_start_id]
dmr_met_adata.uns['iroot'] = root_id
sc.tl.dpt(dmr_met_adata, n_branchings=0, n_dcs=15)

p = sns.lmplot(data=dmr_met_adata.obs, x='dpt_pseudotime', y='EpiBurden')
p.set_xlabels('Diffusion Pseudotime')
p.set_ylabels('Epimutation Burden')

stats.spearmanr(dmr_t.obs['dpt_pseudotime'], dmr_t.obs['EpiBurden'])

leiden_r075_color_dict = dict(zip(sorted(set(dmr_met_adata.obs['leiden_r075'].values)), dmr_met_adata.uns['leiden_r075_colors']))
dmr_met_adata.obs['leiden_r075_colors'] = dmr_met_adata.obs['leiden_r075'].map(lambda x: leiden_r075_color_dict[x])

cimp_color_dict = dict(zip(sorted(set(dmr_met_adata.obs['CIMP'].values)), ['maroon', 'salmon']))
dmr_met_adata.obs['cimp_colors'] = dmr_met_adata.obs['CIMP'].map(lambda x: cimp_color_dict[x])



x = np.array(dmr_met_adata.obs['dpt_pseudotime'])
y = np.array(dmr_met_adata.obs['EpiBurden'])
lowess = sm.nonparametric.lowess(y, x)

graph = sns.lmplot(data=dmr_met_adata.obs, x='dpt_pseudotime', y='EpiBurden', hue='leiden_r075', fit_reg=False)
sns.regplot(data=dmr_met_adata.obs, x='dpt_pseudotime', y='EpiBurden', scatter=False, lowess=True, ax=graph.axes[0, 0], line_kws={"color":"black"})
graph.ax.set_xlabel('Diffusion Pseudotime')
graph.ax.set_ylabel('Epimutation Burden')

graph = sns.lmplot(data=dmr_met_adata.obs, x='dpt_pseudotime', y='EpiBurden', hue='CIMPtype', fit_reg=False, palette={'CIMP(-) tumor (N=51)': 'salmon', 'CIMP(+) tumor (N=33)': 'maroon'})
sns.regplot(data=dmr_met_adata.obs, x='dpt_pseudotime', y='EpiBurden', scatter=False, lowess=True, ax=graph.axes[0, 0], line_kws={"color":"black"})
graph.ax.set_xlabel('Diffusion Pseudotime')
graph.ax.set_ylabel('Epimutation Burden')


lin = tuple(sorted(list(dmr_t.obs['DMR_leiden'].values.unique())))
dmr_t.obs['DMR_leiden'] = dmr_t.obs['DMR_leiden'].cat.reorder_categories(list(lin), ordered=True)





####################################################################################################################################
# Promoter - Gene Expression
# RNA expression processing
## Combat-corrected gene-level VST
gene_vst = pd.read_table("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/STAD_SNUH_vst_ComBat.txt", index_col=0, sep=' ')
#gene_vst = pd.read_table("/home/mhryan/Workspace/02.Projects/02.WC300/02.RNA-seq/STAD_SNUH_vst_ComBat.txt", index_col=0, sep=' ')
gene_vst.columns = list(map(lambda x: 'X'+x, gene_vst.columns))

## DEG results from DESeq2
deg_tn = pd.read_csv("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/Tumor.Normal.compare_ComBat.csv", index_col=0)
#deg_tn = pd.read_csv("/home/mhryan/Workspace/02.Projects/02.WC300/02.RNA-seq/Tumor.Normal.compare_ComBat.csv", index_col=0)
#deg_tn = deg_tn[deg_tn.index.isin(list(map(lambda x: x.split('/')[2], pls.index)))]

deg_genes = deg_tn[(deg_tn['padj'] < 0.01) & (deg_tn['baseMean'] > 10)].index

## Combat-corrected transcript-level counts
#trans_combat = pd.read_table("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/STAD_SNUH_tx_combat_counts.txt", index_col=0, sep=' ')
trans_combat = pd.read_table("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/STAD_SNUH_tx_combat_norm_counts.txt", index_col=0, sep=' ')
#trans_combat = pd.read_table("/home/mhryan/Workspace/02.Projects/02.WC300/02.RNA-seq/STAD_SNUH_tx_combat_norm_counts.txt", index_col=0, sep=' ')
trans_combat.columns = list(map(lambda x: 'X'+x, trans_combat.columns))


## Combat-corrected transcript-level log2(counts+1)
trans_combat_log2 = np.log2(trans_combat + 1)
del trans_combat

## ENSEMBL transcript ID to GeneID table
tx2gene = pd.read_table("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/txdb_geneSymbol.txt")
#tx2gene = pd.read_table("/home/mhryan/Workspace/02.Projects/02.WC300/02.RNA-seq/txdb_geneSymbol.txt")

# Promoter methylation processing

## PLS info
pls_info = pd.read_table("PLS_annotated_table_full_new_non-redundant.txt", index_col=0)
## PLS info with DMR
plsdmr_info = pd.read_table("PLS_annotated_table_full_new_non-redundant_DMR_abs10_woPMD.txt", index_col=0)
## NO-PLS info
nopls_info = pd.read_table("No_PLS_annotated_table_full_new.txt", index_col=0) # Genes w/ multiple overlapped prpomoter don't exist
## NO-PLS info with DMR
noplsdmr_info = pd.read_table("No_PLS_annotated_table_full_new_DMR_abs10_woPMD.txt", index_col=0)

## Gene-biotype used downstream
used_info = ['protein_coding',
             'lincRNA',
             '3prime_overlapping_ncRNA',
             'antisense',
             'bidirectional_promoter_lncRNA',
             'macro_lncRNA',
             'non_coding',
             'processed_transcript',
             'sense_intronic',
             'sense_overlapping']
pls_info = pls_info[pls_info['Type'].isin(used_info)]
nopls_info = nopls_info[nopls_info['Type'].isin(used_info)]
nopls_info = nopls_info[~nopls_info['GeneID'].isin(pls_info['GeneID'])]

plsdmr_info = plsdmr_info[plsdmr_info['Type'].isin(used_info)]
noplsdmr_info = noplsdmr_info[noplsdmr_info['Type'].isin(used_info)]
noplsdmr_info = noplsdmr_info[~noplsdmr_info['GeneID'].isin(plsdmr_info['GeneID'])]

pls = pd.read_table("cCRE_PLS_smooth_modified.txt", index_col=0)
nopls = pd.read_table("NO_PLS_smooth_modified.txt", index_col=0)
pls.columns = list(map(lambda x: 'X'+x, pls.columns))
nopls.columns = list(map(lambda x: 'X'+x, nopls.columns))
pls = pls[pls.index.isin(list(map(lambda x:'/'.join(x.split('/')[:-1]), pls_info.index)))]
nopls = nopls[nopls.index.isin(list(map(lambda x:'/'.join(x.split('/')[:-1]), nopls_info.index)))]

## Merging PLS and No_PLS
comb_pls = pd.concat([pls, nopls])
comb_pls = comb_pls*100
comb_pls_info = pd.concat([pls_info, nopls_info])
comb_pls.index = comb_pls_info.index

comb_plsdmr_info = pd.concat([plsdmr_info, noplsdmr_info])
comb_plsdmr = comb_pls[comb_pls.index.isin(comb_plsdmr_info.index)]
comb_plsdmr = comb_plsdmr
comb_plsdmr_info = comb_plsdmr_info[comb_plsdmr_info.index.isin(comb_plsdmr.index)]

del pls, nopls, pls_info, nopls_info

comb_pls
comb_pls_info
comb_plsdmr
comb_plsdmr_info

a = pd.DataFrame(trans_combat_log2[trans_combat_log2.index.isin(list(comb_pls_info[comb_pls_info['CpGi'] == 'Yes']['ENSTID'].values))].iloc[:,:84].mean(axis=1), columns=['CpGi Promoter'])
b = pd.DataFrame(trans_combat_log2[trans_combat_log2.index.isin(list(comb_pls_info[comb_pls_info['CpGi'] == 'na']['ENSTID'].values))].iloc[:,:84].mean(axis=1), columns=['Non-CpGi Promoter'])
ab = pd.concat([a,b], axis=1)
#p = sns.boxplot(data=ab, palette={'CpGi Promoter':'#00203FFF', 'Non-CpGi Promoter':'#ADEFD1FF'}, width=0.5, showfliers = False)
p = sns.violinplot(data=ab, palette={'CpGi Promoter':'#00203FFF', 'Non-CpGi Promoter':'#ADEFD1FF'}, width=0.5, showfliers = False, scale="count", cut=0)
p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.1)
p.set_ylabel("Gene expression")
p.set_title("CpG islands present in Promoter or not")
plt.tight_layout()
sns.despine()
stats.ttest_ind(a, b)


## Promoter DNA methylation Distribution according to Histone Modifications - DMR
total1 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'Yes') & (comb_plsdmr_info['K27ac'] == 'Yes') & (comb_plsdmr_info['K27me3'] == 'na')].index)]
total2 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'Yes') & (comb_plsdmr_info['K27ac'] == 'Yes') & (comb_plsdmr_info['K27me3'] == 'Yes')].index)]
total3 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'Yes') & (comb_plsdmr_info['K27ac'] == 'na') & (comb_plsdmr_info['K27me3'] == 'Yes')].index)]
total4 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'Yes') & (comb_plsdmr_info['K27ac'] == 'na') & (comb_plsdmr_info['K27me3'] == 'na')].index)]
#total5 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'na') & (comb_plsdmr_info['K27ac'] == 'Yes') & (comb_plsdmr_info['K27me3'] == 'na')].index)]
total6 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'na') & (comb_plsdmr_info['K27ac'] == 'na') & (comb_plsdmr_info['K27me3'] == 'Yes')].index)]
total7 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'na') & (comb_plsdmr_info['K27ac'] == 'na') & (comb_plsdmr_info['K27me3'] == 'na')].index)]

## Only Protein Coding Genes
total1 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'Yes') & (comb_plsdmr_info['K27ac'] == 'Yes') & (comb_plsdmr_info['K27me3'] == 'na') & (comb_plsdmr_info['Type'] == 'protein_coding')].index)]
total2 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'Yes') & (comb_plsdmr_info['K27ac'] == 'Yes') & (comb_plsdmr_info['K27me3'] == 'Yes') & (comb_plsdmr_info['Type'] == 'protein_coding')].index)]
total3 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'Yes') & (comb_plsdmr_info['K27ac'] == 'na') & (comb_plsdmr_info['K27me3'] == 'Yes') & (comb_plsdmr_info['Type'] == 'protein_coding')].index)]
total4 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'Yes') & (comb_plsdmr_info['K27ac'] == 'na') & (comb_plsdmr_info['K27me3'] == 'na') & (comb_plsdmr_info['Type'] == 'protein_coding')].index)]
#total5 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'na') & (comb_plsdmr_info['K27ac'] == 'Yes') & (comb_plsdmr_info['K27me3'] == 'na') & (comb_plsdmr_info['Type'] == 'protein_coding')].index)]
total6 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'na') & (comb_plsdmr_info['K27ac'] == 'na') & (comb_plsdmr_info['K27me3'] == 'Yes') & (comb_plsdmr_info['Type'] == 'protein_coding')].index)]
total7 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'na') & (comb_plsdmr_info['K27ac'] == 'na') & (comb_plsdmr_info['K27me3'] == 'na') & (comb_plsdmr_info['Type'] == 'protein_coding')].index)]

mean_total1 = pd.DataFrame({'AverageMet': total1.iloc[:, :84].mean(axis=1).values, 'Type': ['H3K4me3/H3K27ac']*total1.shape[0]}, index=total1.index)
mean_total2 = pd.DataFrame({'AverageMet': total2.iloc[:, :84].mean(axis=1).values, 'Type': ['H3K4me3/H3K27ac/H3K27me3']*total2.shape[0]}, index=total2.index)
mean_total3 = pd.DataFrame({'AverageMet': total3.iloc[:, :84].mean(axis=1).values, 'Type': ['H3K4me3/H3K27me3']*total3.shape[0]}, index=total3.index)
mean_total4 = pd.DataFrame({'AverageMet': total4.iloc[:, :84].mean(axis=1).values, 'Type': ['H3K4me3']*total4.shape[0]}, index=total4.index)
#mean_total5 = pd.DataFrame({'AverageMet': total5.iloc[:, :84].mean(axis=1).values, 'Type': ['H3K27ac']*total5.shape[0]}, index=total5.index)
mean_total6 = pd.DataFrame({'AverageMet': total6.iloc[:, :84].mean(axis=1).values, 'Type': ['H3K27me3']*total6.shape[0]}, index=total6.index)
mean_total7 = pd.DataFrame({'AverageMet': total7.iloc[:, :84].mean(axis=1).values, 'Type': ['NA']*total7.shape[0]}, index=total7.index)

#mean_total = pd.concat([mean_total1, mean_total2, mean_total3, mean_total4, mean_total5, mean_total6, mean_total7], axis=0)
mean_total = pd.concat([mean_total1, mean_total2, mean_total3, mean_total4, mean_total6, mean_total7], axis=0)

# Ridgeline plots

def label(x, color, label):
    ax = plt.gca()
    ax.text(0.85, 0.18, label, color='black', fontsize=13, ha='left', va='center', transform=ax.transAxes)

sns.set_theme(style="white", font="Arial", rc={"axes.facecolor": (0, 0, 0, 0)})
#sns.color_palette("Accent", 6).as_hex()
grid = sns.FacetGrid(mean_total, row="Type", hue="Type", palette=sns.color_palette("Accent", 6), aspect=5, height=1.15)
grid.map_dataframe(sns.kdeplot, x="AverageMet", fill=True, alpha=1)
grid.map_dataframe(sns.kdeplot, x="AverageMet", color='black')
grid.set(xlim=(0, 100))
grid.fig.subplots_adjust(hspace=-0.83)
grid.set_titles("")
grid.set(yticks=[])
grid.set_axis_labels("", "")
grid.despine(left=True)
#grid.figure.suptitle("Seven Different Combinations of Histone modifications", size=10)
grid.figure.supylabel('Kernel Density', x=0.1, y=0.5)
grid.figure.supxlabel("Average Promoter DNA Methylation (%)", x=0.55, y=0.02)
grid.figure.subplots_adjust(top=1.01)

## Promoter DNA methylation Distribution according to Histone Modifications - ALL including DMR overlapped Promoter
total1 = comb_pls[comb_pls.index.isin(comb_pls_info[(comb_pls_info['K4me3'] == 'Yes') & (comb_pls_info['K27ac'] == 'Yes') & (comb_pls_info['K27me3'] == 'Yes')].index)]
total2 = comb_pls[comb_pls.index.isin(comb_pls_info[(comb_pls_info['K4me3'] == 'Yes') & (comb_pls_info['K27ac'] == 'Yes') & (comb_pls_info['K27me3'] == 'na')].index)]
total3 = comb_pls[comb_pls.index.isin(comb_pls_info[(comb_pls_info['K4me3'] == 'Yes') & (comb_pls_info['K27ac'] == 'na') & (comb_pls_info['K27me3'] == 'Yes')].index)]
total4 = comb_pls[comb_pls.index.isin(comb_pls_info[(comb_pls_info['K4me3'] == 'Yes') & (comb_pls_info['K27ac'] == 'na') & (comb_pls_info['K27me3'] == 'na')].index)]
total5 = comb_pls[comb_pls.index.isin(comb_pls_info[(comb_pls_info['K4me3'] == 'na') & (comb_pls_info['K27ac'] == 'Yes') & (comb_pls_info['K27me3'] == 'na')].index)]
total6 = comb_pls[comb_pls.index.isin(comb_pls_info[(comb_pls_info['K4me3'] == 'na') & (comb_pls_info['K27ac'] == 'na') & (comb_pls_info['K27me3'] == 'Yes')].index)]
total7 = comb_pls[comb_pls.index.isin(comb_pls_info[(comb_pls_info['K4me3'] == 'na') & (comb_pls_info['K27ac'] == 'na') & (comb_pls_info['K27me3'] == 'na')].index)]

total1 = comb_pls[comb_pls.index.isin(comb_pls_info[(comb_pls_info['K4me3'] == 'Yes') & (comb_pls_info['K27ac'] == 'Yes') & (comb_pls_info['K27me3'] == 'Yes') & (comb_pls_info['Type'] == 'protein_coding')].index)]
total2 = comb_pls[comb_pls.index.isin(comb_pls_info[(comb_pls_info['K4me3'] == 'Yes') & (comb_pls_info['K27ac'] == 'Yes') & (comb_pls_info['K27me3'] == 'na') & (comb_pls_info['Type'] == 'protein_coding')].index)]
total3 = comb_pls[comb_pls.index.isin(comb_pls_info[(comb_pls_info['K4me3'] == 'Yes') & (comb_pls_info['K27ac'] == 'na') & (comb_pls_info['K27me3'] == 'Yes') & (comb_pls_info['Type'] == 'protein_coding')].index)]
total4 = comb_pls[comb_pls.index.isin(comb_pls_info[(comb_pls_info['K4me3'] == 'Yes') & (comb_pls_info['K27ac'] == 'na') & (comb_pls_info['K27me3'] == 'na') & (comb_pls_info['Type'] == 'protein_coding')].index)]
total5 = comb_pls[comb_pls.index.isin(comb_pls_info[(comb_pls_info['K4me3'] == 'na') & (comb_pls_info['K27ac'] == 'Yes') & (comb_pls_info['K27me3'] == 'na') & (comb_pls_info['Type'] == 'protein_coding')].index)]
total6 = comb_pls[comb_pls.index.isin(comb_pls_info[(comb_pls_info['K4me3'] == 'na') & (comb_pls_info['K27ac'] == 'na') & (comb_pls_info['K27me3'] == 'Yes') & (comb_pls_info['Type'] == 'protein_coding')].index)]
total7 = comb_pls[comb_pls.index.isin(comb_pls_info[(comb_pls_info['K4me3'] == 'na') & (comb_pls_info['K27ac'] == 'na') & (comb_pls_info['K27me3'] == 'na') & (comb_pls_info['Type'] == 'protein_coding')].index)]

## Normals
mean_total1_N = pd.DataFrame({'AverageMet': total2.iloc[:, :84].mean(axis=1).values, 'Type': ['H3K4me3/H3K27ac']*total2.shape[0]}, index=total2.index)
mean_total2_N = pd.DataFrame({'AverageMet': total1.iloc[:, :84].mean(axis=1).values, 'Type': ['H3K4me3/H3K27ac/H3K27me3']*total1.shape[0]}, index=total1.index)
mean_total3_N = pd.DataFrame({'AverageMet': total3.iloc[:, :84].mean(axis=1).values, 'Type': ['H3K4me3/H3K27me3']*total3.shape[0]}, index=total3.index)
mean_total4_N = pd.DataFrame({'AverageMet': total4.iloc[:, :84].mean(axis=1).values, 'Type': ['H3K4me3']*total4.shape[0]}, index=total4.index)
mean_total5_N = pd.DataFrame({'AverageMet': total5.iloc[:, :84].mean(axis=1).values, 'Type': ['H3K27ac']*total5.shape[0]}, index=total5.index)
mean_total6_N = pd.DataFrame({'AverageMet': total6.iloc[:, :84].mean(axis=1).values, 'Type': ['H3K27me3']*total6.shape[0]}, index=total6.index)
mean_total7_N = pd.DataFrame({'AverageMet': total7.iloc[:, :84].mean(axis=1).values, 'Type': ['NA']*total7.shape[0]}, index=total7.index)

mean_total_N = pd.concat([mean_total1_N, mean_total2_N, mean_total3_N, mean_total4_N, mean_total5_N, mean_total6_N, mean_total7_N], axis=0)

## Tumors
mean_total1_T = pd.DataFrame({'AverageMet': total2.iloc[:, 84:].mean(axis=1).values, 'Type': ['H3K4me3/H3K27ac']*total2.shape[0]}, index=total2.index)
mean_total2_T = pd.DataFrame({'AverageMet': total1.iloc[:, 84:].mean(axis=1).values, 'Type': ['H3K4me3/H3K27ac/H3K27me3']*total1.shape[0]}, index=total1.index)
mean_total3_T = pd.DataFrame({'AverageMet': total3.iloc[:, 84:].mean(axis=1).values, 'Type': ['H3K4me3/H3K27me3']*total3.shape[0]}, index=total3.index)
mean_total4_T = pd.DataFrame({'AverageMet': total4.iloc[:, 84:].mean(axis=1).values, 'Type': ['H3K4me3']*total4.shape[0]}, index=total4.index)
mean_total5_T = pd.DataFrame({'AverageMet': total5.iloc[:, 84:].mean(axis=1).values, 'Type': ['H3K27ac']*total5.shape[0]}, index=total5.index)
mean_total6_T = pd.DataFrame({'AverageMet': total6.iloc[:, 84:].mean(axis=1).values, 'Type': ['H3K27me3']*total6.shape[0]}, index=total6.index)
mean_total7_T = pd.DataFrame({'AverageMet': total7.iloc[:, 84:].mean(axis=1).values, 'Type': ['NA']*total7.shape[0]}, index=total7.index)

mean_total_T = pd.concat([mean_total1_T, mean_total2_T, mean_total3_T, mean_total4_T, mean_total5_T, mean_total6_T, mean_total7_T], axis=0)

# Ridgeline plots

def label(x, color, label):
    ax = plt.gca()
    ax.text(0.85, 0.18, label, color='black', fontsize=13, ha='left', va='center', transform=ax.transAxes)

sns.set_theme(style="white", font="Arial", rc={"axes.facecolor": (0, 0, 0, 0)})
grid = sns.FacetGrid(mean_total_N, row="Type", hue="Type", palette=sns.color_palette("Accent", 7), aspect=5, height=1.15)
grid.map_dataframe(sns.kdeplot, x="AverageMet", fill=True, alpha=1)
grid.map_dataframe(sns.kdeplot, x="AverageMet", color='black')
grid.set(xlim=(0, 100))
grid.fig.subplots_adjust(hspace=-0.83)
grid.set_titles("")
grid.set(yticks=[])
grid.set_axis_labels("", "")
grid.despine(left=True)
#grid.figure.suptitle("Seven Different Combinations of Histone modifications", size=10)
grid.figure.supylabel('Kernel Density', x=0.1, y=0.5, size=12)
grid.figure.supxlabel("Average Promoter DNA Methylation (%) of Normal Samples", x=0.55, y=0.02, size=12)
grid.figure.subplots_adjust(top=1.01)

grid = sns.FacetGrid(mean_total_T, row="Type", hue="Type", palette=sns.color_palette("Accent", 7), aspect=5, height=1.15)
grid.map_dataframe(sns.kdeplot, x="AverageMet", fill=True, alpha=1)
grid.map_dataframe(sns.kdeplot, x="AverageMet", color='black')
grid.set(xlim=(0, 100))
grid.fig.subplots_adjust(hspace=-0.83)
grid.set_titles("")
grid.set(yticks=[])
grid.set_axis_labels("", "")
grid.despine(left=True)
#grid.figure.suptitle("Seven Different Combinations of Histone modifications", size=10)
grid.figure.supylabel('Kernel Density', x=0.1, y=0.5, size=12)
grid.figure.supxlabel("Average Promoter DNA Methylation (%) of Tumor Samples", x=0.55, y=0.02, size=12)
grid.figure.subplots_adjust(top=1.01)


# RNA (transcript) expression distribution plot version 1
p = sns.histplot( trans_combat_log2.iloc[:, :84].sum(axis=1).div(84)[trans_combat_log2.iloc[:, :84].sum(axis=1).div(84).index.isin( list(set(comb_plsdmr_info[comb_plsdmr_info.index.isin(total1.index)]['ENSTID'].values)) )], kde=True, stat='count', color=sns.color_palette("Accent", 7)[0])
p = sns.histplot( trans_combat_log2.iloc[:, :84].sum(axis=1).div(84)[trans_combat_log2.iloc[:, :84].sum(axis=1).div(84).index.isin( list(set(comb_plsdmr_info[comb_plsdmr_info.index.isin(total2.index)]['ENSTID'].values)) )], kde=True, stat='count', color=sns.color_palette("Accent", 7)[1])
p = sns.histplot( trans_combat_log2.iloc[:, :84].sum(axis=1).div(84)[trans_combat_log2.iloc[:, :84].sum(axis=1).div(84).index.isin( list(set(comb_plsdmr_info[comb_plsdmr_info.index.isin(total3.index)]['ENSTID'].values)) )], kde=True, stat='count', color=sns.color_palette("Accent", 7)[2])
p = sns.histplot( trans_combat_log2.iloc[:, :84].sum(axis=1).div(84)[trans_combat_log2.iloc[:, :84].sum(axis=1).div(84).index.isin( list(set(comb_plsdmr_info[comb_plsdmr_info.index.isin(total4.index)]['ENSTID'].values)) )], kde=True, stat='count', color=sns.color_palette("Accent", 7)[3])
p = sns.histplot( trans_combat_log2.iloc[:, :84].sum(axis=1).div(84)[trans_combat_log2.iloc[:, :84].sum(axis=1).div(84).index.isin( list(set(comb_plsdmr_info[comb_plsdmr_info.index.isin(total5.index)]['ENSTID'].values)) )], kde=True, stat='count', color=sns.color_palette("Accent", 7)[4])
p = sns.histplot( trans_combat_log2.iloc[:, :84].sum(axis=1).div(84)[trans_combat_log2.iloc[:, :84].sum(axis=1).div(84).index.isin( list(set(comb_plsdmr_info[comb_plsdmr_info.index.isin(total6.index)]['ENSTID'].values)) )], kde=True, stat='count', color=sns.color_palette("Accent", 7)[5])
p = sns.histplot( trans_combat_log2.iloc[:, :84].sum(axis=1).div(84)[trans_combat_log2.iloc[:, :84].sum(axis=1).div(84).index.isin( list(set(comb_plsdmr_info[comb_plsdmr_info.index.isin(total7.index)]['ENSTID'].values)) )], kde=True, stat='count', color=sns.color_palette("Accent", 7)[6])
p.set_xlabel("Mean Log2 normalized RNA read counts across normal samples")
p.set_ylabel("Number of transcripts")
sns.despine()
plt.xlim((0, 20))
#plt.xlim((0, 12))
plt.tight_layout()
from matplotlib.ticker import StrMethodFormatter
import matplotlib.ticker as ticker
p.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))

# RNA (transcript) expression distribution plot version 2 ==> This is it
fig, ax = plt.subplots(1,1, figsize=(7, 7), constrained_layout=True)
sns.kdeplot(trans_combat_log2.iloc[:, :84].sum(axis=1).div(84)[trans_combat_log2.iloc[:, :84].sum(axis=1).div(84).index.isin( list(set(comb_plsdmr_info[comb_plsdmr_info.index.isin(total1.index)]['ENSTID'].values)) )], color=sns.color_palette("Accent", 7)[0], fill=False, ax=ax)
sns.kdeplot(trans_combat_log2.iloc[:, :84].sum(axis=1).div(84)[trans_combat_log2.iloc[:, :84].sum(axis=1).div(84).index.isin( list(set(comb_plsdmr_info[comb_plsdmr_info.index.isin(total2.index)]['ENSTID'].values)) )], color=sns.color_palette("Accent", 7)[1], fill=False, ax=ax)
sns.kdeplot(trans_combat_log2.iloc[:, :84].sum(axis=1).div(84)[trans_combat_log2.iloc[:, :84].sum(axis=1).div(84).index.isin( list(set(comb_plsdmr_info[comb_plsdmr_info.index.isin(total3.index)]['ENSTID'].values)) )], color=sns.color_palette("Accent", 7)[2], fill=False, ax=ax)
sns.kdeplot(trans_combat_log2.iloc[:, :84].sum(axis=1).div(84)[trans_combat_log2.iloc[:, :84].sum(axis=1).div(84).index.isin( list(set(comb_plsdmr_info[comb_plsdmr_info.index.isin(total4.index)]['ENSTID'].values)) )], color=sns.color_palette("Accent", 7)[3], fill=False, ax=ax)
sns.kdeplot(trans_combat_log2.iloc[:, :84].sum(axis=1).div(84)[trans_combat_log2.iloc[:, :84].sum(axis=1).div(84).index.isin( list(set(comb_plsdmr_info[comb_plsdmr_info.index.isin(total5.index)]['ENSTID'].values)) )], color=sns.color_palette("Accent", 7)[4], fill=False, ax=ax)
sns.kdeplot(trans_combat_log2.iloc[:, :84].sum(axis=1).div(84)[trans_combat_log2.iloc[:, :84].sum(axis=1).div(84).index.isin( list(set(comb_plsdmr_info[comb_plsdmr_info.index.isin(total6.index)]['ENSTID'].values)) )], color=sns.color_palette("Accent", 7)[5], fill=False, ax=ax)
sns.kdeplot(trans_combat_log2.iloc[:, :84].sum(axis=1).div(84)[trans_combat_log2.iloc[:, :84].sum(axis=1).div(84).index.isin( list(set(comb_plsdmr_info[comb_plsdmr_info.index.isin(total7.index)]['ENSTID'].values)) )], color=sns.color_palette("Accent", 7)[6], fill=False, ax=ax)
ax.set_xlim((0, 17.5))
ax.set_xlabel("Mean Log2 normalized RNA read counts across normal samples")
sns.despine(ax=ax)


# Promoter methylation distribution plot version 1
p = sns.histplot(total1.iloc[:, :84].mean(axis=1), kde=True, stat='count', color=sns.color_palette("Accent", 7)[0])
p = sns.histplot(total2.iloc[:, :84].mean(axis=1), kde=True, stat='count', color=sns.color_palette("Accent", 7)[1])
p = sns.histplot(total3.iloc[:, :84].mean(axis=1), kde=True, stat='count', color=sns.color_palette("Accent", 7)[2])
p = sns.histplot(total4.iloc[:, :84].mean(axis=1), kde=True, stat='count', color=sns.color_palette("Accent", 7)[3])
p = sns.histplot(total5.iloc[:, :84].mean(axis=1), kde=True, stat='count', color=sns.color_palette("Accent", 7)[4])
p = sns.histplot(total6.iloc[:, :84].mean(axis=1), kde=True, stat='count', color=sns.color_palette("Accent", 7)[5])
p = sns.histplot(total7.iloc[:, :84].mean(axis=1), kde=True, stat='count', color=sns.color_palette("Accent", 7)[6])
plt.xlim((0, 100))
sns.despine()
p.set_xlabel("Mean Promoter DNA Methylation (%) across normal samples")
plt.tight_layout()

#plt.xlim((0, 12))
plt.tight_layout()
from matplotlib.ticker import StrMethodFormatter
import matplotlib.ticker as ticker
p.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))

# Promoter methylation distribution plot version 2
ax = plt.subplot(1,1,1)
sns.kdeplot(total1.iloc[:, :84].mean(axis=1), color=sns.color_palette("Accent", 7)[0], ax=ax)
sns.kdeplot(total2.iloc[:, :84].mean(axis=1), color=sns.color_palette("Accent", 7)[1], ax=ax)
sns.kdeplot(total3.iloc[:, :84].mean(axis=1), color=sns.color_palette("Accent", 7)[2], ax=ax)
sns.kdeplot(total4.iloc[:, :84].mean(axis=1), color=sns.color_palette("Accent", 7)[3], ax=ax)
sns.kdeplot(total5.iloc[:, :84].mean(axis=1), color=sns.color_palette("Accent", 7)[4], ax=ax)
sns.kdeplot(total6.iloc[:, :84].mean(axis=1), color=sns.color_palette("Accent", 7)[5], ax=ax)
sns.kdeplot(total7.iloc[:, :84].mean(axis=1), color=sns.color_palette("Accent", 7)[6], ax=ax)
plt.xlim((0, 100))
ax.set_xlabel("Mean Promoter DNA Methylation (%) across normal samples")
sns.despine()
plt.tight_layout()




sns.heatmap(pd.concat([pd.DataFrame(comb_plsdmr[comb_plsdmr_info['GeneID'] == 'GPR25'].iloc[:, :84].values), pd.DataFrame(comb_plsdmr[comb_plsdmr_info['GeneID'] ='GPR25'].iloc[:, 84:].values)], axis=0), cmap=cmap, xticklabels=False, yticklabels=False)



# Promoter CpG density distribution plot version 2
a1 = comb_pls_info[comb_pls_info.index.isin(total2.index)]['CpGdensity']
a2 = pd.concat([ a1, pd.DataFrame(["total2"]*len(total2), index=a1.index, columns=["Category"]) ], axis=1)
b1 = comb_pls_info[comb_pls_info.index.isin(total1.index)]['CpGdensity']
b2 = pd.concat([ b1, pd.DataFrame(["total1"]*len(total1), index=b1.index, columns=["Category"]) ], axis=1)
c1 = comb_pls_info[comb_pls_info.index.isin(total3.index)]['CpGdensity']
c2 = pd.concat([ c1, pd.DataFrame(["total3"]*len(total3), index=c1.index, columns=["Category"]) ], axis=1)
d1 = comb_pls_info[comb_pls_info.index.isin(total4.index)]['CpGdensity']
d2 = pd.concat([ d1, pd.DataFrame(["total4"]*len(total4), index=d1.index, columns=["Category"]) ], axis=1)
e1 = comb_pls_info[comb_pls_info.index.isin(total5.index)]['CpGdensity']
e2 = pd.concat([ e1, pd.DataFrame(["total5"]*len(total5), index=e1.index, columns=["Category"]) ], axis=1)
f1 = comb_pls_info[comb_pls_info.index.isin(total6.index)]['CpGdensity']
f2 = pd.concat([ f1, pd.DataFrame(["total6"]*len(total6), index=f1.index, columns=["Category"]) ], axis=1)
g1 = comb_pls_info[comb_pls_info.index.isin(total7.index)]['CpGdensity']
g2 = pd.concat([ g1, pd.DataFrame(["total7"]*len(total7), index=g1.index, columns=["Category"]) ], axis=1)

total = pd.concat([ a2, b2, c2, d2, e2, f2, g2 ], axis=0)
del a1, a2, b1, b2, c1, c2, d1, d2, e1, e2, f1, f2, g1, g2

palette = {"total2": sns.color_palette("Accent", 7)[0],
           "total1": sns.color_palette("Accent", 7)[1],
           "total3": sns.color_palette("Accent", 7)[2],
           "total4": sns.color_palette("Accent", 7)[3],
           "total5": sns.color_palette("Accent", 7)[4],
           "total6": sns.color_palette("Accent", 7)[5],
           "total7": sns.color_palette("Accent", 7)[6]}
p = sns.violinplot(data=total, x='CpGdensity', y='Category', palette=palette, scale='width', cut=0)
p = sns.stripplot(data=total, x='CpGdensity', y='Category', jitter=True, marker='o', color='black', size=1, alpha=0.1)
sns.despine()
p.set_xlabel("CpG Density within DMR Promoter")
p.set_ylabel("")
p.set_yticklabels(["K4me3/K27ac/K27me3","K4me3/K27ac","K4me3/K27me3","K4me3","K27ac","K27me3","None"])
plt.tight_layout()










pls_ttest_tstats = comb_pls.apply(lambda x: stats.ttest_rel(x.iloc[84:], x.iloc[:84])[0], axis=1)
pls_ttest_pvalue = comb_pls.apply(lambda x: stats.ttest_rel(x.iloc[84:], x.iloc[:84])[1], axis=1)
pls_ttest = pd.concat([pls_ttest_tstats, pls_ttest_pvalue], keys=["tstats", "pvalue"], axis=1)
pls_ttest['Bonferroni_pvalue'] = sm.stats.multipletests(pls_ttest['pvalue'], alpha=0.05, method='bonferroni')[1]


g = sns.clustermap(comb_pls[comb_pls.index.isin((comb_pls_info[comb_pls_info.index.isin(pls_ttest[pls_ttest['Bonferroni_pvalue'] < 0.05].index)]['Type'] == 'protein_coding').index)].iloc[:,84:],
                   method='complete',
                   metric='euclidean',
                   row_cluster=True,
                   col_cluster=True,
                   z_score=None,
                   standard_scale=None,
                   cmap=cmap,
                   xticklabels=False,
                   yticklabels=False,
                   col_colors=[dmr_met_adata.obs['leiden_r075_colors'].values, dmr_met_adata.obs['cimp_colors'].values],
                   row_colors=None) # standard scale: 0 (rows) or 1 columns (subtract min and divide by max)

# Hypermethylated Promoters
g1 = sns.clustermap(comb_pls[comb_pls.index.isin((comb_pls_info[comb_pls_info.index.isin(pls_ttest[(pls_ttest['Bonferroni_pvalue'] < 0.05) & (pls_ttest['tstats'] > 0)].index)]['Type'] == 'protein_coding').index)].iloc[:,84:],
                   method='ward',
                   metric='euclidean',
                   row_cluster=True,
                   col_cluster=True,
                   z_score=None,
                   standard_scale=None,
                   cmap=cmap,
                   xticklabels=False,
                   yticklabels=False,
                   col_colors=[dmr_met_adata.obs['leiden_r075_colors'].values, dmr_met_adata.obs['cimp_colors'].values],
                   row_colors=None) # standard scale: 0 (rows) or 1 columns (subtract min and divide by max)
g1.ax_heatmap.set_ylabel('')
# Hypomethylated Promoters
g2 = sns.clustermap(comb_pls[comb_pls.index.isin((comb_pls_info[comb_pls_info.index.isin(pls_ttest[(pls_ttest['Bonferroni_pvalue'] < 0.05) & (pls_ttest['tstats'] < 0)].index)]['Type'] == 'protein_coding').index)].iloc[:,84:],
                   method='ward',
                   metric='euclidean',
                   row_cluster=True,
                   col_cluster=True,
                   z_score=None,
                   standard_scale=None,
                   cmap=cmap,
                   xticklabels=False,
                   yticklabels=False,
                   col_colors=[dmr_met_adata.obs['leiden_r075_colors'].values, dmr_met_adata.obs['cimp_colors'].values],
                   row_colors=None) # standard scale: 0 (rows) or 1 columns (subtract min and divide by max)
g2.ax_heatmap.set_ylabel('')



leiden_r075_color_dict = dict(zip(sorted(set(dmr_met_adata.obs['leiden_r075'].values)), dmr_met_adata.uns['leiden_r075_colors']))

leiden_r075_wNormal_color_dict = dict(zip(sorted(set(dmr_met_adata.obs['leiden_r075'].values)), dmr_met_adata.uns['leiden_r075_colors']))
leiden_r075_wNormal_color_dict['Normal'] = 'grey'

df1 = pd.concat([pd.DataFrame(comb_pls.loc['chr1:217137341-217138084/EH38D4462443,EH38D4462444/ESRRG/ENSG00000196482/ENST00000493603/protein_coding/55'].iloc[84:]), pd.DataFrame(dmr_met_adata.obs['leiden_r075'])], axis=1)
df2 = pd.concat([pd.DataFrame(comb_pls.loc['chr1:217137341-217138084/EH38D4462443,EH38D4462444/ESRRG/ENSG00000196482/ENST00000493603/protein_coding/55'].iloc[:84]), pd.DataFrame(["Normal"]*84, columns=["leiden_r075"], index=comb_pls.iloc[:,:84].columns)], axis=1)
df = pd.concat([df1, df2])
g = sns.boxplot(data=df, x='leiden_r075', y=df.columns[0], order=["Normal", "3", "2", "0", "1"], showfliers=False, palette=leiden_r075_wNormal_color_dict)
g = sns.stripplot(data=df, x='leiden_r075', y=df.columns[0], order=["Normal", "3", "2", "0", "1"], color=".2", size=2, alpha=0.4)
g.set_ylabel("ESRRG Promoter DNA Methylation (%)")
g.set_ylim((0,100))
sns.despine()
plt.tight_layout()


#df1 = pd.concat([pls.loc[pls_info[pls_info['GeneID'] == gene].index].iloc[0, 84:], pd.DataFrame(dmr_met_adata.obs['leiden_r075'])], axis=1)
#df2 = pd.concat([pls.loc[pls_info[pls_info['GeneID'] == gene].index].iloc[0, :84], pd.DataFrame(["Normal"]*84, columns=["DMR Clusters"], index=pls.iloc[::84].columns)], axis=1)
#df = pd.concat([df1, df2])
#sns.boxplot(data=df, x='DMR Clusters', y=df.columns[0], palette=color_dict3, order=["Normal", "leiden_A", "leiden_B", "leiden_C", "leiden_D"],showfliers=False)
#sns.stripplot(data=df, x='DMR Clusters', y=df.columns[0], order=["Normal", "leiden_A", "leiden_B", "leiden_C", "leiden_D"], color=".2", size=2, alpha=0.4)set(ylabel=df.columns[0].split('/')[2] + ' promoter methylation')
#stats.mannwhitneyu(df[df['DMR Clusters'] == 'Normal'].iloc[:,0], df[df['DMR Clusters'] == 'leiden_D'].iloc[:,0])

# RNA expression boxplot
df1 = pd.concat([pd.DataFrame(gene_vst.iloc[:,84:].loc['ESRRG']), pd.DataFrame(dmr_met_adata.obs['leiden_r075'])], axis=1)
df2 = pd.concat([pd.DataFrame(gene_vst.iloc[:,:84].loc['ESRRG']), pd.DataFrame(["Normal"]*84, columns=["leiden_r075"], index=comb_pls.iloc[:,:84].columns)], axis=1)
df = pd.concat([df1, df2])
g = sns.boxplot(data=df, x='leiden_r075', y=df.columns[0], order=["Normal", "3", "2", "0", "1"], showfliers=False, palette=leiden_r075_wNormal_color_dict)
g = sns.stripplot(data=df, x='leiden_r075', y=df.columns[0], order=["Normal", "3", "2", "0", "1"], color=".2", size=2, alpha=0.4)
g.set_ylabel("ESRRG expression (vst)")
sns.despine()
plt.tight_layout()