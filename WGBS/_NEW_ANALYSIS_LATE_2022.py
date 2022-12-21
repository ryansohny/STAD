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
g.ax_heatmap.set_ylabel('')

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
lola_hyper_dmr = lola_dmr[(lola_dmr['userSet'] == 1) & (lola_dmr['collection'] == 'encode_tfbs')]
lola_hypo_dmr = lola_dmr[(lola_dmr['userSet'] == 2) & (lola_dmr['collection'] == 'encode_tfbs')]

lola_hyper_dmr['Antibody from ENCODE'] = lola_hyper_dmr[['antibody', 'cellType']].apply(lambda x: ' from '.join(x), axis=1)
lola_hypo_dmr['Antibody from ENCODE'] = lola_hypo_dmr[['antibody', 'cellType']].apply(lambda x: ' from '.join(x), axis=1)


lola_hypo_dmr[['antibody', 'cellType', 'treatment']].apply(lambda x: ', '.join(x), axis=1)



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
p = sns.violinplot(data=ab, palette={'CpGi Promoter':'#00203FFF', 'Non-CpGi Promoter':'#ADEFD1FF'}, width=0.5, showfliers = False, scale="width", cut=0)
p = sns.stripplot(data=ab, jitter=True, marker='o', color='black', size=1.5, alpha=0.1)
p.set_ylabel("Gene expression")
p.set_title("CpG islands present in Promoter or not")
plt.tight_layout()
sns.despine()
stats.ttest_ind(a, b)


## Transcript Expression Distribution according to Histone Modifications
total1 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'Yes') & (comb_plsdmr_info['K27ac'] == 'Yes') & (comb_plsdmr_info['K27me3'] == 'Yes')].index)]
total2 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'Yes') & (comb_plsdmr_info['K27ac'] == 'Yes') & (comb_plsdmr_info['K27me3'] == 'na')].index)]
total3 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'Yes') & (comb_plsdmr_info['K27ac'] == 'na') & (comb_plsdmr_info['K27me3'] == 'Yes')].index)]
total4 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'Yes') & (comb_plsdmr_info['K27ac'] == 'na') & (comb_plsdmr_info['K27me3'] == 'na')].index)]
total5 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'na') & (comb_plsdmr_info['K27ac'] == 'Yes') & (comb_plsdmr_info['K27me3'] == 'na')].index)]
total6 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'na') & (comb_plsdmr_info['K27ac'] == 'na') & (comb_plsdmr_info['K27me3'] == 'Yes')].index)]
total7 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'na') & (comb_plsdmr_info['K27ac'] == 'na') & (comb_plsdmr_info['K27me3'] == 'na')].index)]

total1 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'Yes') & (comb_plsdmr_info['K27ac'] == 'Yes') & (comb_plsdmr_info['K27me3'] == 'Yes') & (comb_plsdmr_info['Type'] == 'protein_coding')].index)]
total2 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'Yes') & (comb_plsdmr_info['K27ac'] == 'Yes') & (comb_plsdmr_info['K27me3'] == 'na') & (comb_plsdmr_info['Type'] == 'protein_coding')].index)]
total3 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'Yes') & (comb_plsdmr_info['K27ac'] == 'na') & (comb_plsdmr_info['K27me3'] == 'Yes') & (comb_plsdmr_info['Type'] == 'protein_coding')].index)]
total4 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'Yes') & (comb_plsdmr_info['K27ac'] == 'na') & (comb_plsdmr_info['K27me3'] == 'na') & (comb_plsdmr_info['Type'] == 'protein_coding')].index)]
total5 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'na') & (comb_plsdmr_info['K27ac'] == 'Yes') & (comb_plsdmr_info['K27me3'] == 'na') & (comb_plsdmr_info['Type'] == 'protein_coding')].index)]
total6 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'na') & (comb_plsdmr_info['K27ac'] == 'na') & (comb_plsdmr_info['K27me3'] == 'Yes') & (comb_plsdmr_info['Type'] == 'protein_coding')].index)]
total7 = comb_plsdmr[comb_plsdmr.index.isin(comb_plsdmr_info[(comb_plsdmr_info['K4me3'] == 'na') & (comb_plsdmr_info['K27ac'] == 'na') & (comb_plsdmr_info['K27me3'] == 'na') & (comb_plsdmr_info['Type'] == 'protein_coding')].index)]

mean_total1 = pd.DataFrame({'AverageMet': total1.iloc[:, :84].mean(axis=1).values, 'Type': ['H3K4me3/H3K27ac/H3K27me3']*total1.shape[0]}, index=total1.index)
mean_total2 = pd.DataFrame({'AverageMet': total2.iloc[:, :84].mean(axis=1).values, 'Type': ['H3K4me3/H3K27ac']*total2.shape[0]}, index=total2.index)
mean_total3 = pd.DataFrame({'AverageMet': total3.iloc[:, :84].mean(axis=1).values, 'Type': ['H3K4me3/H3K27me3']*total3.shape[0]}, index=total3.index)
mean_total4 = pd.DataFrame({'AverageMet': total4.iloc[:, :84].mean(axis=1).values, 'Type': ['H3K4me3']*total4.shape[0]}, index=total4.index)
mean_total5 = pd.DataFrame({'AverageMet': total5.iloc[:, :84].mean(axis=1).values, 'Type': ['H3K27ac']*total5.shape[0]}, index=total5.index)
mean_total6 = pd.DataFrame({'AverageMet': total6.iloc[:, :84].mean(axis=1).values, 'Type': ['H3K27me3']*total6.shape[0]}, index=total6.index)
mean_total7 = pd.DataFrame({'AverageMet': total7.iloc[:, :84].mean(axis=1).values, 'Type': ['NA']*total7.shape[0]}, index=total7.index)

mean_total = pd.concat([mean_total1, mean_total2, mean_total3, mean_total4, mean_total5, mean_total6, mean_total7], axis=0)


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


