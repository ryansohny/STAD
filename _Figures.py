####################################################################################

import pandas as pd
import seaborn as sns
clinic_info = pd.read_csv('/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/2022_WC300_clinical_information_Xadded_ver2.0.csv', index_col='Sample')

# Bisulfite conversion rate, Non-CpG Methylation Percentage
clinic_info["Conversion_rate"] = clinic_info["Conversion_rate"]*100
clinic_info["Methyl_CHG"] = clinic_info["Methyl_CHG"]*100
clinic_info["Methyl_CHH"] = clinic_info["Methyl_CHH"]*100

ax1 = plt.subplot(1,3,1)
ax2 = plt.subplot(1,3,2)
ax3 = plt.subplot(1,3,3)
sns.boxplot(data=clinic_info, x='Conversion_rate', showfliers = True, color="grey", ax=ax1)
sns.boxplot(data=clinic_info, x='Methyl_CHG', showfliers = True, color="magenta",ax=ax2)
sns.boxplot(data=clinic_info, x='Methyl_CHH', showfliers = True, color="green", ax=ax3)
ax1.set_xlabel("Bisulfite Conversion Rate (%)")
ax2.set_xlabel("Cytosine Methylated in CHG context (%)")
ax3.set_xlabel("Cytosine Methylated in CHH context (%)")
sns.despine()
plt.tight_layout()

"""

"""

# Boxplot : Mean methylation (Tumor vs Normal) ==> deprecated
df = dmr_t.obs[['mean_met_normal', 'mean_met_tumor']]
df.rename(columns={'mean_met_normal': 'Normal', 'mean_met_tumor': 'Tumor'}, inplace=True)

sns.boxplot(data=df, palette={'Normal': '#00008b','Tumor': '#8b0000'})
sns.swarmplot(data=df, color="k").set(ylabel='Mean Methylation Level')

new_df = df['Tumor'] - df['Normal']
sns.boxplot(data=new_df, color='lightsalmon')
sns.swarmplot(data=new_df, color='k').set(ylabel='(Tumor - Normal) Mean methylation')

stats.ttest_rel(df['Normal'], df['Tumor'])

# Violinplot : Mean methylation (Tumor vs Normal) ==> recent version (Figure 1A?)
clinic_info = pd.read_csv('/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/2022_WC300_clinical_information_Xadded_ver2.0.csv', index_col='Sample')
p = sns.violinplot(data=clinic_info.iloc[84:][['PercentMet_COV5_Normal', 'PercentMet_COV5_Tumor']], palette={'PercentMet_COV5_Normal':'navy', 'PercentMet_COV5_Tumor':'darkred'}, cut=0, scale="count")
p = sns.stripplot(data=clinic_info.iloc[84:][['PercentMet_COV5_Normal', 'PercentMet_COV5_Tumor']], color="black")
p.set_xticklabels(['Normal (N=84)', 'Tumor (N=84)'])
#p.set_title('Average Methylation')
p.set_ylabel("Average methylation (%)")

####################################################################################
# Violinplot : PMD counts (Tumor vs Normal)
pmd_counts = pd.read_table("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/03.WGBS/PMD/PMD_counts.txt", sep='\t')
p = sns.violinplot(data=pmd_counts, palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count").set_xticklabels(['Normal', 'Tumor'])
p = sns.stripplot(data=pmd_counts, color="k")
p.set_ylabel("The Number of PMD across genome")
plt.tight_layout()


####################################################################################
# Violinplot : PMD fraction (Tumor vs Normal)
pmdf = pd.read_table("PMD_fraction.txt", index_col=0, sep='\t')
p = sns.violinplot(data=pmdf, x='TN', y='PMD_Fraction', palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count")
p = sns.stripplot(data=pmdf, x='TN', y='PMD_Fraction', color=".3")
p.set_ylabel("Fraction of PMDs in the genome")
p.set_xticklabels(['Normal (N=84)', 'Tumor (N=84)'])
p.set_xlabel("")
plt.tight_layout()
sns.despine()

####################################################################################
# Super Enhancer Overlapped Genes
se_genes = pd.read_table("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/SE_DMR_abs15_hyper_overlapped_genes.txt")
se_genes = list(se_genes.values.flatten()) # 100개

deg_tn[deg_tn.index.isin(se_genes)] # ==> 83개
deg_tn_uplist = deg_tn[(deg_tn['baseMean'] >=10) & (deg_tn['padj'] < 0.005) & (deg_tn['log2FoldChange'] > 0)].index
deg_tn_downlist = deg_tn[(deg_tn['baseMean'] >=10) & (deg_tn['padj'] < 0.005) & (deg_tn['log2FoldChange'] < 0)].index

# deg_tn_uplist[deg_tn_uplist.isin(se_genes)] ==> 12개
# deg_tn_downlist[deg_tn_downlist.isin(se_genes)] ==> 24개
se_deg = list( deg_tn_uplist[deg_tn_uplist.isin(se_genes)] ) + list( deg_tn_downlist[deg_tn_downlist.isin(se_genes)] )

col_colors1 = ['#C0C0C0']*84 + ['#000000']*84
cmap3 = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#0057B8", "#000000", "#ffd700"])
cmap_rkg = matplotlib.colors.LinearSegmentedColormap.from_list("", ["#00FF00", "#000000", "#ff0000"])   #"R"ed blac"k" "G"reen

g = sns.clustermap(gene_vst[gene_vst.index.isin( se_deg )],
                   col_cluster=False,
                    method='ward',
                   metric='euclidean',
                   z_score=0,
                   standard_scale=None,
                   cmap=cmap_rkg,
                    col_colors=[col_colors1],
                   xticklabels=False,
                    yticklabels=True, vmin=-2.5, vmax=2.5)
g.ax_heatmap.set_yticklabels(labels=g.ax_heatmap.get_yticklabels(), fontstyle='italic')

se_deg_rna_roworder = g.dendrogram_row.reordered_ind

#  Array of codes for making SE_DEG_ALL.txt (refer to /mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/DMR/DMR_min55_new/Enrichment/Stomach_SE)
reorder = list(gene_vst[gene_vst.index.isin( se_deg )].iloc[se_deg_rna_roworder].index)

se_deg_met = pd.read_table("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/03.WGBS/NEW/SE_DEG_ALL.txt", index_col=0)
se_deg_met.columns = list(map(lambda x: 'X'+x, se_deg_met.columns))

a = pd.DataFrame(list(map(lambda x: x.split('/')[0], se_deg_met.index)), columns=['Region'], index=se_deg_met.index)
b = pd.DataFrame(list(map(lambda x: x.split('/')[1], se_deg_met.index)), columns=['GeneID'], index=se_deg_met.index)
c = pd.DataFrame(list(map(lambda x: x.split('/')[2], se_deg_met.index)), columns=['CpG'], index=se_deg_met.index)
se_deg_met_info = pd.concat([a,b,c], axis=1)
del a, b, c

reorder_iloc = list()

for i in reorder:
    reorder_iloc.append(list(se_deg_met_info['GeneID'].values).index(i))

spearmanr = list()
pvalues = list()

for i in list(range(36)):
    spearmanr.append(stats.spearmanr(se_deg_met.iloc[reorder_iloc].iloc[i], gene_vst[gene_vst.index.isin( se_deg )].iloc[se_deg_rna_roworder].iloc[i])[0])

for i in list(range(36)):
    pvalues.append(stats.spearmanr(se_deg_met.iloc[reorder_iloc].iloc[i], gene_vst[gene_vst.index.isin( se_deg )].iloc[se_deg_rna_roworder].iloc[i])[1])

from statsmodels.stats.multitest import multipletests
bh_corrected_pval = multipletests(pvals=pvalues, alpha=0.01, method='fdr_bh')[1]

spearmanr_df = pd.DataFrame(spearmanr, index=gene_vst[gene_vst.index.isin( se_deg )].iloc[se_deg_rna_roworder].index, columns=['Spearman_Rho'])
bh_corrected_pval_df = pd.DataFrame(bh_corrected_pval, index=gene_vst[gene_vst.index.isin( se_deg )].iloc[se_deg_rna_roworder].index, columns=['FDR_BH_Pval'])
df_final = pd.concat([spearmanr_df, bh_corrected_pval_df], axis=1)
df_final = df_final[df_final['FDR_BH_Pval'] < 0.01].sort_values(by='Spearman_Rho', ascending=True)

# (Again) Heatmap of SE overlapped gene expression 
g = sns.clustermap(gene_vst[gene_vst.index.isin(df_final.index)].reindex(df_final.index),
                   col_cluster=False,
                   row_cluster=False,
                   z_score=0,
                   standard_scale=None,
                   cmap=cmap_rkg,
                   col_colors=[col_colors1],
                   xticklabels=False,
                   yticklabels=True, vmin=-2.5, vmax=2.5)
g.ax_heatmap.set_yticklabels(labels=g.ax_heatmap.get_yticklabels(), fontstyle='italic')

sorter = list(df_final.index) # https://stackoverflow.com/questions/23482668/sorting-by-a-custom-list-in-pandas
sorterIndex = dict(zip(sorter, range(len(sorter))))
new_se_deg_met_info = se_deg_met_info[se_deg_met_info['GeneID'].isin(df_final.index)]
new_se_deg_met_info['New_order'] = new_se_deg_met_info['GeneID'].map(sorterIndex)
new_se_deg_met_info.sort_values(by='New_order', inplace=True)


# DNA methylation of Super enhancers overlapped with DEG
g = sns.clustermap(se_deg_met[se_deg_met.index.isin(new_se_deg_met_info.index)].reindex(new_se_deg_met_info.index),
                   col_cluster=False,
                   row_cluster=False,
                   cmap=cmap,
                   z_score=None,
                   standard_scale=0,
                   col_colors=[col_colors1],
                   xticklabels=False,
                   yticklabels=False)

sns.clustermap(df_final['Spearman_Rho'], col_cluster=False, row_cluster=False, vmin=-0.65, vmax=0.22)

# Homeobox (HOX) cluster DNA methylation
hox = pd.read_table("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/03.WGBS/NEW/HOX_Clusters_ALL.txt", index_col=0)
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
                   z_score=0,
                   standard_scale=None,
                   cmap=cmap,
                   xticklabels=False,
                   yticklabels=False,
                   col_colors=None,
                   row_colors=None,
                   cbar_kws={'label': 'DNA methylation'},
                   vmin=-3, vmax=2) # Original cmap='gnuplot2'

####################################################################################
## Determining CIMP Tumors

# Call CpGi methylation
#cpgi_met = pd.read_table("CpGi_ALL.txt", index_col=0)
cpgi_met = pd.read_table("CpGi_smooth.txt", index_col=0)
cpgi_met = cpgi_met * 100
# Column name change
cpgi_met.columns = list(map(lambda x: 'X'+x, cpgi_met.columns))
#cpgi_met.index = list(map(lambda x: '/'.join(x.split('/')[:2]), cpgi_met.index))

# Discard missing CpGi DNA methylation rows & pick CpGi sites where Normal DNA methylation < 40
# cpgi_met = cpgi_met[cpgi_met.isna().sum(axis=1) == 0][cpgi_met[cpgi_met.isna().sum(axis=1) == 0].iloc[:, :84].mean(axis=1) < 40]
cpgi_met = cpgi_met[cpgi_met.iloc[:, :84].mean(axis=1) < 40]

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

# Association between PMD methylation and CIMP-CGI DNA methylation
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

####################################################################################
# Violinplot : FOXA2 & CASZ1 (Tumor vs Normal)

foxa2 = pd.concat([pd.DataFrame(gene_vst.loc[['FOXA2']].iloc[:, :84].T.values.flatten(), columns=['Normal']), pd.DataFrame(gene_vst.loc[['FOXA2']].iloc[:, 84:].T.values.flatten(), columns=['Tumor'])], axis=1).set_index(gene_vst.columns[84:])
p = sns.boxplot(data=foxa2, palette={'Normal':'navy', 'Tumor':'darkred'}, width=0.8, showfliers=True)
p = sns.stripplot(data=foxa2, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
p.set_ylabel("Gene Expression")
plt.tight_layout()
sns.despine()

casz1 = pd.concat([pd.DataFrame(gene_vst.loc[['CASZ1']].iloc[:, :84].T.values.flatten(), columns=['Normal']), pd.DataFrame(gene_vst.loc[['CASZ1']].iloc[:, 84:].T.values.flatten(), columns=['Tumor'])], axis=1).set_index(gene_vst.columns[84:])
p = sns.boxplot(data=casz1, palette={'Normal':'navy', 'Tumor':'darkred'}, width=0.8, showfliers=True)
p = sns.stripplot(data=casz1, jitter=True, marker='o', color='black', size=1.5, alpha=0.2)
p.set_ylabel("Gene Expression")
plt.tight_layout()
sns.despine()

# Epimutation Burden-related analyses

dmr_met = pd.read_table("DMR_abs10_smooth.txt", index_col=0)
dmr_met.columns = list(map(lambda x: 'X'+x, dmr_met.columns))
dmr_met = dmr_met*100

dmr = sc.AnnData(dmr_met.iloc[:, 84:].T) # Only Tumor
dmr.raw = dmr
dmr.layers['Percent_met'] = dmr.X

dmr.obs = clinic_info.iloc[84:, :].copy()

sc.pp.scale(dmr)
sc.tl.pca(dmr, n_comps=83, zero_center=True)
#sc.pl.pca(dmr, color='TN', palette={'Normal':'midnightblue', 'Tumor':'darkred'}, annotate_var_explained=True, size=100)
sc.pl.pca(dmr, color=['EBV'], palette={'Negative':'midnightblue', 'Positive':'darkred'}, annotate_var_explained=True, size=300)
sns.despine()
g = sc.pl.pca(dmr, color=['EpiBurden'], color_map=cmap, annotate_var_explained=True, size=300, show=False)
g.set_title("Epimutation Burden")
sns.despine()

PC1 = pd.DataFrame(list(map(lambda x: dmr.obsm['X_pca'][x][0], list(range(0,84)))), columns=['PC1'], index=dmr.obs.index)
df = pd.concat([PC1, dmr.obs['EpiBurden']], axis=1)
g = sns.lmplot(data=df, x="PC1", y="EpiBurden")
g.set_ylabels("Epimutation Burden")
stats.pearsonr(df['EpiBurden'], df['PC1'])

#sc.pl.pca(dmr, color='TN', add_outline=True, size=100, palette={'Normal':'Blue', 'Tumor':'Red'})
#sc.pl.pca_variance_ratio(dmr, log=True)
pca_variance = pd.DataFrame(dmr.uns['pca']['variance_ratio'], index=list(map(lambda x: 'PC' + str(x), list(range(1,101)))), columns=['Variance_ratio'])
np.sum(pca_variance.values.flatten()[:5])
# 0.7084484


####################################################################################
# Violinplot : DNMT (Tumor vs Normal)

dnmt1 = pd.concat([pd.DataFrame(gene_vst.loc[['DNMT1']].iloc[:, :84].T.values.flatten(), columns=['Normal']), pd.DataFrame(gene_vst.loc[['DNMT1']].iloc[:, 84:].T.values.flatten(), columns=['Tumor'])], axis=1).set_index(gene_vst.columns[84:])
sns.violinplot(data=dnmt1, palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count").set_title('DNMT1')

dnmt3a = pd.concat([pd.DataFrame(gene_vst.loc[['DNMT3A']].iloc[:, :84].T.values.flatten(), columns=['Normal']), pd.DataFrame(gene_vst.loc[['DNMT3A']].iloc[:, 84:].T.values.flatten(), columns=['Tumor'])], axis=1).set_index(gene_vst.columns[84:])
sns.violinplot(data=dnmt3a, palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count").set_title('DNMT3A')

dnmt3b = pd.concat([pd.DataFrame(gene_vst.loc[['DNMT3B']].iloc[:, :84].T.values.flatten(), columns=['Normal']), pd.DataFrame(gene_vst.loc[['DNMT3B']].iloc[:, 84:].T.values.flatten(), columns=['Tumor'])], axis=1).set_index(gene_vst.columns[84:])
sns.violinplot(data=dnmt3b, palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count").set_title('DNMT3B')

# Multiple DNMT
dnmt1 = pd.concat([pd.DataFrame(gene_vst.loc['DNMT1'].values, columns=['Expression']), pd.DataFrame(['DNMT1']*168, columns=['DNMT']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)
dnmt3a = pd.concat([pd.DataFrame(gene_vst.loc['DNMT3A'].values, columns=['Expression']), pd.DataFrame(['DNMT3A']*168, columns=['DNMT']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)
dnmt3b = pd.concat([pd.DataFrame(gene_vst.loc['DNMT3B'].values, columns=['Expression']), pd.DataFrame(['DNMT3B']*168, columns=['DNMT']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)
dnmt = pd.concat([dnmt1, dnmt3a, dnmt3b], ignore_index=True)
p = sns.violinplot(data=dnmt, x='DNMT', y='Expression', hue='TN', palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count")
p = sns.stripplot(data=dnmt, x='DNMT', y='Expression', hue='TN', color=".3", size=3, dodge=True)
sns.despine()
p.set_xlabel("")
p.set_ylabel("DNMT Gene Expression")
plt.legend(frameon=False, fontsize=15)
p.set_xticklabels(p.get_xticklabels(), fontstyle='italic')
plt.tight_layout()

# Multiple GATA
gata5 = pd.concat([pd.DataFrame(gene_vst.loc['GATA5'].values, columns=['Expression']), pd.DataFrame(['GATA5']*168, columns=['GATA']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)
gata6 = pd.concat([pd.DataFrame(gene_vst.loc['GATA6'].values, columns=['Expression']), pd.DataFrame(['GATA6']*168, columns=['GATA']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)
gata4 = pd.concat([pd.DataFrame(gene_vst.loc['GATA4'].values, columns=['Expression']), pd.DataFrame(['GATA4']*168, columns=['GATA']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)
gata2 = pd.concat([pd.DataFrame(gene_vst.loc['GATA2'].values, columns=['Expression']), pd.DataFrame(['GATA2']*168, columns=['GATA']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)

gata = pd.concat([gata5, gata6, gata4, gata2], ignore_index=True)
sns.violinplot(data=gata, x='GATA', y='Expression', hue='TN', palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count")
sns.stripplot(data=gata, x='GATA', y='Expression', hue='TN', color=".3", size=3, dodge=True)

# 5 core EMT TFs
snai1 = pd.concat([pd.DataFrame(rna.loc['SNAI1'].values, columns=['Expression']), pd.DataFrame(['SNAI1']*168, columns=['EMT_TFs']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)
twist1 = pd.concat([pd.DataFrame(rna.loc['TWIST1'].values, columns=['Expression']), pd.DataFrame(['TWIST1']*168, columns=['EMT_TFs']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)
snai2 = pd.concat([pd.DataFrame(rna.loc['SNAI2'].values, columns=['Expression']), pd.DataFrame(['SNAI2']*168, columns=['EMT_TFs']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)
zeb1 = pd.concat([pd.DataFrame(rna.loc['ZEB1'].values, columns=['Expression']), pd.DataFrame(['ZEB1']*168, columns=['EMT_TFs']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)
zeb2 = pd.concat([pd.DataFrame(rna.loc['ZEB2'].values, columns=['Expression']), pd.DataFrame(['ZEB2']*168, columns=['EMT_TFs']), pd.DataFrame(['Normal']*84 + ['Tumor']*84, columns=['TN'])], axis=1)

emt_tfs = pd.concat([snai1, twist1, snai2, zeb1, zeb2], ignore_index=True)
sns.violinplot(data=emt_tfs, x='EMT_TFs', y='Expression', hue='TN', palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count")
sns.stripplot(data=emt_tfs, x='EMT_TFs', y='Expression', hue='TN', color=".3", size=3, dodge=True)



df = pd.concat([pd.DataFrame(gb_met.loc['chr7:19020990-19117672/TWIST1/ENSG00000122691/886'][:84].values, columns=['Normal']), pd.DataFrame(gb_met.loc['chr7:19020990-19117672/TWIST1/ENSG00000122691/886'][84:].values, columns=['Tumor'])], axis=1).set_index(gb_met.columns[84:])
p = sns.violinplot(data=df, palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count")
p = sns.stripplot(data=df, color=".3")
p.set_ylabel("Gene Body Methylation (%)")
p.set_title('TWIST1')
stats.ttest_rel(df['Normal'].values, df['Tumor'].values)


gb_met = pd.read_table("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/03.WGBS/NEW/GeneBody_ALL.txt", index_col=0)
gb_met_info = pd.DataFrame(list(zip(list(map(lambda x: x.split('/')[0], gb_met.index)), list(map(lambda x: x.split('/')[1], gb_met.index)))), index=gb_met.index, columns=['Region', 'GeneID'])
hdac9_met = gb_met.loc[gb_met_info[gb_met_info['GeneID'] == 'HDAC9'].index]
hdac9_met_df = pd.concat([pd.DataFrame(hdac9_met.iloc[:, :84].T.values, columns=['Normal']), pd.DataFrame(hdac9_met.iloc[:, 84:].T.values, columns=['Tumor'])], axis=1).set_index(gb_met.columns[84:])
p = sns.violinplot(data=hdac9_met_df, palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count")
p = sns.stripplot(data=hdac9_met_df, color=".3")
p.set_ylabel("Gene Body Methylation (%)")
p.set_title('HDAC9', style='italic')
plt.tight_layout()
stats.ttest_rel(hdac9_met_df['Normal'], hdac9_met_df['Tumor'])



# TCGA subtypes
tcga_info = pd.read_csv('/data/Projects/phenomata/01.Projects/12.Thesis/Part2/tcga_subtypes.csv', index_col='Sample')
dmr.obs['TCGA Subtypes'] = tcga_info['Type']

df = pd.concat([dmr.obs['EpiBurden'], dmr.obs['TCGA Subtypes']], axis=1)
df.rename(columns={"EpiBurden": "Epimutation Burden"}, inplace=True)

p = sns.boxplot(data=df, x='TCGA Subtypes', y='Epimutation Burden', palette={'EBV': '#DC143C', 'MSI': '#4169E1', 'GS/CIN': '#9370DB'}, width=0.8, showfliers=True, order=['EBV', 'MSI', 'GS/CIN'])
p = sns.stripplot(data=df, x='TCGA Subtypes', y='Epimutation Burden', jitter=True, marker='o', color='black', size=2.0, alpha=0.5, order=['EBV', 'MSI', 'GS/CIN'])
p.set_xticklabels(['EBV (N=4)', 'MSI (N=10)', 'GS/CIN (N=70)'])
plt.tight_layout()
sns.despine()

stats.mannwhitneyu(df[df['TCGA Subtypes'] == 'EBV']['Epimutation Burden'], df[df['TCGA Subtypes'] == 'MSI']['Epimutation Burden'])
MannwhitneyuResult(statistic=27.0, pvalue=0.3736263736263737)

stats.mannwhitneyu(df[df['TCGA Subtypes'] == 'EBV']['Epimutation Burden'], df[df['TCGA Subtypes'] == 'GS/CIN']['Epimutation Burden'])
MannwhitneyuResult(statistic=199.0, pvalue=0.16814499237806202)

stats.mannwhitneyu(df[df['TCGA Subtypes'] == 'MSI']['Epimutation Burden'], df[df['TCGA Subtypes'] == 'GS/CIN']['Epimutation Burden'])
MannwhitneyuResult(statistic=356.0, pvalue=0.9362267365375121)
