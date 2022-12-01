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


# CIMP proportional plot
ax = (pd.crosstab(dmr_t.obs['DMR Clusters'], dmr_t.obs['CIMP'], normalize=0)*100).plot.bar(stacked=True, color=['#8b0000ff', '#000080ff'], rot=0)
plt.ylabel("Proportion (%)")
ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1.0))
plt.tight_layout()
sns.despine()




# DMR

dmr_met = pd.read_table("DMR_abs10_smooth.txt", index_col=0)
dmr_met.columns = list(map(lambda x: 'X'+x, dmr_met.columns))
dmr_met = dmr_met*100

