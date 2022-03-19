import seaborn as sns
import pandas as pd
import matplotlib
import scipy.stats as stats

rna = pd.read_table("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/STAD_SNUH_vst_new.txt", index_col=0, sep=' ')
pmd_genes = list(map(lambda x: x.strip('\n'), open("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/PMD_overlapped_genes.txt", 'r').readlines()))
pmd_proteingenes = list(map(lambda x: x.strip('\n'), open("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/PMD_overlapped_protein_genes.txt", 'r').readlines()))
proteingenes = list(map(lambda x: x.strip('\n'), open("/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/GENCODE_V24/protein_genes.txt", 'r').readlines()))

# Expression level of all genes inside PMD and outside PMD
no_pmd_genes_exp = pd.concat([rna.loc[~rna.index.isin(pmd_genes)].iloc[:,:84].mean(axis=1), rna.loc[~rna.index.isin(pmd_genes)].iloc[:,84:].mean(axis=1)], axis=1)
no_pmd_genes_exp.columns = ['Normal', 'Tumor']
p1 = sns.violinplot(data=no_pmd_genes_exp, palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count")
p1.set_title('Outside PMDs')
p1.set_ylabel("Average expression across samples")
plt.tight_layout()

stats.ttest_rel(no_pmd_genes_exp['Normal'], no_pmd_genes_exp['Tumor'])

pmd_genes_exp = pd.concat([rna.loc[rna.index.isin(pmd_genes)].iloc[:,:84].mean(axis=1), rna.loc[rna.index.isin(pmd_genes)].iloc[:,84:].mean(axis=1)], axis=1)
pmd_genes_exp.columns = ['Normal', 'Tumor']
p2 = sns.violinplot(data=pmd_genes_exp, palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count")
p2.set_title('Inside PMDs')
p2.set_ylabel("Average expression across samples")
plt.tight_layout()

stats.ttest_rel(pmd_genes_exp['Normal'], pmd_genes_exp['Tumor'])

# Expression level of "coding genes" inside PMD and outside PMD
rna_protein = rna.loc[rna.index.isin(proteingenes)]

no_pmd_proteingenes_exp = pd.concat([rna_protein.loc[~rna_protein.index.isin(pmd_proteingenes)].iloc[:,:84].mean(axis=1), rna_protein.loc[~rna_protein.index.isin(pmd_proteingenes)].iloc[:,84:].mean(axis=1)], axis=1)
no_pmd_proteingenes_exp.columns = ['Normal', 'Tumor']
p1 = sns.violinplot(data=no_pmd_proteingenes_exp, palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count")
p1.set_title('Outside PMDs')
p1.set_ylabel("Average coding gene expression across samples")
plt.tight_layout()

stats.ttest_rel(no_pmd_proteingenes_exp['Normal'], no_pmd_proteingenes_exp['Tumor'])

pmd_proteingenes_exp = pd.concat([rna_protein.loc[rna_protein.index.isin(pmd_proteingenes)].iloc[:,:84].mean(axis=1), rna_protein.loc[rna_protein.index.isin(pmd_proteingenes)].iloc[:,84:].mean(axis=1)], axis=1)
pmd_proteingenes_exp.columns = ['Normal', 'Tumor']
p2 = sns.violinplot(data=pmd_proteingenes_exp, palette={'Normal':'navy', 'Tumor':'darkred'}, cut=0, scale="count")
p2.set_title('Inside PMDs')
p2.set_ylabel("Average coding gene expression across samples")
plt.tight_layout()

stats.ttest_rel(pmd_proteingenes_exp['Normal'], pmd_proteingenes_exp['Tumor'])

