# %matplotlib
# %autoindent
# /mnt/mone/Project/WC300/02.RNA-seq/08.Kallisto_GencodeV24_Trimmed_Analysis/DESeq2

import gseapy as gp
import pandas as pd
import seaborn as sns

deg = pd.read_csv("Tumor.Normal.compare.csv", index_col=0)
deg_up_list, deg_down_list = list(), list()

for i in deg[(deg.loc[:, 'baseMean'] >= 50) & (deg.loc[:,'log2FoldChange'] > 1) & (deg.loc[:, 'padj'] < 0.05)].index:
	deg_up_list.append(i)


for i in deg[(deg.loc[:, 'baseMean'] >= 50) & (deg.loc[:,'log2FoldChange'] < -1) & (deg.loc[:, 'padj'] < 0.05)].index:
	deg_down_list.append(i)

hallmark_up = gp.enrichr(gene_list=deg_up_list,
			 gene_sets=['MSigDB_Hallmark_2020'],
			 organism='Human',
			 description='Hallmark_deg_up',
			 outdir='./Hallmark_deg_up',
			 cutoff=0.05)

hallmark_down = gp.enrichr(gene_list=deg_down_list,
			 gene_sets=['MSigDB_Hallmark_2020'],
			 organism='Human',
			 description='Hallmark_deg_down',
			 outdir='./Hallmark_deg_down',
			 cutoff=0.05)
