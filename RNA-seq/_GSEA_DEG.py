# %matplotlib
# %autoindent
# /mnt/mone/Project/WC300/02.RNA-seq/08.Kallisto_GencodeV24_Trimmed_Analysis/DESeq2

import pandas as pd
import seaborn as sns
import gseapy as gp
from gseapy.plot import barplot, dotplot


deg = pd.read_csv("Tumor.Normal.compare.csv", index_col=0)
deg_up_list, deg_down_list = list(), list()

for i in deg[(deg.loc[:, 'baseMean'] >= 10) & (deg.loc[:,'log2FoldChange'] > 1) & (deg.loc[:, 'padj'] < 0.05)].index:
	deg_up_list.append(i)


for i in deg[(deg.loc[:, 'baseMean'] >= 10) & (deg.loc[:,'log2FoldChange'] < -1) & (deg.loc[:, 'padj'] < 0.05)].index:
	deg_down_list.append(i)

hallmark_up = gp.enrichr(gene_list=deg_up_list,
			 gene_sets=['MSigDB_Hallmark_2020'],
			 organism='Human',
			 description='Hallmark_deg_up',
			 outdir='./Hallmark_deg_up',
			 cutoff=0.05, 
			 no_plot=True)
dotplot(hallmark_up.res2d, title='HALLMARK DEG UP', column='Adjusted P-value', cutoff=0.05, top_term=10, cmap='viridis_r', ofname='./Hallmark_deg_up/Dotplot_HALLMARK_DEG_UP.pdf')

hallmark_down = gp.enrichr(gene_list=deg_down_list,
			 gene_sets=['MSigDB_Hallmark_2020'],
			 organism='Human',
			 description='Hallmark_deg_down',
			 outdir='./Hallmark_deg_down',
			 cutoff=0.05, no_plot=True)
dotplot(hallmark_down.res2d, title='HALLMARK DEG DOWN', column='Adjusted P-value', cutoff=0.05, top_term=10, cmap='viridis_r', ofname='./Hallmark_deg_down/Dotplot_HALLMARK_DEG_DOWN.pdf')

gobp_up = gp.enrichr(gene_list=deg_up_list,
		     gene_sets=['GO_Biological_Process_2021'],
		     organism='Human',
		     description='GOBP_deg_up',
		     outdir='./GOBP_deg_up',
		     cutoff=0.05, no_plot=True)
dotplot(gobp_up.res2d, title='GOBP DEG UP', column='Adjusted P-value', cutoff=0.05, top_term=10, cmap='viridis_r', ofname='./GOBP_deg_up/Dotplot_GOBP_DEG_UP.pdf')

gobp_down = gp.enrichr(gene_list=deg_down_list,
		       gene_sets=['GO_Biological_Process_2021'],
		       organism='Human',
		       description='GOBP_deg_down',
		       outdir='./GOBP_deg_down',
		       cutoff=0.05, no_plot=True)
dotplot(gobp_down.res2d, title='GOBP DEG DOWN', column='Adjusted P-value', cutoff=0.05, top_term=10, cmap='viridis_r', ofname='./GOBP_deg_down/Dotplot_GOBP_DEG_DOWN.pdf')

kegg_up = gp.enrichr(gene_list=deg_up_list,
		     gene_sets=['KEGG_2021_Human'],
		     organism='Human',
		     description='KEGG_deg_up',
		     outdir='./KEGG_deg_up',
		     cutoff=0.05, no_plot=True)
dotplot(kegg_up.res2d, title='KEGG DEG UP', column='Adjusted P-value', cutoff=0.05, top_term=10, cmap='viridis_r', ofname='./KEGG_deg_up/Dotplot_KEGG_DEG_UP.pdf')

kegg_down = gp.enrichr(gene_list=deg_down_list,
		       gene_sets=['KEGG_2021_Human'],
		       organism='Human',
		       description='KEGG_deg_down',
		       outdir='./KEGG_deg_down',
		       cutoff=0.05, no_plot=True)
dotplot(kegg_down.res2d, title='KEGG DEG DOWN', column='Adjusted P-value', cutoff=0.05, top_term=10, cmap='viridis_r', ofname='./KEGG_deg_down/Dotplot_KEGG_DEG_DOWN.pdf')


