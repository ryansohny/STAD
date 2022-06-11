# se_genes = pd.read_table("/mnt/data/Projects/phenomata/01.Projects/08.StomachCancer_backup/02.RNA-seq/SE_DMR_abs15_hyper_overlapped_genes.txt")
# se_genes = list(se_genes.values.flatten()) # 100개

# deg_tn[deg_tn.index.isin(se_genes)] # ==> 70개

# deg_tn_uplist = deg_tn[(deg_tn['baseMean'] >=10) & (deg_tn['padj'] < 0.005) & (deg_tn['log2FoldChange'] > 0)].index
# deg_tn_downlist = deg_tn[(deg_tn['baseMean'] >=10) & (deg_tn['padj'] < 0.005) & (deg_tn['log2FoldChange'] < 0)].index

# deg_tn_uplist[deg_tn_uplist.isin(se_genes)] ==> 10개
# deg_tn_downlist[deg_tn_downlist.isin(se_genes)] ==> 20개
# se_deg = list( deg_tn_uplist[deg_tn_uplist.isin(se_genes)] ) + list( deg_tn_downlist[deg_tn_downlist.isin(se_genes)] )

# se_deg ==> SE_DMR_abs15_hyper_overlapped_genes_DEG.txt

# HOXA13
# HOXA10
# HOXA11
# TPX2
# HOXA9

import os
se_deg_db = list(map(lambda x: x.strip('\n'), open("SE_DMR_abs15_hyper_overlapped_genes_DEG.txt", 'r').readlines()))

# cat Stomach_SE_DMR_abs15_Hyper.bed | awk '{print $1"\t"$2"\t"$3"\t"$4}' | uniq > Stomach_SE_DMR_abs15_Hyper_redundant_removal.bed

dfh = open("Stomach_SE_DMR_abs15_Hyper_redundant_removal.bed", 'r')
rfh = open("SE_DEG_overlapped.bed", 'w')

for i in dfh:
	line = i.strip('\n').split('\t')
	for gene in se_deg_db:
		if gene in line[3]:
			rfh.write('\t'.join(line[:3]) + '\t' + gene + '\n')
			rfh.flush()