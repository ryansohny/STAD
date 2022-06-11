#import os
#os.system('bedtools intersect -wao \
#-a Stomach_SE_modified_liftover.bed \
#-b DNA_methylation_NT84_minNo55.sorted_qval.0.05.abs15_hyper.bed \
#| grep -vf pattern.txt \
#> Stomach_SE_DMR_abs15_Hyper.bed')

#os.system("/usr/bin/awk '{print $4}' Stomach_SE_DMR_abs15_Hyper.bed > pre_genes.txt")

dfh = open("pre_genes.txt", 'r')
rfh = open("SE_DMR_abs15_hyper_overlapped_genes.txt", 'w')
rfh.write('ID\n')
db = list()
for i in dfh:
	
	line = i.strip('\n').split('/')
	for j in line[1:]:
		if j == '':
			pass
		elif ',' in j:
			j = j.split(',')
			for k in j:
				if k not in db:
					db.append(k)
		else:
			db.append(j)
#print(list(set(sorted(db))))
# ['ABHD4', 'AHNAK', 'AHNAK', 'AHNAK', 'BCAT2', 'BCL2L1', 'BCL2L1', 'BCL2L1', 'BMF', 'BMF', 'BMF', 'C17orf86', 'C7orf50', 'C9ORF167', 'C9orf167', 'C9orf167', 'C9orf169', 'C9orf173', 'CASZ1', 'CASZ1', 'CASZ1', 'CASZ1', 'CASZ1', 'CASZ1', 'CASZ1', 'CERK', 'CHST3', 'CHST3', 'CHST3', 'COBRA1', 'COX4I2', 'CUEDC1', 'CUEDC1', 'CUEDC1', 'DAD1', 'DAD1', 'DAD1', 'EEF1G', 'EGFR', 'EGFR', 'EGFR', 'EPHB3', 'EPHB3', 'EPHB3', 'FAM165B', 'FAM165B', 'FAM165B', 'FAM166A', 'FAM169A', 'FAM46C', 'FAM46C', 'FAM46C', 'FAM46C', 'FAM46C', 'FAM46C', 'FGF21', 'FGF21', 'FGF21', 'FLJ35946', 'FLJ35946', 'FLJ35946', 'FN3K', 'FN3KRP', 'FOXA2', 'FUT1', 'FUT2', 'FUT8', 'FUT8', 'GCNT4', 'GCNT4', 'GDF7', 'GET4', 'GET4', 'GET4', 'GNAQ', 'GNAQ', 'GNAQ', 'GPER', 'GPER', 'GPER', 'GPER', 'GPER', 'GPR146', 'H6PD', 'H6PD', 'H6PD', 'HOXA1', 'HOXA10', 'HOXA11', 'HOXA11-AS1', 'HOXA13', 'HOXA2', 'HOXA3', 'HOXA3', 'HOXA3', 'HOXA4', 'HOXA4', 'HOXA4', 'HOXA5', 'HOXA6', 'HOXA7', 'HOXA9', 'HS1BP3', 'HS1BP3', 'HS1BP3', 'HSD17B14', 'IGF1R', 'IGF1R', 'IGF1R', 'ITPKB', 'ITPKB', 'ITPKB', 'IZUMO1', 'KCNE2', 'LOC148696', 'LOC148696', 'LOC148696', 'LOC643596', 'MAMSTR', 'MIR196B', 'MIR29B2', 'MIR29C', 'MIR4251', 'MIR4251', 'MIR4251', 'MRPS23', 'MTA2', 'NCRNA00261', 'NCRNA00261', 'NCRNA00261', 'NCRNA00261', 'NCRNA00261', 'NDRG1', 'NDRG1', 'NDRG1', 'NEDD4L', 'NEDD4L', 'NEDD4L', 'NRARP', 'PEX14', 'PIWIL2', 'PPP3CC', 'PRDM16', 'PTCH1', 'PTCH1', 'PTCH1', 'PTCH1', 'RAB40B', 'RAB40B', 'RAB40B', 'RASA2', 'RASA2', 'RASIP1', 'RBM47', 'RBM47', 'RBM47', 'RGMB', 'RNF208', 'SCARNA16', 'SEC14L1', 'SEC14L1', 'SEC14L1', 'SEC14L1', 'SEC14L1', 'SEC14L1', 'SEC14L1', 'SEC14L1', 'SEC14L1', 'SEC14L1', 'SEC14L1', 'SEC14L1', 'SLC34A3', 'SLC39A14', 'SLC39A14', 'SLC39A14', 'SLC9A4', 'SLC9A4', 'SLC9A4', 'SLC9A4', 'SLC9A4', 'SLC9A4', 'SPSB1', 'SRP14', 'SUN1', 'SUN1', 'TBC1D22A', 'TBC1D22A', 'TBC1D22A', 'TBCD', 'TNS3', 'TNS3', 'TNS3', 'TPCN2', 'TPCN2', 'TPX2', 'TUBB2C', 'TUT1', 'VEZF1', 'WDR45L', 'ZBTB38', 'ZFAND2A']
db.remove('C9ORF167')
for i in sorted(list(set(db))):
	rfh.write(i + '\n')
dfh.close()
rfh.close()
