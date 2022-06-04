dfh = open("DNA_methylation_NT84_minNo55.sorted_qval.0.05.abs15.bed", 'r')
rfh = open("DMR_abs15_Hyper-Hypo_annotation.txt", 'w')
rfh.write('ID\tType\tCpGnum\tCpGdensity\n')

for i in dfh:
	line = i.strip('\n').split('\t')

	id = line[0] + ':' + line[1] + '-' + line[2]
	region_size = int(line[2]) - int(line[1])
	cpgnum = line[4]
	cpgdensity = str(round((int(cpgnum) / region_size)*1000, 6))

	if float(line[3]) > 0:
		rfh.write(id + '\t' + 'Hyper' + '\t' + cpgnum + '\t' + cpgdensity + '\n')
		rfh.flush()
	else:
		rfh.write(id + '\t' + 'Hypo' + '\t' + cpgnum + '\t' + cpgdensity + '\n')
		rfh.flush()

dfh.close()
rfh.close()

dfh = open("DNA_methylation_NT84_minNo55.sorted_qval.0.05.abs10.bed", 'r')
rfh = open("DMR_abs10_Hyper-Hypo_annotation.txt", 'w')
rfh.write('ID\tType\tCpGnum\tCpGdensity\n')

for i in dfh:
	line = i.strip('\n').split('\t')

	id = line[0] + ':' + line[1] + '-' + line[2]
	region_size = int(line[2]) - int(line[1])
	cpgnum = line[4]
	cpgdensity = str(round((int(cpgnum) / region_size)*1000, 6))

	if float(line[3]) > 0:
		rfh.write(id + '\t' + 'Hyper' + '\t' + cpgnum + '\t' + cpgdensity + '\n')
		rfh.flush()
	else:
		rfh.write(id + '\t' + 'Hypo' + '\t' + cpgnum + '\t' + cpgdensity + '\n')
		rfh.flush()

dfh.close()
rfh.close()