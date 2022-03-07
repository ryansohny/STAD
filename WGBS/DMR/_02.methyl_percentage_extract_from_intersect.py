import sys
import os
dbf = open("/mnt/mone/Project/WC300/03.WGBS/02.DNA_methylation/Re-Analysis/metilene/DMR/sample_20211221.txt", 'r')
samples = sorted(list(map(lambda x: x.strip(), dbf.readlines())))
start = int(sys.argv[1])
end = int(sys.argv[2])
samples = samples[start-1 : end]

os.system('mkdir -p DMR_individual')
for sample in samples:
	dfh = open("/mnt/mone/Project/WC300/03.WGBS/02.DNA_methylation/Re-Analysis/metilene/" + sample + '_WCmerged.txt', 'r')
	rfh = open('./DMR_individual/' + sample + '_WCmerged.bed', 'w')
        
	for line in dfh:
		line = line.strip('\n').split('\t')
		rfh.write(line[1].split(':')[0] + '\t' + str(int(line[1].split(':')[1])-1) + '\t' + line[1].split(':')[1] + '\t' + line[3] + '\t' + line[4] + '\n')
		rfh.flush()
	dfh.close()
	rfh.close()
        
	os.system('/usr/bin/bedtools intersect -wao -a DNA_methylation_NT84_min55.sorted_qval.0.05.abs15.bed -b ./DMR_individual/' + sample + '_WCmerged.bed > ./DMR_individual/Intersect_DMR_' + sample + '.bed')
	os.system('/usr/bin/rm ./DMR_individual/' + sample + '_WCmerged.bed')

	dfh = open("./DMR_individual/Intersect_DMR_" + sample + ".bed", 'r') # Intersect_DMR_15002_T.bed
	rfh = open("./DMR_individual/DMR_" + sample + "_met.txt", 'w')

	line = dfh.readline().strip('\n').split('\t')
	id = line[0] + ':' + line[1] + '-' + line[2]
	if int(line[-3]) + int(line[-2]) >= 5:
		mc = int(line[-3])
		total = int(line[-3]) + int(line[-2])
	else:
		mc = 0
		total = 0

	line = dfh.readline().strip('\n').split('\t')

	while line != ['']:

		if line[0] + ':' + line[1] + '-' + line[2] == id:
			if int(line[-3]) + int(line[-2]) >=5:
				mc += int(line[-3])
				total += (int(line[-3]) + int(line[-2]))
		else:
			if total != 0:
				rfh.write(id + '\t' + str(round((mc/total) * 100, 4)) + '\n')
				rfh.flush()
			else:
				rfh.write(id + '\t' + '' + '\n')
	
			id = line[0] + ':' + line[1] + '-' + line[2]
			if int(line[-3]) + int(line[-2]) >= 5:
				mc = int(line[-3])
				total = int(line[-3]) + int(line[-2])
			else:
				mc = 0
				total = 0
		line = dfh.readline().strip('\n').split('\t')

	# Last DMR
	if total != 0:
		rfh.write(id + '\t' + str(round((mc/total) * 100, 4)) + '\n')
		rfh.flush()
	else:
		rfh.write(id + '\t' + '' + '\n')
	
	dfh.close()
	rfh.close()
	os.system('/usr/bin/rm ./DMR_individual/Intersect_DMR_' + sample + '.bed')
