import sys
import os
dbf = open("/mnt/mone/Project/WC300/03.WGBS/02.DNA_methylation/Re-Analysis/metilene/DMR/sample_20211221.txt", 'r')
samples = sorted(list(map(lambda x: x.strip(), dbf.readlines())))
start = int(sys.argv[1])
end = int(sys.argv[2])
samples = samples[start-1 : end]

os.system('mkdir -p Promoter_up1000down500_individual')
for sample in samples:
	dfh = open("/mnt/mone/Project/WC300/03.WGBS/02.DNA_methylation/Re-Analysis/metilene/" + sample + '_WCmerged.txt', 'r')
	rfh = open('./Promoter_up1000down500_individual/' + sample + '_WCmerged.bed', 'w')
        
	for line in dfh:
		line = line.strip('\n').split('\t')
		rfh.write(line[1].split(':')[0] + '\t' + str(int(line[1].split(':')[1])-1) + '\t' + line[1].split(':')[1] + '\t' + line[3] + '\t' + line[4] + '\n')
		rfh.flush()
	dfh.close()
	rfh.close()
        
	os.system('/usr/bin/bedtools intersect -wao -a TSS_gencode.v38lift37.proteingenes.sorted_up1000down500.bed -b ./Promoter_up1000down500_individual/' + sample + '_WCmerged.bed > ./Promoter_up1000down500_individual/Intersect_DMR_up1000down500_' + sample + '.bed')
	# /usr/bin/bedtools intersect -wao -a TSS_gencode.v38lift37.proteingenes.sorted_up500down500.bed -b ./DMR_individual/14002_N_WCmerged.bed > ./DMR_individual/Intersect_DMR_14002_N.bed
#	os.system('/usr/bin/rm ./Promoter_up1000down500_individual/' + sample + '_WCmerged.bed')
	os.system('grep -vf /mnt/mone/Project/WC300/03.WGBS/02.DNA_methylation/Re-Analysis/metilene/Promoter_methylation/Gencode_V38/pattern.txt ./Promoter_up1000down500_individual/Intersect_DMR_up1000down500_' + sample + '.bed > ./Promoter_up1000down500_individual/test_' + sample + '_intersect.txt')
	os.system('/usr/bin/mv ./Promoter_up1000down500_individual/test_' + sample + '_intersect.txt ./Promoter_up1000down500_individual/Intersect_DMR_up1000down500_' + sample + '.bed')

	dfh = open("./Promoter_up1000down500_individual/Intersect_DMR_up1000down500_" + sample + ".bed", 'r') # Intersect_DMR_15002_T.bed
	rfh = open("./Promoter_up1000down500_individual/Promoter_up1000down500_" + sample + "_met.txt", 'w')

	line = dfh.readline().strip('\n').split('\t')
	id = line[0] + ':' + line[1] + '-' + line[2] + '/' + line[3] + '/' + line[4] 
	c = 1
	if int(line[-3]) + int(line[-2]) >= 5:
		mc = int(line[-3])
		total = int(line[-3]) + int(line[-2])
	else:
		mc = 0
		total = 0

	line = dfh.readline().strip('\n').split('\t')

	while line != ['']:

		if line[0] + ':' + line[1] + '-' + line[2] + '/' + line[3] + '/' + line[4] == id:
			c += 1
			if int(line[-3]) + int(line[-2]) >=5:
				mc += int(line[-3])
				total += (int(line[-3]) + int(line[-2]))
		else:
			if total != 0:
				rfh.write(id + '/' + str(c) + '\t' + str(round((mc/total) * 100, 4)) + '\n')
				rfh.flush()
			else:
				rfh.write(id + '/' + str(c) + '\t' + '' + '\n')
	
			id = line[0] + ':' + line[1] + '-' + line[2] + '/' + line[3] + '/' + line[4]
			c = 1
			if int(line[-3]) + int(line[-2]) >= 5:
				mc = int(line[-3])
				total = int(line[-3]) + int(line[-2])
			else:
				mc = 0
				total = 0
		line = dfh.readline().strip('\n').split('\t')

	# Last DMR
	if total != 0:
		rfh.write(id + '/' + str(c) + '\t' + str(round((mc/total) * 100, 4)) + '\n')
		rfh.flush()
	else:
		rfh.write(id + '/' + str(c) + '\t' + '' + '\n')
	
	dfh.close()
	rfh.close()
	os.system('/usr/bin/rm ./Promoter_up1000down500_individual/Intersect_DMR_up1000down500_' + sample + '.bed')
