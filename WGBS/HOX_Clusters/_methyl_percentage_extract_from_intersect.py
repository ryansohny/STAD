import sys
import os
dbf = open("/mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/CpG_islands/sample_list.txt", 'r')
samples = sorted(list(map(lambda x: x.strip(), dbf.readlines())))
start = int(sys.argv[1])
end = int(sys.argv[2])
samples = samples[start-1 : end]

os.system('mkdir -p HOX_Clusters_individual')
for sample in samples:
	os.system('/usr/bin/bedtools intersect -wao -a HOX_Cluster_1000bp_slop_1kb_windows.bed \
-b /mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/' \
+ sample + '_WCmerged.bed > ./HOX_Clusters_individual/Intersect_HOX_Clusters_' + sample + '.bed')
	os.system('grep -vf /mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/CpG_islands/pattern.txt \
./HOX_Clusters_individual/Intersect_HOX_Clusters_' + sample + '.bed > ./HOX_Clusters_individual/test_' + sample + '_intersect.txt')
	os.system('/usr/bin/mv ./HOX_Clusters_individual/test_' + sample + '_intersect.txt ./HOX_Clusters_individual/Intersect_HOX_Clusters_' + sample + '.bed')

	dfh = open("./HOX_Clusters_individual/Intersect_HOX_Clusters_" + sample + ".bed", 'r') # Intersect_DMR_15002_T.bed
	rfh = open("./HOX_Clusters_individual/HOX_Clusters_" + sample + "_met.txt", 'w')

	line = dfh.readline().strip('\n').split('\t')
	id = line[0] + ':' + line[1] + '-' + line[2]
	c = 1
	if int(line[-3]) + int(line[-2]) >= 5:
		mc = int(line[-3])
		total = int(line[-3]) + int(line[-2])
	else:
		mc = 0
		total = 0

	line = dfh.readline().strip('\n').split('\t')

	while line != ['']:

		if line[0] + ':' + line[1] + '-' + line[2] == id:
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
	
			id = line[0] + ':' + line[1] + '-' + line[2]
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
