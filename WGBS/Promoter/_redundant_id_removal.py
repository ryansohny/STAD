#/mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/Promoter_methylation/Gencode_V24/ENCODE/GENCODE_tag_basic/MERGED_PLS_No_PLS
import os
os.system("awk '{print $6}' PLS_annotated_table_full.txt | sort | uniq -cd > redundant_transcript.txt")

with open("redundant_transcript.txt", 'r') as redf:
	redundancy = list()
	for i in redf:
		redundancy.append(i.strip().split()[1])

with open("TSS_basic_gencode.v24.alltranscripts_sorted.bed", 'r') as tssf:
	tss = dict()
	for i in tssf:
		line = i.strip('\n').split('\t')
		id = line[4].split('/')[1]
		if id in redundancy:
			tss[id] = int(line[2])

with open("PLS_annotated_table_full.txt", 'r') as dfh, open("PLS_annotated_table_full_non-redundant.txt", 'w') as rfh:
	
	header = dfh.readline()
	rfh.write(header)

	redundancy_dict = dict()

	for i in dfh:
		line = i.strip('\n').split('\t')
		if line[5] in redundancy:
			try:
				first_start = int(redundancy_dict[line[5]][1].split(':')[1].split('-')[0])
				first_end = int(redundancy_dict[line[5]][1].split(':')[1].split('-')[1])
				first_center = round( (first_start + first_end) / 2 )

				second_start = int(line[1].split(':')[1].split('-')[0])
				second_end = int(line[1].split(':')[1].split('-')[1])
				second_center = round( (second_start + second_end) / 2 )

				if abs(tss[line[5]] - first_center) == abs(tss[line[5]] - second_center):
					# in the case of ENST00000590949, promoter w/ high CpG Density is chosen
					rfh.write('\t'.join(line) + '\n')

				elif abs(tss[line[5]] - first_center) > abs(tss[line[5]] - second_center):
					rfh.write('\t'.join(line) + '\n')

				elif abs(tss[line[5]] - first_center) < abs(tss[line[5]] - second_center):
					rfh.write('\t'.join(redundancy_dict[line[5]]) + '\n')				
			except KeyError:
				redundancy_dict[line[5]] = line
		
		
		else:
			rfh.write('\t'.join(line) + '\n')
		
		rfh.flush()
