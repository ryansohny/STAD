import sys
try:
	cutoff = sys.argv[1]
except ValueError:
	print('Try typing in correct cutoff (Tumor-Normal)\n')
	sys.exit()
dfh = open("DNA_methylation_NT84_min55.sorted_qval.0.05.out", 'r')
rfh = open("DNA_methylation_NT84_min55.sorted_qval.0.05.abs" + cutoff + '.bed', 'w')

for i in dfh:
	line = i.strip('\n').split('\t')
	if abs(float(line[4])) >= int(cutoff):
		rfh.write('\t'.join(line[:3]) + '\t' + line[4] + '\t' + line[5] + '\n')
		rfh.flush()
