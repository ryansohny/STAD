import gzip as gz
from glob import glob
import os
index = dict()
inf = open("index.txt", 'r')
for i in inf:
	i = i.strip().split()
	index[i[0]] = int(i[1])
inf.close()

dbf = open("sample_list.txt", 'r')
a = [i.strip() for i in dbf]
import sys
start = sys.argv[1]
end = sys.argv[2]
d = a[int(start)-1:int(end)]


for i in d:
	dfh = gz.open('../' + i + '_1_bismark_bt2_pe.deduplicated.CpG_report.txt.gz', 'r')
	rfh = open("test" + start + '_' + end + '.txt', 'w')
	for j in dfh:
		line = j.decode('utf-8').strip('\n').split()
		try:
			rfh.write(str((index[line[0]] * 1000000000) + int(line[1])) + '\t' + line[0] + ':' + line[1] + '\t' + '\t'.join(line[2:]) + '\n')
			rfh.flush()
		except KeyError:
			pass
	dfh.close()
	rfh.close()
	os.system("sort -n -k 1 test" + start + '_' + end + '.txt > ../' + i + '_bismark_bt2_pe.deduplicated.CpG_report.2.txt')
  os.system("rm test" + start + '_' + end + '.txt')
