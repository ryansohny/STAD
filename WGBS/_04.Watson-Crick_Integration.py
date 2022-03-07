dbf = open("sample_list.txt", 'r')
a = [i.strip() for i in dbf]
import sys
start = sys.argv[1]
end = sys.argv[2]
d = a[int(start)-1:int(end)]

for i in d:
        dfh = open('../' + i + '_val_1_bismark_bt2_pe.deduplicated.CpG_report.2.txt', 'r')
        rfh = open(i + '_WCmerged.txt', 'w')
        line1 = dfh.readline().strip().split()
        line2 = dfh.readline().strip().split()
        while line1 != []:
                rfh.write(line1[0] + '\t' + line1[1] + '\t' + line1[2] + '\t' + str(int(line1[3]) + int(line2[3])) + '\t' + str(int(line1[4]) + int(line2[4])) + '\t' + line1[5] + '\t' + line1[6] + '\n')
                rfh.flush()
                line1 = dfh.readline().strip().split()
                line2 = dfh.readline().strip().split()
        dfh.close()
        rfh.close()
