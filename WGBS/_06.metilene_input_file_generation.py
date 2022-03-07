dbf = open("sample_list.txt", 'r')
a = [i.strip() for i in dbf]
import sys
start = sys.argv[1]
end = sys.argv[2]
d = a[int(start)-1:int(end)]

for i in d:
        dfh = open(i + '_WCmerged.txt', 'r')
        rfh = open(i + '.bg', 'w')
        for j in dfh:
                line = j.strip('\n').split('\t')
                if int(line[3]) + int(line[4]) >= 5:
                        chrom = line[1].split(':')[0]
                        start = str(int(line[1].split(':')[1]) -1)
                        end = line[1].split(':')[1]
                        met = round(int(line[3]) / (int(line[3]) + int(line[4])) * 100, 4)
                        # round(a / (a+b) * 100, 4)
                        rfh.write(chrom + '\t' + start + '\t' + end + '\t' + str(met) + '\n')
                        rfh.flush()
        dfh.close()
        rfh.close()
