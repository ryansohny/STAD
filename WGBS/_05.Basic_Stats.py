import os
a = os.listdir('.')
d = []
for i in a:
        if i[-12:] == 'WCmerged.txt':
                d.append(i)
d = sorted(d)
#for i in d:
#       print(i)

rfh = open("Basic_stats_GC_methyl_integrated.txt", 'w')
rfh.write("Sample\tTotalCpG\tCoveredLessThan5\tCoveredAbove5\tCoveredAbove10\tCoveredAbove15\n")
for i in d:
        dfh = open(i, 'r')
        #t', '20087_T_WCmerged.txt''20095_TD1_WCmerged.txt', 
        id = i[:-13]
#       id = i.split('_')[0] + '_' + i.split('_')[1]
        total = int()
        less5 = int()
        above5 = int()
        above10 = int()
        above15 = int()
        for j in dfh:
                total += 1
                line = j.strip().split()
                if int(line[3]) + int(line[4]) < 5:
                        less5 += 1
                elif int(line[3]) + int(line[4]) >= 5:
                        above5 += 1
                        if int(line[3]) + int(line[4]) >= 10:
                                above10 += 1
                                if int(line[3]) + int(line[4]) >= 15:
                                        above15 += 1
        rfh.write(id + '\t' + str(total) + '\t' + str(less5) + '\t' + str(above5) + '\t' + str(above10) + '\t' + str(above15) + '\n')
        rfh.flush()
        dfh.close()
