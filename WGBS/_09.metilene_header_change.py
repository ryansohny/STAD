import sys
dfh = open(sys.argv[1], 'r') # DNA_methylation_NT84.bg
rfh = open(sys.argv[1][:-3] + '_met.txt', 'w')

line = dfh.readline().strip('\n').split('\t')
rfh.write('chr\tpos\t' + '\t'.join(list(map(lambda x: 'n_' + x, line[3:87]))) + '\t' + '\t'.join(list(map(lambda x: 't_' + x, line[87:]))) + '\n')

for i in dfh:
        line = i.strip('\n').split('\t')
        rfh.write(line[0] + '\t' + line[2] + '\t' + '\t'.join(line[3:]) + '\n')
        rfh.flush()
