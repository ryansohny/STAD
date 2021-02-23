import sys
## input file foramt ==> sampleid on the first column 

dbf = open(sys.argv[1], 'r')
for i in dbf:
        rfh = open(i.strip() + '.fofn', 'w')
        rfh.write(i.strip() + ".sorted.dp.ir.recal.bam")
        rfh.close()
