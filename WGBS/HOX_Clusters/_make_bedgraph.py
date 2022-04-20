import os
import sys
from glob import glob
dbf = open("/mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/CpG_islands/sample_list.txt", 'r')
samples = sorted(list(map(lambda x: x.strip(), dbf.readlines())))
dbf.close()

#chr2:176091890-176092890/71    18.2857
#chr2:176092890-176093890/94    12.5448
#chr2:176093890-176094890/22    15.0
#chr2:176094890-176095890/8     23.8462
#chr2:176095890-176096890/13    22.3684

for sample in samples:
        dfh = open("./HOX_Clusters_individual/HOX_Clusters_" + sample + "_met.txt", 'r')
        rfh = open("./HOX_Clusters_individual/HOX_Clusters_" + sample + ".bedgraph", 'w')
        rfh.write('track type=bedGraph name=' + sample + ' description=HOX_Clusters\n')
        for i in dfh:
                line = i.strip('\n').split('\t')
                rfh.write(line[0].split(':')[0] + '\t' + line[0].split(':')[1].split('-')[0] + '\t' + line[0].split(':')[1].split('-')[1].split('/')[0] + '\t' + line[1] + '\n')
                rfh.flush()
        dfh.close()
        rfh.close()


os.system('/usr/bin/tar -zvcf HOX_Clusters_bedGrpah.tar.gz HOX_Clusters_individual/*.bedgraph')
