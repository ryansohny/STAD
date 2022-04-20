import os
from glob import glob
dbf = open("/mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/CpG_islands/sample_list.txt", 'r')
samples = sorted(list(map(lambda x: x.strip(), dbf.readlines())))

dbf.close()
for sample in samples:
        os.system('/usr/bin/cut -f2 ./HOX_Clusters_individual/HOX_Clusters_' + sample + '_met.txt > ./HOX_Clusters_individual/HOX_Clusters_' + sample + '_met_cut.txt')

# HOX_Clusters_individual/HOX_Clusters_14002_N_met.txt

normal_head_files = ['./HOX_Clusters_individual/HOX_Clusters_14002_N_met.txt']
normal_files_rest = sorted(glob('./HOX_Clusters_individual/HOX_Clusters_*_N*_met_cut.txt'))[1:]

all = ' '.join(normal_head_files + normal_files_rest + sorted(glob('./HOX_Clusters_individual/HOX_Clusters_*_T*_met_cut.txt')))

os.system('/usr/bin/paste ' + all + ' > HOX_Clusters_ALL_pre.txt')
os.system('/usr/bin/cat header_ALL.txt HOX_Clusters_ALL_pre.txt > HOX_Clusters_ALL.txt')
os.system('/usr/bin/rm HOX_Clusters_ALL_pre.txt')
