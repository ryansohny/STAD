import os
from glob import glob
dbf = open("/mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/DMR/DMR_min55_new/sample_list.txt", 'r')
samples = sorted(list(map(lambda x: x.strip(), dbf.readlines())))
dbf.close()
for sample in samples:
        os.system('/usr/bin/cut -f2 ./SE_DEG_individual/SE_' + sample + '_met.txt > ./SE_DEG_individual/SE_' + sample + '_met_cut.txt')

normal_head_files = ['./SE_DEG_individual/SE_14002_N_met.txt']
normal_files_rest = sorted(glob('./SE_DEG_individual/SE_*_N*_met_cut.txt'))[1:]

tumor_files_rest = sorted(glob('./SE_DEG_individual/SE_*_T*_met_cut.txt'))[1:]

all = ' '.join(normal_head_files + normal_files_rest + sorted(glob('./SE_DEG_individual/SE_*_T*_met_cut.txt')))

os.system('/usr/bin/paste ' + all + ' > SE_DEG_ALL_pre.txt')
os.system('/usr/bin/cat header_ALL.txt SE_DEG_ALL_pre.txt > SE_DEG_ALL.txt')
os.system('/usr/bin/rm SE_DEG_ALL_pre.txt')