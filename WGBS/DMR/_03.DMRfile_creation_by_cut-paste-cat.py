import os
from glob import glob
dbf = open("/mnt/mone/Project/WC300/03.WGBS/02.DNA_methylation/Re-Analysis/metilene/DMR/sample_20211221.txt", 'r')
samples = sorted(list(map(lambda x: x.strip(), dbf.readlines())))
dbf.close()
for sample in samples:
        os.system('/usr/bin/cut -f2 ./DMR_individual/DMR_' + sample + '_met.txt > ./DMR_individual/DMR_' + sample + '_met_cut.txt')

normal_head_files = ['./DMR_individual/DMR_14002_N_met.txt']
normal_files_rest = sorted(glob('./DMR_individual/DMR_*_N_met_cut.txt'))[1:]
normal = ' '.join(normal_head_files + normal_files_rest)
os.system('/usr/bin/paste ' + normal + ' > DMR_abs15_normal_pre.txt')
os.system('/usr/bin/cat header_Normal.txt DMR_abs15_normal_pre.txt > DMR_abs15_Normal.txt')
os.system('/usr/bin/rm DMR_abs15_normal_pre.txt')

tumor_head_files = ['./DMR_individual/DMR_14002_T_met.txt']
tumor_files_rest = sorted(glob('./DMR_individual/DMR_*_T_met_cut.txt'))[1:]
tumor = ' '.join(tumor_head_files + tumor_files_rest)
os.system('/usr/bin/paste ' + tumor + ' > DMR_abs15_tumor_pre.txt')
os.system('/usr/bin/cat header_Tumor.txt DMR_abs15_tumor_pre.txt > DMR_abs15_Tumor.txt')
os.system('/usr/bin/rm DMR_abs15_tumor_pre.txt')

all = ' '.join(normal_head_files + normal_files_rest + sorted(glob('./DMR_individual/DMR_*_T_met_cut.txt')))

os.system('/usr/bin/paste ' + all + ' > DMR_abs15_ALL_pre.txt')
os.system('/usr/bin/cat header_ALL.txt DMR_abs15_ALL_pre.txt > DMR_abs15_ALL.txt')
os.system('/usr/bin/rm DMR_abs15_ALL_pre.txt')
