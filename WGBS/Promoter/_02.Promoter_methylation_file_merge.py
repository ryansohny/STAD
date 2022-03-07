import os
from glob import glob
dbf = open("/mnt/mone/Project/WC300/03.WGBS/02.DNA_methylation/Re-Analysis/metilene/DMR/sample_20211221.txt", 'r')
samples = sorted(list(map(lambda x: x.strip(), dbf.readlines())))
dbf.close()
for sample in samples:
	os.system('/usr/bin/cut -f2 ./Promoter_individual/Promoter_up500down500_' + sample + '_met.txt > ./Promoter_individual/Promoter_up500down500_' + sample + '_met_cut.txt')

normal_head_files = ['./Promoter_individual/Promoter_up500down500_14002_N_met.txt']
normal_files_rest = sorted(glob('./Promoter_individual/Promoter_up500down500_*_N_met_cut.txt'))[1:]
normal = ' '.join(normal_head_files + normal_files_rest)
os.system('/usr/bin/paste ' + normal + ' > Promoter_up500down500_normal_pre.txt')
os.system('/usr/bin/cat header_Normal.txt Promoter_up500down500_normal_pre.txt > Promoter_up500down500_Normal.txt')
os.system('/usr/bin/rm Promoter_up500down500_normal_pre.txt')

tumor_head_files = ['./Promoter_individual/Promoter_up500down500_14002_T_met.txt']
tumor_files_rest = sorted(glob('./Promoter_individual/Promoter_up500down500_*_T_met_cut.txt'))[1:]
tumor = ' '.join(tumor_head_files + tumor_files_rest)
os.system('/usr/bin/paste ' + tumor + ' > Promoter_up500down500_tumor_pre.txt')
os.system('/usr/bin/cat header_Tumor.txt Promoter_up500down500_tumor_pre.txt > Promoter_up500down500_Tumor.txt')
os.system('/usr/bin/rm Promoter_up500down500_tumor_pre.txt')

all = ' '.join(normal_head_files + normal_files_rest + sorted(glob('./Promoter_individual/Promoter_up500down500_*_T_met_cut.txt')))

os.system('/usr/bin/paste ' + all + ' > Promoter_up500down500_ALL_pre.txt')
os.system('/usr/bin/cat header_ALL.txt Promoter_up500down500_ALL_pre.txt > Promoter_up500down500_ALL.txt')
os.system('/usr/bin/rm Promoter_up500down500_ALL_pre.txt')
