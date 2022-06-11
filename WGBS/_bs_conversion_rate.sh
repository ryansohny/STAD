#!/bin/bash


echo Sample$'\t'Conversion_rate > WC300_2022_BS-conversion-rate.txt

for((i=$2;i<=$3;i++))
do
        sed -n ${i}p $1 > tmp${i}
        sample=$(awk '{print $1}' tmp${i})
        rm tmp${i}
#/clinix1/Analysis/mongol/phenomata/04.GC_CC/Whole_new_2020/03.WGBS/02.DNA_methylation/sample_list.txt
conversion=`zgrep 'chrL' ${sample}_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt.gz | awk '{sum+=$4} {SUM+=$5} END {print 1-(sum/SUM)}'`
echo ${sample}$'\t'${conversion} >> WC300_2022_BS-conversion-rate.txt
done
