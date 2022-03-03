#!/bin/bash

export PATH=/mnt/mone/Project/WC300/Tools/Anaconda3/bin:$PATH
source activate bismark

ref=/mnt/mone/Project/WC300/03.WGBS_New/Reference
fasta=/mnt/mone/PMI/WC300_WGBS/Sequence_Reads
outdir1=/mnt/mone/Project/WC300/03.WGBS_New/01.Alignment
outdir2=/mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation

mkdir -p logs
for((i=$2;i<=$3;i++))
do

        sed -n ${i}p $1 > tmp${i}
        sample=$(awk '{print $1}' tmp${i})
        rm tmp${i}
/usr/bin/time -f "%E" -o ./logs/trim_${sample}.out \
trim_galore \
-o /mnt/mone/Project/WC300/03.WGBS_New/00.Sequence_Reads \
--path_to_cutadapt /mnt/mone/Project/WC300/Tools/Anaconda3/envs/bismark/bin/cutadapt \
--clip_R1 10 \
--clip_R2 15 \
--three_prime_clip_R1 12 \
--three_prime_clip_R2 13 \
--illumina \
--paired \
--basename ${sample} \
--cores 8 \
${fasta}/${sample}_1.fastq.gz \
${fasta}/${sample}_2.fastq.gz

done
