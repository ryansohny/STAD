export PATH=/mnt/mone/Project/WC300/Tools/Anaconda3/bin:$PATH
source activate igvtools

for((i=$2;i<=$3;i++))
do
        sed -n ${i}p $1 > tmp${i}
        sample=$(awk '{print $1}' tmp${i})
        rm tmp${i}
mkdir -p Coverage_file_cov5

echo 'bedGraph sorting'
sort -k1,1 -k2,2n \
/mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/${sample}.bg \
> /mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/${sample}.sorted.bg

echo 'BigWig Generation'
bedGraphToBigWig \
/mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/${sample}.sorted.bg \
/mnt/mone/Project/WC300/03.WGBS_New/Reference/hg38_primary_chrX.chrom.sizes \
/mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/Coverage_file_cov5/${sample}.bw

echo 'Wig file generation'
bigWigToWig \
/mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/Coverage_file_cov5/${sample}.bw \
/mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/Coverage_file_cov5/${sample}.wig

echo 'TDF file generation'
igvtools toTDF \
/mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/Coverage_file_cov5/${sample}.wig \
/mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/Coverage_file_cov5/${sample}.tdf \
/mnt/mone/Project/WC300/03.WGBS_New/Reference/hg38_primary_chrX.chrom.sizes

echo 'Removing BigWig and Wig'
/bin/rm \
/mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/${sample}.sorted.bg \
/mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/Coverage_file_cov5/${sample}.bw \
/mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/Coverage_file_cov5/${sample}.wig

done
