#!/bin/bash
#export PATH=/mnt/mone/Project/WC300/Tools/Anaconda3/bin:$PATH
#source activate deeptools

input_dir=/mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/Coverage_file_cov5
cpgi=${input_dir}/deepTools/CpGi
se=${input_dir}/deepTools/Super_Enhancers
genebody=${input_dir}/deepTools/Gene_Body

blacklist=/mnt/mone/Project/AK1_PacBio/01.DNA/Public_Data/GRCh38_unified_blacklist.bed

# Gene Body
## Normal
computeMatrix scale-regions \
--regionsFileName ${genebody}/GeneBody_gencode.v24.allgenes.sorted.strand.bed \
--scoreFileName ${input_dir}/*_N.bw \
--outFileName ${genebody}/Normal_GeneBody_allgenes.mat.gz \
--outFileNameMatrix ${genebody}/Normal_GeneBody_allgenes.tab \
--outFileSortedRegions ${genebody}/Normal_GeneBody_sortedRegions.bed \
--startLabel "TSS" \
--endLabel "TES" \
--upstream 1000 \
--downstream 1000 \
--regionBodyLength 1000 \
--binSize 10 \
--skipZeros \
--blackListFileName ${blacklist} \
--verbose \
--numberOfProcessors 50
## Tumor
computeMatrix scale-regions \
--regionsFileName ${genebody}/GeneBody_gencode.v24.allgenes.sorted.strand.bed \
--scoreFileName ${input_dir}/*_T.bw \
--outFileName ${genebody}/Tumor_GeneBody_allgenes.mat.gz \
--outFileNameMatrix ${genebody}/Tumor_GeneBody_allgenes.tab \
--outFileSortedRegions ${genebody}/Tumor_GeneBody_sortedRegions.bed \
--startLabel "TSS" \
--endLabel "TES" \
--upstream 1000 \
--downstream 1000 \
--regionBodyLength 1000 \
--binSize 10 \
--skipZeros \
--blackListFileName ${blacklist} \
--verbose \
--numberOfProcessors 50
