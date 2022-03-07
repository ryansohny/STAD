# "hg38.analysisSet.fa" downloaded from http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/analysisSet/
# "chrL.fa" downloaded from https://www.ncbi.nlm.nih.gov/nuccore/J02459.1?report=fasta
# cat hg38.analysisSet.fa chrL.fa > hg38.analysisSet_chrL.fa

#!/bin/bash

export PATH=/mnt/mone/Project/WC300/Tools/Anaconda3/bin:$PATH
source activate bismark

bismark_genome_preparation \
--path_to_aligner /mnt/mone/Project/WC300/Tools/Anaconda3/envs/bismark/bin/ \
--parallel 55 \
--verbose \
/mnt/mone/Project/WC300/03.WGBS_New/Reference

# This creates CT_converted (CT_conversion) and GA_coverted genome (GA_conversion) under newly created Bisulfite_Genome directory
