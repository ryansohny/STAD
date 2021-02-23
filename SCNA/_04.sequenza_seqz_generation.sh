source activate sequenza

sequenza=/clinix1/Analysis/mongol/phenomata/Tools/anaconda3/envs/sequenza/bin/sequenza-utils

## input file format ==> Normalid on the first column tab separated wih Tumorid on the second column

for((i=$2;i<=$3;i++))
do
        sed -n ${i}p $1 > tmp${i}
        normal=$(awk '{print $1}' tmp${i})
        tumor=$(awk '{print $2}' tmp${i})
        rm tmp${i}
mkdir -p ${tumor}
echo 'Generation of seqz file : ' ${tumor}
${sequenza} \
bam2seqz \
--pileup \
-n ${normal}.pileup.gz \
-t ${tumor}.pileup.gz \
--fasta genome.fa \
--parallel 25 \
-C 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y \
--samtools samtools_executable \
--tabix tabix_executable \
-gc genome_primary.gc50Base.wig.gz \
-o ${tumor}/${tumor}.seqz.gz
done
