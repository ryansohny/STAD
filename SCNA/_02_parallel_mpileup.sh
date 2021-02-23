source activate parallel

## input file format ==> sampleid on the first column
## usage: /usr/bin/sh _parallel_mpileup.sh input_file start end

referenceIndex=genome.fa.fai

for((i=$2;i<=$3;i++))
do
        sed -n ${i}p $1 > tmp${i}
        sample=$(awk '{print $1}' tmp${i})
        rm tmp${i}
parallel \
--gnu \
--compress \
--tmpdir tmp \
--colsep '\t' \
--keep-order \
--jobs 40 \
samtools mpileup -b ${sample}.fofn \
-q 20 \
-Q 20 \
-r {1} \
:::: \
${referenceIndex} \
> ${sample}.pileup
done
