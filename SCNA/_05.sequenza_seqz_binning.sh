source activate sequenza

sequenza=sequenza-utils

for((i=$2;i<=$3;i++))
do
        sed -n ${i}p $1 > tmp${i}
        normal=$(awk '{print $1}' tmp${i})
        tumor=$(awk '{print $2}' tmp${i})
        rm tmp${i}
mkdir -p ${tumor}
echo 'Post-process by binning the original seqz file : ' ${tumor}
zcat \
./${tumor}/${tumor}_1.seqz.gz \
./${tumor}/${tumor}_2.seqz.gz \
./${tumor}/${tumor}_3.seqz.gz \
./${tumor}/${tumor}_4.seqz.gz \
./${tumor}/${tumor}_5.seqz.gz \
./${tumor}/${tumor}_6.seqz.gz \
./${tumor}/${tumor}_7.seqz.gz \
./${tumor}/${tumor}_8.seqz.gz \
./${tumor}/${tumor}_9.seqz.gz \
./${tumor}/${tumor}_10.seqz.gz \
./${tumor}/${tumor}_11.seqz.gz \
./${tumor}/${tumor}_12.seqz.gz \
./${tumor}/${tumor}_13.seqz.gz \
./${tumor}/${tumor}_14.seqz.gz \
./${tumor}/${tumor}_15.seqz.gz \
./${tumor}/${tumor}_16.seqz.gz \
./${tumor}/${tumor}_17.seqz.gz \
./${tumor}/${tumor}_18.seqz.gz \
./${tumor}/${tumor}_19.seqz.gz \
./${tumor}/${tumor}_20.seqz.gz \
./${tumor}/${tumor}_21.seqz.gz \
./${tumor}/${tumor}_22.seqz.gz \
./${tumor}/${tumor}_X.seqz.gz \
./${tumor}/${tumor}_Y.seqz.gz | \
awk '{if (NR==1) {print $0} else {if ($1!="chromosome"){print $0}}}' | \
${sequenza} seqz_binning \
-s - \
--window 50 \
--tabix tabix_executable \
-o ./${tumor}/${tumor}_bin50.seqz.gz
done
