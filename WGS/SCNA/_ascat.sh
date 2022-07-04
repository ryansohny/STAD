export PATH=/mnt/mone/Project/WC300/Tools/Anaconda3/bin:$PATH
source activate ASCAT

for((i=$2;i<=$3;i++))
do
        sed -n ${i}p $1 > tmp${i}
        tumor=$(awk '{print $1}' tmp${i})
        normal=$(awk '{print $2}' tmp${i})
        gender=$(awk '{print $3}' tmp${i})
        rm tmp${i}

mkdir -p ${tumor}

/usr/bin/cp _ascat.R ${tumor}

cd ${tumor}

Rscript _ascat.R ${tumor} ${normal} ${gender}

cd ..

done