#!/bin/bash

for((i=$2;i<=$3;i++))
do
        sed -n ${i}p $1 > tmp${i}
        sample=$(awk '{print $1}' tmp${i})
        rm tmp${i}
tabix -s 1 -b 2 -e 2 ${sample}.pileup.gz
done
