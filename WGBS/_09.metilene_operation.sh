#!/bin/bash

metilene=/mnt/mone/Project/WC300/Tools/metilene_v0.2-8/metilene_linux64
metilene_output=/mnt/mone/Project/WC300/Tools/metilene_v0.2-8/metilene_output.pl

$metilene \
--groupA t \
-groupB n \
--maxdist 500 \
-mincpgs 5 \
--minMethDiff 0.1 \
--threads 56 \
--mode 1 \
--minNoA 55 \
--minNoB 55 \
--valley 0.700000 \
DNA_methylation_NT84_new_met.txt | sort -V -k1,1 -k2,2n > DNA_methylation_NT84_min55.sorted.output

$metilene_output \
-q ./DNA_methylation_NT84_min55.sorted.output \
-o ./DNA_methylation_NT84_min55.sorted \
-p 0.05 \
-a t \
-b n
