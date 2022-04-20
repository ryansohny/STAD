#!/bin/bash

bedtools slop -i HOX_Cluster.bed \
-g hg38.analysisSet_chrL.fa.fai \
-b 1000 \
> HOX_Cluster_1000bp_slop.bed


bedtools makewindows \
-b HOX_Cluster_1000bp_slop.bed \
-w 1000 \
> HOX_Cluster_1000bp_slop_1kb_windows.bed
