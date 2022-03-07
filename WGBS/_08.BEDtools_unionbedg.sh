bedtools=/usr/bin/bedtools

normal_files=$(ls *_N.bg *_TD4.bg)
tumor_files=$(ls *_T.bg *_TD1.bg)
normal_id=$(echo $normal_files | tr -d '.bg')
tumor_id=$(echo $tumor_files | tr -d '.bg')

$bedtools unionbedg -i $normal_files $tumor_files -filler . -header -names $normal_id $tumor_id > DNA_methylation_NT84.bg
