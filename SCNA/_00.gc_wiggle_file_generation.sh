source activate sequenza
reference=genome.fa # PCAWG reference (hg19)
sequenza-utils \
gc_wiggle \
-w 50 \
--fasta ${refernece} \
-o genome.gc50Base.wig.gz
