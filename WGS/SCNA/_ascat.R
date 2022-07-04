library(ASCAT)

args = commandArgs(trailingOnly=TRUE)
tumorid <- list(args[1])
normalid <- list(args[2])
gender <- list(args[3])
bamDir <- "/mnt/mone/Project/WC300/06.WGS/01.Aligned"

# Prepare BAF, LogR file for ASCAT (High-Throughput Sequencing (HTS) mode)
tumorBam <- paste(paste(bamDir, tumorid, sep="/"), ".recal.bam", sep='')
normalBam <- paste(paste(bamDir, normalid, sep="/"), ".recal.bam", sep='')
ascat.prepareHTS(tumorBam,
                 normalBam,
                 tumorid,
                 normalid,
                 "/mnt/mone/Project/WC300/Tools/Anaconda3/envs/ASCAT/bin/alleleCounter",
                 "/mnt/mone/Project/WC300/06.WGS/06.ASCAT/Reference_Files/alleles/G1000_alleles_hg38_chr",
                 "/mnt/mone/Project/WC300/06.WGS/06.ASCAT/Reference_Files/loci/G1000_loci_hg38_chr",
                 gender,
                 "hg38",
                 nthreads=10,
                 minCounts=10,
                 chrom_names=c(1:22, 'X'),
                 min_base_qual=20,
                 min_map_qual=35)

# Load BAF, LogR Data
workdir <- paste("/mnt/mone/Project/WC300/06.WGS/06.ASCAT/", tumorid, '/', sep='')
tumorLogR = paste(workdir, tumorid, "_tumourLogR.txt", sep='')
tumorBAF = paste(workdir, tumorid, "_tumourBAF.txt", sep='')
normalLogR = paste(workdir, tumorid, "_normalLogR.txt", sep='')
normalBAF = paste(workdir, tumorid, "_normalBAF.txt", sep='')
ascat.bc = ascat.loadData(Tumor_LogR_file = tumorLogR,
                          Tumor_BAF_file = tumorBAF,
                          Germline_LogR_file = normalLogR,
                          Germline_BAF_file = normalBAF,
                          gender = gender,
                          chrs = c(1:22, "X"),
                          genomeVersion = "hg38")

# Plotting Data before correction (GC content, Replication timing)
ascat.plotRawData(ascat.bc, img.prefix = "BC_")

# Correction (GC content, Replication timing)
ascat.bc = ascat.correctLogR(ascat.bc, 
                             GCcontentfile = "/mnt/mone/Project/WC300/06.WGS/06.ASCAT/Reference_Files/GC_G1000_hg38.txt",
                             replictimingfile = "/mnt/mone/Project/WC300/06.WGS/06.ASCAT/Reference_Files/RT_G1000_hg38.txt")

# Plotting Data After correction
ascat.plotRawData(ascat.bc, img.prefix = "AC_")

# Segmentation using allele-specific piecewise constant fitting (ASPCF)
ascat.bc = ascat.aspcf(ascat.bc,
                       penalty = 70,
                       out.prefix = paste(tumorid, "_", sep=''))

# Plotting segmentation Data
ascat.plotSegmentedData(ascat.bc)

# ASCAT RUN (Gamma should be 1 for HTS)
ascat.output = ascat.runAscat(ascat.bc, 
                              gamma=1,
                              pdfPlot=TRUE,
                              write_segments=TRUE)

# ASCAT Quality check
QC = ascat.metrics(ascat.bc,ascat.output)
write.csv(QC, 
          paste(workdir, "Results_table_", tumorid, ".csv", sep = ''),
          row.names=TRUE)
