suppressMessages(library(LOLA))
suppressMessages(library(GenomicRanges))
suppressMessages(library(qvalue))

# Working Directory : /mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/DSS/metilene/DMR/DMR_min55_new/LOLA

# Loading up the Database
db_path <- "/mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/metilene/DMR/DMR_min55_new/Motif_Analysis_wHOMER_wLOLA/LOLACoreCaches/t1/resources/regions/LOLACore/hg38"
input_path <- "/mnt/mone/Project/WC300/03.WGBS_New/02.DNA_methylation/DSS/metilene/DMR/DMR_min55_new/LOLA"
regionDB <- loadRegionDB( dbLocation=db_path, useCache=TRUE )

# Loading up the background regions (Total tested Regions)
Background <- readBed( file.path( input_path, "DMR_background.bed" ) )

# (i)
# Loading up the query (DMR w/ abs(methyl) > 10) 
HyperDMR_10 <- readBed( file.path( input_path, "DMR_abs10_hyper.bed" ) )
HypoDMR_10 <- readBed( file.path( input_path, "DMR_abs10_hypo.bed" ) )

# Combining the query regions into a single GRangesList (For the sake of simultaneous/faster calculation)
DMR_10_sets <- GRangesList( HyperDMR_10, HypoDMR_10 )

# Running Enrichment test 
locResults = runLOLA(DMR_10_sets, Background, regionDB, cores=20)
#> class(locResults)
#[1] "data.table" "data.frame"

# Slicing the results based on various criteria
#locResults[userSet==1][order(maxRnk, decreasing=FALSE),]

# Write the enrichment result to a designated folder
writeCombinedEnrichment(locResults, outFolder= "LOLA_Results_1", includeSplits=TRUE)


# (ii)
# Loading up the query (DMR w/ abs(methyl) > 10, For Hypo-DMR we exclude the regions overlapped with PMD) 
HyperDMR_10 <- readBed( file.path( input_path, "DMR_abs10_hyper.bed" ) )
HypoDMR_10_woPMD <- readBed( file.path( input_path, "DMR_abs10_hypo_woPMD.bed" ) )

# Combining the query regions into a single GRangesList (For the sake of simultaneous/faster calculation)
DMR_10_sets <- GRangesList( HyperDMR_10, HypoDMR_10_woPMD )

# Running Enrichment test 
locResults = runLOLA(DMR_10_sets, Background, regionDB, cores=20)
#> class(locResults)
#[1] "data.table" "data.frame"

# Slicing the results based on various criteria
#locResults[userSet==1][order(maxRnk, decreasing=FALSE),]

# Write the enrichment result to a designated folder
writeCombinedEnrichment(locResults, outFolder= "LOLA_Results_2", includeSplits=TRUE)
