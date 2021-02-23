library(sequenza)

sample_list <- list(...) # list of sample id
dir <- "..." # directory where *_bin50.seqz.gz files exist
for(sample in sample_list) {
        datadir <- paste(dir, sample, sample, sep='/')
        data.file <- paste(datadir, "_bin50.seqz.gz", sep='')
        test <- sequenza.extract(data.file, verbose = TRUE)
        CP <- sequenza.fit(test)
        sequenza.results(sequenza.extract = test, cp.table = CP, sample.id = sample, out.dir=sample)
}
