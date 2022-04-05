library(dada2)

fnFs <- unique(snakemake@input$for_reads)
fnRs <- unique(snakemake@input$rev_reads)

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filtFs <- unique(snakemake@output$for_reads)
filtRs <- unique(snakemake@output$rev_reads)

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
		     truncLen=0,
		     maxN=0, maxEE=c(2, 2),
		     truncQ=2, rm.phix=TRUE, compress=TRUE,
		     multithread=TRUE)

saveRDS(out, file = snakemake@output$out_rds)
