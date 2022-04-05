library(dada2)

fnFs <- unique(snakemake@input$for_read)

print(fnFs)

# sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
# 
# pdf(
# plotQualityProfile(fnFs[1:2])
# 
# 
# # filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
# # filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
# 
# filtFs <- unique(snakemake@output$for_reads)
# filtRs <- unique(snakemake@output$rev_reads)
# names(filtFs) <- sample.names
# names(filtRs) <- sample.names
# 
# out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
#               maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
#               compress=TRUE, multithread=TRUE)
# 
# saveRDS(out, file = snakemake@output$out_rds)
