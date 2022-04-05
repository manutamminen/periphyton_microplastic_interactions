library(dada2)

filtFs <- unique(snakemake@input$for_reads)
filtRs <- unique(snakemake@input$rev_reads)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

saveRDS(errF, snakemake@output$for_errors)
saveRDS(errR, snakemake@output$rev_errors)

pdf(snakemake@output$for_error_plots)
plotErrors(errF, nominalQ=TRUE)
dev.off()

pdf(snakemake@output$rev_error_plots)
plotErrors(errR, nominalQ=TRUE)
dev.off()
