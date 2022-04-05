library(dada2)

filtFs <- unique(snakemake@input$for_reads)
filtRs <- unique(snakemake@input$rev_reads)

errF <- readRDS(snakemake@input$for_errors)
errR <- readRDS(snakemake@input$rev_errors)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool=TRUE)

saveRDS(dadaFs, snakemake@output$for_samples)
saveRDS(dadaRs, snakemake@output$rev_samples)

