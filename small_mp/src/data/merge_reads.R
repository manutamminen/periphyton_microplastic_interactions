library(dada2)

filtFs <- unique(snakemake@input$for_reads)
filtRs <- unique(snakemake@input$rev_reads)

dadaFs <- readRDS(snakemake@input$for_samples)
dadaRs <- readRDS(snakemake@input$rev_samples)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

saveRDS(mergers, snakemake@output$merged_samples)


