library(dada2)
library(tidyverse)
library(DECIPHER)

out <- readRDS(snakemake@input$out_rds)
dadaFs <- readRDS(snakemake@input$for_samples)
dadaRs <- readRDS(snakemake@input$rev_samples)
mergers <- readRDS(snakemake@input$merged_samples)

fnFs <- unique(snakemake@input$for_reads)
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

seqtab <- makeSequenceTable(mergers)

len_distr <- data.frame(table(nchar(getSequences(seqtab))))

seqtab.nochim <- removeBimeraDenovo(seqtab,
				    method="consensus",
				    multithread=TRUE,
				    verbose=TRUE)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN),
	       sapply(dadaRs, getN),
	       sapply(mergers, getN),
	       rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered",
		     "denoisedF", "denoisedR",
		     "merged", "nonchim")
rownames(track) <- sample.names

# taxa <- assignTaxonomy(seqtab.nochim,
#            snakemake@params[['tax_db']],
#            multithread=TRUE)
# 
write_csv(len_distr, snakemake@output$len_distr)
write_csv(as.data.frame(track), snakemake@output$sample_stats)
# saveRDS(seqtab.nochim, snakemake@output$seqtab_nochim)
# saveRDS(taxa, snakemake@output$taxa)
# 


dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
# load("~/tax/IDTaxa/SILVA_SSU_r132_March2018.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
load(snakemake@params$tax_db) # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)


saveRDS(seqtab.nochim, snakemake@output$seqtab_nochim)
saveRDS(taxid, snakemake@output$taxa)
