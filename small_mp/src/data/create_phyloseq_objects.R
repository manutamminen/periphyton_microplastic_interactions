
library(phyloseq)
library(tidyverse)
library(Biostrings)
library(readxl)
library(readxl)


bact_otu_tbl <-
  readRDS(snakemake@input$bact_seqtab_nochim)
  # readRDS("../../data/results/16S_seqtab_nochim.RDS")


rownames(bact_otu_tbl) <-
  rownames(bact_otu_tbl) %>%
  str_split("_") %>%
  map_chr(1) %>%
  str_replace_all("-", "_")


euk_otu_tbl <-
  readRDS(snakemake@input$euk_seqtab_nochim)
  # readRDS("../../data/results/18S_seqtab_nochim.RDS")


rownames(euk_otu_tbl) <-
  rownames(euk_otu_tbl) %>%
  str_split("_") %>%
  map_chr(1) %>%
  str_replace_all("-", "_")


bact_taxa_tbl <-
  readRDS(snakemake@input$bact_taxa)
  # readRDS("../../data/results/16S_taxa.RDS")


euk_taxa_tbl <-
  readRDS(snakemake@input$euk_taxa)
  # readRDS("../../data/results/18S_taxa.RDS")


bact_smp_tbl <-
  read_csv(snakemake@input$samples) %>%
  # read_csv("../../big_mp_samples.csv") %>%
    filter(str_detect(Sample, "16s")) %>%
    mutate(Sample = str_replace_all(Sample, "-", "_")) %>%
    data.frame(row.names = "Sample") %>%
    sample_data


euk_smp_tbl <-
  read_csv(snakemake@input$samples) %>%
  # read_csv("../../big_mp_samples.csv") %>%
    filter(str_detect(Sample, "18s")) %>%
    mutate(Sample = str_replace_all(Sample, "-", "_")) %>%
    data.frame(row.names = "Sample") %>%
    sample_data

bact_ps <-
  phyloseq(otu_table(bact_otu_tbl, taxa_are_rows=FALSE),
           bact_smp_tbl,
           tax_table(bact_taxa_tbl)) %>%
  subset_taxa(domain == "Bacteria")


euk_ps <-
  phyloseq(otu_table(euk_otu_tbl, taxa_are_rows=FALSE),
           euk_smp_tbl,
           tax_table(euk_taxa_tbl)) %>%
  subset_taxa(domain == "Eukaryota")


bact_dna <- Biostrings::DNAStringSet(taxa_names(bact_ps))
names(bact_dna) <- taxa_names(bact_ps)
bact_ps <- merge_phyloseq(bact_ps, bact_dna)
taxa_names(bact_ps) <- paste0("ASV", seq(ntaxa(bact_ps)))


euk_dna <- Biostrings::DNAStringSet(taxa_names(euk_ps))
names(euk_dna) <- taxa_names(euk_ps)
euk_ps <- merge_phyloseq(euk_ps, euk_dna)
taxa_names(euk_ps) <- paste0("ASV", seq(ntaxa(euk_ps)))


saveRDS(bact_ps, snakemake@output$bact_phyloseq)
saveRDS(euk_ps, snakemake@output$euk_phyloseq)
