
library(phyloseq)
library(tidyverse)
library(Biostrings)

bact_ps <-
  readRDS(snakemake@input$bact_phyloseq) %>%
  # readRDS("../../data/results/16S_phyloseq.RDS") %>%
  subset_taxa(family  != "Mitochondria" &
              class   != "Chloroplast") %>%
  filter_taxa(function(x) sum(x) > 3, TRUE)


euk_ps <-
  readRDS(snakemake@input$euk_phyloseq) %>%
  # readRDS("../../data/results/18S_phyloseq.RDS") %>%
  filter_taxa(function(x) sum(x) > 3, TRUE)


refseq(bact_ps) %>%
  as.character %>%
  tibble(Id = names(.), Seq = .) %>%
  mutate(Out = paste0(">", Id, "\n", Seq)) %>%
  .$Out %>%
  write(snakemake@output[[1]])


refseq(euk_ps) %>%
  as.character %>%
  tibble(Id = names(.), Seq = .) %>%
  mutate(Out = paste0(">", Id, "\n", Seq)) %>%
  .$Out %>%
  write(snakemake@output[[2]])

