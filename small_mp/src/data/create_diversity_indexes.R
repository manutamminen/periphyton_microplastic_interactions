
library(phyloseq)
library(tidyverse)

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

bact_div <-
  bact_ps %>%
  estimate_richness(measures = c("Chao1", "Shannon")) %>%
  rownames_to_column(var = "Sample") %>%
  mutate(Treatment = case_when(str_detect(Sample, "SC") ~ "control",
                               str_detect(Sample, "SMP") ~ "MP",
                               str_detect(Sample, "SA") ~ "A"))

euk_div <-
  euk_ps %>%
  estimate_richness(measures = c("Chao1", "Shannon")) %>%
  rownames_to_column(var = "Sample") %>%
  mutate(Treatment = case_when(str_detect(Sample, "SC") ~ "control",
                               str_detect(Sample, "SMP") ~ "MP",
                               str_detect(Sample, "SA") ~ "A"))

bind_rows(
  tibble(Type = "Chao1, bacteria",
       C_vs_A = wilcox.test(filter(bact_div, Treatment == "control")$Chao1,
                            filter(bact_div, Treatment == "A")$Chao1)$p.value,
       C_vs_MP = wilcox.test(filter(bact_div, Treatment == "control")$Chao1,
                             filter(bact_div, Treatment == "MP")$Chao1)$p.value,
       A_vs_MP = wilcox.test(filter(bact_div, Treatment == "A")$Chao1,
                             filter(bact_div, Treatment == "MP")$Chao1)$p.value),
  tibble(Type = "Shannon, bacteria",
       C_vs_A = wilcox.test(filter(bact_div, Treatment == "control")$Shannon,
                            filter(bact_div, Treatment == "A")$Shannon)$p.value,
       C_vs_MP = wilcox.test(filter(bact_div, Treatment == "control")$Shannon,
                             filter(bact_div, Treatment == "MP")$Shannon)$p.value,
       A_vs_MP = wilcox.test(filter(bact_div, Treatment == "A")$Shannon,
                             filter(bact_div, Treatment == "MP")$Shannon)$p.value),
  tibble(Type = "Chao1, eukaryotes",
       C_vs_A = wilcox.test(filter(euk_div, Treatment == "control")$Chao1,
                            filter(euk_div, Treatment == "A")$Chao1)$p.value,
       C_vs_MP = wilcox.test(filter(euk_div, Treatment == "control")$Chao1,
                             filter(euk_div, Treatment == "MP")$Chao1)$p.value,
       A_vs_MP = wilcox.test(filter(euk_div, Treatment == "A")$Chao1,
                             filter(euk_div, Treatment == "MP")$Chao1)$p.value),
  tibble(Type = "Shannon, eukaryotes",
       C_vs_A = wilcox.test(filter(euk_div, Treatment == "control")$Shannon,
                            filter(euk_div, Treatment == "A")$Shannon)$p.value,
       C_vs_MP = wilcox.test(filter(euk_div, Treatment == "control")$Shannon,
                             filter(euk_div, Treatment == "MP")$Shannon)$p.value,
       A_vs_MP = wilcox.test(filter(euk_div, Treatment == "A")$Shannon,
                             filter(euk_div, Treatment == "MP")$Shannon)$p.value)) %>%
  mutate_if(is.numeric, ~.*3) %>%
  mutate_if(is.numeric, ~ifelse(. > 1, 1, .)) %>%
  write_csv(snakemake@output$pval_table)



pdf(snakemake@output$bact_div, useDingbats=FALSE)
plot_richness(bact_ps, x="Treatment", measures=c("Chao1", "Shannon"))
dev.off()


pdf(snakemake@output$euk_div, useDingbats=FALSE)
plot_richness(euk_ps, x="Treatment", measures=c("Chao1", "Shannon"))
dev.off()

