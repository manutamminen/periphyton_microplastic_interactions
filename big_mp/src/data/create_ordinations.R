
library(phyloseq)
library(tidyverse)
library(readxl)
library(vegan)
library(ggrepel)
library(broom)


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

calc_rda <- function(ps_obj, formula) {
    form <- as.formula(formula)

    rda_otut <-
        data.frame(otu_table(ps_obj)) %>%
        mutate(Otu=rownames(.)) %>%
        pivot_longer(cols = -"Otu",
                     names_to = "Sample_ID",
                     values_to = "Count") %>%
        filter(complete.cases(.))

    smp_tbl <-
        data.frame(sample_data(ps_obj)) %>%
        mutate(Sample_ID=rownames(.))

    rda_table <-
        left_join(rda_otut, smp_tbl, by=c("Otu" = "Sample_ID")) %>%
        select(Otu, Sample_ID, Count, Treatment, Ix) %>%
        pivot_wider(id_cols = c("Otu", "Treatment", "Ix"),
                    names_from = "Sample_ID",
                    values_from = "Count")



    spe <- select(rda_table, -(1:3))


    FULL.cap <- capscale(form, data=rda_table)

    basplot <- plot(FULL.cap)

    taxonomy <- data.frame(tax_table(ps_obj)) %>%
        mutate(Otu=rownames(.))

    species <- data.frame(basplot$species) %>%
        mutate(Otu=rownames(.),
               dist=sqrt(CAP1^2 + CAP2^2)) %>%
        arrange(desc(dist)) %>%
        left_join(taxonomy, by="Otu")

    sites <- data.frame(basplot$sites) %>%
        mutate(Treatment = rda_table$Treatment,
               Ix = rda_table$Ix,
               Sample = rda_table$Otu)

    list(sites=sites,
         species=species,
         cap=FULL.cap)
}


plot_rda <- function(rda_obj) {
  eigvals <-
    eigenvals(rda_obj$cap)/sum(eigenvals(rda_obj$cap))*100
  xtitle <- paste0("CAP1 (", round(eigvals[1], 2), " %)")
  ytitle <- paste0("CAP2 (", round(eigvals[2], 2), " %)")
  species <- head(rda_obj$species, 20)
  pval <-
    rda_obj$cap %>%
    anova %>%
    tidy %>%
    .$p.value %>%
    .[1]
  plot_title <- paste0("Permutational ANOVA p = ", pval)
  ggplot() +
    geom_point(data=rda_obj$sites,
               aes(x=CAP1, y=CAP2, color=Treatment), size=5) +
    # geom_text_repel(data = rda_obj$sites,
    #                 aes(x= CAP1, y = CAP2,
    #                     label = Sample),
    #                 size = 2,
    #                 hjust = 1,
    #                 max.overlaps=30) +
    # geom_point(data=species, aes(x=CAP1, y=CAP2)) +
    # geom_text_repel(data = species,
    #                 aes(x= CAP1, y = CAP2,
    #                     label = genus),
    #                 size = 2,
    #                 hjust = 1,
    #                 max.overlaps=30) +
    xlab(xtitle) +
    ylab(ytitle) +
    ggtitle(plot_title) +
    theme(text = element_text(size = 20))
}


bact_rda <- calc_rda(bact_ps, "spe ~ Treatment")
plot_rda(bact_rda)
ggsave(snakemake@output$bact_rda)


euk_rda <- calc_rda(euk_ps, "spe ~ Treatment")
plot_rda(euk_rda)
ggsave(snakemake@output$euk_rda)

bact_rda <-
  bact_ps %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  calc_rda("spe ~ Treatment")
plot_rda(bact_rda)
ggsave(snakemake@output$bact_rda_norm)


euk_rda <-
  euk_ps %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  calc_rda("spe ~ Treatment")
plot_rda(euk_rda)
ggsave(snakemake@output$euk_rda_norm)


bact_rda$species %>%
  select(-CAP1, -CAP2, -dist) %>%
  write_csv(snakemake@output$bact_species)


euk_rda$species %>%
  select(-CAP1, -CAP2, -dist) %>%
  write_csv(snakemake@output$euk_species)


bact_otus <-
  otu_table(bact_ps)
rows <-
  rownames(bact_otus)
bact_otus %>%
  data.frame %>%
  as_tibble %>%
  mutate(Sample = rows) %>%
  write_csv(snakemake@output$bact_otus)

euk_otus <-
  otu_table(euk_ps)
rows <-
  rownames(euk_otus)
euk_otus %>%
  data.frame %>%
  as_tibble %>%
  mutate(Sample = rows) %>%
  write_csv(snakemake@output$euk_otus)


