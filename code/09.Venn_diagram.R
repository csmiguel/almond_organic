# Venn Diagram
# Miguel Camacho-Sanchez https://orcid.org/0000-0002-6385-7963
# September 2022

library(ggVennDiagram)
library(plyr) # always attach before dplyr to avoid conflicts
library(phyloseq)
library(tidyverse)
library(cowplot)

# create list of genera per treatment for Bacteria
bact <-
  readRDS("data/intermediate/ps_bact.rds") %>%
  tax_glom("genus") %>%
  psmelt() %>%
  group_by(treatment, genus) %>%
  mutate(sumreps = sum(Abundance),
         count = n()) %>%
  dplyr::select(genus, treatment, sumreps, count) %>%
  dplyr::filter(sumreps > 0 & count == 3) %>%
  plyr::dlply(~treatment, function(x) {
    x$genus
  })

# plot Venn bact
pbact <-
  ggVennDiagram(bact) +
  scale_fill_gradient(low="white",high = "grey20")

# create list of genera per treatment for Fungi
fung <-
  readRDS("data/intermediate/ps_fung.rds") %>%
  tax_glom("genus") %>%
  psmelt() %>%
  group_by(treatment, genus) %>%
  mutate(sumreps = sum(Abundance),
         count = n()) %>%
  dplyr::select(genus, treatment, sumreps, count) %>%
  dplyr::filter(sumreps > 0 & count == 3) %>%
  plyr::dlply(~treatment, function(x) {
    x$genus
  })

# plot Venn fungi
pfun <-
  ggVennDiagram(fung) +
  scale_fill_gradient(low="white",high = "grey20")

# plot all and save
pall <-
  cowplot::plot_grid(pbact, pfun, labels = c("A. Bacteria", "B. Fungi"))

ggsave(filename = "output/Venn_genera.pdf",
       plot = pall,
       width = 8,
       height = 5.2)
