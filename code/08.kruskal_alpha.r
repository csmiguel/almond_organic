# alpha diversity models
# Miguel Camacho-Sanchez https://orcid.org/0000-0002-6385-7963
# September 2022

library(phyloseq)
library(tidyverse)

# test Kruskal Wallis to evaluate effect of management of alpha diversity.

# 16S
ps_bact <- readRDS("data/intermediate/ps_bact.rds")

psbact_rar <- phyloseq::rarefy_even_depth(physeq = ps_bact, sample.size = 67200)

alpha_bact <-
  phyloseq::estimate_richness(psbact_rar) %>%
  tibble::rownames_to_column("sample") %>%
  mutate(manag = str_remove(sample, "[0-9]"))

# ITS
ps_its <- readRDS("data/intermediate/ps_fung.rds")

psfung_rar <- phyloseq::rarefy_even_depth(physeq = ps_its, sample.size = 34800)

alpha_fung <-
  phyloseq::estimate_richness(psfung_rar) %>%
  tibble::rownames_to_column("sample") %>%
  mutate(manag = str_remove(sample, "[0-9]"))

# save results to text file
sink(file = "output/kruskal-alpha.txt")

cat("\nBacteria\n")
kruskal.test(Observed ~ manag, data = alpha_bact)
kruskal.test(Shannon ~ manag, data = alpha_bact)

cat("\nFungi\n")
kruskal.test(Observed ~ manag, data = alpha_fung)
kruskal.test(Shannon ~ manag, data = alpha_fung)

sink()
