# plot pathogenic fungi. Fig 6 in publication
# Miguel Camacho-Sanchez https://orcid.org/0000-0002-6385-7963
# November 2022.
library(phyloseq)
library(tidyverse)
library(ggpubr)

psfung <-
  readRDS("data/intermediate/ps_fung.rds") %>%
  tax_glom("genus", NArm = FALSE) %>%
  phyloseq::transform_sample_counts(function(x) x / sum(x)) %>%
  psmelt()

pathogenic <-
  c("Fusarium",
    "Cylindrocarpon")
dfpath <-
  psfung %>%
  dplyr::filter(genus %in% pathogenic)

p1 <-
  ggbarplot(dfpath, x = "treatment",
          y = "Abundance", color = "genus",
          add = c("mean_se", "jitter"),
          label = TRUE, lab.nb.digits = 2,
          palette = c("#00AFBB", "#E7B800"),
          position = position_dodge(),
          xlab = "Management",
          ylab = "Relative abundance")
ggsave("output/path_fung.pdf",
       p1,
       width = 5,
       height = 4)
