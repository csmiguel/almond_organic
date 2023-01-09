# create relative abundance plots
# Miguel Camacho-Sanchez https://orcid.org/0000-0002-6385-7963
# September 2022

library(phyloseq)
library(dplyr)
library(ggplot2)
library(cowplot)

# function for plotting
source("code/functions/plot_relative_abundance.r")

# write phyloseq objects
ps_bact <- readRDS("data/intermediate/ps_bact.rds")
ps_fung <- readRDS("data/intermediate/ps_fung.rds")

# plot rel abundances bact
pbact_p <-
  plot_rel_abundace(
  ps = ps_bact,
  ntop = 12,
  glomtax = "phylum",
  pal_col = "Paired"
)
pbact_o <-
  plot_rel_abundace(
  ps = ps_bact,
  ntop = 12,
  glomtax = "order",
  pal_col = "Paired"
)
pbact_f <-
  plot_rel_abundace(
  ps = ps_bact,
  ntop = 12,
  glomtax = "family",
  pal_col = "Paired"
)

# plot rel abundances fungi
pfungi_p <-
  plot_rel_abundace(
    ps = ps_fung,
    ntop = 12,
    glomtax = "phylum",
    pal_col = "Paired"
  )
pfungi_o <-
  plot_rel_abundace(
    ps = ps_fung,
    ntop = 12,
    glomtax = "order",
    pal_col = "Paired"
  )
pfungi_f <-
  plot_rel_abundace(
    ps = ps_fung,
    ntop = 12,
    glomtax = "family",
    pal_col = "Paired"
  )

# compose plot
p <- cowplot::plot_grid(pfungi_p, pfungi_o, pfungi_f,
                        pbact_p, pbact_o, pbact_f,
                        labels = LETTERS[1:6],
                        ncol = 3,
                        nrow = 2,
                        byrow = TRUE)
# save plot
ggsave("output/rel_ab.pdf",
       p,
       width = 30,
       height = 20,
       units = "cm")

# merged by treatments
source("code/functions/merge_ps_treatment.r")
merged_ps_bact <- merge_ps_treatment(ps2merge = ps_bact)
merged_ps_fung <- merge_ps_treatment(ps2merge = ps_fung)

# plot rel abundances bact
mpbact_p <-
  plot_rel_abundace(
  ps = merged_ps_bact,
  ntop = 12,
  glomtax = "phylum",
  pal_col = "Paired"
)
mpbact_o <-
  plot_rel_abundace(
  ps = merged_ps_bact,
  ntop = 12,
  glomtax = "order",
  pal_col = "Paired"
)
mpbact_f <-
  plot_rel_abundace(
  ps = merged_ps_bact,
  ntop = 12,
  glomtax = "family",
  pal_col = "Paired"
)

# plot rel abundances fungi
mpfungi_p <-
  plot_rel_abundace(
    ps = merged_ps_fung,
    ntop = 12,
    glomtax = "phylum",
    pal_col = "Paired"
  )
mpfungi_o <-
  plot_rel_abundace(
    ps = merged_ps_fung,
    ntop = 12,
    glomtax = "order",
    pal_col = "Paired"
  )
mpfungi_f <-
  plot_rel_abundace(
    ps = merged_ps_fung,
    ntop = 12,
    glomtax = "family",
    pal_col = "Paired"
  )

# compose plot
p2 <- cowplot::plot_grid(mpbact_p, mpbact_o, mpbact_f,
                        mpfungi_p, mpfungi_o, mpfungi_f,
                        labels = LETTERS[1:6],
                        ncol = 3,
                        nrow = 2,
                        byrow = TRUE)
# save plot
ggsave("output/rel_ab_merged.pdf",
       p2,
       width = 27,
       height = 20,
       units = "cm")
