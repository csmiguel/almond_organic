# create relative abundanc tables from phyloseq objects.
# Miguel Camacho-Sanchez https://orcid.org/0000-0002-6385-7963
# September 2022

library(phyloseq)
library(phyloseq)
library(dplyr)
library(reshape2)
library(xlsx)

# load function
source("code/functions/table_rel_abundaces.r")

# write phyloseq objects
ps_bact <- readRDS("data/intermediate/ps_bact.rds")
ps_fung <- readRDS("data/intermediate/ps_fung.rds")

## 16S
# create Excel
xlsx::write.xlsx("Relative abundance data by taxonomic rank",
                 "output/rel_ab_bact.xlsx")

# write tables to Excel
rank_names(ps_bact)[1:6] %>%
  lapply(function(x) {
    h <- table_rel_anbundace_ranks(ps_bact, glomtax = x)
    xlsx::write.xlsx(x = h,
                     file = "output/rel_ab_bact.xlsx",
                     sheetName =  x,
                     append = T)
  })


## ITS
# create Excel
xlsx::write.xlsx("Relative abundance data by taxonomic rank",
                 "output/rel_ab_fungi.xlsx")

# write tables to Excel
rank_names(ps_fung)[1:6] %>%
  lapply(function(x) {
    h <- table_rel_anbundace_ranks(ps_fung, glomtax = x)
    xlsx::write.xlsx(x = h,
                     file = "output/rel_ab_fungi.xlsx",
                     sheetName =  x,
                     append = T)
  })

