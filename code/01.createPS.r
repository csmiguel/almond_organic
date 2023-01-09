# create phyloseq object from BIOM data from AllGenetics
# Miguel Camacho-Sanchez https://orcid.org/0000-0002-6385-7963
# September 2022
library(phyloseq)
library(dplyr)

source("code/functions/clean_taxTable.r")
# read biom and export phyloseq object
# set ranks
myranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

#### 16S
# read biom files from filtered asvs
temp_bact <-
  phyloseq::import_biom("data/raw/16S_table_after_filtering.biom")
colnames(tax_table(temp_bact)) <- myranks
# edit sample names
sn1 <-
    gsub("CON", "CM", sample_names(temp_bact)) %>%
    gsub(pattern = "ECO", replacement = "OM")
stopifnot(sample_names(temp_bact) == dimnames(temp_bact@otu_table@.Data)[[2]])
dimnames(temp_bact@otu_table@.Data)[[2]] <- sn1

# create table with sample data
samp_data <-
  data.frame(
    treatment = gsub(pattern = "[0-9]",
                     replacement = "",
                     sn1),
    row.names = sn1)

# add sample data to phyloseq
ps_bact <-
  phyloseq::phyloseq(
    otu_table(temp_bact),
    tax_table(clean_taxTable(tax_table(temp_bact), "unidentifi|uncultu|metagenome")),
    sample_data(samp_data))

### ITS
temp_fung <-
  phyloseq::import_biom("data/raw/ITS_table_after_filtering.biom")
colnames(tax_table(temp_fung)) <- myranks

# add sample data to phyloseq
stopifnot(sample_names(temp_fung) == dimnames(temp_fung@otu_table@.Data)[[2]])
dimnames(temp_fung@otu_table@.Data)[[2]] <- sn1

ps_fung <-
  phyloseq::phyloseq(
    otu_table(temp_fung),
    tax_table(clean_taxTable(tax_table(temp_fung), "unidentifi|uncultu|metagenome")),
    sample_data(samp_data))

# write phyloseq objects
saveRDS(ps_bact, "data/intermediate/ps_bact.rds")
saveRDS(ps_fung, "data/intermediate/ps_fung.rds")
