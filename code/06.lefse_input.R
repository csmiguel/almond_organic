# create input for LEfSe from BIOM tables.
# I create phyloseq object. I clean the taxa_table with a custom function called clean_taxTable2
# in which lower taxonomic ranks are filled with the highest full tax rank. Then, it resturns a text file
# which is the input for LEfSe in a dedicated Galaxy Server.http://huttenhower.sph.harvard.edu/galaxy/
# Miguel Camacho-Sanchez https://orcid.org/0000-0002-6385-7963
# September 2022
library(phyloseq)
library(dplyr)
library(phyloseqCompanion)
# load important function
source("code/functions/clean_taxTable.r")

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
    tax_table(clean_taxTable(tax_table(temp_bact), "unidentifi|uncultu|metagenome|ineage|Subgroup")),
    sample_data(samp_data)) %>%
  tax_glom("genus", NArm = F)
# format data for lefse
phyloseqCompanion::phyloseq2lefse(ps_bact,
                                  file.name = "data/intermediate/16s_lefse.txt",
                                  covars = "treatment",
                                  taxa.levels = rank_names(ps_bact),
                                  transpose.otus = taxa_are_rows(ps_bact))

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
    tax_table(clean_taxTable(tax_table(temp_fung), "unidentifi|uncultu|metagenome|ineage|Subgroup")),
    sample_data(samp_data)) %>%
  tax_glom("genus", NArm = F)

phyloseqCompanion::phyloseq2lefse(ps_fung,
                                  file.name = "data/intermediate/its_lefse.txt",
                                  covars = "treatment",
                                  taxa.levels = rank_names(ps_fung),
                                  transpose.otus = taxa_are_rows(ps_fung))

# next steps:
# 2. import in Galaxy as "tabular" http://huttenhower.sph.harvard.edu/galaxy/
# 3. run "Format data for LEfSe" and "LDA Effect Size".
# 4. download output from analysis to "data/intermediate" and remove those rows containing ".NA"
# cat data/intermediate/its_int_res | grep  -v '.NA' > data/intermediate/its_int_res_noNAs
# 5. import its_int_res_noNAs to galaxy  as "lefse_internal_res".
# 6. do plots.