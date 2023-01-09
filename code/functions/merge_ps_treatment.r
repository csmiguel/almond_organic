# accessory functions to merge phyloseq by factor
# Miguel Camacho-Sanchez https://orcid.org/0000-0002-6385-7963
# September 2022
# since phyloseq::merge_samples is not working for my phyloseq objects (I think
# the issue is related to the fact that I have only variable in sample_data),
# I created my own function to merge ps objects by treatment.
merge_ps_treatment <- function(ps2merge = NULL) {
  ps_rel  <- phyloseq::transform_sample_counts(ps2merge, function(x) x / sum(x))
  merged_otu <-
    otu_table(ps_rel) %>%
    as.data.frame() %>%
    apply(1, function(x){
      c(sum(x[1:3]), sum(x[4:6]), sum(x[7:9]))
    }) %>%
    t %>%
    as.data.frame %>%
    setNames(c("CM", "OM", "NM"))
  merged_sample <-
    data.frame(treatment = c("CM", "OM", "NM"), row.names = c("CM", "OM", "NM"))
  phyloseq(
    otu_table(merged_otu, taxa_are_rows = T),
    sample_data(merged_sample),
    tax_table(ps2merge)
    )
  }
