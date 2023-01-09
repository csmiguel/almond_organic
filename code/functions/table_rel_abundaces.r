# accessory function to create relative abundance tables
# Miguel Camacho-Sanchez https://orcid.org/0000-0002-6385-7963
# September 2022
table_rel_anbundace_ranks <- function(ps = NULL,
                                      glomtax = NULL) {
  # create table with relative abundances: cols are treatments and rows are taxa.
  # transform counts to relative
  psx <-
    phyloseq::tax_glom(ps,
                       taxrank = glomtax,
                       NArm = FALSE) %>%
    phyloseq::transform_sample_counts(function(x) { # transform to proportions
      x / sum(x) })

  # create table
  phyloseq::psmelt(psx) %>%
    mutate(rankn = eval(parse(text = glomtax))) %>%
    group_by(OTU, treatment, rankn) %>%
    summarise(av_abundance = mean(Abundance)) %>%
    reshape2::dcast(formula = rankn ~ treatment,
                    value.var = "av_abundance",
                    fun.aggregate = mean) %>%
    mutate(average = (CM + OM + NM)/3) %>%
    arrange(desc(average))
}
