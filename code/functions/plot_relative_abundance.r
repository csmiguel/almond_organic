# accessory functions to plot relative abundance by tax rank
# Miguel Camacho-Sanchez https://orcid.org/0000-0002-6385-7963
# September 2022
plot_rel_abundace <- function(ps = NULL, ntop = NULL, glomtax = "phylum", pal_col = "Paired") {
  # ps, phyloseq object
  # ntop, number of taxa with highest abundances to plot
  # glomtax, rank to agglomerate by
  # pal_col, name of the palette color to use
  # transform counts and agglomerate
  psx <-
    phyloseq::tax_glom(ps,
                       taxrank = glomtax,
                       NArm = FALSE) %>%
    phyloseq::transform_sample_counts(function(x) { # transform to proportions
      x / sum(x) })
  #

  if (phyloseq::taxa_are_rows(psx)) {
    h <- rowSums(phyloseq::otu_table(psx))
  } else if (!phyloseq::taxa_are_rows(psx)) {
    h <- colSums(phyloseq::otu_table(psx))
  }
  names_ntop <-
    h %>%
    sort(decreasing = T) %>%
    names() %>%
    .[1:ntop]
  #tax_table
  tax_table_rank <-
    phyloseq::tax_table(psx)@.Data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("asv")
  #vector with most abundant phyla
  names_legend <-
    c(tax_table_rank[[glomtax]][match(names_ntop, tax_table_rank$asv)],
      "other")
  # vector sample order
  sample_order <- c("OM1", "OM2", "OM3", "CM1", "CM2", "CM3", "NM1", "NM2", "NM3")
  # agglomerate low frequency taxa
  tidy_counts <-
    phyloseq::psmelt(psx) %>%
    as_tibble() %>%
    dplyr::select(OTU, Sample, Abundance, treatment, all_of(glomtax)) %>%
    dplyr::mutate(labels_ntop = ifelse(OTU %in% names_ntop,
                                eval(parse(text = glomtax)),
                                "other")) %>%
    group_by(Sample, labels_ntop) %>%
    summarise(sum_tax = sum(Abundance)) %>%
    dplyr::mutate(tax_col = factor(labels_ntop,
                            levels = names_legend)) %>%
    dplyr::select(-labels_ntop) %>%
    dplyr::mutate(Sample = factor(Sample,
                           levels = sample_order))

  # relative abundance plots
  ggplot(tidy_counts) +
    geom_bar(aes(x = Sample, y = sum_tax, fill = tax_col),
             stat = "identity",
             position = "stack",
             color = "black",
             size = 0.1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = -90, hjust = 0),
          strip.background = element_blank(),
          strip.text.x = element_text(size = 10, face = "italic"),
          legend.text = element_text(size = 8)) +
    scale_fill_brewer(palette = pal_col, name = glomtax) +
    ylab("Relative abundance") +
    xlab(NULL)
}
