# plots and calculations for relative amount of bacterial and fungal DNA in soil
# Miguel Camacho-Sanchez https://orcid.org/0000-0002-6385-7963
# September 2022

library(plyr)
library(tidyverse)
library(cowplot)

# relative values
# the input data are relative copies of 16s and ITS based on standards used in the qPCR. Therefore, they can only be # expressed in relative amounts, unless genome size and no of copies are known for the standards used.

qpcr <-
  read.delim("data/raw/qpcr.txt") %>% # raw data
  # compute foldchange of bact respect to fungi
  plyr::ddply(~sample, function(x) {
    qbact <-
      filter(x, org == "bact") %>% pull(quantity)
    qfung <-
      filter(x, org == "fung") %>% pull(quantity)
    x %>%
      mutate(
        foldchange_fb = qfung / qbact)
  }) %>%
# make all values relative to 1 per factor
  plyr::ddply(~org, function(x) { # compute values relative the the maximum values observed.
    x %>%
      dplyr::mutate(
        rel_qt = if (unique(org) == "bact"){
        quantity / max(quantity)
      } else {
        quantity / max(quantity)
      })
  }) %>%
  mutate(
    numsamples = as.numeric(as.factor(sample)),
    treat = gsub("[0-9]", "", sample) %>%
            factor(levels = c("OM", "CM", "NM")),
    org = factor(org, levels = c("fung", "bact")))


p1 <-
  ggplot(qpcr, aes(x = treat,
                 y = rel_qt,
                 fill = org)) +
  geom_boxplot() +
  geom_point(aes(fill = org), size = 5,
  shape = 21, position = position_jitterdodge(jitter.width = .1)) +
  ylab("Relative number of copies from ITS or 16S") +
  xlab("Management") +
  scale_fill_manual(name = "",
                      labels = c("ITS (Fungi)", "16S (Bacteria)"),
                      values = c("grey30", "grey70")) +
  theme_classic() +
  theme(
    legend.position = c(.7, .9)
  )

# plot on relative change of bact respect to fungi
p2 <-
filter(qpcr, org == "bact") %>%
  ggplot(aes(x = treat,
             y = foldchange_fb)) +
  geom_boxplot() +
  geom_jitter(width = .1, size = 3, color = "black", alpha = 0.6) +
  ylab("Copies of ITS/16S") +
  xlab("Management") +
  theme_classic()

# write data frame to table
write.csv(qpcr, file = "output/qpcr_calculations.txt")
# save plot
pall <-
  cowplot::plot_grid(p1, p2,
                     labels = c("A", "B"),
                     rel_widths = c(1, .6))
ggsave(filename = "output/qpcr_plots.pdf",
       plot = pall,
       height = 4,
       width = 8)
