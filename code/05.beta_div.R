# beta diversity
# Miguel Camacho-Sanchez https://orcid.org/0000-0002-6385-7963
# September 2022

library(tidyverse)
library(flextable)
library(phyloseq)
library(ggrepel)
library(cowplot)
library(grid)

# read phyloseq
ps16s <-
  readRDS("data/intermediate/ps_bact.rds") %>%
  # transformation of the abundances to natural logarithmic scale
  transform_sample_counts(function(x) 1E6 * x/sum(x)) %>%
  phyloseq::transform_sample_counts(function(x) log(1 + x))

psITS <-
  readRDS("data/intermediate/ps_fung.rds") %>%
  # transformation of the abundances to natural logarithmic scale
  transform_sample_counts(function(x) 1E6 * x/sum(x)) %>%
  phyloseq::transform_sample_counts(function(x) log(1 + x))

# 1. MDS with Bray distances
bray_16s <-
  phyloseq::ordinate(ps16s,
                     method = "MDS",
                     distance = "bray")
bray_ITS <-
  phyloseq::ordinate(psITS,
                     method = "MDS",
                     distance = "bray")

#get eigenvalues
evals16s <- bray_16s$values$Eigenvalues
evalsITS <- bray_ITS$values$Eigenvalues

# col management
colomanag <- c('#7fc97f','#beaed4','#fdc086')

#plot MDS
# PC1 vs PC2
# plot_ordination returns error whenever I am using sample_data(ps) with one column only. Therefore, I create extra artificial variables to bypass this error.
# 16s
h <- ps16s
sample_data(h) <-
  as(sample_data(ps16s), "data.frame") %>%
  mutate(samples = sample_names(h), random = rep(c("red", "blue", "yellow"),3)) %>%
  rename(management = treatment)

# replace label "Axis" with "Dimension
colnames(bray_16s$vectors) <-
  colnames(bray_16s$vectors) %>% stringr::str_replace("Axis.", "Dimension_")

# plot
plot16s <-
  plot_ordination(h,
                  bray_16s,
                  type = "samples",
                  axes = c(1, 2),
                  color = "management") +
  labs("Management") +
  geom_point(size = 2) +
  coord_fixed(sqrt(evals16s[2] / evals16s[1])) +
  ggrepel::geom_text_repel(
    aes(label = samples),
    size = 3,
    min.segment.length = 1,
    point.padding = 0.5) +
  scale_colour_manual(
    values = colomanag,
    breaks=c("OM", "CM", "NM")) +
  theme_cowplot() +
  theme(legend.position="none")

#ITS
hh <- psITS
sample_data(hh) <-
  as(sample_data(psITS), "data.frame") %>%
  mutate(samples = sample_names(hh), random = rep(c("red", "blue", "yellow"),3)) %>%
  rename(management = treatment)

# replace label "Axis" with "Dimension
colnames(bray_ITS$vectors) <-
  colnames(bray_ITS$vectors) %>% stringr::str_replace("Axis.", "Dimension_")


plotITS <-
  plot_ordination(hh,
                  bray_ITS,
                  type = "samples",
                  axes = c(1, 2),
                  color = "management") +
  labs("Management") +
  geom_point(size = 2) +
  coord_fixed(sqrt(evalsITS[2] / evalsITS[1])) +
  ggrepel::geom_text_repel(
    aes(label = samples),
    size = 3,
    min.segment.length = 1,
    point.padding = 0.5) +
  scale_colour_manual(
    values = colomanag,
    breaks=c("OM", "CM", "NM")) +
  theme_cowplot() +
  theme(legend.position = "top")

pmds <- cowplot::plot_grid(plotITS,
                   plot16s, ncol = 2,
                   labels = c("A. Fungi", "B. Bacteria"),
                   vjust = 4)

ggsave("output/mds.pdf",
       pmds,
       width = 10,
       height = 6)
# plot dimensions
#16s
pcs16s <-
  data.frame(value = 100 * bray_16s$values$Relative_eig,
             Dimension = seq_along(bray_16s$values$Relative_eig)) %>%
  head(n = 8)

p_dims16s <-
  ggplot2::ggplot(aes(x = Dimension, y = value),
                  data = pcs16s) +
  geom_bar(stat = "identity") +
  scale_x_continuous(labels = pcs16s$Dimension, breaks = pcs16s$Dimension) +
  xlab("Dimension") +
  ylab("Variation %") +
  theme_cowplot()
#ITS
pcsITS <-
  data.frame(value = 100 * bray_ITS$values$Relative_eig,
             Dimension = seq_along(bray_ITS$values$Relative_eig)) %>%
  head(n = 8)

p_dimsITS <-
  ggplot2::ggplot(aes(x = Dimension, y = value),
                  data = pcsITS) +
  geom_bar(stat = "identity") +
  scale_x_continuous(labels = pcsITS$Dimension, breaks = pcsITS$Dimension) +
  xlab("Dimension") +
  ylab("Variation %") +
  theme_cowplot()

###########
# 2. ADONIS
#permutations
nperm <- 9999
set.seed(123)

#16s
#data frame with metadata
metadf16s <- data.frame(sample_data(ps16s))
# weighted unifrac distances
braydist16s <- phyloseq::distance(ps16s, method = "bray")
t16s <-
  vegan::adonis(
    braydist16s ~ treatment,
    data = metadf16s,
    permutations = nperm)

# ITS
#data frame with metadata
metadfITS <- data.frame(sample_data(psITS))
# weighted unifrac distances
braydistITS <- phyloseq::distance(psITS, method = "bray")
tITS <-
  vegan::adonis(
    braydistITS ~ treatment,
    data = metadfITS,
    permutations = nperm)

#summary table for adonis results
table_adonis <-
  list(t16s, tITS) %>%
  lapply(function(x) {
    x$aov.tab %>%
      as("data.frame") %>%
      cbind(formula  = x$call %>%
              as.character() %>% .[2]) %>%
      dplyr::select(formula, Df, `F.Model`, R2, `Pr(>F)`)  %>%
      .[1, ]
  }) %>%
  do.call(what = rbind) %>%
  dplyr::mutate(`F.Model` = round(`F.Model`, 2),
                R2 = round(R2, 2),
                `Pr(>F)` = round(`Pr(>F)`, 4))
ft <- flextable::flextable(table_adonis)
ft %>%
  flextable::save_as_docx(path = "output/summary_table_adonis.docx")

### plot "everything" related to beta diversity
# composite plot
# top row: mods
top_row <-
  plot_grid(plot16s + theme(legend.position = "none"),
            plotITS,
            labels = c("A", "B"))

# middle row PCs
middle_row <-
  plot_grid(p_dims16s,
            p_dimsITS,
            labels = c("C", "D"))
# bottom row
# flextable to ggplot
ft_raster <- as_raster(ft)
ftgg <-
  ggplot() + 
  theme_void() + 
  annotation_custom(rasterGrob(ft_raster), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
bottom_row <-
  plot_grid(ftgg,
            labels = c("E"))
pall <-
  plot_grid(top_row, 
            middle_row,
            bottom_row,
            ncol = 1, rel_heights = c(3,2,1))


#save plots
ggsave("output/beta_div.pdf",
       pall,
       width = 9,
       height = 7,
       units = "in")
