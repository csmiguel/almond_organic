# funguilds
# Miguel Camacho-Sanchez https://orcid.org/0000-0002-6385-7963
# September 2022

library(tidyverse)
library(xlsx)
library(reshape2)

fung <-
  xlsx::read.xlsx("data/raw/Fungal_functional_properties_predictions.xlsx", 1) %>%
  filter(Confidence.Ranking %in% c("Probable", "Highly Probable"),
         Guild != "NULL") %>% # filter confident guilds
  mutate( # rename variables and transform read counts to relative abundances
    CM1 = CON1 / sum(CON1),
    CM2 = CON2 / sum(CON2),
    CM3 = CON3 / sum(CON3),
    OM1 = ECO1 / sum(ECO1),
    OM2 = ECO2 / sum(ECO2),
    OM3 = ECO3 / sum(ECO3),
    NM1 = NM1 / sum(NM1),
    NM2 = NM2 / sum(NM2),
    NM3 = NM3 / sum(NM3)) %>%
  select(CM1, CM2, CM3, OM1, OM2, OM3, NM1, NM2, NM3, Guild) %>%
  melt(value.name = "abundance", variable.name = "management", id.vars = "Guild") %>%
  mutate(guild = gsub(pattern = ".*aprotroph.*", replacement = "Saprotroph", x = Guild) %>% # merge guilds
           gsub(pattern = ".*Animal.*", replacement = "Animal related"))

# determine order of guilds by relative abundance
levels_guild <-
  fung %>%
  ddply(~guild, function(x){
    sum(x$abundance)
  }) %>%
  arrange(desc(V1)) %>%
  pull(guild)
h <- # order levels for plotting purposes
  fung %>%
  mutate(guild = factor(guild, levels = levels_guild),
    management = factor(management,
    levels = c("OM1", "OM2", "OM3", "CM1", "CM2", "CM3", "NM1", "NM2", "NM3")))

p_guild <-
  ggplot(h) +
  geom_bar(aes(x = management, y = abundance, fill = guild), stat = "identity") +
  scale_fill_manual(values = c('#7fc97f','#bf5b17','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f')) +
  ylab("Relative abundance") +
  theme(legend.position="none") +
  theme_classic()
ggsave("output/guilds_fungi.pdf", p_guild, width = 7, height = 5)

# save table
saveRDS(fung, "data/intermediate/fung_guilds.rds")

# make table with averages
fung_av <-
  fung %>%
  mutate(manag = as.factor(
    stringr::str_remove(as.character(management), "[0-9]"))) %>%
  group_by(manag, guild) %>%
  dplyr::summarise(avab = round(sum(abundance) / 3, 2)) %>%
  dcast(guild ~ manag, value.var = "avab")
