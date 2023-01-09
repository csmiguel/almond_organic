# accessory functions to clean taxTable
# Miguel Camacho-Sanchez https://orcid.org/0000-0002-6385-7963
# September 2022
# this function takes the transformed output tax table from idTaxa and removes instances of "unidentified" ASVs matching "patterns"
clean_taxTable <- function(taxTable = NULL, patterns = NULL) {
# remove assigned taxonomy from first instance of "unidentified" to lowest taxonomic rank
    hh <-
    gsub("^.__", "", taxTable) %>%
    apply(1, function(asv) {
      first_instance_unidentified <- grep(patterns, asv)[1]
      if(!is.na(first_instance_unidentified)) {
        asv[first_instance_unidentified:length(asv)] <- NA
      }
    # remove duplicated tanomic names in lower taxonomic ranks
    asv[duplicated(asv)] <- NA
    asv
    }) %>% t()
    # place NAs whenever a species name has the format "Bacillus_sp."
    hh[, colnames(hh) != "species"]
    }

# as above but fills lower rank NA values with highest kown taxonomic rank.
clean_taxTable2 <- function(taxTable = NULL, patterns = NULL) {
  # remove assigned taxonomy from first instance of "unidentified" to lowest taxonomic rank
  hh <-
    gsub("^.__", "", taxTable) %>%
    apply(1, function(asv) {
      first_instance_unidentified <- grep(patterns, asv)[1]
      if(!is.na(first_instance_unidentified)) {
        asv[first_instance_unidentified:length(asv)] <- NA
      }
      # remove duplicated tanomic names in lower taxonomic ranks
      asv[duplicated(asv)] <- NA
      first_instance_NA <- which(is.na(asv))[1]
      if(!identical(first_instance_NA, integer(0)))
       lowestRankNamed <- asv[first_instance_NA - 1]
        asv[is.na(asv)] <- lowestRankNamed
        asv
    }) %>% t()
  hh[, colnames(hh) != "species"]
  }
