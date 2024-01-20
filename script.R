library(dplyr)
library(stringr)
library(purrr)
library(worrms)
library(furrr)

dna_files <- list.files("output", "*DNADerivedData*", full.names = TRUE)
occurrence_files <- list.files("output", "*Occurrence*", full.names = TRUE)

dna <- map(dna_files, read.table, sep = "\t", quote = "", header = TRUE) %>%
  bind_rows() %>%
  mutate_if(is.character, na_if, "")

occurrence <- map(occurrence_files, read.table, sep = "\t", quote = "", header = TRUE) %>%
  bind_rows() %>%
  mutate_if(is.character, na_if, "") %>%
  mutate(
    species = ifelse(taxonRank == "species", scientificName, NA),
    aphiaid = as.numeric(str_replace(scientificNameID, "urn:lsid:marinespecies.org:taxname:", ""))
  ) %>%
  left_join(dna, by = "occurrenceID")

# resolve to accepted

resolve_to_accepted <- function(occurrence) {
  aphiaids <- unique(occurrence$aphiaid)
  aphiaid_batches <- split(aphiaids, as.integer((seq_along(aphiaids) - 1) / 50))
  plan(multisession, workers = 4)
  aphiaid_mapping <- future_map(aphiaid_batches, wm_record) %>%
    bind_rows() %>%
    select(aphiaid = AphiaID, valid_aphiaid = valid_AphiaID) %>%
    distinct() %>%
    filter(aphiaid != valid_aphiaid)
  
  valid_aphiaids <- unique(aphiaid_mapping$valid_aphiaid)
  valid_aphiaid_batches <- split(valid_aphiaids, as.integer((seq_along(valid_aphiaids) - 1) / 50))
  valid_taxa <- map(valid_aphiaid_batches, wm_record) %>%
    bind_rows() %>%
    select(valid_aphiaid = AphiaID, scientificName = scientificname, scientificNameID = lsid, taxonRank = rank, kingdom, phylum, class, order, family, genus) %>%
    mutate(taxonRank = tolower(taxonRank))
  occurrence %>%
    mutate(verbatimScientificName = scientificName) %>%
    left_join(aphiaid_mapping, by = "aphiaid") %>%
    rows_update(valid_taxa, by = "valid_aphiaid") %>%
    mutate(
      species = ifelse(taxonRank == "species", scientificName, NA),
      aphiaid = as.numeric(str_replace(scientificNameID, "urn:lsid:marinespecies.org:taxname:", ""))
    ) %>%
    select(-valid_aphiaid)
}

occurrence <- resolve_to_accepted(occurrence)

# stats

stats <- occurrence %>%
  filter(!is.na(higherGeography)) %>%
  group_by(higherGeography) %>%
  summarize(species = n_distinct(species), asvs = n_distinct(DNA_sequence), reads = sum(organismQuantity))

marker_stats <- occurrence %>%
  filter(!is.na(higherGeography)) %>%
  group_by(higherGeography, pcr_primer_name_forward) %>%
  summarize(species = n_distinct(species), asvs = n_distinct(DNA_sequence), reads = sum(organismQuantity))

family_stats <- occurrence %>%
  filter(!is.na(higherGeography)) %>%
  group_by(higherGeography, pcr_primer_name_forward, phylum, family) %>%
  summarize(species = n_distinct(species), asvs = n_distinct(DNA_sequence), reads = sum(organismQuantity))

write.csv(stats, "results/stats.csv", row.names = FALSE, quote = FALSE, na = "")
write.csv(marker_stats, "results/marker_stats.csv", row.names = FALSE, quote = FALSE, na = "")
write.csv(family_stats, "results/family_stats.csv", row.names = FALSE, quote = FALSE, na = "")
