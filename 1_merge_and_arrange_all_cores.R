# TAXALOSS Project - Import and merge all OBITools 3 processed data
# Script from Jeremy Courtin
# Script last update - 16.03.2023

# Layout of the script: 
# 0 - GENERAL OPTION - SETTING UP SCRIPT
# 1 - Set variables
# 2 - Read in the raw data from different cores + plots
# 3 - Merge all cores
# 4 - Get only ASVs present with a minimum of 9 reads / in 9 samples
# 5 - Save the data

############################################################################
# 0 - GENERAL OPTION - SETTING UP SCRIPT
############################################################################
rm(list = ls()) # Remove all the objects we created so far.

# Set up options for the rest of the script
options(stringsAsFactors=FALSE) # states that everything is reads as character 

#Load in R packages
library("tidyverse")
library("readxl")
library("vegan")
library("scales")
library("rioja")
library("tidypaleo")
theme_set(theme_paleo(15))


############################################################################
# 1 - Set variables
############################################################################
#btoko
seq1 <- "HUA-9"
core1 <- "PG2133"

#ilerney16kp
seq2 <- "ALRK-5"
core2 <- "16-KP-01-L02-Long3"

#Bilyakh
seq3 <- "ALRK-8"
core3 <- "PG1755"

#Levinson Lessing
seq4 <- "ALRK-11"
core4 <- "Co1401"

#Ilerney
seq5 <- "ALRK-10"
core5 <- "EN18208"

# E5
seq6 <- "ALRK-13"
core6 <- "E5"

# Rauchuagytgyn
seq7 <- "APMG-30"
core7 <- "EN18218"

# Emanda
seq8 <- "APMG-42"
core8 <- "Co1412"

############################################################################
# 2 - Read in the raw data from different cores + plots
############################################################################
######################
# Bolshoe Toko
######################
#tl_btoko_raw <- read_delim(paste("cores/", seq1, "_", core1, "_file02.2_identitylevel_90_merged_replicates.csv", sep = ""), delim = ";", col_names = T)
tl_btoko_raw <- read_delim(paste("cores/", seq1, "_", core1, "_file02.2_identitylevel_90_merged_replicates_new.csv", sep = ""), delim = ";", col_names = T)
tl_btoko_raw %>% select(best_identity, NUC_SEQ_arc) %>% filter(best_identity == 1) %>% distinct()
tl_btoko_raw %>% select(best_identity, NUC_SEQ_arc) %>% filter(best_identity < 1) %>% distinct()

btoko_asv <- tl_btoko_raw %>% select(NUC_SEQ_arc, best_identity, best_family, best_genus, best_species, scientific_name)
tl_btoko <- tl_btoko_raw %>% mutate(taxaid = paste(NUC_SEQ_arc, best_identity, best_family, best_genus, best_species, scientific_name, sep = "_")) %>% select(-NUC_SEQ_arc, -best_identity, -best_family, -best_genus, -best_species, -scientific_name)
btoko_metadata <- read_excel("metadata/cores/HUA-9_BToko_agefile.xlsx", col_names = T)
btoko_samp <- paste(btoko_metadata$Extraction_number, btoko_metadata$Mean_collection_depth_round, btoko_metadata$Section_age_round, sep = "_")

tl_btoko_long <- tl_btoko %>% pivot_longer(!taxaid, names_to = "sampleid", values_to = "count") %>% filter(count > 0) %>% group_by(taxaid) %>% mutate(total = sum(count)) %>% ungroup() %>% mutate(sampleid = paste("btoko", sampleid, sep = "_"))
tl_btoko_long_sup10 <- tl_btoko_long %>% filter(total > 9) %>% arrange(total)

#####
# Check and plots
#####
tl_btoko_long_sup10 %>% filter(grepl("_1_",taxaid)) %>% filter(grepl("Rosace",taxaid)) %>% select(taxaid) %>% distinct()

wide <- tl_btoko_long_sup10 %>% filter(grepl("_1_",taxaid)) %>% pivot_wider(names_from = taxaid, values_from = count) # make the data from long to wide
wide[is.na(wide)] <- 0 
t_plot <- pivot_longer(wide, cols = names(wide[,-c(1,2)]), names_to = "taxaid", values_to = "count")

# Plot on genus level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!genus == "NA") %>%
  dplyr::group_by(genus, age) %>% dplyr::mutate(genus_count = sum(count)) %>% ungroup() %>% select(-count, -species) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(genus_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = genus_count/tot_sample*100) %>% arrange(desc(rel_abund)) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, genus) %>% dplyr::mutate(keep = as.numeric(sum(rel_abund))) %>% ungroup() %>% 
  dplyr::group_by(genus) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 4.99) %>%
  arrange(family) %>% select(genus) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!genus == "NA") %>%
  dplyr::group_by(genus, age) %>% dplyr::mutate(genus_count = sum(count)) %>% ungroup() %>% select(-count, -species) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(genus_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = genus_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, genus) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(genus) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>% mutate(max = ifelse(is.na(max), 0, max)) %>%
  subset(max > 4.99) %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(genus)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq1, "_", core1, "_stratplot_100precent_genus_old.png"), width = 60, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq1, "_", core1, "_stratplot_100precent_genus_new.png"), width = 60, height = 20, units = "cm", dpi = 300)

# Plot on family level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!family == "NA") %>%
  dplyr::group_by(family, age) %>% dplyr::mutate(family_count = sum(count)) %>% ungroup() %>% select(-count, -species, -genus) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(family_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = family_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, family) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(family) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 4.99) %>%
  arrange(family) %>% select(family) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!family == "NA") %>%
  dplyr::group_by(family, age) %>% dplyr::mutate(family_count = sum(count)) %>% ungroup() %>% select(-count, -species, -genus) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(family_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = family_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, family) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(family) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 4.99)  %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(family)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq1, "_", core1, "_stratplot_100precent_family_old.png"), width = 60, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq1, "_", core1, "_stratplot_100precent_family_new.png"), width = 60, height = 20, units = "cm", dpi = 300)

# Plot on species level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!species == "NA") %>%
  dplyr::group_by(species, age) %>% dplyr::mutate(species_count = sum(count)) %>% ungroup() %>% select(-count) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(species_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = species_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, species) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(species) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 1) %>%
  arrange(family) %>% select(species) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!species == "NA") %>%
  dplyr::group_by(species, age) %>% dplyr::mutate(species_count = sum(count)) %>% ungroup() %>% select(-count) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(species_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = species_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, species) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(species) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 1) %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(species)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq1, "_", core1, "_stratplot_100precent_species_old.png"), width = 100, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq1, "_", core1, "_stratplot_100precent_species_new.png"), width = 100, height = 20, units = "cm", dpi = 300)

######################
# Ilerney 16KP
######################
#tl_ilerney16kp_raw <- read_delim(paste("cores/", seq2, "_", core2, "_file02.2_identitylevel_90_merged_replicates.csv", sep = ""), delim = ";", col_names = T)
tl_ilerney16kp_raw <- read_delim(paste("cores/", seq2, "_", core2, "_file02.2_identitylevel_90_merged_replicates_new.csv", sep = ""), delim = ";", col_names = T)
ilerney16kp_asv <- tl_ilerney16kp_raw %>% select(NUC_SEQ_arc, best_identity, best_family, best_genus, best_species, scientific_name)
tl_ilerney16kp <- tl_ilerney16kp_raw %>% mutate(taxaid = paste(NUC_SEQ_arc, best_identity, best_family, best_genus, best_species, scientific_name, sep = "_")) %>% select(-NUC_SEQ_arc, -best_identity, -best_family, -best_genus, -best_species, -scientific_name)
ilerney16kp_metadata <- read_excel("metadata/cores/ALRK-5_Ilirney_16KP_agefile.xlsx", col_names = T)
ilerney16kp_samp <- paste(ilerney16kp_metadata$Extraction_number, ilerney16kp_metadata$Mean_collection_depth_round, ilerney16kp_metadata$Section_age_round, sep = "_")

tl_ilerney16kp_long <- tl_ilerney16kp %>% pivot_longer(!taxaid, names_to = "sampleid", values_to = "count") %>% filter(count > 0) %>% group_by(taxaid) %>% mutate(total = sum(count)) %>% ungroup() %>% mutate(sampleid = paste("ilerney16kp", sampleid, sep = "_"))
tl_ilerney16kp_long_sup10 <- tl_ilerney16kp_long %>% filter(total > 9) %>% arrange(total)

#####
# Check and plots
#####
tl_ilerney16kp_long_sup10 %>% filter(grepl("_1_",taxaid)) %>% filter(grepl("Rosace",taxaid)) %>% select(taxaid) %>% distinct()

wide <- tl_ilerney16kp_long_sup10 %>% filter(grepl("_1_",taxaid)) %>% pivot_wider(names_from = taxaid, values_from = count) # make the data from long to wide
wide[is.na(wide)] <- 0 
t_plot <- pivot_longer(wide, cols = names(wide[,-c(1,2)]), names_to = "taxaid", values_to = "count")

# Plot on genus level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!genus == "NA") %>%
  dplyr::group_by(genus, age) %>% dplyr::mutate(genus_count = sum(count)) %>% ungroup() %>% select(-count, -species) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(genus_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = genus_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, genus) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(genus) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>% 
  subset(max > 4.99) %>%
  arrange(family) %>% select(genus) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!genus == "NA") %>%
  dplyr::group_by(genus, age) %>% dplyr::mutate(genus_count = sum(count)) %>% ungroup() %>% select(-count, -species) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(genus_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = genus_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, genus) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(genus) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>% 
  subset(max > 4.99) %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(genus)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq2, "_", core2, "_stratplot_100precent_genus_old.png"), width = 60, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq2, "_", core2, "_stratplot_100precent_genus_new.png"), width = 60, height = 20, units = "cm", dpi = 300)

# Plot on family level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!family == "NA") %>%
  dplyr::group_by(family, age) %>% dplyr::mutate(family_count = sum(count)) %>% ungroup() %>% select(-count, -species, -genus) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(family_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = family_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, family) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(family) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>%
  subset(max > 4.99) %>%
  arrange(family) %>% select(family) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!family == "NA") %>%
  dplyr::group_by(family, age) %>% dplyr::mutate(family_count = sum(count)) %>% ungroup() %>% select(-count, -species, -genus) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(family_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = family_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, family) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(family) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>%
  subset(max > 4.99)  %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(family)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq2, "_", core2, "_stratplot_100precent_family_old.png"), width = 60, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq2, "_", core2, "_stratplot_100precent_family_new.png"), width = 60, height = 20, units = "cm", dpi = 300)

# Plot on species level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!species == "NA") %>%
  dplyr::group_by(species, age) %>% dplyr::mutate(species_count = sum(count)) %>% ungroup() %>% select(-count) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(species_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = species_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, species) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(species) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>%
  subset(max > 1) %>%
  arrange(family) %>% select(species) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!species == "NA") %>%
  dplyr::group_by(species, age) %>% dplyr::mutate(species_count = sum(count)) %>% ungroup() %>% select(-count) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(species_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = species_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, species) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(species) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>%
  subset(max > 1) %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(species)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq2, "_", core2, "_stratplot_100precent_species_old.png"), width = 100, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq2, "_", core2, "_stratplot_100precent_species_new.png"), width = 100, height = 20, units = "cm", dpi = 300)

######################
# Bilyakh
######################
#tl_bilyakh_raw <- read_delim(paste("cores/", seq3, "_", core3, "_file02.2_identitylevel_90_merged_replicates.csv", sep = ""), delim = ";", col_names = T)
tl_bilyakh_raw <- read_delim(paste("cores/", seq3, "_", core3, "_file02.2_identitylevel_90_merged_replicates_new.csv", sep = ""), delim = ";", col_names = T)
bilyakh_asv <- tl_bilyakh_raw %>% select(NUC_SEQ_arc, best_identity, best_family, best_genus, best_species, scientific_name)
tl_bilyakh <- tl_bilyakh_raw %>% mutate(taxaid = paste(NUC_SEQ_arc, best_identity, best_family, best_genus, best_species, scientific_name, sep = "_")) %>% select(-NUC_SEQ_arc, -best_identity, -best_family, -best_genus, -best_species, -scientific_name)
bilyakh_metadata <- read_excel("metadata/cores/ALRK-8_Bilyakh_agefile.xlsx", col_names = T)
bilyakh_samp <- paste(bilyakh_metadata$Extraction_number, bilyakh_metadata$Mean_collection_depth_round, bilyakh_metadata$Section_age_round, sep = "_")

tl_bilyakh_long <- tl_bilyakh %>% pivot_longer(!taxaid, names_to = "sampleid", values_to = "count") %>% filter(count > 0) %>% group_by(taxaid) %>% mutate(total = sum(count)) %>% ungroup() %>% mutate(sampleid = paste("bilyakh", sampleid, sep = "_"))
tl_bilyakh_long_sup10 <- tl_bilyakh_long %>% filter(total > 9) %>% arrange(total)

#####
# Check and plots
#####
tl_bilyakh_long_sup10 %>% filter(grepl("_1_",taxaid)) %>% filter(grepl("Rosace",taxaid)) %>% select(taxaid) %>% distinct()

wide <- tl_bilyakh_long_sup10 %>% filter(grepl("_1_",taxaid)) %>% pivot_wider(names_from = taxaid, values_from = count) # make the data from long to wide
wide[is.na(wide)] <- 0 
t_plot <- pivot_longer(wide, cols = names(wide[,-c(1,2)]), names_to = "taxaid", values_to = "count")

# Plot on genus level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!genus == "NA") %>%
  dplyr::group_by(genus, age) %>% dplyr::mutate(genus_count = sum(count)) %>% ungroup() %>% select(-count, -species) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(genus_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = genus_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, genus) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(genus) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 4.99) %>%
  arrange(family) %>% select(genus) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!genus == "NA") %>%
  dplyr::group_by(genus, age) %>% dplyr::mutate(genus_count = sum(count)) %>% ungroup() %>% select(-count, -species) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(genus_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = genus_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, genus) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(genus) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 4.99) %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(genus)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq3, "_", core3, "_stratplot_100precent_genus_old.png"), width = 60, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq3, "_", core3, "_stratplot_100precent_genus_new.png"), width = 60, height = 20, units = "cm", dpi = 300)

# Plot on family level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!family == "NA") %>%
  dplyr::group_by(family, age) %>% dplyr::mutate(family_count = sum(count)) %>% ungroup() %>% select(-count, -species, -genus) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(family_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = family_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, family) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(family) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 4.99) %>%
  arrange(family) %>% select(family) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!family == "NA") %>%
  dplyr::group_by(family, age) %>% dplyr::mutate(family_count = sum(count)) %>% ungroup() %>% select(-count, -species, -genus) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(family_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = family_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, family) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(family) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 4.99)  %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(family)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq3, "_", core3, "_stratplot_100precent_family_old.png"), width = 60, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq3, "_", core3, "_stratplot_100precent_family_new.png"), width = 60, height = 20, units = "cm", dpi = 300)

# Plot on species level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!species == "NA") %>%
  dplyr::group_by(species, age) %>% dplyr::mutate(species_count = sum(count)) %>% ungroup() %>% select(-count) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(species_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = species_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, species) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(species) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 1) %>%
  arrange(family) %>% select(species) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!species == "NA") %>%
  dplyr::group_by(species, age) %>% dplyr::mutate(species_count = sum(count)) %>% ungroup() %>% select(-count) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(species_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = species_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, species) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(species) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 1) %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(species)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq3, "_", core3, "_stratplot_100precent_species_old.png"), width = 100, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq3, "_", core3, "_stratplot_100precent_species_new.png"), width = 100, height = 20, units = "cm", dpi = 300)

######################
# Levinson Lessing
######################
#tl_lele_raw <- read_delim(paste("cores/", seq4, "_", core4, "_file02.2_identitylevel_90_merged_replicates.csv", sep = ""), delim = ";", col_names = T)
tl_lele_raw <- read_delim(paste("cores/", seq4, "_", core4, "_file02.2_identitylevel_90_merged_replicates_new.csv", sep = ""), delim = ";", col_names = T)
lele_asv <- tl_lele_raw %>% select(NUC_SEQ_arc, best_identity, best_family, best_genus, best_species, scientific_name)
tl_lele <- tl_lele_raw %>% mutate(taxaid = paste(NUC_SEQ_arc, best_identity, best_family, best_genus, best_species, scientific_name, sep = "_")) %>% select(-NUC_SEQ_arc, -best_identity, -best_family, -best_genus, -best_species, -scientific_name)
lele_metadata <- read_excel("metadata/cores/ALRK-11_LeLe_agefile.xlsx", col_names = T)
lele_samp <- paste(lele_metadata$Extraction_number, lele_metadata$Mean_collection_depth_round, lele_metadata$Section_age_round, sep = "_")

tl_lele_long <- tl_lele %>% pivot_longer(!taxaid, names_to = "sampleid", values_to = "count") %>% filter(count > 0) %>% subset(sampleid %in% all_of(lele_samp), drop = F) %>% group_by(taxaid) %>% mutate(total = sum(count)) %>% ungroup() %>% mutate(sampleid = paste("lele", sampleid, sep = "_")) %>% filter(!grepl("_NA", sampleid))
tl_lele_long_sup10 <- tl_lele_long %>% filter(total > 9) %>% arrange(total)

#####
# Check and plots
#####
tl_lele_long_sup10 %>% filter(grepl("_1_",taxaid)) %>% filter(grepl("Rosace",taxaid)) %>% select(taxaid) %>% distinct()

wide <- tl_lele_long_sup10 %>% filter(grepl("_1_",taxaid)) %>% pivot_wider(names_from = taxaid, values_from = count) # make the data from long to wide
wide[is.na(wide)] <- 0 
t_plot <- pivot_longer(wide, cols = names(wide[,-c(1,2)]), names_to = "taxaid", values_to = "count")

# Plot on genus level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!genus == "NA") %>%
  dplyr::group_by(genus, age) %>% dplyr::mutate(genus_count = sum(count)) %>% ungroup() %>% select(-count, -species) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(genus_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = genus_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, genus) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(genus) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 4.99) %>%
  arrange(family) %>% select(genus) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!genus == "NA") %>%
  dplyr::group_by(genus, age) %>% dplyr::mutate(genus_count = sum(count)) %>% ungroup() %>% select(-count, -species) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(genus_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = genus_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, genus) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(genus) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 4.99) %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(genus)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq4, "_", core4, "_stratplot_100precent_genus_old.png"), width = 60, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq4, "_", core4, "_stratplot_100precent_genus_new.png"), width = 60, height = 20, units = "cm", dpi = 300)

# Plot on family level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!family == "NA") %>%
  dplyr::group_by(family, age) %>% dplyr::mutate(family_count = sum(count)) %>% ungroup() %>% select(-count, -species, -genus) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(family_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = family_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, family) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(family) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 4.99) %>%
  arrange(family) %>% select(family) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!family == "NA") %>%
  dplyr::group_by(family, age) %>% dplyr::mutate(family_count = sum(count)) %>% ungroup() %>% select(-count, -species, -genus) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(family_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = family_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, family) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(family) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 4.99)  %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(family)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq4, "_", core4, "_stratplot_100precent_family_old.png"), width = 60, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq4, "_", core4, "_stratplot_100precent_family_new.png"), width = 60, height = 20, units = "cm", dpi = 300)

# Plot on species level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!species == "NA") %>%
  dplyr::group_by(species, age) %>% dplyr::mutate(species_count = sum(count)) %>% ungroup() %>% select(-count) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(species_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = species_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, species) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(species) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 1) %>%
  arrange(family) %>% select(species) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!species == "NA") %>%
  dplyr::group_by(species, age) %>% dplyr::mutate(species_count = sum(count)) %>% ungroup() %>% select(-count) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(species_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = species_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, species) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(species) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 1) %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(species)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq4, "_", core4, "_stratplot_100precent_species_old.png"), width = 100, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq4, "_", core4, "_stratplot_100precent_species_new.png"), width = 100, height = 20, units = "cm", dpi = 300)

######################
# Ilerney EN
######################
#tl_ilerneyEN_raw <- read_delim(paste("cores/", seq5, "_", core5, "_file02.2_identitylevel_90_merged_replicates.csv", sep = ""), delim = ";", col_names = T)
tl_ilerneyEN_raw <- read_delim(paste("cores/", seq5, "_", core5, "_file02.2_identitylevel_90_merged_replicates_new.csv", sep = ""), delim = ";", col_names = T)
ilerneyEN_asv <- tl_ilerneyEN_raw %>% select(NUC_SEQ_arc, best_identity, best_family, best_genus, best_species, scientific_name)
tl_ilerneyEN <- tl_ilerneyEN_raw %>% mutate(taxaid = paste(NUC_SEQ_arc, best_identity, best_family, best_genus, best_species, scientific_name, sep = "_")) %>% select(-NUC_SEQ_arc, -best_identity, -best_family, -best_genus, -best_species, -scientific_name)
ilerneyEN_metadata <- read_excel("metadata/cores/ALRK-10_Ilirney_agefile.xlsx", col_names = T)
ilerneyEN_samp <- paste(ilerneyEN_metadata$Extraction_number, ilerneyEN_metadata$Mean_collection_depth_round, ilerneyEN_metadata$Section_age_round, sep = "_")

tl_ilerneyEN_long <- tl_ilerneyEN %>% pivot_longer(!taxaid, names_to = "sampleid", values_to = "count") %>% filter(count > 0) %>% group_by(taxaid) %>% mutate(total = sum(count)) %>% ungroup() %>% mutate(sampleid = paste("ilerneyEN", sampleid, sep = "_"))
tl_ilerneyEN_long_sup10 <- tl_ilerneyEN_long %>% filter(total > 9) %>% arrange(total)

#####
# Check and plots
#####
tl_ilerneyEN_long_sup10 %>% filter(grepl("_1_",taxaid)) %>% filter(grepl("Rosace",taxaid)) %>% select(taxaid) %>% distinct()

wide <- tl_ilerneyEN_long_sup10 %>% filter(grepl("_1_",taxaid)) %>% pivot_wider(names_from = taxaid, values_from = count) # make the data from long to wide
wide[is.na(wide)] <- 0 
t_plot <- pivot_longer(wide, cols = names(wide[,-c(1,2)]), names_to = "taxaid", values_to = "count")

# Plot on genus level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!genus == "NA") %>%
  dplyr::group_by(genus, age) %>% dplyr::mutate(genus_count = sum(count)) %>% ungroup() %>% select(-count, -species) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(genus_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = genus_count/tot_sample*100) %>% arrange(desc(rel_abund)) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, genus) %>% dplyr::mutate(keep = as.numeric(sum(rel_abund))) %>% ungroup() %>% 
  dplyr::group_by(genus) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 4.99) %>%
  arrange(family) %>% select(genus) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!genus == "NA") %>%
  dplyr::group_by(genus, age) %>% dplyr::mutate(genus_count = sum(count)) %>% ungroup() %>% select(-count, -species) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(genus_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = genus_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, genus) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(genus) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>% mutate(max = ifelse(is.na(max), 0, max)) %>%
  subset(max > 4.99) %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(genus)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq5, "_", core5, "_stratplot_100precent_genus_old.png"), width = 60, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq5, "_", core5, "_stratplot_100precent_genus_new.png"), width = 60, height = 20, units = "cm", dpi = 300)

# Plot on family level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!family == "NA") %>%
  dplyr::group_by(family, age) %>% dplyr::mutate(family_count = sum(count)) %>% ungroup() %>% select(-count, -species, -genus) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(family_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = family_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, family) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(family) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 4.99) %>%
  arrange(family) %>% select(family) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!family == "NA") %>%
  dplyr::group_by(family, age) %>% dplyr::mutate(family_count = sum(count)) %>% ungroup() %>% select(-count, -species, -genus) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(family_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = family_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, family) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(family) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 4.99)  %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(family)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq5, "_", core5, "_stratplot_100precent_family_old.png"), width = 60, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq5, "_", core5, "_stratplot_100precent_family_new.png"), width = 60, height = 20, units = "cm", dpi = 300)

# Plot on species level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!species == "NA") %>%
  dplyr::group_by(species, age) %>% dplyr::mutate(species_count = sum(count)) %>% ungroup() %>% select(-count) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(species_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = species_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, species) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(species) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 1) %>%
  arrange(family) %>% select(species) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!species == "NA") %>%
  dplyr::group_by(species, age) %>% dplyr::mutate(species_count = sum(count)) %>% ungroup() %>% select(-count) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(species_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = species_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, species) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(species) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 1) %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(species)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq5, "_", core5, "_stratplot_100precent_species_old.png"), width = 100, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq5, "_", core5, "_stratplot_100precent_species_new.png"), width = 100, height = 20, units = "cm", dpi = 300)

######################
# E5
######################
#tl_E5_raw <- read_delim(paste("cores/", seq6, "_", core6, "_file02.2_identitylevel_90_merged_replicates.csv", sep = ""), delim = ";", col_names = T)
tl_E5_raw <- read_delim(paste("cores/", seq6, "_", core6, "_file02.2_identitylevel_90_merged_replicates_new.csv", sep = ""), delim = ";", col_names = T)
E5_asv <- tl_E5_raw %>% select(NUC_SEQ_arc, best_identity, best_family, best_genus, best_species, scientific_name)
tl_E5 <- tl_E5_raw %>% mutate(taxaid = paste(NUC_SEQ_arc, best_identity, best_family, best_genus, best_species, scientific_name, sep = "_")) %>% select(-NUC_SEQ_arc, -best_identity, -best_family, -best_genus, -best_species, -scientific_name)
E5_metadata <- read_excel("metadata/cores/ALRK-13_E5_agefile.xlsx", col_names = T)
E5_samp <- paste(E5_metadata$Extraction_number, E5_metadata$Mean_collection_depth_round, E5_metadata$Section_age_round, sep = "_")

tl_E5_long <- tl_E5 %>% pivot_longer(!taxaid, names_to = "sampleid", values_to = "count") %>% filter(count > 0) %>% group_by(taxaid) %>% mutate(total = sum(count)) %>% ungroup() %>% mutate(sampleid = paste("E5", sampleid, sep = "_"))
tl_E5_long_sup10 <- tl_E5_long %>% filter(total > 9) %>% arrange(total)

#####
# Check and plots
#####
tl_E5_long_sup10 %>% filter(grepl("_1_",taxaid)) %>% filter(grepl("Rosace",taxaid)) %>% select(taxaid) %>% distinct()

wide <- tl_E5_long_sup10 %>% filter(grepl("_1_",taxaid)) %>% pivot_wider(names_from = taxaid, values_from = count) # make the data from long to wide
wide[is.na(wide)] <- 0 
t_plot <- pivot_longer(wide, cols = names(wide[,-c(1,2)]), names_to = "taxaid", values_to = "count")

# Plot on genus level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!genus == "NA") %>%
  dplyr::group_by(genus, age) %>% dplyr::mutate(genus_count = sum(count)) %>% ungroup() %>% select(-count, -species) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(genus_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = genus_count/tot_sample*100) %>% arrange(desc(rel_abund)) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, genus) %>% dplyr::mutate(keep = as.numeric(sum(rel_abund))) %>% ungroup() %>% 
  dplyr::group_by(genus) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 4.99) %>%
  arrange(family) %>% select(genus) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!genus == "NA") %>%
  dplyr::group_by(genus, age) %>% dplyr::mutate(genus_count = sum(count)) %>% ungroup() %>% select(-count, -species) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(genus_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = genus_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, genus) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(genus) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>% mutate(max = ifelse(is.na(max), 0, max)) %>%
  subset(max > 4.99) %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(genus)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq6, "_", core6, "_stratplot_100precent_genus_old.png"), width = 60, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq6, "_", core6, "_stratplot_100precent_genus_new.png"), width = 60, height = 20, units = "cm", dpi = 300)

# Plot on family level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!family == "NA") %>%
  dplyr::group_by(family, age) %>% dplyr::mutate(family_count = sum(count)) %>% ungroup() %>% select(-count, -species, -genus) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(family_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = family_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, family) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(family) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 4.99) %>%
  arrange(family) %>% select(family) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!family == "NA") %>%
  dplyr::group_by(family, age) %>% dplyr::mutate(family_count = sum(count)) %>% ungroup() %>% select(-count, -species, -genus) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(family_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = family_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, family) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(family) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 4.99)  %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(family)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq6, "_", core6, "_stratplot_100precent_family_old.png"), width = 60, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq6, "_", core6, "_stratplot_100precent_family_new.png"), width = 60, height = 20, units = "cm", dpi = 300)

# Plot on species level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!species == "NA") %>%
  dplyr::group_by(species, age) %>% dplyr::mutate(species_count = sum(count)) %>% ungroup() %>% select(-count) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(species_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = species_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, species) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(species) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 1) %>%
  arrange(family) %>% select(species) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!species == "NA") %>%
  dplyr::group_by(species, age) %>% dplyr::mutate(species_count = sum(count)) %>% ungroup() %>% select(-count) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(species_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = species_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, species) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(species) %>% dplyr::mutate(max = max(keep)) %>%
  subset(max > 1) %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(species)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq6, "_", core6, "_stratplot_100precent_species_old.png"), width = 100, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq6, "_", core6, "_stratplot_100precent_species_new.png"), width = 100, height = 20, units = "cm", dpi = 300)

######################
# Rauchuagytgyn
######################
#tl_rauchuagytgyn_raw <- read_delim(paste("cores/", seq7, "_", core7, "_file02.2_identitylevel_90_merged_replicates.csv", sep = ""), delim = ";", col_names = T)
tl_rauchuagytgyn_raw <- read_delim(paste("cores/", seq7, "_", core7, "_file02.2_identitylevel_90_merged_replicates_new.csv", sep = ""), delim = ";", col_names = T)
rauchuagytgyn_asv <- tl_rauchuagytgyn_raw %>% select(NUC_SEQ_arc, best_identity, best_family, best_genus, best_species, scientific_name)
tl_rauchuagytgyn <- tl_rauchuagytgyn_raw %>% mutate(taxaid = paste(NUC_SEQ_arc, best_identity, best_family, best_genus, best_species, scientific_name, sep = "_")) %>% select(-NUC_SEQ_arc, -best_identity, -best_family, -best_genus, -best_species, -scientific_name)
rauchuagytgyn_metadata <- read_excel("metadata/cores/APMG-30_Rauch_agefile.xlsx", col_names = T)
rauchuagytgyn_samp <- paste(rauchuagytgyn_metadata$Extraction_number, rauchuagytgyn_metadata$Mean_collection_depth_round, rauchuagytgyn_metadata$Section_age_round, sep = "_")

tl_rauchuagytgyn_long <- tl_rauchuagytgyn %>% pivot_longer(!taxaid, names_to = "sampleid", values_to = "count") %>% filter(count > 0) %>% group_by(taxaid) %>% mutate(total = sum(count)) %>% ungroup() %>% mutate(sampleid = paste("rauchuagytgyn", sampleid, sep = "_"))
tl_rauchuagytgyn_long_sup10 <- tl_rauchuagytgyn_long %>% filter(total > 9) %>% arrange(total)

#####
# Check and plots
#####
tl_rauchuagytgyn_long_sup10 %>% filter(grepl("_1_",taxaid)) %>% filter(grepl("Rosace",taxaid)) %>% select(taxaid) %>% distinct()

wide <- tl_rauchuagytgyn_long_sup10 %>% filter(grepl("_1_",taxaid)) %>% pivot_wider(names_from = taxaid, values_from = count) # make the data from long to wide
wide[is.na(wide)] <- 0 
t_plot <- pivot_longer(wide, cols = names(wide[,-c(1,2)]), names_to = "taxaid", values_to = "count")

# Plot on genus level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!genus == "NA") %>%
  dplyr::group_by(genus, age) %>% dplyr::mutate(genus_count = sum(count)) %>% ungroup() %>% select(-count, -species) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(genus_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = genus_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, genus) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(genus) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>% 
  subset(max > 4.99) %>%
  arrange(family) %>% select(genus) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!genus == "NA") %>%
  dplyr::group_by(genus, age) %>% dplyr::mutate(genus_count = sum(count)) %>% ungroup() %>% select(-count, -species) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(genus_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = genus_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, genus) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(genus) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>% 
  subset(max > 4.99) %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(genus)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq7, "_", core7, "_stratplot_100precent_genus_old.png"), width = 60, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq7, "_", core7, "_stratplot_100precent_genus_new.png"), width = 60, height = 20, units = "cm", dpi = 300)

# Plot on family level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!family == "NA") %>%
  dplyr::group_by(family, age) %>% dplyr::mutate(family_count = sum(count)) %>% ungroup() %>% select(-count, -species, -genus) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(family_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = family_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, family) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(family) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>%
  subset(max > 4.99) %>%
  arrange(family) %>% select(family) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!family == "NA") %>%
  dplyr::group_by(family, age) %>% dplyr::mutate(family_count = sum(count)) %>% ungroup() %>% select(-count, -species, -genus) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(family_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = family_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, family) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(family) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>%
  subset(max > 4.99)  %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(family)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq7, "_", core7, "_stratplot_100precent_family_old.png"), width = 60, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq7, "_", core7, "_stratplot_100precent_family_new.png"), width = 60, height = 20, units = "cm", dpi = 300)

# Plot on species level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!species == "NA") %>%
  dplyr::group_by(species, age) %>% dplyr::mutate(species_count = sum(count)) %>% ungroup() %>% select(-count) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(species_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = species_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, species) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(species) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>%
  subset(max > 1) %>%
  arrange(family) %>% select(species) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!species == "NA") %>%
  dplyr::group_by(species, age) %>% dplyr::mutate(species_count = sum(count)) %>% ungroup() %>% select(-count) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(species_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = species_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, species) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(species) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>%
  subset(max > 1) %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(species)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq7, "_", core7, "_stratplot_100precent_species_old.png"), width = 100, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq7, "_", core7, "_stratplot_100precent_species_new.png"), width = 100, height = 20, units = "cm", dpi = 300)

######################
# Emanda
######################
#tl_emanda_raw <- read_delim(paste("cores/", seq8, "_", core8, "_file02.2_identitylevel_90_merged_replicates.csv", sep = ""), delim = ";", col_names = T)
tl_emanda_raw <- read_delim(paste("cores/", seq8, "_", core8, "_file02.2_identitylevel_90_merged_replicates_new.csv", sep = ""), delim = ";", col_names = T)
emanda_asv <- tl_emanda_raw %>% select(NUC_SEQ_arc, best_identity, best_family, best_genus, best_species, scientific_name)
tl_emanda <- tl_emanda_raw %>% mutate(taxaid = paste(NUC_SEQ_arc, best_identity, best_family, best_genus, best_species, scientific_name, sep = "_")) %>% select(-NUC_SEQ_arc, -best_identity, -best_family, -best_genus, -best_species, -scientific_name)
emanda_metadata <- read_excel("metadata/cores/APMG-42_Emanda_agefile.xlsx", col_names = T)
emanda_samp <- paste(emanda_metadata$Extraction_number, emanda_metadata$Mean_collection_depth_round, emanda_metadata$Section_age_round, sep = "_")

tl_emanda_long <- tl_emanda %>% pivot_longer(!taxaid, names_to = "sampleid", values_to = "count") %>% filter(count > 0) %>% group_by(taxaid) %>% mutate(total = sum(count)) %>% ungroup() %>% mutate(sampleid = paste("emanda", sampleid, sep = "_"))
tl_emanda_long_sup10 <- tl_emanda_long %>% filter(total > 9) %>% arrange(total)

#####
# Check and plots
#####
tl_emanda_long_sup10 %>% filter(grepl("_1_",taxaid)) %>% filter(grepl("Rosace",taxaid)) %>% select(taxaid) %>% distinct()

wide <- tl_emanda_long_sup10 %>% filter(grepl("_1_",taxaid)) %>% pivot_wider(names_from = taxaid, values_from = count) # make the data from long to wide
wide[is.na(wide)] <- 0 
t_plot <- pivot_longer(wide, cols = names(wide[,-c(1,2)]), names_to = "taxaid", values_to = "count")

# Plot on genus level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!genus == "NA") %>%
  dplyr::group_by(genus, age) %>% dplyr::mutate(genus_count = sum(count)) %>% ungroup() %>% select(-count, -species) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(genus_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = genus_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, genus) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(genus) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>% 
  subset(max > 4.99) %>%
  arrange(family) %>% select(genus) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!genus == "NA") %>%
  dplyr::group_by(genus, age) %>% dplyr::mutate(genus_count = sum(count)) %>% ungroup() %>% select(-count, -species) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(genus_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = genus_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, genus) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(genus) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>% 
  subset(max > 4.99) %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(genus)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq8, "_", core8, "_stratplot_100precent_genus_old.png"), width = 60, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq8, "_", core8, "_stratplot_100precent_genus_new.png"), width = 60, height = 20, units = "cm", dpi = 300)

# Plot on family level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!family == "NA") %>%
  dplyr::group_by(family, age) %>% dplyr::mutate(family_count = sum(count)) %>% ungroup() %>% select(-count, -species, -genus) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(family_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = family_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, family) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(family) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>%
  subset(max > 4.99) %>%
  arrange(family) %>% select(family) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!family == "NA") %>%
  dplyr::group_by(family, age) %>% dplyr::mutate(family_count = sum(count)) %>% ungroup() %>% select(-count, -species, -genus) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(family_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = family_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, family) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(family) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>%
  subset(max > 4.99)  %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(family)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq8, "_", core8, "_stratplot_100precent_family_old.png"), width = 60, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq8, "_", core8, "_stratplot_100precent_family_new.png"), width = 60, height = 20, units = "cm", dpi = 300)

# Plot on species level (100%)
reorder <- t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!species == "NA") %>%
  dplyr::group_by(species, age) %>% dplyr::mutate(species_count = sum(count)) %>% ungroup() %>% select(-count) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(species_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = species_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, species) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(species) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>%
  subset(max > 1) %>%
  arrange(family) %>% select(species) %>% distinct()

t_plot %>% 
  separate(taxaid, into = c("seq", "id", "family", "genus", "species"), sep = "_") %>% filter(id == "1") %>% select(-seq, -id) %>% 
  separate(sampleid, into = c("core", "sample", "depth", "age")) %>% filter(!is.na(age)) %>% select(-core, -sample, -depth, -total) %>% distinct() %>%
  dplyr::filter(!species == "NA") %>%
  dplyr::group_by(species, age) %>% dplyr::mutate(species_count = sum(count)) %>% ungroup() %>% select(-count) %>% distinct() %>%
  dplyr::group_by(age) %>% dplyr::mutate(tot_sample = sum(species_count)) %>% ungroup() %>%
  dplyr::mutate(rel_abund = species_count/tot_sample*100) %>%
  dplyr::mutate(age = as.numeric(age), rel_abund = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  dplyr::group_by(age, species) %>% dplyr::mutate(keep = sum(rel_abund)) %>% ungroup() %>% 
  dplyr::group_by(species) %>% dplyr::mutate(max = max(keep)) %>% ungroup() %>%
  subset(max > 1) %>%
  ggplot(aes(x = rel_abund, y = age)) +
  geom_areah(aes(fill = family)) +
  geom_col_segsh() + 
  #geom_lineh() +
  scale_y_reverse() +
  facet_abundanceh(vars(species)) +
  labs(x = "Relative abundance (%)", y = "Age (kyrs BP)")
#ggsave(paste0("cores/", seq8, "_", core8, "_stratplot_100precent_species_old.png"), width = 100, height = 20, units = "cm", dpi = 300)
ggsave(paste0("cores/", seq8, "_", core8, "_stratplot_100precent_species_new.png"), width = 100, height = 20, units = "cm", dpi = 300)

############################################################################
# 3 - Merge all cores
############################################################################
tl_merge_long <- full_join(tl_bilyakh_long, tl_btoko_long)
tl_merge_long <- full_join(tl_merge_long, tl_E5_long)
tl_merge_long <- full_join(tl_merge_long, tl_emanda_long)
tl_merge_long <- full_join(tl_merge_long, tl_ilerney16kp_long)
tl_merge_long <- full_join(tl_merge_long, tl_ilerneyEN_long)
tl_merge_long <- full_join(tl_merge_long, tl_lele_long)
tl_merge_long <- full_join(tl_merge_long, tl_rauchuagytgyn_long)

tl_merge_long_raw <- tl_merge_long

# Separate sample_id into metadata info
tl_merge_long <- separate(data=tl_merge_long_raw,
                          col = "sampleid",
                          into = c("site", "EXT_no", "depth", "age"),
                          sep = "_",
                          remove = FALSE,
                          extra = "drop")

# Check if all samples have an age
tl_merge_long %>% select(site, depth, age) %>% distinct() %>% arrange(desc(as.numeric(age)))
tl_merge_long %>% select(site, depth, age) %>% distinct() %>% arrange(as.numeric(age))
tl_merge_long %>% select(site, depth, age) %>% distinct() %>% filter(age == "NA")

# Measure some metrics: total reads count per raxa, number of sample occurence, number of cores occurence and total number of reads per sample
tl_merge_long1 <- tl_merge_long %>% filter(!is.na(age)) %>%
  mutate(sample_id = paste(site, age, sep = "_")) %>% 
  group_by(taxaid) %>% mutate(reads_count = sum(count), nsamples = n_distinct(sampleid), ncores = n_distinct(site)) %>% ungroup() %>% 
  group_by(taxaid, sample_id) %>% mutate(reads_age = sum(count)) %>% ungroup()

tl_merge_long1 %>% select(taxaid, reads_count) %>% distinct()
unique(tl_merge_long1$sampleid)

############################################################################
# 4 - Get only ASVs present with a minimum of 9 reads / in 9 samples
############################################################################
tl_merge_clean <- tl_merge_long1 %>% filter(reads_count > 9) %>% filter(nsamples > 9) %>% arrange(age) %>% select(taxaid, sample_id, reads_age, nsamples, ncores, reads_count) %>% distinct() 
tl_clean_wide <- tl_merge_clean %>% pivot_wider(names_from = sample_id, values_from = reads_age) #set as wide dataset
tl_clean_wide[is.na(tl_clean_wide)] <- 0 # replace NAs per 0
tl_clean_wide <- tl_clean_wide %>% mutate(id = 1:nrow(tl_clean_wide)) # set an id for each row of the data
tl_clean_wide <- separate(data=tl_clean_wide,
                          col = "taxaid",
                          into = c("NUC_SEQ", "identity", "family", "genus", "species", "scientific_name"),
                          sep = "_",
                          remove = FALSE,
                          extra = "drop") 
tl_clean_wide <- tl_clean_wide %>% mutate(uniq_id = paste(id, scientific_name, identity, sep = "_")) # set a uniq taxid for each ASV
tl_clean_wide$identity <- str_replace_all(tl_clean_wide$identity, ",", ".") %>% as.numeric() # change indeitty values to numeric

############################################################################
# 5 - Save the data
############################################################################
write_delim(tl_merge_clean, file = "output/merged_sites_clean_10-reads_10-samples_long.csv", delim = ",")
write_delim(tl_clean_wide, file = "output/merged_sites_clean_10-reads_10-samples_wide.csv", delim = ",")
