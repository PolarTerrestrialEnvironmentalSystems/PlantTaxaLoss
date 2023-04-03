# TAXALOSS Project - script to see if the candidate taxa are present in modern databases
# Script from Jeremy Courtin
# Script last update - 03.04.2023

# Layout of the script: 
# 0 - GENERAL OPTION - SETTING UP SCRIPT
# 1 - Assessment if the candidates lost are absent from EMBL, arctic and phylonorway databases as well
# 1.1 - EMBL
# 1.2 - Arctic database
# 1.3 - Phylonorway database

###############################################################################
# 0 - GENERAL OPTION - SETTING UP SCRIPT
###############################################################################
rm(list = ls()) # Remove all the objects we created so far.

# Set up options for the rest of the script
options(stringsAsFactors=FALSE) # state that everything is reads as character 

# Load all needed packages
library(tidyverse)
library(Biostrings)

# Import the resampled data
cand_ASVs <- read_delim("Intermediate_results/ASVs_taxa_after_louvain_community.csv", delim = ",", col_names = T)
median_read_resampled <- read_delim("median_reads_resampled_1000_final.csv", delim = ",", col_names = T)

################################################################################
# 1 - Assessment if the candidates lost are absent from EMBL, arctic and phylonorway databases as well
################################################################################
################################################################################
# 1.1 - EMBL
################################################################################
fastaFile = readDNAStringSet("check_against_databases/gh_embl143_db_97.fasta")
seq_name = names(fastaFile)
NUC_SEQ = paste(fastaFile)
EMBL <- tibble(NUC_SEQ, seq_name)

cand_absent_EMBL <- median_read_resampled %>% dplyr::select(uniq_label, `13`) %>% 
  separate(col = "uniq_label", into = c("g", "s", "t", "f", "n"), sep = "_", remove = F, extra = "drop") %>%
  mutate(uniq_label = paste(g, s, f, n, sep = "_")) %>% subset(`13` == 0 & s == "cand") %>% dplyr::select(uniq_label) %>% 
  left_join(dplyr::select(cand_ASVs, uniq_label, identity, NUC_SEQ), by = "uniq_label") %>% 
  mutate(NUC_SEQ = toupper(NUC_SEQ)) %>% left_join(EMBL) %>% arrange(uniq_label) %>% mutate(db = "embl")
write_delim(cand_absent_EMBL, "check_against_databases/cand_absent_EMBL.csv", delim = ",", col_names = T)

n_nucseq_embl <- cand_absent_EMBL %>% group_by(uniq_label) %>% summarise(n = n_distinct(NUC_SEQ))
write_delim(n_nucseq_embl, "check_against_databases/n_nuc_seq_EMBL.csv", delim = ",", col_names = T)

cand_absent_EMBL_2 <- median_read_resampled %>% dplyr::select(uniq_label, `13`) %>% 
  separate(col = "uniq_label", into = c("g", "s", "t", "f", "n"), sep = "_", remove = F, extra = "drop") %>%
  mutate(uniq_label = paste(g, s, f, n, sep = "_")) %>% subset(s == "cand") %>% dplyr::select(uniq_label) %>% 
  left_join(dplyr::select(cand_ASVs, uniq_label, identity, NUC_SEQ), by = "uniq_label") %>% 
  mutate(NUC_SEQ = toupper(NUC_SEQ)) %>% left_join(EMBL) %>% filter(!is.na(seq_name)) %>% arrange(uniq_label) %>% mutate(db = "embl")
write_delim(cand_absent_EMBL_2, "check_against_databases/cand_EMBL.csv", delim = ",", col_names = T)

################################################################################
# 1.2 - Arctic database
################################################################################
fastaFile = readDNAStringSet("check_against_databases/arctborbryo_gh.fasta")
seq_name = names(fastaFile)
NUC_SEQ = paste(fastaFile)
arc <- tibble(NUC_SEQ, seq_name)

cand_absent_arc <- median_read_resampled %>% dplyr::select(uniq_label, `13`) %>% 
  separate(col = "uniq_label", into = c("g", "s", "t", "f", "n"), sep = "_", remove = F, extra = "drop") %>%
  mutate(uniq_label = paste(g, s, f, n, sep = "_")) %>% subset(`13` == 0 & s == "cand") %>% dplyr::select(uniq_label) %>% 
  left_join(dplyr::select(cand_ASVs, uniq_label, identity, NUC_SEQ), by = "uniq_label") %>% 
  mutate(NUC_SEQ = toupper(NUC_SEQ)) %>% left_join(arc) %>% arrange(uniq_label) %>% mutate(db = "arc")
write_delim(cand_absent_arc, "check_against_databases/cand_absent_arcborbryo.csv", delim = ",", col_names = T)

n_nucseq_arc <- cand_absent_arc %>% group_by(uniq_label) %>% summarise(n = n_distinct(NUC_SEQ))
write_delim(n_nucseq_arc, "check_against_databases/n_nuc_seq_arcborbryo.csv", delim = ",", col_names = T)

cand_absent_arc_2 <- median_read_resampled %>% dplyr::select(uniq_label, `13`) %>% 
  separate(col = "uniq_label", into = c("g", "s", "t", "f", "n"), sep = "_", remove = F, extra = "drop") %>%
  mutate(uniq_label = paste(g, s, f, n, sep = "_")) %>% subset(s == "cand") %>% dplyr::select(uniq_label) %>% 
  left_join(dplyr::select(cand_ASVs, uniq_label, identity, NUC_SEQ), by = "uniq_label") %>% 
  mutate(NUC_SEQ = toupper(NUC_SEQ)) %>% left_join(arc) %>% filter(!is.na(seq_name)) %>% arrange(uniq_label) %>% mutate(db = "arc")
write_delim(cand_absent_arc_2, "check_against_databases/cand_arcborbry.csv", delim = ",", col_names = T)

################################################################################
# 1.3 - Phylonorway database
################################################################################
fastaFile = readDNAStringSet("check_against_databases/PhyloNorway_GH_database.fasta")
seq_name = names(fastaFile)
NUC_SEQ = paste(fastaFile)
phylo <- tibble(NUC_SEQ, seq_name)

cand_absent_phylo <- median_read_resampled %>% dplyr::select(uniq_label, `13`) %>% 
  separate(col = "uniq_label", into = c("g", "s", "t", "f", "n"), sep = "_", remove = F, extra = "drop") %>%
  mutate(uniq_label = paste(g, s, f, n, sep = "_")) %>% subset(`13` == 0 & s == "cand") %>% dplyr::select(uniq_label) %>% 
  left_join(dplyr::select(cand_ASVs, uniq_label, identity, NUC_SEQ), by = "uniq_label") %>% 
  mutate(NUC_SEQ = toupper(NUC_SEQ)) %>% left_join(phylo) %>% arrange(seq_name) %>% dplyr::mutate(db = "phylo")
write_delim(cand_absent_phylo, "check_against_databases/cand_absent_phylo.csv", delim = ",", col_names = T)

n_nucseq_phylo <- cand_absent_phylo %>% group_by(uniq_label) %>% summarise(n = n_distinct(NUC_SEQ))
write_delim(n_nucseq_phylo, "check_against_databases/n_nuc_seq_phylo.csv", delim = ",", col_names = T)

cand_absent_phylo_2 <- median_read_resampled %>% dplyr::select(uniq_label, `13`) %>% 
  separate(col = "uniq_label", into = c("g", "s", "t", "f", "n"), sep = "_", remove = F, extra = "drop") %>%
  mutate(uniq_label = paste(g, s, f, n, sep = "_")) %>% subset(s == "cand") %>% dplyr::select(uniq_label) %>% 
  left_join(dplyr::select(cand_ASVs, uniq_label, identity, NUC_SEQ), by = "uniq_label") %>% 
  mutate(NUC_SEQ = toupper(NUC_SEQ)) %>% left_join(phylo) %>% filter(!is.na(seq_name)) %>% arrange(uniq_label) %>% dplyr::mutate(db = "phylo")
write_delim(cand_absent_phylo_2, "check_against_databases/cand_phylo.csv", delim = ",", col_names = T)

cand_all <- rbind(cand_absent_EMBL_2, cand_absent_arc_2, cand_absent_phylo_2) %>% arrange(uniq_label)
cand_all %>% select(uniq_label) %>% distinct()
cand_all %>% select(uniq_label, NUC_SEQ) %>% distinct()
cand_absent_all <- rbind(cand_absent_EMBL, cand_absent_arc, cand_absent_phylo) %>% arrange(seq_name)
cand_absent_all %>% filter(!is.na(seq_name)) %>% select(uniq_label) %>% distinct()
cand_absent_all %>% filter(!is.na(seq_name)) %>% select(uniq_label, NUC_SEQ) %>% distinct()

write_delim(cand_all, "check_against_databases/candidate_alldatabase_check.csv", delim = ",", col_names = T)
write_delim(cand_absent_all, "check_against_databases/candidate_absent_alldatabase_check.csv", delim = ",", col_names = T)

