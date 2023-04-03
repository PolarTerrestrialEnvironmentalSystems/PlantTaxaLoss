# TAXALOSS Project - transform ASVs info to taxa info
# Co-occurence / Network / Community
# Script from Jeremy Courtin
# Script last update - 03.04.2023

# Layout of the script: 
# 0 - GENERAL OPTION - SETTING UP SCRIPT
# 01 - Import your raw data
# 1 - Set parameters and save general information
# 2 - Co-occurence and community building -> TCorr.test and iGraph
# 3 - Prepare for network analysis
# 4 - Keep only the positive and sup to corr_score correlation scores
# 5 - Make network analysis using  graph_from_data_frame
# 6 - Detect communities with the cluster label prop method
# 7 - Create a dataframe with membership info
# /!\ Check for homopolymers within detected communities /!\
# 7.1 - Add number of ASV per community as metadata
# 8 - Filter the ASVs into final potential taxa
# 8.1 - Merge ASVs with same taxa names - genus level (unused in manuscript)
# 8.2 - Merge ASVs with same taxa names - family level (used in manuscript)

############################################################################
# 0 - GENERAL OPTION - SETTING UP SCRIPT
############################################################################
rm(list = ls()) # Remove all the objects we created so far.

# Set up options for the rest of the script
options(stringsAsFactors=FALSE) # states that everything is reads as character 

# Load all needed packages
library(tidyverse)
library(vegan)
library(igraph)
library(psych)
library(picante)

############################################################################
# 01 - Import your raw data
############################################################################
df_raw <- read_delim("output/merged_sites_clean_10-reads_10-samples_wide.csv", delim = ",", col_names = T) 

# get the total number of reads we handle
sum(df_raw$reads_count)

# add uniq - short labels
df <- df_raw %>% mutate(short_label = paste(id, identity, scientific_name, sep = "_"))

# Check
df
df %>% subset(identity == 1) # 367 ASVs assigned at 100%
df %>% subset(reads_count > 99) # 5421 ASVs with minimum 100 reads
df %>% subset(reads_count > 99) %>% subset(identity == 1) # 332 ASVs assigned at 100% and with minimum 100 reads
df %>% subset(reads_count > 99) %>% subset(identity == 1) %>% select(family, genus, species, scientific_name) %>% distinct()#

# make a subset with 100% identity ASVS
df_NULL <- df %>% subset(identity == 1)
NULL_fam <- data.frame(table(df_NULL$family)) %>% arrange(desc(Freq)) # get the number of database ASVs per family

# make a subset with candidates ASVS (identity <100%)
df_candidates <- df %>% subset(!identity == 1)
candidates_fam <- data.frame(table(df_candidates$family)) %>% arrange(desc(Freq)) # get the number of candidate ASVs per family

############################################################################
# 1 - Set parameters and save general information
############################################################################


names(df)
# create metadata taxa
metadata_taxa <- df %>% select(uniq_id, short_label, taxaid, NUC_SEQ, identity, family, genus, species, scientific_name, nsamples, ncores, reads_count)

# subset for minimum 100 reads
df_100 <- df %>% filter(reads_count >= 100) %>% arrange(reads_count) %>% select(-id, -uniq_id, -taxaid, -NUC_SEQ, -identity, -family, -genus, -species, -scientific_name, -nsamples, -ncores, -reads_count)

# Put taxa_id as row_names
df_rdy_com <- as.data.frame(df_100)
rownames(df_rdy_com) <- df_rdy_com$short_label
df_rdy_com$short_label <- NULL

# get the homopolymers
homopol <- df %>% subset(identity == 1 & reads_count > 99) %>% subset(!is.na(species)) %>% select(family, genus, species, scientific_name) %>% arrange(species) %>% group_by(species) %>% mutate(n_occur = length(species)) %>% subset(n_occur > 1) %>% select(scientific_name) %>% distinct()
# df <- df_NULL # Check for homopolymers - OPTIONAL

############################################################################
# 2 - Co-occurence and community building -> TCorr.test and iGraph
############################################################################
# Here, I do not resample but apply the fourth root transformation:
df_trans <- sqrt(sqrt(sqrt(sqrt(df_rdy_com))))
save(df_trans, file = "output/df_for_co-occurence_test.rda")

# Do correlation of occurence matrix - very time consuming
corr_df <- corr.test(t(df_trans), # a matrix or dataframe
                      use = "pairwise", # pairwise" is the default value and will do pairwise deletion of cases. use="complete" will select just complete cases
                      method="spearman", # default is Pearson. spearman is slower especially for large datasets
                      adjust="holm", # What adjustment for multiple tests should be used? ("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none").
                      alpha=.05, # alpha level of confidence intervals
                      ci=FALSE, # default is TRUE but set false. For very large matrices (>200x200), noticeable speed improvement if confidence interval are not found
                      minlength=5, # What is the minimum length for abbreviations. Defaults to 5.
                      normal=F) # normal = FALSE -> we do not have normal distribution. But it is much slower to process

# Save just in case of crash
save(corr_df, file = "output/Co-occurence_raw_results.rda")

############################################################################
# 3 - Prepare for network analysis
############################################################################
# Make p-values as long
corr_p <- as.data.frame(corr_df$p)
corr_p$sp1 <- rownames(corr_p)
corr_p_long <- gather(corr_p, sp2, prob, all_of(colnames(as.data.frame(corr_df$p))))
#as_tibble(corr_p_long)

# Make uniq_ids
corr_p_long$uniqid <- paste(corr_p_long$sp1, corr_p_long$sp2, sep = "_")
corr_p_long <- corr_p_long %>% select(uniqid, prob)

# Re-format correlation matrix and merge p-value
corr_r <- as.data.frame(corr_df$r)
corr_r$sp1 <- rownames(corr_r)
corr_r_long <- gather(corr_r, sp2, corr, all_of(colnames(as.data.frame(corr_df$r))))
#as_tibble(corr_r_long)
corr_r_long$uniqid <- paste(corr_r_long$sp1, corr_r_long$sp2, sep = "_")

corr_net <- left_join(corr_r_long, corr_p_long, by = "uniqid")
#as_tibble(corr_net)
corr_net <- corr_net %>% select(sp1, sp2, corr, prob)

# remove corr = 1 (correlation to same asv)
corr_net <- corr_net %>% subset(!corr == 1) 
#as_tibble(corr_net)
#unique(corr_net$sp1)

# save nodes information
nodes <- data.frame(label = rownames(df_trans))
nodes <- left_join(nodes, metadata_taxa, by = c("label" = "short_label"))
#as_tibble(nodes)

############################################################################
# 4 - Keep only the positive and sup to corr_score correlation scores
############################################################################
corr_score = 0.4 # minimum correlation score we use in the manuscript
#corr_score = 0.3
corr_net_pos <- corr_net %>% subset(corr > 0) %>% subset(corr >= corr_score)

#Save the correlation results (just in case and time saving)
links_pos<- data.frame(from = corr_net_pos$sp1, 
                       to = corr_net_pos$sp2, 
                       corr = corr_net_pos$corr)
write_delim(links_pos, "output/links_pos_new.csv", delim = ",", col_names = TRUE)
write_delim(nodes, "output/nodes_new.csv", delim = ",", col_names = TRUE)

#Re-read (avoid re-running the calculation of co-occurrence)
links_pos <- read_delim("output/links_pos_new.csv", delim = ",", col_names = TRUE) %>% as.data.frame()
nodes <- read_delim("output/nodes_new.csv", delim = ",", col_names = TRUE) %>% as.data.frame()

############################################################################
# 5 - Make network analysis using  graph_from_data_frame 
############################################################################
net_pos <- graph_from_data_frame(d=links_pos, vertices=nodes, directed=F)
#plot(net_pos)

# change nodes info
V(net_pos)$label <- paste(nodes$family, strtrim(nodes$identity, 4), sep = ".")

############################################################################
# 6 - Detect communities with the cluster louvain method to optimise modularity
############################################################################
#ceb <- cluster_edge_betweenness(net_pos) # ideally would use ceb from replicability but cannot run (matrix too big)
#clp <- cluster_label_prop(net_pos) # use clp instead but will create new communities everytime not the most robust approach but ceb is too cpu demanding
cl <- cluster_louvain(net_pos) # alternative used in the manuscript

# check number of communities
#length(ceb)
#length(clp)
length(cl) # 275 communities detected in manuscript

#check membership
#membership(ceb)
#membership(clp)
membership(cl)

############################################################################
# 7 - Create a dataframe with membership info
############################################################################
#memberceb <- data.frame(group = ceb$membership, label = ceb$names)
#memberclp <- data.frame(group = clp$membership, label = clp$names)
membercl <- data.frame(group = cl$membership, label = cl$names)

############################################################################
# 7.1 - Add number of ASV per community as metadata
############################################################################
#memberceb <- memberceb %>% group_by(group) %>% mutate(nb_in_group = n_distinct(label))
#memberclp <- memberclp %>% group_by(group) %>% mutate(nb_in_group = n_distinct(label))
membercl <- membercl %>% group_by(group) %>% mutate(nb_in_group = n_distinct(label))

# created dataset
#memberceb <- left_join(memberceb, taxa_metadata, by = c("label" = "short_label"))
#memberclp <- left_join(memberclp, metadata_taxa, by = c("label" = "short_label"))
membercl <- left_join(membercl, metadata_taxa, by = c("label" = "short_label"))

############################################################################
# /!\ Check for homopolymers within detected communities /!\
############################################################################
# Check if all ASVs assigned at 100% to same taxon name (homopolymers) are all present in the same communties -> indicate that communities make sense somehow
homo_member <- membercl %>% subset(species %in% homopol$species & identity == 1) %>% arrange(scientific_name) %>% select(group, species)
write_delim(homo_member, "homopolymers_100percent_communities_louvain_only100%.csv", delim = ",", col_names = T)

# save the raw communities tables
# write_delim(memberclp, paste("output/Communities_clp_corrsup", corr_score, "more_100_reads_new.csv", sep = "_"), 
#             delim = ",", col_names = TRUE)
# write_delim(df, paste("output/DF_Communities_clp_corrsup", corr_score, "more_100_reads_new.csv", sep = "_"), 
#             delim = ",", col_names = TRUE)

write_delim(membercl, paste("output/Communities_louvain_corrsup", corr_score, "more_100_reads_new.csv", sep = "_"),
            delim = ",", col_names = TRUE)
write_delim(df, paste("output/DF_Communities_louvain_corrsup", corr_score, "more_100_reads_new.csv", sep = "_"),
            delim = ",", col_names = TRUE)

# Re-read again, avoid re-running the community building process
memberclp <- read_delim("output/Communities_louvain_corrsup_0.4_more_100_reads_new.csv",
            delim = ",", col_names = TRUE)

memberclp %>% filter(identity == 1) %>% select(family, genus, species, scientific_name) %>% distinct() %>%
  filter(!is.na(species))


df <- read_delim("output/DF_Communities_louvain_corrsup_0.4_more_100_reads_new.csv", 
                         delim = ",", col_names = TRUE)

############################################################################
# 8 - Filter the ASVs into final potential taxa
############################################################################
# keep reads  > 100
df_join <- df %>% filter(reads_count >= 100) %>% arrange(reads_count) %>% select(-id, -taxaid, -NUC_SEQ, -identity, -family, -genus, -species, -scientific_name, -nsamples, -ncores, -reads_count)

# Get only ASVs part of communities with minimum 5 ASVs
member <- left_join(memberclp, df_join, by = "uniq_id") %>% ungroup()
member_sup5_group <- member %>% filter(nb_in_group >= 5) %>% arrange(group)

member_sequence_sup5 <- member_sup5_group %>% select(group, uniq_id, NUC_SEQ, identity, family, genus, species, scientific_name) %>% mutate(status = ifelse(identity == 1, 1, "cand")) %>% mutate(uniq_label = paste(group, status, family, scientific_name, sep = "_"))
member_sequence_sup5 %>% filter(identity == 1)

suptable1 <- member_sup5_group %>% mutate(type = ifelse(identity == 1, "dbASV", "candidate ASV")) %>%
  select(NUC_SEQ, type, identity, family, genus, species, scientific_name, nsamples, ncores, reads_count) %>% arrange(family) %>% arrange(desc(type))

# Save for supplementary information
write_delim(member_sequence_sup5, "output/ASVs_taxa_after_louvain_community.csv", delim = ",", col_names = T)
write_delim(suptable1, "output/Supplementary_table1.csv", delim = ",", col_names = T)

samples_name <- member %>% 
  select(-group, -label, -nb_in_group, -uniq_id, -taxaid, -NUC_SEQ, -identity, -family, -genus, 
         -species, -scientific_name, -nsamples, -ncores, -reads_count, -short_label) %>%
  names()

write_delim(member_sup5_group, "output/Supplementary_info_after_louvain_community_detection.csv", delim = ",", col_names = T)
unique(member_sup5_group$group)

############################################################################
# 8.1 - Merge ASVs with same taxa names - genus level (unused in manuscript)
############################################################################
merge_assignments_prep <- member_sup5_group %>% select(group, nb_in_group, identity, family, genus, species, scientific_name, all_of(samples_name)) %>%
  pivot_longer(cols = all_of(samples_name), names_to = "sample_id", values_to = "reads_count")

merge_assignments_100percent <- merge_assignments_prep %>% filter(identity == 1 & reads_count > 0) %>% arrange(sample_id) %>%
  group_by(group, scientific_name, sample_id) %>% mutate(reads_count = sum(reads_count)) %>% ungroup() %>% distinct() 

merge_assignments_90percent <- merge_assignments_prep %>% filter(identity < 1 & reads_count > 0) %>% arrange(sample_id) %>% mutate(identity = 0.9) %>%
  group_by(group, scientific_name, sample_id) %>% mutate(reads_count = sum(reads_count)) %>% ungroup() %>% distinct() 

merge_assignments <- full_join(merge_assignments_100percent, merge_assignments_90percent) %>%
  group_by(group, scientific_name) %>% mutate(identity = max(identity)) %>% ungroup() %>% 
  group_by(group, identity, scientific_name, sample_id) %>% mutate(reads_count_new = sum(reads_count)) %>% ungroup() %>% select(-reads_count) %>% distinct()

merge_assignments <- merge_assignments %>% group_by(group, nb_in_group, family, genus) %>% mutate(minid = min(identity), maxid = max(identity)) %>% 
  arrange(group) %>% ungroup() %>%
  group_by(group, nb_in_group, family, genus, minid, sample_id) %>% mutate(reads_count_final = sum(reads_count_new)) %>% ungroup()

merge_assignments$keep <- ifelse(merge_assignments$identity == 1, "keep", 
                                 ifelse(merge_assignments$maxid < 1, "keep", "out"))

merge_assignments <- merge_assignments %>% filter(keep == "keep")

#
merge_assignments_clean <- merge_assignments %>% select(-reads_count_new, -minid, -maxid, -keep) %>% distinct() %>% 
  group_by(group, identity, scientific_name) %>% mutate(tot_reads = sum(reads_count_final), nsamples = n_distinct(sample_id)) %>% ungroup() %>% 
  filter(tot_reads >= 100) %>% filter(nsamples >= 10) 

merge_assignments_clean %>% select(identity, scientific_name, tot_reads, nsamples) %>% distinct() %>% filter(identity == 1)
merge_assignments_clean %>% select(identity, scientific_name, tot_reads, nsamples) %>% distinct() %>% filter(!identity == 1)
merge_assignments_wide <- merge_assignments_clean %>% pivot_wider(names_from = sample_id, values_from = reads_count_final) %>% arrange(tot_reads) %>% 
  group_by(group) %>% mutate(nb_in_group = n_distinct(scientific_name)) %>% arrange(group) %>% ungroup()
merge_assignments_wide %>% subset(!is.na(species)) %>% select(scientific_name) %>% distinct()
merge_assignments_wide %>% filter(!identity == 1) 
long_to_wide <- merge_assignments_wide %>% pivot_longer(cols = all_of(samples_name), names_to = "sample_id", values_to = "reads_count")
long_to_wide$reads_count[is.na(long_to_wide$reads_count)] <- 0
merge_assignments_wide_clean <- long_to_wide %>% pivot_wider(names_from = sample_id, values_from = reads_count) 

write_delim(merge_assignments_wide_clean, "output/louvain_merged_assignment_rdy_for_resampling.csv", delim = ",")

############################################################################
# 8.2 - Merge ASVs with same taxa names - family level (used in manuscript)
############################################################################
merge_assignments1 <- full_join(merge_assignments_100percent, merge_assignments_90percent) %>%
  group_by(group, nb_in_group, scientific_name) %>% mutate(identity = max(identity)) %>% 
  arrange(scientific_name) %>% arrange(group) %>% ungroup() %>%
  group_by(group, identity, scientific_name, sample_id) %>% mutate(reads_count_new = sum(reads_count)) %>% ungroup() %>% select(-reads_count) %>% distinct()

check <- merge_assignments1 %>% pivot_wider(names_from = sample_id, values_from = reads_count_new)
merge_assignments1 %>% filter(group == 615 & genus == "Andromeda")
merge_assignments1 %>% filter(group == 1 & genus == "Hedysarum")
merge_assignments1 %>% filter(group == 5 & genus == "Populus") %>% arrange(desc(sample_id))
full_join(merge_assignments_100percent, merge_assignments_90percent) %>% filter(group == 5 & genus == "Populus") %>% arrange(desc(sample_id))

merge_assignments1 <- merge_assignments1 %>% group_by(group, nb_in_group, family) %>% mutate(minid = min(identity), maxid = max(identity)) %>% 
  arrange(group) %>% ungroup() 

merge_assignments_100percent1 <- merge_assignments1 %>% subset(minid == 1) %>% mutate(reads_count_final = reads_count_new)
merge_assignments_90percent1 <- merge_assignments1 %>% subset(maxid <1) %>% mutate(reads_count_final = reads_count_new)
merge_assignments_90percent2 <- merge_assignments1 %>% subset(minid <1 & maxid == 1) %>% 
  mutate(keep = ifelse(identity == 1, "keep", "out")) %>%
  group_by(group, nb_in_group, family, keep, sample_id) %>% mutate(new = ifelse(identity == 1, 0, sum(reads_count_new))) %>% ungroup() %>%
  group_by(group, nb_in_group, family, sample_id) %>% mutate(tosum = max(new)) %>% ungroup() %>%
  mutate(reads_count_final = ifelse(identity == 1, reads_count_new+tosum, tosum)) %>% arrange(keep) %>% filter(keep == "keep") %>% select(-keep, -new, -tosum) 

merge_assignments2 <- rbind(merge_assignments_100percent1, merge_assignments_90percent1, merge_assignments_90percent2) 
# /!\ by merging the same families, it creates artifacts
# As within one community, candidate ASVs can be part of the same family of different 100% assigned taxa
# I do not think it matter too much but just as a side note
# + already discussed with Ulrike and she agreed as they are all impacted in the same way and it do not affect extinction detection

merge_assignments_clean1 <- merge_assignments2 %>% select(-reads_count_new, -minid, -maxid) %>% distinct() %>% 
  group_by(group, identity, scientific_name) %>% mutate(tot_reads = sum(reads_count_final), nsamples = n_distinct(sample_id)) %>% ungroup() %>% 
  filter(tot_reads >= 100) %>% filter(nsamples >= 10)

merge_assignments_wide1 <- merge_assignments_clean1 %>% pivot_wider(names_from = sample_id, values_from = reads_count_final) %>% arrange(tot_reads) %>% 
  group_by(group) %>% mutate(nb_in_group = n_distinct(scientific_name)) %>% arrange(group) %>% ungroup()
merge_assignments_wide1 %>% filter(!identity == 1)
long_to_wide1 <- merge_assignments_wide1 %>% pivot_longer(cols = all_of(samples_name), names_to = "sample_id", values_to = "reads_count")
long_to_wide1$reads_count[is.na(long_to_wide1$reads_count)] <- 0
merge_assignments_wide_clean1 <- long_to_wide1 %>% pivot_wider(names_from = sample_id, values_from = reads_count)

merge_assignments_wide_clean1

merge_assignments_wide_clean1 %>% subset(identity == 1)
merge_assignments_wide_clean1 %>% subset(!identity == 1)

write_delim(merge_assignments_wide_clean1, "output/louvain_merged_assignment_rdy_for_resampling_merge_family_level.csv", delim = ",")
