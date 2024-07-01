# TAXALOSS Project - script to prepare most Figures for the manuscript
# Script from Jeremy Courtin
# Script last update - 03.04.2023

# Layout of the script: 
# 0 - GENERAL OPTION - SETTING UP SCRIPT
# 1 - BASIC statistics
# 2 - Mean and Median reads count after resampling
# 3 - Set the candidates with similar proportion
# 4 - Measure Richness (just as indicator, not really reported in the manuscript)
# 5.1 - Beta diversity - PODANI
# 5.2 - Beta diversity - beta_div other beta diveristy estimate
# 6.1 - Plant types proportions - number of taxa
# 6.2 - Plant types proportions - read counts (not plotted but alternative plot)
# 7 - Family proportions
# 8.1 - Proportion of plant types and phylogenetic tree
# 8.2 - Proportion of plant types
# 8.3 - Plot tree
# 9 - Hypotheses testing
# 9.1 - cophenetic distances
# 9.2 - SCBD
# 9.3 - Community size - Specialist
# 9.4 - Metrics - plants type - measure at each iteration
# 9.5 - Metrics - families - measure at each iteration

# vocabulary used: candidate ASVs / taxa = non-dbASVs / taxa

###############################################################################
# 0 - GENERAL OPTION - SETTING UP SCRIPT
###############################################################################
rm(list = ls()) # Remove all the objects we created so far.

# Set up options for the rest of the script
options(stringsAsFactors=FALSE) # state that everything is reads as character 

# Load all needed packages
library(tidyverse)
library(ggplot2)
library(adespatial)
library(viridis)
library(hrbrthemes)

# Import the resampled data
load("resampled_data_final_1000iterations.rda")

# Rename the raw data
resampleAll <- smpl_raw
resampleAll[[1]][[1]]

# Set age timeslice: 
age <- tibble(timeslice = c(13:0), age = c(1000, 3000, 5000, 7000, 9000, 11000, 13000, 15000, 17000, 19000, 21000, 23000, 25000, 27000))

# Check the file before resampling / useful later
df_long_for_rar <- read_delim("output/df_long_for_rar.csv", delim = ",", col_names = T) %>%
  subset(read_counts > 0)

# Here you have the number of samples 
df_long_for_rar %>% #subset(mean_age < 28001) %>% 
  dplyr::select(mean_age, sample_id) %>%
  #summarise(sum = sum(read_counts)) %>%
  distinct() 

# Here you have general info about the taxa
taxa_info <- read_delim("output/final_dataset_before_resampling_family_merge_new.csv", delim = ",", col_names = T) %>%
  dplyr::select(identity, family, genus, species, scientific_name, group, nb_in_group) %>% mutate(uniq = 1:nrow(.)) %>% 
  mutate(level = ifelse(is.na(species), ifelse(is.na(genus), "family", "genus"), "species")) %>%
  mutate(genus = ifelse(is.na(genus), family, genus)) %>% mutate(species = ifelse(is.na(species), paste(genus, uniq, sep = " "), paste(species, uniq, sep = " "))) %>% mutate(identity = ifelse(identity == 0.9, "cand", "1")) %>%
  mutate(key = paste(nb_in_group, identity, family, scientific_name, sep = "_"))

taxa_info %>% filter(identity == 1) %>% filter(!level == "species")

# make the abundance table
abund_table <- lapply(resampleAll, function(x) x[[1]])
abund_table <- lapply(abund_table, function(x) replace(x, is.na(x), 0))

################################################################################
# 1 - BASIC statistics
################################################################################
metrics <- lapply(abund_table, function(x) pivot_longer(x, cols = `13`:`0`, names_to = "timeslice", values_to = "reads", values_drop_na = TRUE) %>%
                    mutate(timeslice = as.numeric(timeslice)) %>% mutate(timeslice = max(timeslice)-timeslice) %>%
                    separate(col = "uniq_label", into = c("group", "status", "type", "family"), sep = "_", remove = F, extra = "drop") %>%
                    group_by(status, timeslice) %>% subset(reads > 0) %>% mutate(n_status = n_distinct(uniq_label)) %>% dplyr::select(type, timeslice, n_status) %>% ungroup() %>% distinct())

metrics

################################################################################
# 2 - Mean and Median reads count after resampling
################################################################################
# mean_read_resamp <- lapply(resampleAll, function(x) x[[1]])
# mean_read_resamp <- lapply(mean_read_resamp, function(x) replace(x, is.na(x), 0))
# 
# # MEAN
# mean_read_resampled <- mean_read_resamp[[1]][,1] %>%
#   bind_cols(apply(abind::abind(lapply(mean_read_resamp, function(x) {x[,-1]}), along = 3), 1:2,
#                   function(x) mean(x))) %>%  mutate(sumVar = rowSums(.[-1])) %>%
#   filter(sumVar>0) %>% dplyr::select(-sumVar)
# 
# #plot mean richness
# data.frame(NonZeroCounts=colSums(mean_read_resampled!=0)[colSums(mean_read_resampled!=0)!=0]) %>% mutate(name = rownames(.)) %>% setNames(c("richness", "name")) %>%
#   subset(!name == "uniq_label") %>% as_tibble() %>% mutate(name = as.numeric(name)) %>%
#   ggplot(aes(x=name, y=richness)) +
#   geom_point() +
#   geom_line()
# 
# write_delim(mean_read_resampled, "mean_reads_resampled_1000_final.csv", delim = ",", col_names = T)
# 
# # MEDIAN
# median_read_resampled <- mean_read_resamp[[1]][,1] %>%
#   bind_cols(apply(abind::abind(lapply(mean_read_resamp, function(x) {x[,-1]}), along = 3), 1:2,
#                   function(x) median(x))) %>%  mutate(sumVar = rowSums(.[-1])) %>%
#   filter(sumVar>0) %>% dplyr::select(-sumVar)
# 
# #plot median richness
# data.frame(NonZeroCounts=colSums(median_read_resampled!=0)[colSums(median_read_resampled!=0)!=0]) %>% mutate(name = rownames(.)) %>% setNames(c("richness", "name")) %>%
#   subset(!name == "uniq_label") %>% as_tibble() %>% mutate(name = as.numeric(name)) %>%
#   ggplot(aes(x=name, y=richness)) +
#   geom_point() +
#   geom_line()
# 
# write_delim(median_read_resampled, "median_reads_resampled_1000_final.csv", delim = ",", col_names = T)

median_read_resampled <- read_delim("median_reads_resampled_1000_final.csv", delim = ",", col_names = T)

hist(median_read_resampled$`0`, breaks = 200)

################################################################################
# 3 - Set the candidates with similar proportion
################################################################################
# use the median data
candidates_absent_modern <- median_read_resampled %>% filter(grepl('_cand_', uniq_label)) %>% subset(`13` == 0) %>% dplyr::select(uniq_label)

# Check from the raw table
norm_raw <- df_long_for_rar %>% group_by(uniq_label) %>% mutate(nsample = n_distinct(sample_id), tot_read = sum(read_counts)) %>% 
  dplyr::select(uniq_label, type, nsample, tot_read) %>% distinct() %>% ungroup()

# modified to MAX
norm_stat_cand <- norm_raw %>% subset(uniq_label %in% candidates_absent_modern$uniq_label) %>% group_by(type) %>% 
  mutate(maxnsamples = max(nsample), maxreads = max(tot_read)) %>% ungroup() %>% 
  dplyr::select(type, maxnsamples, maxreads) %>% distinct()

cand_abs <- norm_raw %>% subset(uniq_label %in% candidates_absent_modern$uniq_label) %>% group_by(type) %>% 
  mutate(maxnsamples = max(nsample), maxreads = max(tot_read))

# General statistics to report in Supplementary material and also to check if we want to use the maximum or the quantile information to use for the modern candidates.
mean(cand_abs$tot_read)
mean(cand_abs$nsample)
quantile(cand_abs$nsample, 0.95)
quantile(cand_abs$tot_read, 0.95)

# when using the maximum values
norm_raw$keep <- ifelse(norm_raw$nsample < (norm_stat_cand$maxnsamples + norm_stat_cand$maxnsamples*0.1),
                        ifelse(norm_raw$tot_read < (norm_stat_cand$maxreads + norm_stat_cand$maxreads*0.1), "yes", "no"), "no")

# # when using the 95% quantile instead of the maximum values
# norm_raw$keep <- ifelse(norm_raw$nsample < 107 ,
#                          ifelse(norm_raw$tot_read < 10032.5, "yes", "no"), "no")

norm <- norm_raw %>% subset(keep == "yes" & type == 1)
norm_taxa <- norm$uniq_label 

#
cand_norm <- norm_raw %>% subset(keep == "yes" & type == "cand")
cand_norm_taxa <- cand_norm %>% subset(!uniq_label %in% all_of(candidates_absent_modern$uniq_label), drop = F)

# Report in supplementary material
mean(cand_abs$tot_read)
max(cand_abs$tot_read)
mean(cand_norm_taxa$tot_read)
max(cand_norm_taxa$tot_read)

mean(cand_abs$nsample)
max(cand_abs$nsample)
mean(cand_norm_taxa$nsample)
max(cand_norm_taxa$nsample)

# Plot for supplementary material
plot_cand_abs <- cand_abs %>% ungroup() %>% dplyr::select(nsample, tot_read) %>% mutate(type = "lost_candidate")
plot_cand_norm_taxa <- cand_norm_taxa %>% dplyr::select(nsample, tot_read) %>% mutate(type = "modern_candidate")

rbind(plot_cand_abs, plot_cand_norm_taxa) %>%
  #mutate(day = fct_reorder(test, value)) %>%
  ggplot(aes(fill=type, y=tot_read, x=type)) + 
  geom_boxplot(position="dodge", width=0.1, color="black", alpha=0.2) +
  xlab("") +
  ylab("reads count")

rbind(plot_cand_abs, plot_cand_norm_taxa) %>%
  #mutate(day = fct_reorder(test, value)) %>%
  ggplot(aes(fill=type, y=nsample, x=type)) + 
  geom_boxplot(position="dodge", width=0.5, color="black", alpha=0.2) +
  xlab("") +
  ylab("number of samples")

################################################################################
# 4 - Measure Richness (just as indicator, not really reported in the manuscript)
################################################################################
richness <- lapply(abund_table, function(x) x %>% setNames(c("uniq_label", "1000", "3000", "5000", "7000", "9000", "11000", "13000", "15000", "17000", "19000", "21000", "23000", "25000", "27000")) %>%
                     pivot_longer(cols = "1000":"27000") %>% group_by(name) %>% mutate(richness = sum(value>0)) %>% dplyr::select(name, richness) %>% distinct()) 
richness <- do.call(rbind, richness)

# median richness
med_richness <- richness %>% dplyr::group_by(name) %>% dplyr::summarise(median_rich = median(richness)) %>% setNames(c("age", "median_rich")) %>% mutate(age = as.numeric(age))
med_richness %>% arrange(desc(age))

# Plot richness
richness %>% 
  mutate(name = factor(name, levels=c("1000", "3000", "5000", "7000", "9000", "11000", "13000", "15000", "17000", "19000", "21000", "23000", "25000", "27000"))) %>%
  ggplot(aes(y=richness, x=name)) + 
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
  geom_boxplot(position="dodge", width=0.1, color="black", alpha=0.2) +
  xlab("") +
  ylab("Richness")

################################################################################
# 5.1 - Beta diversity - PODANI
################################################################################
# Here measure Podani with Jaccard (not used in the manuscript anymore)
beta_div_POD <- lapply(abund_table, function(x) x %>% setNames(c("uniq_label", "1000", "3000", "5000", "7000", "9000", "11000", "13000", "15000", "17000", "19000", "21000", "23000", "25000", "27000")) %>% 
                 as.data.frame() %>% `rownames<-`(.[,1]) %>% dplyr::select(-uniq_label) %>% t(.) %>% 
                 beta.div.comp(coef="J", quant=FALSE))

beta_div_POD[[1]]

beta_div_POD_rich <- lapply(beta_div_POD, function(x) x$rich)
beta_div_POD_rich1 <- lapply(beta_div_POD_rich, function(x) x %>% as.matrix() %>% as.data.frame() %>% mutate(row = as.numeric(rownames(.))) %>% 
         pivot_longer(cols = '1000':'27000') %>% mutate(name = as.numeric(name), keep = ifelse(name == (row - 2000), "yes", "no")) %>%
         filter(keep == "yes") %>% dplyr::select(name, value))
beta_div_POD_D <- lapply(beta_div_POD, function(x) x$D)
beta_div_POD_D1 <- lapply(beta_div_POD_D, function(x) x %>% as.matrix() %>% as.data.frame() %>% mutate(row = as.numeric(rownames(.))) %>% 
                               pivot_longer(cols = '1000':'27000') %>% mutate(name = as.numeric(name), keep = ifelse(name == (row - 2000), "yes", "no")) %>%
                               filter(keep == "yes") %>% dplyr::select(name, value))
podJ <-  lapply(beta_div_POD, function(x) x$repl)
podJ2 <- lapply(podJ, function(x) x %>% as.matrix() %>% as.data.frame() %>% mutate(row = as.numeric(rownames(.))) %>% 
         pivot_longer(cols = '1000':'27000') %>% mutate(name = as.numeric(name), keep = ifelse(name == (row - 2000), "yes", "no")) %>%
         filter(keep == "yes") %>% dplyr::select(name, value))

pod_plot1 <- do.call(rbind, beta_div_POD_rich1) %>% setNames(c("age", "Beta_rich"))

pod_plot2 <- do.call(rbind, beta_div_POD_D1) %>% setNames(c("age1", "Beta_D"))

pod_plot3 <- do.call(rbind, podJ2) %>% setNames(c("age2", "Beta_repl"))

pod_plot <- cbind(pod_plot1, pod_plot2, pod_plot3) %>% dplyr::select(age, Beta_rich, Beta_D, Beta_repl) %>% as_tibble()

# Extratc info on richness changes, replacement rates and dissimilarity
med_pod <- pod_plot %>% group_by(age) %>% summarise(median_Beta_repl = median(Beta_repl), 
                                                    median_Beta_D = median(Beta_D),
                                                    median_Beta_rich = median(Beta_rich))

# Plot it (Replacement rate)
pod_plot %>%
  mutate(age = factor(age, levels=c("1000", "3000", "5000", "7000", "9000", "11000", "13000", "15000", "17000", "19000", "21000", "23000", "25000"))) %>%
  ggplot(aes(y=Beta_repl, x=age)) + 
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
  geom_boxplot(position="dodge", width=0.1, color="black", alpha=0.2) +
  xlab("") +
  ylab("Podani")

# Save it to re-use later
write_delim(pod_plot, "beta_diversity_1000iterations.csv", delim = ",")
pod_plot <- read_delim("beta_diversity_1000iterations.csv", delim = ",")

# Plot it again - Here use the replacement rate curve in Figure 2
pod_plot %>% 
  group_by(age) %>% mutate(medianD = median(Beta_D), meanD = mean(Beta_D), q95D = quantile(Beta_D, 0.95), q5D = quantile(Beta_D, 0.05),
                           medianR = median(Beta_rich), meanR = mean(Beta_rich), q95R = quantile(Beta_rich, 0.95), q5R = quantile(Beta_rich, 0.05),
                           medianREP = median(Beta_repl), meanREP = mean(Beta_repl), q95REP = quantile(Beta_repl, 0.95), q5REP = quantile(Beta_repl, 0.05)) %>% 
  ungroup() %>% subset(age > 1000) %>%
  ggplot(aes(x = age)) +
  geom_ribbon(aes(ymin = q5D, ymax = q95D), fill = "grey90", alpha=0.5) +
  #geom_ribbon(aes(ymin = q5R, ymax = q95R), fill = "red", alpha=0.5) +
  geom_ribbon(aes(ymin = q5REP, ymax = q95REP), fill = "red", alpha=0.05) +
  geom_point(aes(y = medianD), colour = "black") +
  #geom_point(aes(y = medianR, color = "richness difference"), colour = "blue") +
  geom_point(aes(y = medianREP), colour = "red") +
  geom_line(aes(y = medianD, color = "dissimilarity"), colour = "black") +
  #geom_line(aes(y = meanR, color = "richness difference"), colour = "blue") +
  geom_line(aes(y = medianREP, color = "replacement"), colour = "red") +
  scale_x_reverse(breaks=seq(2000, 25000, 2000)) +
  xlab("Age (cal. yrs BP)") +
  ylab("") +
  labs(title = "black: dissimilarity; red: replacement") +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

################################################################################
# 5.2 - Beta diversity - beta_div other beta diveristy estimate
################################################################################
# Measure changes of betadiversity between samples
beta_div_samples <- lapply(abund_table, function(x) x %>% as.data.frame() %>% `rownames<-`(.[,1]) %>% dplyr::select (-uniq_label) %>% beta.div("hellinger", nperm=999))
lapply(beta_div_samples, function(x) x[[2]])

# Measure beta diversity of taxa (contribution to functionality of the system)
beta_div_taxa <- do.call(rbind, lapply(beta_div_samples, function(x) x$LCBD %>% as.matrix() %>% as.data.frame() %>% summarise(taxon = rownames(.), value = V1))) %>% as_tibble() %>% group_by(taxon) %>% summarise(median = median(value), mean = mean(value))

# Save it
write_delim(beta_div_taxa, "beta_div_taxa.csv", delim = ",")

# SCBD per sample
beta_div_samples1 <- do.call(rbind, lapply(beta_div_samples, function(x) x[[2]])) %>% as.data.frame()
beta_div_samples1 %>% pivot_longer(cols = 13:0) %>% group_by(name) %>% summarise(median_SCBD = median(value)) %>% mutate(timeslice = as.numeric(name)) %>% 
  left_join(age) %>%
  ggplot(aes(x = age, y = median_SCBD)) +
  geom_line()

# merge the betadiversity estimates
beta_div_merge <- beta_div_samples1 %>% pivot_longer(cols = 13:0) %>% group_by(name) %>% summarise(median_SCBD = median(value)) %>% mutate(timeslice = as.numeric(name)) %>% 
  left_join(age) %>% dplyr::select(age, median_SCBD)
do.call(rbind, lapply(beta_div_samples, function(x) x[[2]])) %>% as.data.frame()

# make beta-diversity dataset
med_div <- left_join(med_pod, med_richness)
med_div <- left_join(med_div, beta_div_merge)

write_delim(med_div, "Median_Diversity_Estimates.csv", delim = ",", col_names = T)

# Check how many taxa are present in all timeslices 
countall <- lapply(abund_table, function(x) pivot_longer(x, cols = `13`:`0`, names_to = "timeslice", values_to = "reads", values_drop_na = TRUE) %>% 
                     mutate(count0 = if_else(reads == 0,  1, 0)) %>% dplyr::select(uniq_label, count0) %>% distinct() %>% group_by(uniq_label) %>% mutate(count = sum(count0)) %>%
                     ungroup() %>% subset(count == 0) %>% summarise(count = n_distinct(uniq_label)))
countall <- do.call(rbind, countall)
median(countall$count) # median = 90
sd(countall$count) # sd = 10.24554

################################################################################
# 6.1 - Plant types proportions - number of taxa
################################################################################
# Candidates taxa count
plant_type_abund_rawc_cand <- lapply(abund_table, function(x) pivot_longer(x, cols = `13`:`0`, names_to = "timeslice", values_to = "reads", values_drop_na = TRUE) %>% 
                                       mutate(timeslice = as.numeric(timeslice)) %>%
                                       separate(col = "uniq_label", into = c("group", "status", "type", "family"), sep = "_", remove = F, extra = "drop") %>% 
                                       filter(status == "cand") %>% group_by(type, timeslice) %>% subset(reads > 0) %>% mutate(n_type = n_distinct(uniq_label)) %>% dplyr::select(type, timeslice, n_type) %>% ungroup() %>% distinct())

plant_type_abundc_cand <- do.call(rbind, plant_type_abund_rawc_cand) %>% group_by(type, timeslice) %>% mutate(mean_type = median(n_type), sd_type = sd(n_type)) %>% 
  dplyr::select(type, timeslice, mean_type, sd_type) %>% ungroup() %>% distinct() %>% pivot_wider(names_from = type, values_from = c(mean_type, sd_type)) %>% left_join(age)
plant_type_abundc_cand <-  plant_type_abundc_cand %>% mutate(type = "cand") %>% 
  dplyr::select(age, mean_type_grass, mean_type_herb, mean_type_shrub, mean_type_tree, type) %>% 
  set_names(c("age", "grass", "herb", "shrub", "tree", "type")) %>% 
  pivot_longer(cols = c("grass", "herb", "shrub", "tree"))

# 100% taxa count
plant_type_abund_rawc_1 <- lapply(abund_table, function(x) pivot_longer(x, cols = `13`:`0`, names_to = "timeslice", values_to = "reads", values_drop_na = TRUE) %>% 
                                    mutate(timeslice = as.numeric(timeslice)) %>%
                                    separate(col = "uniq_label", into = c("group", "status", "type", "family"), sep = "_", remove = F, extra = "drop") %>% 
                                    filter(status == "1") %>% group_by(type, timeslice) %>% subset(reads > 0) %>% mutate(n_type = n_distinct(uniq_label)) %>% dplyr::select(type, timeslice, n_type) %>% ungroup() %>% distinct())

plant_type_abundc_1 <- do.call(rbind, plant_type_abund_rawc_1) %>% group_by(type, timeslice) %>% mutate(mean_type = median(n_type), sd_type = sd(n_type)) %>% 
  dplyr::select(type, timeslice, mean_type, sd_type) %>% ungroup() %>% distinct() %>% pivot_wider(names_from = type, values_from = c(mean_type, sd_type)) %>% left_join(age)
plant_type_abundc_1 <- plant_type_abundc_1 %>% mutate(type = "100perc") %>% 
  dplyr::select(age, mean_type_grass, mean_type_herb, mean_type_shrub, mean_type_tree, type) %>% 
  set_names(c("age", "grass", "herb", "shrub", "tree", "type")) %>% 
  pivot_longer(cols = c("grass", "herb", "shrub", "tree"))
plant_type_abundc_1

# plot plant types 
p_types <- rbind(plant_type_abundc_cand, plant_type_abundc_1) %>% mutate(label = paste(name, type, sep = " "))
p_types_p <- p_types %>% group_by(age) %>%
  mutate(n = sum(value)) %>% ungroup() %>%
  group_by(label, age) %>%
  mutate(percentage = value / n) %>%
  ungroup()

# Figure 1 - For illustrator
p_types_p %>% mutate(label1 = paste(type, name)) %>% dplyr::select(age, percentage, label, label1) %>%
ggplot(aes(x=age, y=percentage, fill=label1)) + 
  geom_area(alpha=0.6, size=.25, colour="black") +
  scale_fill_viridis(discrete = T) +
  scale_x_reverse(n.breaks = 6) +
#  theme_ipsum(axis_text_size = 20, subtitle_size = 20, axis_title_size = 25) +
  theme(legend.text = element_text(size=20))

#get some metrics to report in the manuscript
save_ptypes <- p_types_p %>%dplyr::mutate(big_type = ifelse(name == "grass", "herb", ifelse(name == "herb", "herb", "woody"))) %>%
  select(age, n, type, name, value, percentage, big_type) %>% 
  dplyr::group_by(age, big_type) %>% dplyr::mutate(v_bt = sum(value)) %>% dplyr::mutate(perc_bt = v_bt/n) %>% ungroup() %>%
  dplyr::group_by(age, type) %>% dplyr::mutate(v_type = sum(value)) %>% dplyr::mutate(perc_type = v_type/n) %>% ungroup()
  arrange(desc(age))

# save
write_delim(save_ptypes, "Plant_types_Fig1_metrics.csv", delim = ",")

# 
p_types_p %>% filter(age > 15000) %>% 
  group_by(type, name) %>% mutate(mean = mean(percentage)*100) %>% ungroup() %>%
  group_by(type, age) %>% mutate(nt = sum(value)) %>% ungroup() %>%
  group_by(age, type, name) %>% mutate(percnt = value/nt) %>% ungroup() %>%
  group_by(type, name) %>% mutate(meannt = mean(percnt)*100) %>% ungroup() %>%
  select(type, name, mean, meannt) %>% distinct() %>%
  group_by(name) %>% mutate(all = sum(mean)) %>% ungroup()

p_types_p %>% filter(age < 9000) %>% 
  group_by(type, name) %>% mutate(mean = mean(percentage)*100) %>% ungroup() %>%
  group_by(type, age) %>% mutate(nt = sum(value)) %>% ungroup() %>%
  group_by(age, type, name) %>% mutate(percnt = value/nt) %>% ungroup() %>%
  group_by(type, name) %>% mutate(meannt = mean(percnt)*100) %>% ungroup() %>%
  select(type, name, mean, meannt) %>% distinct() %>%
  group_by(name) %>% mutate(all = sum(mean)) %>% ungroup()

p_types_p %>% filter(age == c(15000, 13000, 11000, 9000)) %>% 
  group_by(type, name) %>% mutate(mean = mean(percentage)*100) %>% ungroup() %>%
  group_by(type, age) %>% mutate(nt = sum(value)) %>% ungroup() %>%
  group_by(age, type, name) %>% mutate(percnt = value/nt) %>% ungroup() %>%
  group_by(type, name) %>% mutate(meannt = mean(percnt)*100) %>% ungroup() %>%
  select(type, name, mean, meannt) %>% distinct() %>%
  group_by(name) %>% mutate(all = sum(mean)) %>% ungroup()
  
################################################################################
# 6.2 - Plant types proportions - read counts (not plotted but alternative plot)
################################################################################
# Candidates reads count
plant_type_abund_raw_cand <- lapply(abund_table, function(x) pivot_longer(x, cols = `13`:`0`, names_to = "timeslice", values_to = "reads", values_drop_na = TRUE) %>% 
                                       mutate(timeslice = as.numeric(timeslice)) %>%
                                       separate(col = "uniq_label", into = c("group", "status", "type", "family"), sep = "_", remove = F, extra = "drop") %>% 
                                       filter(status == "cand") %>% group_by(type, timeslice) %>% subset(reads > 0) %>% mutate(n_type = sum(reads)) %>% dplyr::select(type, timeslice, n_type) %>% ungroup() %>% distinct())

plant_type_abund_cand <- do.call(rbind, plant_type_abund_raw_cand) %>% group_by(type, timeslice) %>% mutate(mean_type = median(n_type), sd_type = sd(n_type)) %>% 
  dplyr::select(type, timeslice, mean_type, sd_type) %>% ungroup() %>% distinct() %>% pivot_wider(names_from = type, values_from = c(mean_type, sd_type)) %>% left_join(age)
plant_type_abund_cand <-  plant_type_abund_cand %>% mutate(type = "cand") %>% 
  dplyr::select(age, mean_type_grass, mean_type_herb, mean_type_shrub, mean_type_tree, type) %>% 
  set_names(c("age", "grass", "herb", "shrub", "tree", "type")) %>% 
  pivot_longer(cols = c("grass", "herb", "shrub", "tree"))

# 100% reads count
plant_type_abund_raw_1 <- lapply(abund_table, function(x) pivot_longer(x, cols = `13`:`0`, names_to = "timeslice", values_to = "reads", values_drop_na = TRUE) %>% 
                                    mutate(timeslice = as.numeric(timeslice)) %>% 
                                    separate(col = "uniq_label", into = c("group", "status", "type", "family"), sep = "_", remove = F, extra = "drop") %>% 
                                    filter(status == "1") %>% group_by(type, timeslice) %>% subset(reads > 0) %>% mutate(n_type = n_distinct(uniq_label)) %>% dplyr::select(type, timeslice, n_type) %>% ungroup() %>% distinct())

plant_type_abund_1 <- do.call(rbind, plant_type_abund_raw_1) %>% group_by(type, timeslice) %>% mutate(mean_type = median(n_type), sd_type = sd(n_type)) %>% 
  dplyr::select(type, timeslice, mean_type, sd_type) %>% ungroup() %>% distinct() %>% pivot_wider(names_from = type, values_from = c(mean_type, sd_type)) %>% left_join(age)
plant_type_abund_1 <- plant_type_abund_1 %>% mutate(type = "100perc") %>% 
  dplyr::select(age, mean_type_grass, mean_type_herb, mean_type_shrub, mean_type_tree, type) %>% 
  set_names(c("age", "grass", "herb", "shrub", "tree", "type")) %>% 
  pivot_longer(cols = c("grass", "herb", "shrub", "tree"))
plant_type_abund_1

# plot plant types
p_types_reads <- rbind(plant_type_abund_cand, plant_type_abund_1) %>% mutate(label = paste(name, type, sep = " "))
p_types_p_reads <- p_types_reads %>% group_by(age) %>%
  mutate(n = sum(value)) %>% ungroup() %>%
  group_by(label, age) %>%
  mutate(percentage = value / n) %>%
  ungroup()

p_types_p %>% dplyr::select(age, percentage, label) %>%
  ggplot(aes(x=age, y=percentage, fill=label)) + 
  geom_area(alpha=0.6, size=.25, colour="black") +
  scale_fill_viridis(discrete = T) +
  # scale_x_reverse(n.breaks = 6) +
  #  theme_ipsum(axis_text_size = 20, subtitle_size = 20, axis_title_size = 25) +
  theme(legend.text = element_text(size=20))

# To report in the supplementary information
p_types_p %>% filter(age > 15000) %>% 
  group_by(type, name) %>% mutate(mean = mean(percentage)*100) %>% ungroup() %>%
  group_by(type, age) %>% mutate(nt = sum(value)) %>% ungroup() %>%
  group_by(age, type, name) %>% mutate(percnt = value/nt) %>% ungroup() %>%
  group_by(type, name) %>% mutate(meannt = mean(percnt)*100) %>% ungroup() %>%
  select(type, name, mean, meannt) %>% distinct() %>%
  group_by(name) %>% mutate(all = sum(mean)) %>% ungroup()

p_types_p %>% filter(age < 9000) %>% 
  group_by(type, name) %>% mutate(mean = mean(percentage)*100) %>% ungroup() %>%
  group_by(type, age) %>% mutate(nt = sum(value)) %>% ungroup() %>%
  group_by(age, type, name) %>% mutate(percnt = value/nt) %>% ungroup() %>%
  group_by(type, name) %>% mutate(meannt = mean(percnt)*100) %>% ungroup() %>%
  select(type, name, mean, meannt) %>% distinct() %>%
  group_by(name) %>% mutate(all = sum(mean)) %>% ungroup()

p_types_p %>% filter(age == c(15000, 13000, 11000, 9000)) %>% 
  group_by(type, name) %>% mutate(mean = mean(percentage)*100) %>% ungroup() %>%
  group_by(type, age) %>% mutate(nt = sum(value)) %>% ungroup() %>%
  group_by(age, type, name) %>% mutate(percnt = value/nt) %>% ungroup() %>%
  group_by(type, name) %>% mutate(meannt = mean(percnt)*100) %>% ungroup() %>%
  select(type, name, mean, meannt) %>% distinct() %>%
  group_by(name) %>% mutate(all = sum(mean)) %>% ungroup()

################################################################################
# 7 - Family proportions
################################################################################
fam_sumamry <- lapply(abund_table, function(x) x %>% pivot_longer(cols = `13`:`0`, names_to = "timeslice", values_to = "reads", values_drop_na = TRUE) %>%
                        mutate(timeslice = as.numeric(timeslice)) %>% left_join(age) %>%
                        separate(col = "uniq_label", into = c("group", "status", "type", "family"), sep = "_", remove = F, extra = "drop") %>%
                        group_by(family, timeslice) %>% subset(reads > 0) %>% mutate(n_family = n_distinct(uniq_label), reads_fam = sum(reads)) %>% dplyr::select(family, timeslice, age, n_family, reads_fam) %>% ungroup() %>% distinct())

fam_summary <- do.call(rbind, fam_sumamry) %>% group_by(family, timeslice) %>% mutate(median_fam_n = median(n_family), sd_family_n = sd(n_family), mean_fam_reads = mean(reads_fam), sd_family_reads = sd(reads_fam)) %>% 
  dplyr::select(family, timeslice, age, median_fam_n, sd_family_n, mean_fam_reads, sd_family_reads) %>% ungroup() %>% distinct()

uniqfam <- fam_summary %>% select(family) %>% distinct()

write_delim(uniqfam, "family.csv", delim = ",")

# Prepare and read the family type file
type_fam <- read_delim("type_family.csv")

# add plant type info + relative abundance
fam_p <- fam_summary %>% group_by(age) %>%
  mutate(n = sum(median_fam_n), n1 = sum(mean_fam_reads)) %>% ungroup() %>%
  group_by(family, age) %>%
  mutate(percentage_n = median_fam_n / n*100, percentage_reads = mean_fam_reads / n1*100) %>%
  ungroup() %>%
  left_join(type_fam)

order <- fam_p %>% subset(age == 1000) %>% arrange(desc(percentage_reads)) %>% dplyr::select(family)
order1 <- fam_p %>% group_by(family) %>% mutate(max = max(percentage_reads)) %>% dplyr::select(family, max, type) %>% distinct() %>% arrange(desc(max)) %>% arrange(desc(type))
                          
library(tidypaleo)
theme_set(theme_paleo(8))
# Nice looking stratigraphic plot - confirm composition turnover all reads (candidates and dbtaxa) - not used in manuscript
fam_p %>% dplyr::select(age, percentage_n, percentage_reads, family, type) %>%
  group_by(family) %>% mutate(max = max(percentage_reads)) %>%
  mutate(family = factor(family, levels= order1$family)) %>%
  #subset(max > 0.5) %>%
  ggplot(aes(x = percentage_reads, y = age, fill = type)) +
  geom_areah(alpha=0.6, size=.25, colour="black") +
  scale_y_reverse() +
  facet_geochem_gridh(vars(family)) +
  labs(x = "Relative proportion of reads (%)", y = "Age (cal. yrs BP")

################################################################################
# 8.1 - Proportion of plant types and phylogenetic tree
################################################################################
# Set metrics for phylogenetic tree
metrics <- norm_raw
metrics$status <- ifelse(metrics$type == "cand", "cand", ifelse(metrics$keep == "yes", "sim100", "100"))
metrics$status <- factor(metrics$status, levels = c("100", "sim100", "cand"))
metrics <- metrics %>% separate(uniq_label, into = c("group", "type", "plant", "family", "scientific"), sep =  "_", remove = F) %>% mutate(key = paste(group, type, family, scientific, sep = "_"))

#List of taxa + info
taxalist <- df_long_for_rar %>% dplyr::select(uniq_label, type, family_name, scientific_name, group) %>% distinct() %>% arrange(family_name) %>%
  mutate(key = paste(group, type, family_name, scientific_name, sep = "_"))

# Merge the list of taxa to extra info
plot_x <- left_join(taxalist, dplyr::select(taxa_info, key, family, genus, species, level), by = "key") %>% distinct()
plot_x %>% group_by(key) %>% dplyr::mutate(n_d = sum(group > 0)) %>% arrange(desc(n_d))
plot_x %>% select(uniq_label, key) %>% distinct()

################################################################################
# 8.2 - Proportion of plant types
################################################################################
plot_x %>% mutate(plot = paste(type, level, sep = "_")) %>% 
  group_by(plot) %>% mutate(n = n_distinct(key)) %>% dplyr::select(plot, n) %>% distinct() %>%
  ggplot(aes(x="", y=n, fill=plot)) + 
  geom_bar(stat="identity", width=1, color="white") +
 # scale_fill_viridis(discrete = T) +
  theme(legend.text = element_text(size=20))

# Number of taxa per family
df_long_for_rar %>% dplyr::select(uniq_label, type, family_name) %>% distinct() %>%
  #group_by(type2, type) %>% 
  group_by(family_name) %>% 
  summarise(n = n_distinct(uniq_label)) %>% arrange(desc(n))

# number of taxa per plant type
df_long_for_rar %>% dplyr::select(uniq_label, type, type2) %>% distinct() %>%
  #group_by(type2, type) %>% 
  group_by(type2) %>% 
  summarise(n = n_distinct(uniq_label))

# summary with the value to report in the manuscript
plot_x %>% mutate(plot = paste(type, level, sep = "_")) %>% 
  group_by(plot) %>% mutate(n = n_distinct(key)) %>% dplyr::select(plot, n, type, level) %>% distinct() %>%
  group_by(type) %>% mutate(sum = sum(n)) %>% ungroup() %>% mutate(rel_abund = n/sum*100) %>% arrange(plot) %>%
  mutate(tot_sum = sum(n)) %>%
  group_by(level) %>% mutate(rel_abund_level = sum(n)/tot_sum*100)

# Phylogenetic tree
taxalist <- left_join(taxalist, dplyr::select(taxa_info, key, family, genus, species), by = "key")
taxalist <- left_join(taxalist, dplyr::select(metrics, uniq_label, status), by = "uniq_label")

# load libraries
library("V.PhyloMaker2")
library("ape")
library("Biostrings")
library("ggplot2")
library("ggtree")

### make the example file
taxa_tree <- taxalist %>% dplyr::select(species, genus, family, status, uniq_label)
cand_tax <-  taxa_tree %>% filter(status == "cand") %>% filter(!uniq_label %in% cand_abs$uniq_label) %>% dplyr::select(species)
cand_tax$species <- gsub(" ", "_", cand_tax$species)
lost_cand <- taxa_tree %>% filter(uniq_label %in% cand_abs$uniq_label) %>% dplyr::select(species)
lost_cand$species <- gsub(" ", "_", lost_cand$species)
tax_100 <- taxa_tree %>% filter(status == "100") %>% dplyr::select(species)
tax_100$species <- gsub(" ", "_", tax_100$species)
tax_sim100 <- taxa_tree %>% filter(status == "sim100") %>% dplyr::select(species)
tax_sim100$species <- gsub(" ", "_", tax_sim100$species)

group_tree <- list(cand = cand_tax$species, "100" = tax_100$species, sim100 = tax_sim100$species, lost_taxa = lost_cand$species) 

### run the function
tree <- phylo.maker(taxa_tree, tree = GBOTB.extended.TPL, scenarios=c("S1","S2","S3"))

tree1 <- groupOTU(tree$scenario.1, group_tree)
#tree1 <- groupOTU(tree$scenario.2, group_tree)
#tree1 <- groupOTU(tree$scenario.3, group_tree)

# usually, used tree$scenario.3. Maybe tree$scenario.1 is better -> length of branches are good and shows the family and genus level more accurately.
# But what is best? I want to highlight that each tip is still a taxon.

tree_tib <- as_tibble(tree1) %>% separate(col = "label", into = "genus", sep = "_", remove = F) 
fam_info <- taxa_tree %>% dplyr::select(family, genus) %>% distinct()
tree_tib <- left_join(tree_tib, fam_info, by = "genus")

fam_label <- tree_tib %>% dplyr::select(parent, family) %>% distinct() %>% 
  group_by(family) %>% summarise(parent = min(parent)) %>% arrange(parent) %>% subset(!is.na(family))

################################################################################
# 8.3 - Plot tree
################################################################################
# ggtree(tree$scenario.3, branch.length='none', layout='circular') + 
#   geom_tiplab

# Use this tree in addition to the second one to modify in illustrator and create the final plot.
ggtree(tree1, aes(color=group), branch.length='none', layout='circular') +
  scale_color_manual(values=c("red3", "black", "blue", "black")) +
  geom_tiplab() +
  theme(
    legend.position = "none",
    plot.margin = grid::unit(c(-15, -15, -15, -15), "mm")
  )

######
# Plot tree but need manual finish with ADOBE illustrator
ggtree(tree1, aes(color=group), branch.length='none', layout='circular') +
  scale_color_manual(values=c("red3", "black", "blue", "black")) +
  #  geom_tiplab() +
  geom_cladelab(
    node = 315,
    label = "",
    align = TRUE,
    angle = -35,
    offset.text = 0.05,
    hjust = "center",
    fontsize = 2,
    offset = .2,
    barsize = 1.5
  ) +
  theme(
    legend.position = "none",
    plot.margin = grid::unit(c(-15, -15, -15, -15), "mm")
  )

################################################################################
# 9 - Hypotheses testing
################################################################################
################################################################################
# 9.1 - cophenetic distances
################################################################################
# cophenetic.phylo computes the pairwise distances between the pairs of tips from a phylogenetic tree using its branch lengths.
# dist.nodes does the same but between all nodes, internal and terminal, of the tree. 

#Manually change the column names
phydist <- cophenetic.phylo(tree1) %>% as.matrix() %>% as.data.frame() %>% mutate(name1 = row.names(.)) %>% 
  pivot_longer(cols = "Asteraceae_328":"Equisetum_25") %>% mutate(keep = ifelse(name1==name, 0, 1)) %>% subset(keep == 1) %>% dplyr::select(-keep)
tree_merge <- tree_tib %>% dplyr::select(family, genus, group, label)

phydist <- left_join(phydist, tree_merge, by = c("name" = "label")) %>% 
  setNames(c("label", "label1", "value", "family1", "genus1", "group1")) %>%
  left_join(tree_merge, by ="label")

lostcand_100 <- phydist %>% subset(group1 == "lost_taxa" & !group == "cand") %>% summarise(type = "lostcand_vs_100", value = value)
lostcand_rest <- phydist %>% subset(group1 == "lost_taxa" & !group == "lost_taxa") %>% summarise(type = "lostcand_vs_rest", value = value)
lostcand_lostcand <- phydist %>% subset(group1 == "lost_taxa" & group == "lost_taxa") %>% summarise(type = "lostcand_vs_lostcand", value = value) 
d100_100 <- phydist %>% subset(!group1 == "cand" & !group == "cand") %>% summarise(type = "100_vs_100", value = value) 

rbind(lostcand_rest, lostcand_lostcand) %>%
  ggplot(aes( y=value, x=type)) + 
  geom_boxplot(position="dodge", width=0.5, color="black", alpha=0.2) +
  geom_jitter(width=0.05, alpha = 0.1)+
  xlab("") +
  ylab("Phylogenetic distances (cophenetic)") +
  ylim(150, 300) +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

wilcox.test(lostcand_rest$value, lostcand_lostcand$value) # p-value = 0.2743
wilcox.test(lostcand_rest$value, lostcand_lostcand$value, alternative = "greater") #  p-value = 0.8629
mean(lostcand_rest$value) # 253.0174
mean(lostcand_lostcand$value) # 241.8972

################################################################################
# 9.2 - SCBD
################################################################################
library("iNEXT")
library("adespatial")
library("adegraphics")

beta_div <- lapply(abund_table, function(x) x %>% as.data.frame() %>% mutate(uniq_label = paste(uniq_label, `13`, sep = "_")) %>% `rownames<-`(.[,1]) %>% dplyr::select (-uniq_label) %>% t(.) %>% beta.div("hellinger", nperm=999))

# Prepare for supplementary document
beta_list <- lapply(beta_div, function(x) x$SCBD %>% as.data.frame())

comm_beta_list <- do.call(rbind, beta_list) %>% setNames("scbd") %>% mutate(uniq = row.names(.)) %>% as_tibble() %>%
  select(uniq, scbd) %>% separate(uniq, into = c("community", "status", "type", "family", "scientificname"), sep = "_", extra = "drop") %>%
  dplyr::group_by(community, status, type, family, scientificname) %>% dplyr::summarise(meanscbd = mean(scbd)) %>% ungroup() %>% distinct() %>% arrange(desc(meanscbd)) %>% arrange(community)

# add abundance
med_reads <- lapply(abund_table, function(x) x %>% as.data.frame() %>% dplyr::group_by(uniq_label) %>% dplyr::summarise(reads_sum = sum(`13`:`0`)))
med_reads1 <- do.call(rbind, med_reads) %>% dplyr::group_by(uniq_label) %>% dplyr::summarise(reads_mean = mean(reads_sum)) %>% separate(uniq_label, into = c("community", "status", "type", "family", "scientificname"), sep = "_", extra = "drop") 

comm_beta_list1 <- comm_beta_list %>% left_join(med_reads1) %>% dplyr::group_by(community) %>% dplyr::mutate(n_taxa = n_distinct(scientificname), n_candidate = sum(status == "cand")) %>% ungroup()

commASV <- read_delim("Intermediate_results/ASVs_taxa_after_louvain_community.csv", delim = ",", col_names = T) %>% dplyr::group_by(group) %>% dplyr::summarise(n_ASVs = n_distinct(uniq_id), n_candidatesASVs = sum(identity < 1)) %>% 
  mutate(n_dbASVs = n_ASVs-n_candidatesASVs, community = as.character(group)) %>% select(-group)

comm_beta_list2 <- comm_beta_list1 %>% left_join(commASV)

write_delim(comm_beta_list2, "Supplementary_community_compo.csv", delim = ",")

# lost candidates
beta_div_mean_cand <- lapply(beta_div, function(x) as.data.frame(x$SCBD) %>% mutate(uniq_label = rownames(.)) %>% as_tibble() %>% setNames(c("SCBD", "uniq_label")) %>%
                                 separate(uniq_label, into = c("group", "status", "type", "family", "scientific_name", "pres_modern"), sep = "_", remove = F) %>%
                                 mutate(id = paste(group, status, type, family, scientific_name, sep = "_")) %>% subset(pres_modern == "0") %>% subset(status == "cand") %>%
                                 group_by(status) %>% summarise(mean_scbd = mean(SCBD)))
SCBD_cand <- do.call(rbind, beta_div_mean_cand) %>% mutate(status = "lost_candidate")

# candidates present modern with similar parameters as candidates absent modern
beta_div_mean_cand_pres <- lapply(beta_div, function(x) as.data.frame(x$SCBD) %>% mutate(uniq_label = rownames(.)) %>% as_tibble() %>% setNames(c("SCBD", "uniq_label")) %>%
                                    separate(uniq_label, into = c("group", "status", "type", "family", "scientific_name", "pres_modern"), sep = "_", remove = F) %>%
                                    mutate(id = paste(group, status, type, family, scientific_name, sep = "_")) %>% subset(!pres_modern == "0") %>% subset(id %in% cand_norm_taxa$uniq_label, drop = F) %>%
                                    group_by(status) %>% summarise(mean_scbd = mean(SCBD)))
SCBD_cand_pres <- do.call(rbind, beta_div_mean_cand_pres) %>% mutate(status = "present_candidate")

# 100% taxa
beta_div_mean_100 <- lapply(beta_div, function(x) as.data.frame(x$SCBD) %>% mutate(uniq_label = rownames(.)) %>% as_tibble() %>% setNames(c("SCBD", "uniq_label")) %>%
                                    separate(uniq_label, into = c("group", "status", "type", "family", "scientific_name", "pres_modern"), sep = "_", remove = F) %>%
                                    mutate(id = paste(group, status, type, family, scientific_name, sep = "_")) %>% subset(status == "1") %>%
                                    group_by(status) %>% summarise(mean_scbd = mean(SCBD)))
SCBD_100 <- do.call(rbind, beta_div_mean_100) %>% mutate(status = "all_100%")

#
mean(SCBD_cand$mean_scbd) 
mean(SCBD_cand_pres$mean_scbd)
mean(SCBD_100$mean_scbd) 
wilcox.test(SCBD_cand$mean_scbd, SCBD_cand_pres$mean_scbd) 
wilcox.test(SCBD_cand$mean_scbd, SCBD_100$mean_scbd) 

SCBD_pres <- rbind(SCBD_cand_pres, SCBD_cand)
SCBD_pres <- rbind(SCBD_pres, SCBD_100)

# Merge both Figures
SCBD_pres %>%
  # mutate(day = fct_reorder(test, value)) %>%
  mutate(status = factor(status, levels=c("all_100%", "present_candidate", "lost_candidate"))) %>%
  ggplot(aes(fill=status, y=mean_scbd, x=status)) + 
  #geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
  geom_boxplot(position="dodge", width=0.5, color="black", alpha=0.2) +
  geom_jitter(width=0.05, alpha = 0.1)+
  xlab("") +
  ylab("SCBD score") +
  ylim(0, 0.0003) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

SCBD_pres %>%
  # mutate(day = fct_reorder(test, value)) %>%
  mutate(status = factor(status, levels=c("all_100%", "present_candidate", "lost_candidate"))) %>%
  ggplot(aes(fill=status, y=mean_scbd, x=status)) + 
  #geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
  geom_jitter(width=0.05, alpha = 0.1)+
  geom_boxplot(position="dodge", width=0.5, color="black", alpha=0.2) +
  xlab("") +
  ylab("SCBD score") +
  ylim(0.0043, 0.0046) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

################################################################################
# 9.3 - Community size - Specialist
################################################################################
# candidates present modern similar to candidates absent modern
specialist_prescand <- lapply(abund_table, function(x) separate(x, uniq_label, into = c("group", "status", "type", "family", "scientific_name"), sep = "_", remove = F) %>%
                                group_by(uniq_label) %>% mutate(sum = sum(`13`:`0`)) %>% subset(sum > 0) %>% ungroup() %>%
                              group_by(group) %>% mutate(community_size = n_distinct(uniq_label)) %>% subset(!`13`==0) %>% 
                              subset(uniq_label %in% cand_norm_taxa$uniq_label, drop = F) %>% mutate(group = as.numeric(group)) %>% group_by(status) %>% summarise(community_size = mean(community_size)))
specialist_prescand <- do.call(rbind, specialist_prescand) %>% mutate(status = "present_candidate")

# candidates absent modern
specialist_cand <- lapply(abund_table, function(x) separate(x, uniq_label, into = c("group", "status", "type", "family", "scientific_name"), sep = "_", remove = F) %>%
                            group_by(uniq_label) %>% mutate(sum = sum(`13`:`0`)) %>% subset(sum > 0) %>% ungroup() %>%
                              group_by(group) %>% mutate(community_size = n_distinct(uniq_label)) %>% subset(`13`==0) %>% 
                              subset(status == "cand") %>% mutate(group = as.numeric(group)) %>% group_by(status) %>% summarise(community_size = mean(community_size)))
specialist_cand <- do.call(rbind, specialist_cand) %>% mutate(status = "lost_candidate")

# candidates present modern similar to candidates absent modern
specialist_100 <- lapply(abund_table, function(x) separate(x, uniq_label, into = c("group", "status", "type", "family", "scientific_name"), sep = "_", remove = F) %>%
                           group_by(uniq_label) %>% mutate(sum = sum(`13`:`0`)) %>% subset(sum > 0) %>% ungroup() %>%
                           group_by(group) %>% mutate(community_size = n_distinct(uniq_label)) %>%
                           subset(status == "1") %>% mutate(group = as.numeric(group)) %>% group_by(status) %>% summarise(community_size = mean(community_size)))
specialist_100 <- do.call(rbind, specialist_100) %>% mutate(status = "all_100%")

mean(specialist_cand$community_size)
mean(specialist_prescand$community_size)
mean(specialist_100$community_size)

wilcox.test(specialist_cand$community_size, specialist_prescand$community_size)
wilcox.test(specialist_cand$community_size, specialist_100$community_size)

specialist <- rbind(specialist_cand, specialist_prescand)
specialist <- rbind(specialist, specialist_100)

# Figure -> modify in illustrator
specialist %>%
  # mutate(day = fct_reorder(test, value)) %>%
  mutate(status = factor(status, levels=c("all_100%", "present_candidate", "lost_candidate"))) %>%
  ggplot(aes(fill=status, y=community_size, x=status)) + 
 # geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
  geom_boxplot(position="dodge", width=0.5, color="black", alpha=0.2) +
  geom_jitter(width=0.05, alpha = 0.1)+
  xlab("") +
  ylab("Community size")+
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

################################################################################
# 9.4 - Metrics - plants type - measure at each iteration
################################################################################
# make table for each group
type_work1 <- lapply(abund_table, function(x) separate(x, col = "uniq_label", into = c("group", "status", "type", "family"), sep = "_", remove = F, extra = "drop") %>%
                      group_by(type, status) %>% summarise(n_absent = sum(`13` == 0), n_recent = sum(`13` > 0)) %>% ungroup() %>%
                      group_by(type) %>% mutate(sumfam = sum(n_absent) + sum(n_recent)) %>% ungroup() %>%
                      group_by(type, status) %>% mutate(sumfam_group = sum(n_absent) + sum(n_recent)) )
type_work <- do.call(rbind, type_work1) %>% subset(status == "cand" & sumfam > 4) %>% mutate (prop_abs = n_absent/sumfam_group) %>% dplyr::select(type, prop_abs)

#
type_work_cand <- lapply(abund_table, function(x) separate(x, col = "uniq_label", into = c("group", "status", "type", "family"), sep = "_", remove = F, extra = "drop") %>% 
                           subset(status == "cand") %>% group_by(type) %>% summarise(n_absent = sum(`13` == 0), n_recent = sum(`13` > 0)) %>% ungroup())
type_work_cand <- do.call(rbind, type_work_cand)

type_work_cand_metrics <- type_work_cand %>% group_by(type) %>% 
  mutate(meanabsent_cand = mean(n_absent), sdabsent_cand = sd(n_absent), meanrecent_cand = mean(n_recent), sdrecent_cand = sd(n_recent)) %>%
  dplyr::select(type, meanabsent_cand, sdabsent_cand, meanrecent_cand, sdrecent_cand) %>% distinct()

sum(type_work_cand_metrics$meanabsent_cand)

#
type_work_sim100 <- lapply(abund_table, function(x) separate(x, col = "uniq_label", into = c("group", "status", "type", "family"), sep = "_", remove = F, extra = "drop") %>% 
                             subset(status == "1") %>% subset(uniq_label %in% norm_taxa, drop = F) %>% group_by(type) %>%  summarise(n_absent = sum(`13` == 0), n_recent = sum(`13` > 0)) %>% ungroup())
type_work_sim100 <- do.call(rbind, type_work_sim100)

type_work_sim100_metrics <- type_work_sim100 %>% group_by(type) %>% 
  mutate(meanabsent_sim100 = mean(n_absent), sdabsent_sim100 = sd(n_absent), meanrecent_sim100 = mean(n_recent), sdrecent_sim100 = sd(n_recent)) %>%
  dplyr::select(type, meanabsent_sim100, sdabsent_sim100, meanrecent_sim100, sdrecent_sim100) %>% distinct()

#
type_work_sim100_plot <- type_work_sim100 %>% mutate(prop = n_absent/(n_absent+n_recent)) %>% dplyr::select(type, prop) %>% setNames(c("type", "100%")) %>% pivot_longer(cols = "100%")
type_work_cand_plot <- type_work_cand %>% mutate(prop = n_absent/(n_absent+n_recent)) %>% dplyr::select(type, prop) %>% setNames(c("type", "candidate")) %>% pivot_longer(cols = "candidate")

type_plot <- rbind(type_work_sim100_plot, type_work_cand_plot)

# statistic test
##### Custom t.test #####
pairwise.t.test.with.t.and.df <- function (x, g, p.adjust.method = p.adjust.methods, pool.sd = !paired, 
                                           paired = FALSE, alternative = c("two.sided", "less", "greater"), 
                                           ...) 
{
  if (paired & pool.sd) 
    stop("pooling of SD is incompatible with paired tests")
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
  g <- factor(g)
  p.adjust.method <- match.arg(p.adjust.method)
  alternative <- match.arg(alternative)
  if (pool.sd) {
    METHOD <- "t tests with pooled SD"
    xbar <- tapply(x, g, mean, na.rm = TRUE)
    s <- tapply(x, g, sd, na.rm = TRUE)
    n <- tapply(!is.na(x), g, sum)
    degf <- n - 1
    total.degf <- sum(degf)
    pooled.sd <- sqrt(sum(s^2 * degf)/total.degf)
    compare.levels <- function(i, j) {
      dif <- xbar[i] - xbar[j]
      se.dif <- pooled.sd * sqrt(1/n[i] + 1/n[j])
      t.val <- dif/se.dif
      if (alternative == "two.sided") 
        2 * pt(-abs(t.val), total.degf)
      else pt(t.val, total.degf, lower.tail = (alternative == 
                                                 "less"))
    }
    compare.levels.t <- function(i, j) {
      dif <- xbar[i] - xbar[j]
      se.dif <- pooled.sd * sqrt(1/n[i] + 1/n[j])
      t.val = dif/se.dif 
      t.val
    }       
  }
  else {
    METHOD <- if (paired) 
      "paired t tests"
    else "t tests with non-pooled SD"
    compare.levels <- function(i, j) {
      xi <- x[as.integer(g) == i]
      xj <- x[as.integer(g) == j]
      t.test(xi, xj, paired = paired, alternative = alternative, 
             ...)$p.value
    }
    compare.levels.t <- function(i, j) {
      xi <- x[as.integer(g) == i]
      xj <- x[as.integer(g) == j]
      t.test(xi, xj, paired = paired, alternative = alternative, 
             ...)$statistic
    }
    compare.levels.df <- function(i, j) {
      xi <- x[as.integer(g) == i]
      xj <- x[as.integer(g) == j]
      t.test(xi, xj, paired = paired, alternative = alternative, 
             ...)$parameter
    }
  }
  PVAL <- pairwise.table(compare.levels, levels(g), p.adjust.method)
  TVAL <- pairwise.table.t(compare.levels.t, levels(g), p.adjust.method)
  if (pool.sd) 
    DF <- total.degf
  else
    DF <- pairwise.table.t(compare.levels.df, levels(g), p.adjust.method)           
  ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL, 
              p.adjust.method = p.adjust.method, t.value = TVAL, dfs = DF)
  class(ans) <- "pairwise.htest"
  ans
}
pairwise.table.t <- function (compare.levels.t, level.names, p.adjust.method) 
{
  ix <- setNames(seq_along(level.names), level.names)
  pp <- outer(ix[-1L], ix[-length(ix)], function(ivec, jvec) sapply(seq_along(ivec), 
                                                                    function(k) {
                                                                      i <- ivec[k]
                                                                      j <- jvec[k]
                                                                      if (i > j)
                                                                        compare.levels.t(i, j)               
                                                                      else NA
                                                                    }))
  pp[lower.tri(pp, TRUE)] <- pp[lower.tri(pp, TRUE)]
  pp
}
#####
#fit the one-way ANOVA model
type_stat <- type_plot %>% mutate(name_test = paste(type, name, sep = "_"))
model <- aov(value ~ name_test, data = type_stat)

#view model output
summary(model)

type_test <- pairwise.t.test.with.t.and.df(type_stat$value, type_stat$name_test, p.adjust.method="bonferroni")
type_testpv <- as.data.frame(type_test$p.value) %>% mutate(param1 = row.names(.)) %>% 
  pivot_longer(cols =  'grass_100%':'shrub_candidate', names_to = "param2") %>% 
  separate(col = "param1", into = c("type1", "status1"), sep = "_", remove = F, extra = "drop") %>%
  separate(col = "param2", into = c("type2", "status2"), sep = "_", remove = F, extra = "drop") %>%
  mutate(keep = ifelse(type1 == type2, "yes", "no")) %>% subset(keep == "yes") %>% filter(!is.na(value)) %>% dplyr::select(-keep)
type_testpv %>% dplyr::select(type1, status1, status2, value)

# Plot
type_plot %>%
 #mutate(day = fct_reorder(name, value)) %>%
  #  mutate(day = factor(test, levels=c("Thur", "Fri", "Sat", "Sun"))) %>%
  ggplot(aes(fill=name, y=value, x=type)) + 
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
#  geom_boxplot(position="dodge", width=0.5, color="black", alpha=0.2) +
  xlab("") +
  ylab("Proportion of taxa absent from modern")

# Compare proportion between types
#type_test <- pairwise.t.test(type_stat$value, type_stat$name_test, p.adjust.method="bonferroni")
type_test <- pairwise.t.test.with.t.and.df(type_stat$value, type_stat$name_test, p.adjust.method="bonferroni")
type_test$t.value
type_test$p.value

type_testpv <- as.data.frame(type_test$p.value) %>% mutate(param1 = row.names(.)) %>% 
  pivot_longer(cols =  'grass_100%':'shrub_candidate', names_to = "param2") %>% 
  separate(col = "param1", into = c("type1", "status1"), sep = "_", remove = F, extra = "drop") %>%
  separate(col = "param2", into = c("type2", "status2"), sep = "_", remove = F, extra = "drop") %>%
  mutate(keep = ifelse(type1 == type2, "yes", "no")) %>% subset(keep == "no" & status1 == "candidate" & status2 == "candidate") %>% filter(!is.na(value)) %>% dplyr::select(-keep)
type_testpv %>% dplyr::select(type1, status1, status2, value)

type_stat %>% subset(name == "candidate") %>% group_by(type) %>% summarise(mean = mean(value))

# Plot
order_type <- type_plot %>% subset(name == "candidate") %>% group_by(type) %>% mutate(mean = mean(value)) %>% dplyr::select(type, mean) %>% distinct() %>% arrange(desc(mean))
type_plot %>% subset(name == "candidate") %>%
  # mutate(day = fct_reorder(name, value)) %>%
  mutate(type = factor(type, levels=order_type$type)) %>%
  ggplot(aes(fill=type, y=value, x=type)) + 
 # geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
  geom_boxplot(position="dodge", width=0.5, color="black", alpha=0.5) +
  geom_jitter(width=0.05, alpha = 0.1)+
  xlab("") +
  ylab("Proportion of candidate taxa absent from modern timeslice") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

################################################################################
# 9.5 - Metrics - families - measure at each iteration
################################################################################
# make table for each group
fam_work <- lapply(abund_table, function(x) separate(x, col = "uniq_label", into = c("group", "status", "type", "family"), sep = "_", remove = F, extra = "drop") %>% 
                     group_by(family, status) %>% summarise(n_absent = sum(`13` == 0), n_recent = sum(`13` > 0)) %>% ungroup() %>%
                     group_by(family) %>% mutate(sumfam = sum(n_absent) + sum(n_recent)) %>% ungroup() %>%
                     group_by(family, status) %>% mutate(sumfam_group = sum(n_absent) + sum(n_recent)) )
fam_work_plot <- do.call(rbind, fam_work) %>% subset(status == "cand" & sumfam > 4) %>% mutate (prop_abs = n_absent/sumfam) %>% dplyr::select(family, prop_abs)
fam_work_plot1 <- do.call(rbind, fam_work) %>% subset(status == "cand" & sumfam_group > 4) %>% mutate (prop_abs = n_absent/sumfam_group) %>% dplyr::select(family, prop_abs)
fam_work_plot2 <- do.call(rbind, fam_work) %>% subset(status == "cand" & sumfam > 4) %>% mutate (prop_abs = n_absent/sumfam_group) %>% dplyr::select(family, prop_abs)

fam_work_plot1 %>% dplyr::select(family) %>% distinct()

#
fam_work_cand <- lapply(abund_table, function(x) separate(x, col = "uniq_label", into = c("group", "status", "type", "family"), sep = "_", remove = F, extra = "drop") %>% 
                          subset(status == "cand") %>% group_by(family) %>% summarise(n_absent = sum(`13` == 0), n_recent = sum(`13` > 0)) %>% ungroup())
fam_work_cand <- do.call(rbind, fam_work_cand)

#
fam_work_sim100 <- lapply(abund_table, function(x) separate(x, col = "uniq_label", into = c("group", "status", "type", "family"), sep = "_", remove = F, extra = "drop") %>% 
                            subset(status == "1") %>% subset(uniq_label %in% norm_taxa, drop = F) %>% group_by(family) %>%  summarise(n_absent = sum(`13` == 0), n_recent = sum(`13` > 0)) %>% ungroup())
fam_work_sim100 <- do.call(rbind, fam_work_sim100)

#
fam_work_sim100_plot <- fam_work_sim100 %>% mutate(prop = n_absent/(n_absent+n_recent)) %>% dplyr::select(family, prop) %>% setNames(c("family", "100%")) %>% pivot_longer(cols = "100%")
fam_work_cand_plot <- fam_work_cand %>% mutate(prop = n_absent/(n_absent+n_recent)) %>% dplyr::select(family, prop) %>% setNames(c("family", "candidate")) %>% pivot_longer(cols = "candidate")

fam_plot <- rbind(fam_work_sim100_plot, fam_work_cand_plot)

# statistic test
#fit the one-way ANOVA model
fam_stat <- fam_plot %>% mutate(name_test = paste(family, name, sep = "_"))
fam_stat %>% group_by(name_test) %>% mutate(mean = mean(value)) %>% ungroup() %>% 
  dplyr::select(family, name, mean) %>% distinct() %>% arrange(family)
model <- aov(value ~ name_test, data = type_stat)

#view model output
summary(model)

fam_test <- pairwise.t.test(fam_stat$value, fam_stat$name_test, p.adjust.method="bonferroni")
as.data.frame(fam_test$p.value) %>% names()
fam_testpv <- as.data.frame(fam_test$p.value) %>% mutate(param1 = row.names(.)) %>% 
  pivot_longer(cols =  'Adoxaceae_candidate':'Saxifragaceae_candidate', names_to = "param2") %>% 
  separate(col = "param1", into = c("family1", "status1"), sep = "_", remove = F, extra = "drop") %>%
  separate(col = "param2", into = c("family2", "status2"), sep = "_", remove = F, extra = "drop") %>%
  mutate(keep = ifelse(family1 == family2, "yes", "no")) %>% subset(keep == "yes") %>% filter(!is.na(value)) %>% dplyr::select(-keep) %>% dplyr::select(family1, value)
fam_sign <- fam_testpv %>% subset(value < 0.05)
fam_sign_n <- fam_sign$family1

#
fam_stat <- fam_work_plot %>% mutate(name_test = paste(family, status, sep = "_"))
fam_stat %>% group_by(name_test) %>% mutate(mean = mean(prop_abs)) %>% ungroup() %>% 
  dplyr::select(family, status, mean) %>% distinct() %>% arrange(family)
model <- aov(value ~ name_test, data = type_stat)

#view model output
summary(model)

fam_test <- pairwise.t.test(fam_stat$prop_abs, fam_stat$name_test, p.adjust.method="bonferroni")
as.data.frame(fam_test$p.value) %>% names()
fam_testpv <- as.data.frame(fam_test$p.value) %>% mutate(param1 = row.names(.)) %>% 
  pivot_longer(cols =  'Asparagaceae_cand':'Salicaceae_cand', names_to = "param2") %>% 
  separate(col = "param1", into = c("family1", "status1"), sep = "_", remove = F, extra = "drop") %>%
  separate(col = "param2", into = c("family2", "status2"), sep = "_", remove = F, extra = "drop") %>%
  mutate(keep = ifelse(family1 == family2, "yes", "no")) %>% subset(keep == "no") %>% filter(!is.na(value)) %>% dplyr::select(-keep) %>% dplyr::select(family1, family2, value)
fam_sign <- fam_testpv %>% subset(value > 0.05) %>% setNames(c("Family_1", "Family_2", "p-value"))
fam_sign_n <- fam_sign$Family_1

#
fam_stat <- fam_work_plot1 %>% mutate(name_test = paste(family, status, sep = "_"))
fam_stat %>% group_by(name_test) %>% mutate(mean = mean(prop_abs)) %>% ungroup() %>% 
  dplyr::select(family, status, mean) %>% distinct() %>% arrange(family)
model <- aov(value ~ name_test, data = type_stat)

fam_stat %>% group_by(family) %>% summarise(mean = mean(prop_abs)) %>% arrange(desc(mean))
#view model output
summary(model)

#fam_test <- pairwise.t.test(fam_stat$prop_abs, fam_stat$name_test, p.adjust.method="bonferroni")
fam_test <- pairwise.t.test.with.t.and.df(fam_stat$prop_abs, fam_stat$name_test, p.adjust.method="bonferroni")
fam_test$p.value

fam_test$t.value

as.data.frame(fam_test$p.value) %>% names()
fam_testpv <- as.data.frame(fam_test$p.value) %>% mutate(param1 = row.names(.)) %>% 
  pivot_longer(cols =  'Asparagaceae_cand':'Salicaceae_cand', names_to = "param2") %>% 
  separate(col = "param1", into = c("family1", "status1"), sep = "_", remove = F, extra = "drop") %>%
  separate(col = "param2", into = c("family2", "status2"), sep = "_", remove = F, extra = "drop") %>%
  mutate(keep = ifelse(family1 == family2, "yes", "no")) %>% subset(keep == "no") %>% filter(!is.na(value)) %>% dplyr::select(-keep) %>% dplyr::select(family1, family2, value)
fam_testpv %>% subset(value > 0) %>% as.data.frame()

fam_sign <- fam_testpv %>% subset(value > 0.05)
fam_sign_n <- fam_sign$family1

# Plot
fam_work_plot %>%
  ggplot(aes(fill=family, y=prop_abs, x=family)) + 
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
  geom_boxplot(position="dodge", width=0.1, color="black", alpha=0.2) +
  xlab("") +
  ylab("proportion of candidate absent from modern - relative to total number of taxa per family")
fam_work_plot %>% group_by(family) %>% mutate(mean = mean(prop_abs)) %>% dplyr::select(family, mean) %>% distinct() %>% arrange(desc(mean))

order_fam_1 <- fam_work_plot1 %>% group_by(family) %>% mutate(mean = mean(prop_abs)) %>% dplyr::select(family, mean) %>% distinct() %>% arrange(desc(mean))
fam_work_plot1 %>%
  mutate(family = factor(family, levels=order_fam_1$family)) %>%
  ggplot(aes(fill=family, y=prop_abs, x=family)) + 
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
  geom_boxplot(position="dodge", width=0.1, color="black", alpha=0.2) +
  xlab("") +
  ylab("proportion of candidate absent from modern")
fam_work_plot1 %>% group_by(family) %>% mutate(mean = mean(prop_abs)) %>% dplyr::select(family, mean) %>% distinct() %>% arrange(desc(mean))

fam_work_plot1 %>%
  mutate(family = factor(family, levels=order_fam_1$family)) %>%
  ggplot(aes(y=prop_abs, x=family)) + 
  #geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
  geom_boxplot(position="dodge", width=0.5, color="black", alpha=0.5) +
  geom_jitter(width=0.05, alpha = 0.1)+
  xlab("") +
  ylab("proportion of candidate absent from modern") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

fam_work_plot2 %>%
  ggplot(aes(fill=family, y=prop_abs, x=family)) + 
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
  geom_boxplot(position="dodge", width=0.1, color="black", alpha=0.2) +
  xlab("") +
  ylab("proportion of candidate absent from modern")

fam_work_plot2 %>% group_by(family) %>% mutate(mean = mean(prop_abs)) %>% dplyr::select(family, mean) %>% distinct() %>% arrange(desc(mean))

