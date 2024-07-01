# TAXALOSS Project - script to confirm if the candidate taxa hold potential to be extinct
# Script from Jeremy Courtin
# Script last update - 03.04.2023

# Layout of the script: 
# 0 - GENERAL OPTION - SETTING UP SCRIPT
# 1 - Import the row data
# 2 - Check simulated coverage in modern timeslice
# 3 - Check simulated coverage in the old timeslices

###############################################################################
# 0 - GENERAL OPTION - SETTING UP SCRIPT
###############################################################################
rm(list = ls()) # Remove all the objects we created so far.

# Set up options for the rest of the script
options(stringsAsFactors=FALSE) # state that everything is reads as character 

# Load all needed packages
library(tidyverse)

###############################################################################
# 1 - Import the row data
###############################################################################
coverage <- read_delim("Simulate_sampling_GBIF/GBIF_DB_coverage.csv")

coverage %>% select(ncbi_family_acc, ncbi_genus_acc, ncbi_species_acc, GBIF_occurence, coverage_ext) %>% 
  distinct() %>% arrange(coverage_ext)

data <- coverage %>% mutate(taxlevel = "species") %>% select(ncbi_species_acc, taxlevel, GBIF_occurence, coverage_ext) %>% distinct() %>%
  setNames(c("final_name", "taxlevel", "occurence", "coverage")) %>% 
  arrange(desc(occurence)) %>% arrange(desc(final_name)) %>% dplyr::group_by(final_name) %>% fill(occurence) %>% ungroup() %>% distinct() %>%
  mutate(cov = ifelse(coverage == 0, NA, coverage)) %>% 
  arrange(desc(coverage)) %>% arrange(desc(final_name)) %>% dplyr::group_by(final_name) %>% mutate(occ = max(occurence)) %>% mutate(occurence = occ) %>%
  fill(cov) %>% ungroup() %>%
  mutate(coverage = ifelse(is.na(cov), 0, 1)) %>% select(-cov) %>% distinct() %>%
  group_by(final_name) %>% mutate(occ = length(final_name)) %>% arrange(desc(occ)) %>% ungroup() %>% select(-occ) %>%
  mutate(occurence = ifelse(is.na(occurence), 10, occurence))

# import the median data
median_read_resampled <- read_delim("median_reads_resampled_1000_final.csv", delim = ",", col_names = T)

median_read_resampled %>% subset(grepl("_cand_", uniq_label)) %>% summarise(sum(`13` > 0))
median_read_resampled %>% subset(grepl("_1_", uniq_label)) %>% summarise(sum(`13` > 0))
hist(median_read_resampled$`13`, breaks = 200)

cand <- median_read_resampled %>% subset(grepl("_cand_", uniq_label)) %>% pivot_longer(cols = `13`:`0`) %>% mutate(type = "cand")
db <- median_read_resampled %>% subset(grepl("_1_", uniq_label)) %>% pivot_longer(cols = `13`:`0`) %>% mutate(type = "db")
all <- rbind(cand, db) %>% mutate(name = as.numeric(name))

data %>% summarise(covered = sum(coverage == 1), uncovered = sum(coverage == 0)) %>%
  mutate(covered = covered/sum(covered, uncovered)*100) %>% mutate(uncovered = 100-covered) %>%
  pivot_longer(cols = covered:uncovered) %>%
  ggplot(aes(x = name, y = value)) +
  geom_point() +
  ylim(0, 100) +
  theme_light()

# gbif coverage by the database
ggplot(data) +
  geom_histogram(mapping = aes(x = occurence), position = "stack") +
  facet_wrap(vars(as.factor(coverage))) +
  theme_light()

ggplot(data) +
  geom_density(mapping = aes(x = occurence), position = "stack") +
  facet_wrap(vars(as.factor(coverage))) +
  theme_light()

# all taxa
ggplot(all) +
  geom_histogram(mapping = aes(x = value), position = "stack") +
  facet_wrap(vars(as.factor(name))) +
  theme_light()

# candidate
all %>% subset(type == "cand") %>%
  ggplot() +
  geom_histogram(mapping = aes(x = value), position = "stack") +
  facet_wrap(vars(as.factor(name))) +
  theme_light()

#database
all %>% subset(type == "db") %>%
  ggplot() +
  geom_histogram(mapping = aes(x = value), position = "stack") +
  facet_wrap(vars(as.factor(name))) +
  theme_light()

###############################################################################
# 2 - Check simulated coverage in modern timeslice
###############################################################################
data <- data %>% as.data.frame()

smpl <- lapply(1:999, function(x) {
  t <- data[sample(1:nrow(data), size = 232, replace = T),] %>% # Here set size = number of taxa detected in the modern timeslice
    group_by(coverage) %>% summarise(rel = n()) %>% pull(rel)
  (t/sum(t))*100
}) %>% do.call("rbind", .)

smpl %>% as_tibble() %>% setNames(c("uncovered", "covered")) %>% pivot_longer(cols = "uncovered":"covered") %>%
  ggplot(aes(x = name, y = value)) +
  geom_boxplot() +
  ylim(0, 100) +
  theme_light()

as_tibble(apply(smpl, 2, median))[1,]/100*232

# Plot
plot(0:1, apply(smpl, 2, median), ylim = c(0, 100)) 
segments(0:1, apply(smpl, 2, quantile, probs = 0.025), 0:1, apply(smpl, 2, quantile, probs = 0.975)) 
(apply(smpl, 2, quantile, probs = 0.025)/sum(apply(smpl, 2, quantile, probs = 0.025)))*100

###############################################################################
# 3 - Check simulated coverage in the old timeslices
###############################################################################
time_slices <- 13
sample <- c(278, 286, 290, 285, 282, 265, 269, 279, 272, 260, 259, 259, 257) # Here list the number of taxa detected in each of the 13 other timeslices (but modern)
data$occurence

# Measure expected loss
expect_loss1 <- lapply(1:999, function(z) {
  
  for(ts in 1:time_slices) {
    
    tab <- lapply(1:5, function(i) {
      data[sample(1:nrow(data), size = sample[ts], prob = data$occurence, replace = T),] %>% 
        filter(!duplicated(final_name))}) %>% 
      do.call("rbind",.) %>% 
      filter(!duplicated(final_name))
    
    if(ts == 1) {
      out <- tibble(taxa = data$final_name, db = data$coverage) %>%
        bind_cols(tibble(sapply(data$final_name, function(y) y %in% tab$final_name)) %>% setNames(glue::glue("ts_{ts}")))
    } else {
      out <- out %>%
        bind_cols(tibble(sapply(data$final_name, function(y) y %in% tab$final_name)) %>% setNames(glue::glue("ts_{ts}")))
    }
  }
  
  ext <- out %>% mutate(extinct = apply(out[,-c(1:2)], 1, function(x) x[length(x)] == 0),
                        ext_ts  = apply(out[,-c(1:2)], 1, function(x) which.min(rev(!x)))) %>% filter(extinct)
  
  tmp <- ext %>% group_by(db) %>% summarise(count = n()) 
  c((tmp$count/sum(tmp$count))*100, tmp$count) 
  
}) %>% do.call("rbind", .)

# Plot
expect_loss1 %>% as_tibble() %>% setNames(c("uncovered_prop", "covered_prop", "uncovered_count", "covered_count")) %>% pivot_longer(cols = "uncovered_prop":"covered_count") %>%
  subset(grepl("_prop", name)) %>%
  ggplot(aes(x = name, y = value)) +
  geom_boxplot() +
  ylim(0, 100) +
  theme_light()

bp <- barplot(height = apply(expect_loss1, 2, median)[2:1], ylim = c(0, 100), col = rev(c("orange", "cornflowerblue")))
segments(bp, apply(expect_loss1, 2, quantile, probs = 0.025)[2:1], bp, apply(expect_loss1, 2, quantile, probs = 0.975)[2:1], col = "black", lwd = 3) 
text(bp, 95, apply(expect_loss1, 2, median)[4:3]) 
axis(1, bp, labels = c("db_exist", "db_absent"))






