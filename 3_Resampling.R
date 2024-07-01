# TAXALOSS Project - resampling script for timeslices
# Script from Simeon Lisovski & Jeremy Courtin
# Script last update - 03.04.2023

# Layout of the script: 
# 0 - GENERAL OPTION - SETTING UP SCRIPT
# 1 - Import the raw data
# 2 - Import plant types data and manual check it!
# 3 - Prepare the similar100% taxa list
# 4 - Resampling
# 4.1 - Set variables for the resmapling
# 4.2 - Resample -> to run on the linux server (using parallel package)
# 4.3 - Save resampled data!

############################################################################
# 0 - GENERAL OPTION - SETTING UP SCRIPT
############################################################################
rm(list = ls()) # Remove all the objects we created so far.

# Set up options for the rest of the script
options(stringsAsFactors=FALSE)

# Load all needed packages
library(tidyverse)
library(data.table)

############################################################################
# 1 - Import the raw data
############################################################################
# Start resampling script
df_raw <- read_delim("output/louvain_merged_assignment_rdy_for_resampling_merge_family_level.csv", delim = ",", col_names = T)

# Get sample names
samples_name <- df_raw %>%
  select(-group, -nb_in_group, -identity, -family, -genus, -species, -scientific_name, -tot_reads, -nsamples) %>%
  names()

############################################################################
# 2 - Import plant types data and manual check it!
############################################################################
# Import the plant types data
types <- read_delim("plant_types_list.csv", delim = ",", col_names = T) %>% select(scientific_name, type2, type3) %>% distinct()

# Create the dataset for resampling - check the plant - types if everything is in order
df_for_rar <- left_join(df_raw, types, by = "scientific_name")
write_delim(df_for_rar, "check_plant_types_final_manually_family_merge_new.csv", delim = ",")
df_for_rar %>% filter(is.na(type2)) # check if all taxa have a type

# Manually check that every entry have a plant type assigned and save new file before reimporting it:
# df_for_rar <- read_delim("plant_types/done_check_plant_types_final_manually_family_merge_new.csv", delim = ",", col_names = T)

# Exclude fern, mosses and aquatic taxa
df_for_rar <- df_for_rar %>% filter(!type2 == "fern") %>% filter(!type2 == "moss") %>% filter(type3 == "terrestrial")
df_for_rar %>% filter(identity == 1)
df_for_rar %>% filter(!identity == 1)
write_delim(df_for_rar, "output/final_dataset_before_resampling_family_merge_new.csv", delim = ",")

# Finalize layout for resampling script
df_long <- df_for_rar %>% pivot_longer(cols = all_of(samples_name), names_to = "sample_id", values_to = "reads_count") %>%
  select(-genus, -species) %>% subset(reads_count > 0)
df_long$type <- ifelse(df_long$identity == 1, 1, "cand")
df_long <- separate(data = df_long,
                    col = "sample_id",
                    into = c("core", "mean_age"),
                    sep = "_",
                    remove = FALSE,
                    extra = "drop")
df_long_for_rar <- df_long %>% mutate(uniq_label = paste(group, type, type2, family, scientific_name, sep = "_")) %>%
  select(uniq_label, group, type, type2, family, scientific_name, type3, core, mean_age, sample_id, reads_count, nb_in_group) %>%
  setNames(c("uniq_label", "ngroup", "type", "type2", "family_name", "scientific_name", "type3", "core", "mean_age", "sample_id", "read_counts", "group"))
df_long_for_rar$group <- as.integer(df_long_for_rar$group)
df_long_for_rar$mean_age <- as.numeric(df_long_for_rar$mean_age)
df_long_for_rar %>% arrange(desc(mean_age))

df_long %>% group_by(core) %>% summarise(sum = n_distinct(sample_id)) %>% summarise(sum = sum(sum))

write_delim(df_long_for_rar, "output/df_long_for_rar.csv", delim = ",")

############################################################################
# 3 - Prepare the similar100% taxa list
############################################################################
cand_absent_modern <- df_long_for_rar %>% filter(type == "cand" & mean_age < 2001 & read_counts > 0) %>% select(uniq_label) %>% distinct()
cand_absent_modern <- cand_absent_modern$uniq_label
cand_absent_stat <- df_long_for_rar %>% filter(type == "cand") %>% subset(!uniq_label %in% all_of(cand_absent_modern), drop = F) %>%
  group_by(uniq_label) %>% mutate(tot_reads = sum(read_counts), nsamples = n_distinct(sample_id)) %>% ungroup() %>%
  select(uniq_label, tot_reads, nsamples) %>% mutate(avgreads = mean(tot_reads), minreads = min(tot_reads), maxreads = max(tot_reads), avgsamples = mean(nsamples), minsamples = min(nsamples), maxsamples = max(nsamples)) %>%
  select(-uniq_label, -tot_reads, -nsamples) %>% distinct()

norm <- df_long_for_rar %>% group_by(uniq_label) %>% mutate(tot_reads = sum(read_counts), nsamples = n_distinct(sample_id)) %>% ungroup() %>%
  filter(type == "1" & tot_reads <= cand_absent_stat$maxreads + (0.1*cand_absent_stat$maxreads) & nsamples <= cand_absent_stat$maxsamples + (0.1*cand_absent_stat$maxsamples)) %>%
  select(uniq_label, nsamples) %>% distinct()
norm_taxa <- df_long_for_rar %>% subset(uniq_label %in% all_of(cand_absent_modern), drop = F) %>% select(uniq_label) %>% distinct()
norm_taxa <- full_join(norm_taxa, norm, by = "uniq_label")
norm_taxa <- unique(norm_taxa$uniq_label)

df_long_for_rar <- as.data.table(df_long_for_rar)

############################################################################
# 4 - Resampling
############################################################################
############################################################################
# 4.1 - Set variables for the resmapling
############################################################################
# Set time slices
t_slices <- seq(0, 32000, 2000)

# Set names of time slices
m_slices <- t_slices + median(diff(t_slices))/2
df_long_for_rar[,time_slice := as.integer(cut(mean_age, breaks = 40))]

# number of resampling iterations
nSample <- 1000

# Get only timeslices as old as 28000
df_long_for_rar[mean_age<=28000,
][,time_slice := as.integer(cut(mean_age, t_slices))
][,time_slice := abs(time_slice - max(.SD[,"time_slice"]))][order(time_slice)] %>% as_tibble() %>% subset(mean_age <= 28000) %>% group_by(time_slice) %>%
  summarise(nreads = sum(read_counts), ntaxa = n_distinct(uniq_label), ncore = n_distinct(core), nsamples = n_distinct(sample_id))

## Plotting
tbHist <- df_long_for_rar[, .(.N, age = mean(mean_age)), by = time_slice]
ggplot(tbHist, mapping = aes(x = age, y = N)) +
  geom_rect(data = tibble(slice1 = t_slices, slice2 = c(t_slices[-1], max(df_long_for_rar$mean_age)), l = 1:length(t_slices)),
            inherit.aes=FALSE,
            aes(xmin = slice1, xmax = slice2, ymin = 0, ymax = max(tbHist$N),
                fill = factor(l)), alpha=0.5) +
  scale_fill_manual(values = viridis::magma(length(t_slices))) +
  geom_bar(stat = "identity") +
  theme_bw()

############################################################################
# 4.2 - Resample -> to run on the linux server (using parallel package)
############################################################################
# resampling the full 1000 times do not always work so now I repeat the following script 10 times with 100 time iterations and save the results.
# Then we merge the 10 saved files into one!

###############################################################################
smpl <- parallel::mclapply(1:100, function(y) {
  
  resampleList <- df_long_for_rar[mean_age<=28000,
  ][,time_slice := as.integer(cut(mean_age, t_slices))
  ][,time_slice := abs(time_slice - max(.SD[,"time_slice"]))][order(time_slice)]
  #resampleList %>% as_tibble(.) %>% select(time_slice, mean_age) %>% group_by(time_slice) %>% mutate(min_age = min(mean_age), max_age = max(mean_age), mean_age = mean(mean_age))%>% distinct()
  
  resampleCores   <- resampleList[, .SD[core %in% unique(core)[
    sample(1:length(unique(core)), min(resampleList[, length(unique(core)), by = time_slice]$V1))
  ]] , by = time_slice]
  
  resampleSamples <- resampleCores[,.SD[sample_id %in% unique(sample_id)[
    sample(1:length(unique(sample_id)), min(resampleCores[, length(unique(sample_id)), by = time_slice]$V1))
  ]], by = time_slice]
  
  resampleReads   <- resampleSamples[,num:=1:.N
  ][, cbind(.SD, dup = 1:read_counts), by = "num"
  ][,read_counts := 1,][,c("uniq_label", "time_slice")][,
                                                        .SD[sample(1:.N, min(resampleSamples[, sum(read_counts), by = time_slice]$V1))],
                                                        by = time_slice][,.(.N), by = c("uniq_label", "time_slice")]
  
  NcountTab       <- dcast(resampleReads, uniq_label ~ time_slice, value.var = "N")  %>%
    mutate(across(-uniq_label, ~if_else(is.na(.), 0, 1))) %>%
    right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label)), by = "uniq_label")
  
  NcountTab2 <- as_tibble(resampleReads) %>% group_by(uniq_label, time_slice) %>% summarise(n_count = sum(N)) %>% arrange(desc(time_slice)) %>%
    pivot_wider(values_from = n_count, names_from = time_slice) %>%
    right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label))) %>%
    arrange(uniq_label) %>% ungroup()
  
  reapTab         <- rbind(
    
    apply(t(apply(NcountTab[,-1], 1, function(x) {
      diff(as.numeric(x))
    })), 2, function(y) sum(y<0, na.rm = T)),
    
    apply(t(apply(NcountTab[,-1], 1, function(x) {
      sapply(1:length(diff(as.numeric(x))), function(y) {
        diff(as.numeric(x))[y]<0 & all(diff(as.numeric(x))[-c(1:y)] == 0)
      })})), 2, function(z) sum(z, na.rm = T)),
    
    do.call("rbind", lapply((ncol(NcountTab[,-1])-1):1, function(x) {                     ## reappear after x time slices
      apply(
        t(
          apply(NcountTab[,-1], 1, function(y) {                                          ## specific sample y
            sapply(1:length(diff(as.numeric(y))), function(z) {                           ## specific time slice z
              ifelse(diff(as.numeric(y))[z]<0 && any(diff(as.numeric(y))[-c(1:z)]==1) &&
                       min(which(diff(as.numeric(y))[-c(1:z)]==1))==x, TRUE, FALSE)
            })
          })
        ),
        2, sum, na.rm = T)
    }))
    
  )
  
  obs <- apply(reapTab, 2, function(x) 1-(x[2]/x[1]))
  # opar <- par(mfrow = c(2,1))
  # plot(rev(m_slices[1:length(obs)]), obs, type= 'o', xlim = rev(range(m_slices[1:length(obs)])), xaxt = "n", xlab = "")
  # axis(1, at = rev(m_slices[1:length(obs)]))
  # abline(v = rev(m_slices[1:length(obs)])[c(TRUE, FALSE)], lty = 3)
  
  exp <- apply(sapply(2:ncol(reapTab), function(x)
    sapply(1:ncol(reapTab), function(y) sum(reapTab[y+2,1:x])/sum(reapTab[1,1:x]))
  ),
  2, function(z) cumsum(rev(z)))[ncol(reapTab):1,]
  
  # invisible(apply(exp[-1,], 2, function(x) lines(rev(m_slices[2:(length(obs))]), x, col = "grey70")))
  # lines(rev(m_slices[2:length(obs)]), exp[-1,4], lwd = 2)
  # plot(rev(m_slices[2:length(obs)]), exp[-1,4]- obs[-length(obs)], type= "o")
  # abline(h = c(0, 0.05))
  
  out <- apply(exp[-1,], 2, function(x) x - obs[-length(obs)])
  # plot(rev(m_slices[2:length(obs)]), out[,4], type= 'o', xlim = rev(range(m_slices[1:length(obs)])), xaxt = "n", xlab = "")  
  # axis(1, at = rev(m_slices[1:length(obs)]))
  # abline(h = 0)
  # abline(v = rev(m_slices[1:length(obs)])[c(TRUE, FALSE)], lty = 3)
  
  list(
    NcountTab2
    ,
    out)
}, mc.cores = 5)

#
save(smpl, file = "output/smpl1_corrected.rda")

###############################################################################
smpl <- parallel::mclapply(1:100, function(y) {
  
  resampleList <- df_long_for_rar[mean_age<=28000,
  ][,time_slice := as.integer(cut(mean_age, t_slices))
  ][,time_slice := abs(time_slice - max(.SD[,"time_slice"]))][order(time_slice)]
  
  resampleCores   <- resampleList[, .SD[core %in% unique(core)[
    sample(1:length(unique(core)), min(resampleList[, length(unique(core)), by = time_slice]$V1))
  ]] , by = time_slice]
  
  resampleSamples <- resampleCores[,.SD[sample_id %in% unique(sample_id)[
    sample(1:length(unique(sample_id)), min(resampleCores[, length(unique(sample_id)), by = time_slice]$V1))
  ]], by = time_slice]
  
  resampleReads   <- resampleSamples[,num:=1:.N
  ][, cbind(.SD, dup = 1:read_counts), by = "num"
  ][,read_counts := 1,][,c("uniq_label", "time_slice")][,
                                                        .SD[sample(1:.N, min(resampleSamples[, sum(read_counts), by = time_slice]$V1))],
                                                        by = time_slice][,.(.N), by = c("uniq_label", "time_slice")]
  
  NcountTab       <- dcast(resampleReads, uniq_label ~ time_slice, value.var = "N")  %>%
    mutate(across(-uniq_label, ~if_else(is.na(.), 0, 1))) %>%
    right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label)), by = "uniq_label")
  
  NcountTab2 <- as_tibble(resampleReads) %>% group_by(uniq_label, time_slice) %>% summarise(n_count = sum(N)) %>% arrange(desc(time_slice)) %>%
    pivot_wider(values_from = n_count, names_from = time_slice) %>%
    right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label))) %>%
    arrange(uniq_label) %>% ungroup()
  
  reapTab         <- rbind(
    
    apply(t(apply(NcountTab[,-1], 1, function(x) {
      diff(as.numeric(x))
    })), 2, function(y) sum(y<0, na.rm = T)),
    
    apply(t(apply(NcountTab[,-1], 1, function(x) {
      sapply(1:length(diff(as.numeric(x))), function(y) {
        diff(as.numeric(x))[y]<0 & all(diff(as.numeric(x))[-c(1:y)] == 0)
      })})), 2, function(z) sum(z, na.rm = T)),
    
    do.call("rbind", lapply((ncol(NcountTab[,-1])-1):1, function(x) {                     ## reappear after x time slices
      apply(
        t(
          apply(NcountTab[,-1], 1, function(y) {                                          ## specific sample y
            sapply(1:length(diff(as.numeric(y))), function(z) {                           ## specific time slice z
              ifelse(diff(as.numeric(y))[z]<0 && any(diff(as.numeric(y))[-c(1:z)]==1) &&
                       min(which(diff(as.numeric(y))[-c(1:z)]==1))==x, TRUE, FALSE)
            })
          })
        ),
        2, sum, na.rm = T)
    }))
    
  )
  
  obs <- apply(reapTab, 2, function(x) 1-(x[2]/x[1]))
  # opar <- par(mfrow = c(2,1))
  # plot(rev(m_slices[1:length(obs)]), obs, type= 'o', xlim = rev(range(m_slices[1:length(obs)])), xaxt = "n", xlab = "")
  # axis(1, at = rev(m_slices[1:length(obs)]))
  # abline(v = rev(m_slices[1:length(obs)])[c(TRUE, FALSE)], lty = 3)
  
  exp <- apply(sapply(2:ncol(reapTab), function(x)
    sapply(1:ncol(reapTab), function(y) sum(reapTab[y+2,1:x])/sum(reapTab[1,1:x]))
  ),
  2, function(z) cumsum(rev(z)))[ncol(reapTab):1,]
  
  # invisible(apply(exp[-1,], 2, function(x) lines(rev(m_slices[2:(length(obs))]), x, col = "grey70")))
  # lines(rev(m_slices[2:length(obs)]), exp[-1,4], lwd = 2)
  # plot(rev(m_slices[2:length(obs)]), exp[-1,4]- obs[-length(obs)], type= "o")
  # abline(h = c(0, 0.05))
  
  out <- apply(exp[-1,], 2, function(x) x - obs[-length(obs)])
  # plot(rev(m_slices[2:length(obs)]), out[,4], type= 'o', xlim = rev(range(m_slices[1:length(obs)])), xaxt = "n", xlab = "")  
  # axis(1, at = rev(m_slices[1:length(obs)]))
  # abline(h = 0)
  # abline(v = rev(m_slices[1:length(obs)])[c(TRUE, FALSE)], lty = 3)
  
  list(
    NcountTab2
    ,
    out)
}, mc.cores = 5)

#
save(smpl, file = "output/smpl2_corrected.rda")

###############################################################################
smpl <- parallel::mclapply(1:100, function(y) {
  
  resampleList <- df_long_for_rar[mean_age<=28000,
  ][,time_slice := as.integer(cut(mean_age, t_slices))
  ][,time_slice := abs(time_slice - max(.SD[,"time_slice"]))][order(time_slice)]
  
  resampleCores   <- resampleList[, .SD[core %in% unique(core)[
    sample(1:length(unique(core)), min(resampleList[, length(unique(core)), by = time_slice]$V1))
  ]] , by = time_slice]
  
  resampleSamples <- resampleCores[,.SD[sample_id %in% unique(sample_id)[
    sample(1:length(unique(sample_id)), min(resampleCores[, length(unique(sample_id)), by = time_slice]$V1))
  ]], by = time_slice]
  
  resampleReads   <- resampleSamples[,num:=1:.N
  ][, cbind(.SD, dup = 1:read_counts), by = "num"
  ][,read_counts := 1,][,c("uniq_label", "time_slice")][,
                                                        .SD[sample(1:.N, min(resampleSamples[, sum(read_counts), by = time_slice]$V1))],
                                                        by = time_slice][,.(.N), by = c("uniq_label", "time_slice")]
  
  NcountTab       <- dcast(resampleReads, uniq_label ~ time_slice, value.var = "N")  %>%
    mutate(across(-uniq_label, ~if_else(is.na(.), 0, 1))) %>%
    right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label)), by = "uniq_label")
  
  NcountTab2 <- as_tibble(resampleReads) %>% group_by(uniq_label, time_slice) %>% summarise(n_count = sum(N)) %>% arrange(desc(time_slice)) %>%
    pivot_wider(values_from = n_count, names_from = time_slice) %>%
    right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label))) %>%
    arrange(uniq_label) %>% ungroup()
  
  reapTab         <- rbind(
    
    apply(t(apply(NcountTab[,-1], 1, function(x) {
      diff(as.numeric(x))
    })), 2, function(y) sum(y<0, na.rm = T)),
    
    apply(t(apply(NcountTab[,-1], 1, function(x) {
      sapply(1:length(diff(as.numeric(x))), function(y) {
        diff(as.numeric(x))[y]<0 & all(diff(as.numeric(x))[-c(1:y)] == 0)
      })})), 2, function(z) sum(z, na.rm = T)),
    
    do.call("rbind", lapply((ncol(NcountTab[,-1])-1):1, function(x) {                     ## reappear after x time slices
      apply(
        t(
          apply(NcountTab[,-1], 1, function(y) {                                          ## specific sample y
            sapply(1:length(diff(as.numeric(y))), function(z) {                           ## specific time slice z
              ifelse(diff(as.numeric(y))[z]<0 && any(diff(as.numeric(y))[-c(1:z)]==1) &&
                       min(which(diff(as.numeric(y))[-c(1:z)]==1))==x, TRUE, FALSE)
            })
          })
        ),
        2, sum, na.rm = T)
    }))
    
  )
  
  obs <- apply(reapTab, 2, function(x) 1-(x[2]/x[1]))
  # opar <- par(mfrow = c(2,1))
  # plot(rev(m_slices[1:length(obs)]), obs, type= 'o', xlim = rev(range(m_slices[1:length(obs)])), xaxt = "n", xlab = "")
  # axis(1, at = rev(m_slices[1:length(obs)]))
  # abline(v = rev(m_slices[1:length(obs)])[c(TRUE, FALSE)], lty = 3)
  
  exp <- apply(sapply(2:ncol(reapTab), function(x)
    sapply(1:ncol(reapTab), function(y) sum(reapTab[y+2,1:x])/sum(reapTab[1,1:x]))
  ),
  2, function(z) cumsum(rev(z)))[ncol(reapTab):1,]
  
  # invisible(apply(exp[-1,], 2, function(x) lines(rev(m_slices[2:(length(obs))]), x, col = "grey70")))
  # lines(rev(m_slices[2:length(obs)]), exp[-1,4], lwd = 2)
  # plot(rev(m_slices[2:length(obs)]), exp[-1,4]- obs[-length(obs)], type= "o")
  # abline(h = c(0, 0.05))
  
  out <- apply(exp[-1,], 2, function(x) x - obs[-length(obs)])
  # plot(rev(m_slices[2:length(obs)]), out[,4], type= 'o', xlim = rev(range(m_slices[1:length(obs)])), xaxt = "n", xlab = "")  
  # axis(1, at = rev(m_slices[1:length(obs)]))
  # abline(h = 0)
  # abline(v = rev(m_slices[1:length(obs)])[c(TRUE, FALSE)], lty = 3)
  
  list(
    NcountTab2
    ,
    out)
}, mc.cores = 5)

#
save(smpl, file = "output/smpl3_corrected.rda")

###############################################################################
smpl <- parallel::mclapply(1:100, function(y) {
  
  resampleList <- df_long_for_rar[mean_age<=28000,
  ][,time_slice := as.integer(cut(mean_age, t_slices))
  ][,time_slice := abs(time_slice - max(.SD[,"time_slice"]))][order(time_slice)]
  
  resampleCores   <- resampleList[, .SD[core %in% unique(core)[
    sample(1:length(unique(core)), min(resampleList[, length(unique(core)), by = time_slice]$V1))
  ]] , by = time_slice]
  
  resampleSamples <- resampleCores[,.SD[sample_id %in% unique(sample_id)[
    sample(1:length(unique(sample_id)), min(resampleCores[, length(unique(sample_id)), by = time_slice]$V1))
  ]], by = time_slice]
  
  resampleReads   <- resampleSamples[,num:=1:.N
  ][, cbind(.SD, dup = 1:read_counts), by = "num"
  ][,read_counts := 1,][,c("uniq_label", "time_slice")][,
                                                        .SD[sample(1:.N, min(resampleSamples[, sum(read_counts), by = time_slice]$V1))],
                                                        by = time_slice][,.(.N), by = c("uniq_label", "time_slice")]
  
  NcountTab       <- dcast(resampleReads, uniq_label ~ time_slice, value.var = "N")  %>%
    mutate(across(-uniq_label, ~if_else(is.na(.), 0, 1))) %>%
    right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label)), by = "uniq_label")
  
  NcountTab2 <- as_tibble(resampleReads) %>% group_by(uniq_label, time_slice) %>% summarise(n_count = sum(N)) %>% arrange(desc(time_slice)) %>%
    pivot_wider(values_from = n_count, names_from = time_slice) %>%
    right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label))) %>%
    arrange(uniq_label) %>% ungroup()
  
  reapTab         <- rbind(
    
    apply(t(apply(NcountTab[,-1], 1, function(x) {
      diff(as.numeric(x))
    })), 2, function(y) sum(y<0, na.rm = T)),
    
    apply(t(apply(NcountTab[,-1], 1, function(x) {
      sapply(1:length(diff(as.numeric(x))), function(y) {
        diff(as.numeric(x))[y]<0 & all(diff(as.numeric(x))[-c(1:y)] == 0)
      })})), 2, function(z) sum(z, na.rm = T)),
    
    do.call("rbind", lapply((ncol(NcountTab[,-1])-1):1, function(x) {                     ## reappear after x time slices
      apply(
        t(
          apply(NcountTab[,-1], 1, function(y) {                                          ## specific sample y
            sapply(1:length(diff(as.numeric(y))), function(z) {                           ## specific time slice z
              ifelse(diff(as.numeric(y))[z]<0 && any(diff(as.numeric(y))[-c(1:z)]==1) &&
                       min(which(diff(as.numeric(y))[-c(1:z)]==1))==x, TRUE, FALSE)
            })
          })
        ),
        2, sum, na.rm = T)
    }))
    
  )
  
  obs <- apply(reapTab, 2, function(x) 1-(x[2]/x[1]))
  # opar <- par(mfrow = c(2,1))
  # plot(rev(m_slices[1:length(obs)]), obs, type= 'o', xlim = rev(range(m_slices[1:length(obs)])), xaxt = "n", xlab = "")
  # axis(1, at = rev(m_slices[1:length(obs)]))
  # abline(v = rev(m_slices[1:length(obs)])[c(TRUE, FALSE)], lty = 3)
  
  exp <- apply(sapply(2:ncol(reapTab), function(x)
    sapply(1:ncol(reapTab), function(y) sum(reapTab[y+2,1:x])/sum(reapTab[1,1:x]))
  ),
  2, function(z) cumsum(rev(z)))[ncol(reapTab):1,]
  
  # invisible(apply(exp[-1,], 2, function(x) lines(rev(m_slices[2:(length(obs))]), x, col = "grey70")))
  # lines(rev(m_slices[2:length(obs)]), exp[-1,4], lwd = 2)
  # plot(rev(m_slices[2:length(obs)]), exp[-1,4]- obs[-length(obs)], type= "o")
  # abline(h = c(0, 0.05))
  
  out <- apply(exp[-1,], 2, function(x) x - obs[-length(obs)])
  # plot(rev(m_slices[2:length(obs)]), out[,4], type= 'o', xlim = rev(range(m_slices[1:length(obs)])), xaxt = "n", xlab = "")  
  # axis(1, at = rev(m_slices[1:length(obs)]))
  # abline(h = 0)
  # abline(v = rev(m_slices[1:length(obs)])[c(TRUE, FALSE)], lty = 3)
  
  list(
    NcountTab2
    ,
    out)
}, mc.cores = 5)

#
save(smpl, file = "output/smpl4_corrected.rda")

###############################################################################
smpl <- parallel::mclapply(1:100, function(y) {
  
  resampleList <- df_long_for_rar[mean_age<=28000,
  ][,time_slice := as.integer(cut(mean_age, t_slices))
  ][,time_slice := abs(time_slice - max(.SD[,"time_slice"]))][order(time_slice)]
  
  resampleCores   <- resampleList[, .SD[core %in% unique(core)[
    sample(1:length(unique(core)), min(resampleList[, length(unique(core)), by = time_slice]$V1))
  ]] , by = time_slice]
  
  resampleSamples <- resampleCores[,.SD[sample_id %in% unique(sample_id)[
    sample(1:length(unique(sample_id)), min(resampleCores[, length(unique(sample_id)), by = time_slice]$V1))
  ]], by = time_slice]
  
  resampleReads   <- resampleSamples[,num:=1:.N
  ][, cbind(.SD, dup = 1:read_counts), by = "num"
  ][,read_counts := 1,][,c("uniq_label", "time_slice")][,
                                                        .SD[sample(1:.N, min(resampleSamples[, sum(read_counts), by = time_slice]$V1))],
                                                        by = time_slice][,.(.N), by = c("uniq_label", "time_slice")]
  
  NcountTab       <- dcast(resampleReads, uniq_label ~ time_slice, value.var = "N")  %>%
    mutate(across(-uniq_label, ~if_else(is.na(.), 0, 1))) %>%
    right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label)), by = "uniq_label")
  
  NcountTab2 <- as_tibble(resampleReads) %>% group_by(uniq_label, time_slice) %>% summarise(n_count = sum(N)) %>% arrange(desc(time_slice)) %>%
    pivot_wider(values_from = n_count, names_from = time_slice) %>%
    right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label))) %>%
    arrange(uniq_label) %>% ungroup()
  
  reapTab         <- rbind(
    
    apply(t(apply(NcountTab[,-1], 1, function(x) {
      diff(as.numeric(x))
    })), 2, function(y) sum(y<0, na.rm = T)),
    
    apply(t(apply(NcountTab[,-1], 1, function(x) {
      sapply(1:length(diff(as.numeric(x))), function(y) {
        diff(as.numeric(x))[y]<0 & all(diff(as.numeric(x))[-c(1:y)] == 0)
      })})), 2, function(z) sum(z, na.rm = T)),
    
    do.call("rbind", lapply((ncol(NcountTab[,-1])-1):1, function(x) {                     ## reappear after x time slices
      apply(
        t(
          apply(NcountTab[,-1], 1, function(y) {                                          ## specific sample y
            sapply(1:length(diff(as.numeric(y))), function(z) {                           ## specific time slice z
              ifelse(diff(as.numeric(y))[z]<0 && any(diff(as.numeric(y))[-c(1:z)]==1) &&
                       min(which(diff(as.numeric(y))[-c(1:z)]==1))==x, TRUE, FALSE)
            })
          })
        ),
        2, sum, na.rm = T)
    }))
    
  )
  
  obs <- apply(reapTab, 2, function(x) 1-(x[2]/x[1]))
  # opar <- par(mfrow = c(2,1))
  # plot(rev(m_slices[1:length(obs)]), obs, type= 'o', xlim = rev(range(m_slices[1:length(obs)])), xaxt = "n", xlab = "")
  # axis(1, at = rev(m_slices[1:length(obs)]))
  # abline(v = rev(m_slices[1:length(obs)])[c(TRUE, FALSE)], lty = 3)
  
  exp <- apply(sapply(2:ncol(reapTab), function(x)
    sapply(1:ncol(reapTab), function(y) sum(reapTab[y+2,1:x])/sum(reapTab[1,1:x]))
  ),
  2, function(z) cumsum(rev(z)))[ncol(reapTab):1,]
  
  # invisible(apply(exp[-1,], 2, function(x) lines(rev(m_slices[2:(length(obs))]), x, col = "grey70")))
  # lines(rev(m_slices[2:length(obs)]), exp[-1,4], lwd = 2)
  # plot(rev(m_slices[2:length(obs)]), exp[-1,4]- obs[-length(obs)], type= "o")
  # abline(h = c(0, 0.05))
  
  out <- apply(exp[-1,], 2, function(x) x - obs[-length(obs)])
  # plot(rev(m_slices[2:length(obs)]), out[,4], type= 'o', xlim = rev(range(m_slices[1:length(obs)])), xaxt = "n", xlab = "")  
  # axis(1, at = rev(m_slices[1:length(obs)]))
  # abline(h = 0)
  # abline(v = rev(m_slices[1:length(obs)])[c(TRUE, FALSE)], lty = 3)
  
  list(
    NcountTab2
    ,
    out)
}, mc.cores = 5)

#
save(smpl, file = "output/smpl5_corrected.rda")

###############################################################################
smpl <- parallel::mclapply(1:100, function(y) {
  
  resampleList <- df_long_for_rar[mean_age<=28000,
  ][,time_slice := as.integer(cut(mean_age, t_slices))
  ][,time_slice := abs(time_slice - max(.SD[,"time_slice"]))][order(time_slice)]
  
  resampleCores   <- resampleList[, .SD[core %in% unique(core)[
    sample(1:length(unique(core)), min(resampleList[, length(unique(core)), by = time_slice]$V1))
  ]] , by = time_slice]
  
  resampleSamples <- resampleCores[,.SD[sample_id %in% unique(sample_id)[
    sample(1:length(unique(sample_id)), min(resampleCores[, length(unique(sample_id)), by = time_slice]$V1))
  ]], by = time_slice]
  
  resampleReads   <- resampleSamples[,num:=1:.N
  ][, cbind(.SD, dup = 1:read_counts), by = "num"
  ][,read_counts := 1,][,c("uniq_label", "time_slice")][,
                                                        .SD[sample(1:.N, min(resampleSamples[, sum(read_counts), by = time_slice]$V1))],
                                                        by = time_slice][,.(.N), by = c("uniq_label", "time_slice")]
  
  NcountTab       <- dcast(resampleReads, uniq_label ~ time_slice, value.var = "N")  %>%
    mutate(across(-uniq_label, ~if_else(is.na(.), 0, 1))) %>%
    right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label)), by = "uniq_label")
  
  NcountTab2 <- as_tibble(resampleReads) %>% group_by(uniq_label, time_slice) %>% summarise(n_count = sum(N)) %>% arrange(desc(time_slice)) %>%
    pivot_wider(values_from = n_count, names_from = time_slice) %>%
    right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label))) %>%
    arrange(uniq_label) %>% ungroup()
  
  reapTab         <- rbind(
    
    apply(t(apply(NcountTab[,-1], 1, function(x) {
      diff(as.numeric(x))
    })), 2, function(y) sum(y<0, na.rm = T)),
    
    apply(t(apply(NcountTab[,-1], 1, function(x) {
      sapply(1:length(diff(as.numeric(x))), function(y) {
        diff(as.numeric(x))[y]<0 & all(diff(as.numeric(x))[-c(1:y)] == 0)
      })})), 2, function(z) sum(z, na.rm = T)),
    
    do.call("rbind", lapply((ncol(NcountTab[,-1])-1):1, function(x) {                     ## reappear after x time slices
      apply(
        t(
          apply(NcountTab[,-1], 1, function(y) {                                          ## specific sample y
            sapply(1:length(diff(as.numeric(y))), function(z) {                           ## specific time slice z
              ifelse(diff(as.numeric(y))[z]<0 && any(diff(as.numeric(y))[-c(1:z)]==1) &&
                       min(which(diff(as.numeric(y))[-c(1:z)]==1))==x, TRUE, FALSE)
            })
          })
        ),
        2, sum, na.rm = T)
    }))
    
  )
  
  obs <- apply(reapTab, 2, function(x) 1-(x[2]/x[1]))
  # opar <- par(mfrow = c(2,1))
  # plot(rev(m_slices[1:length(obs)]), obs, type= 'o', xlim = rev(range(m_slices[1:length(obs)])), xaxt = "n", xlab = "")
  # axis(1, at = rev(m_slices[1:length(obs)]))
  # abline(v = rev(m_slices[1:length(obs)])[c(TRUE, FALSE)], lty = 3)
  
  exp <- apply(sapply(2:ncol(reapTab), function(x)
    sapply(1:ncol(reapTab), function(y) sum(reapTab[y+2,1:x])/sum(reapTab[1,1:x]))
  ),
  2, function(z) cumsum(rev(z)))[ncol(reapTab):1,]
  
  # invisible(apply(exp[-1,], 2, function(x) lines(rev(m_slices[2:(length(obs))]), x, col = "grey70")))
  # lines(rev(m_slices[2:length(obs)]), exp[-1,4], lwd = 2)
  # plot(rev(m_slices[2:length(obs)]), exp[-1,4]- obs[-length(obs)], type= "o")
  # abline(h = c(0, 0.05))
  
  out <- apply(exp[-1,], 2, function(x) x - obs[-length(obs)])
  # plot(rev(m_slices[2:length(obs)]), out[,4], type= 'o', xlim = rev(range(m_slices[1:length(obs)])), xaxt = "n", xlab = "")  
  # axis(1, at = rev(m_slices[1:length(obs)]))
  # abline(h = 0)
  # abline(v = rev(m_slices[1:length(obs)])[c(TRUE, FALSE)], lty = 3)
  
  list(
    NcountTab2
    ,
    out)
}, mc.cores = 5)

#
save(smpl, file = "output/smpl6_corrected.rda")

###############################################################################
smpl <- parallel::mclapply(1:100, function(y) {
  
  resampleList <- df_long_for_rar[mean_age<=28000,
  ][,time_slice := as.integer(cut(mean_age, t_slices))
  ][,time_slice := abs(time_slice - max(.SD[,"time_slice"]))][order(time_slice)]
  
  resampleCores   <- resampleList[, .SD[core %in% unique(core)[
    sample(1:length(unique(core)), min(resampleList[, length(unique(core)), by = time_slice]$V1))
  ]] , by = time_slice]
  
  resampleSamples <- resampleCores[,.SD[sample_id %in% unique(sample_id)[
    sample(1:length(unique(sample_id)), min(resampleCores[, length(unique(sample_id)), by = time_slice]$V1))
  ]], by = time_slice]
  
  resampleReads   <- resampleSamples[,num:=1:.N
  ][, cbind(.SD, dup = 1:read_counts), by = "num"
  ][,read_counts := 1,][,c("uniq_label", "time_slice")][,
                                                        .SD[sample(1:.N, min(resampleSamples[, sum(read_counts), by = time_slice]$V1))],
                                                        by = time_slice][,.(.N), by = c("uniq_label", "time_slice")]
  
  NcountTab       <- dcast(resampleReads, uniq_label ~ time_slice, value.var = "N")  %>%
    mutate(across(-uniq_label, ~if_else(is.na(.), 0, 1))) %>%
    right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label)), by = "uniq_label")
  
  NcountTab2 <- as_tibble(resampleReads) %>% group_by(uniq_label, time_slice) %>% summarise(n_count = sum(N)) %>% arrange(desc(time_slice)) %>%
    pivot_wider(values_from = n_count, names_from = time_slice) %>%
    right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label))) %>%
    arrange(uniq_label) %>% ungroup()
  
  reapTab         <- rbind(
    
    apply(t(apply(NcountTab[,-1], 1, function(x) {
      diff(as.numeric(x))
    })), 2, function(y) sum(y<0, na.rm = T)),
    
    apply(t(apply(NcountTab[,-1], 1, function(x) {
      sapply(1:length(diff(as.numeric(x))), function(y) {
        diff(as.numeric(x))[y]<0 & all(diff(as.numeric(x))[-c(1:y)] == 0)
      })})), 2, function(z) sum(z, na.rm = T)),
    
    do.call("rbind", lapply((ncol(NcountTab[,-1])-1):1, function(x) {                     ## reappear after x time slices
      apply(
        t(
          apply(NcountTab[,-1], 1, function(y) {                                          ## specific sample y
            sapply(1:length(diff(as.numeric(y))), function(z) {                           ## specific time slice z
              ifelse(diff(as.numeric(y))[z]<0 && any(diff(as.numeric(y))[-c(1:z)]==1) &&
                       min(which(diff(as.numeric(y))[-c(1:z)]==1))==x, TRUE, FALSE)
            })
          })
        ),
        2, sum, na.rm = T)
    }))
    
  )
  
  obs <- apply(reapTab, 2, function(x) 1-(x[2]/x[1]))
  # opar <- par(mfrow = c(2,1))
  # plot(rev(m_slices[1:length(obs)]), obs, type= 'o', xlim = rev(range(m_slices[1:length(obs)])), xaxt = "n", xlab = "")
  # axis(1, at = rev(m_slices[1:length(obs)]))
  # abline(v = rev(m_slices[1:length(obs)])[c(TRUE, FALSE)], lty = 3)
  
  exp <- apply(sapply(2:ncol(reapTab), function(x)
    sapply(1:ncol(reapTab), function(y) sum(reapTab[y+2,1:x])/sum(reapTab[1,1:x]))
  ),
  2, function(z) cumsum(rev(z)))[ncol(reapTab):1,]
  
  # invisible(apply(exp[-1,], 2, function(x) lines(rev(m_slices[2:(length(obs))]), x, col = "grey70")))
  # lines(rev(m_slices[2:length(obs)]), exp[-1,4], lwd = 2)
  # plot(rev(m_slices[2:length(obs)]), exp[-1,4]- obs[-length(obs)], type= "o")
  # abline(h = c(0, 0.05))
  
  out <- apply(exp[-1,], 2, function(x) x - obs[-length(obs)])
  # plot(rev(m_slices[2:length(obs)]), out[,4], type= 'o', xlim = rev(range(m_slices[1:length(obs)])), xaxt = "n", xlab = "")  
  # axis(1, at = rev(m_slices[1:length(obs)]))
  # abline(h = 0)
  # abline(v = rev(m_slices[1:length(obs)])[c(TRUE, FALSE)], lty = 3)
  
  list(
    NcountTab2
    ,
    out)
}, mc.cores = 5)

#
save(smpl, file = "output/smpl7_corrected.rda")

###############################################################################
smpl <- parallel::mclapply(1:100, function(y) {
  
  resampleList <- df_long_for_rar[mean_age<=28000,
  ][,time_slice := as.integer(cut(mean_age, t_slices))
  ][,time_slice := abs(time_slice - max(.SD[,"time_slice"]))][order(time_slice)]
  
  resampleCores   <- resampleList[, .SD[core %in% unique(core)[
    sample(1:length(unique(core)), min(resampleList[, length(unique(core)), by = time_slice]$V1))
  ]] , by = time_slice]
  
  resampleSamples <- resampleCores[,.SD[sample_id %in% unique(sample_id)[
    sample(1:length(unique(sample_id)), min(resampleCores[, length(unique(sample_id)), by = time_slice]$V1))
  ]], by = time_slice]
  
  resampleReads   <- resampleSamples[,num:=1:.N
  ][, cbind(.SD, dup = 1:read_counts), by = "num"
  ][,read_counts := 1,][,c("uniq_label", "time_slice")][,
                                                        .SD[sample(1:.N, min(resampleSamples[, sum(read_counts), by = time_slice]$V1))],
                                                        by = time_slice][,.(.N), by = c("uniq_label", "time_slice")]
  
  NcountTab       <- dcast(resampleReads, uniq_label ~ time_slice, value.var = "N")  %>%
    mutate(across(-uniq_label, ~if_else(is.na(.), 0, 1))) %>%
    right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label)), by = "uniq_label")
  
  NcountTab2 <- as_tibble(resampleReads) %>% group_by(uniq_label, time_slice) %>% summarise(n_count = sum(N)) %>% arrange(desc(time_slice)) %>%
    pivot_wider(values_from = n_count, names_from = time_slice) %>%
    right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label))) %>%
    arrange(uniq_label) %>% ungroup()
  
  reapTab         <- rbind(
    
    apply(t(apply(NcountTab[,-1], 1, function(x) {
      diff(as.numeric(x))
    })), 2, function(y) sum(y<0, na.rm = T)),
    
    apply(t(apply(NcountTab[,-1], 1, function(x) {
      sapply(1:length(diff(as.numeric(x))), function(y) {
        diff(as.numeric(x))[y]<0 & all(diff(as.numeric(x))[-c(1:y)] == 0)
      })})), 2, function(z) sum(z, na.rm = T)),
    
    do.call("rbind", lapply((ncol(NcountTab[,-1])-1):1, function(x) {                     ## reappear after x time slices
      apply(
        t(
          apply(NcountTab[,-1], 1, function(y) {                                          ## specific sample y
            sapply(1:length(diff(as.numeric(y))), function(z) {                           ## specific time slice z
              ifelse(diff(as.numeric(y))[z]<0 && any(diff(as.numeric(y))[-c(1:z)]==1) &&
                       min(which(diff(as.numeric(y))[-c(1:z)]==1))==x, TRUE, FALSE)
            })
          })
        ),
        2, sum, na.rm = T)
    }))
    
  )
  
  obs <- apply(reapTab, 2, function(x) 1-(x[2]/x[1]))
  # opar <- par(mfrow = c(2,1))
  # plot(rev(m_slices[1:length(obs)]), obs, type= 'o', xlim = rev(range(m_slices[1:length(obs)])), xaxt = "n", xlab = "")
  # axis(1, at = rev(m_slices[1:length(obs)]))
  # abline(v = rev(m_slices[1:length(obs)])[c(TRUE, FALSE)], lty = 3)
  
  exp <- apply(sapply(2:ncol(reapTab), function(x)
    sapply(1:ncol(reapTab), function(y) sum(reapTab[y+2,1:x])/sum(reapTab[1,1:x]))
  ),
  2, function(z) cumsum(rev(z)))[ncol(reapTab):1,]
  
  # invisible(apply(exp[-1,], 2, function(x) lines(rev(m_slices[2:(length(obs))]), x, col = "grey70")))
  # lines(rev(m_slices[2:length(obs)]), exp[-1,4], lwd = 2)
  # plot(rev(m_slices[2:length(obs)]), exp[-1,4]- obs[-length(obs)], type= "o")
  # abline(h = c(0, 0.05))
  
  out <- apply(exp[-1,], 2, function(x) x - obs[-length(obs)])
  # plot(rev(m_slices[2:length(obs)]), out[,4], type= 'o', xlim = rev(range(m_slices[1:length(obs)])), xaxt = "n", xlab = "")  
  # axis(1, at = rev(m_slices[1:length(obs)]))
  # abline(h = 0)
  # abline(v = rev(m_slices[1:length(obs)])[c(TRUE, FALSE)], lty = 3)
  
  list(
    NcountTab2
    ,
    out)
}, mc.cores = 5)

#
save(smpl, file = "output/smpl8_corrected.rda")

###############################################################################
smpl <- parallel::mclapply(1:100, function(y) {
  
  resampleList <- df_long_for_rar[mean_age<=28000,
  ][,time_slice := as.integer(cut(mean_age, t_slices))
  ][,time_slice := abs(time_slice - max(.SD[,"time_slice"]))][order(time_slice)]
  
  resampleCores   <- resampleList[, .SD[core %in% unique(core)[
    sample(1:length(unique(core)), min(resampleList[, length(unique(core)), by = time_slice]$V1))
  ]] , by = time_slice]
  
  resampleSamples <- resampleCores[,.SD[sample_id %in% unique(sample_id)[
    sample(1:length(unique(sample_id)), min(resampleCores[, length(unique(sample_id)), by = time_slice]$V1))
  ]], by = time_slice]
  
  resampleReads   <- resampleSamples[,num:=1:.N
  ][, cbind(.SD, dup = 1:read_counts), by = "num"
  ][,read_counts := 1,][,c("uniq_label", "time_slice")][,
                                                        .SD[sample(1:.N, min(resampleSamples[, sum(read_counts), by = time_slice]$V1))],
                                                        by = time_slice][,.(.N), by = c("uniq_label", "time_slice")]
  
  NcountTab       <- dcast(resampleReads, uniq_label ~ time_slice, value.var = "N")  %>%
    mutate(across(-uniq_label, ~if_else(is.na(.), 0, 1))) %>%
    right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label)), by = "uniq_label")
  
  NcountTab2 <- as_tibble(resampleReads) %>% group_by(uniq_label, time_slice) %>% summarise(n_count = sum(N)) %>% arrange(desc(time_slice)) %>%
    pivot_wider(values_from = n_count, names_from = time_slice) %>%
    right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label))) %>%
    arrange(uniq_label) %>% ungroup()
  
  reapTab         <- rbind(
    
    apply(t(apply(NcountTab[,-1], 1, function(x) {
      diff(as.numeric(x))
    })), 2, function(y) sum(y<0, na.rm = T)),
    
    apply(t(apply(NcountTab[,-1], 1, function(x) {
      sapply(1:length(diff(as.numeric(x))), function(y) {
        diff(as.numeric(x))[y]<0 & all(diff(as.numeric(x))[-c(1:y)] == 0)
      })})), 2, function(z) sum(z, na.rm = T)),
    
    do.call("rbind", lapply((ncol(NcountTab[,-1])-1):1, function(x) {                     ## reappear after x time slices
      apply(
        t(
          apply(NcountTab[,-1], 1, function(y) {                                          ## specific sample y
            sapply(1:length(diff(as.numeric(y))), function(z) {                           ## specific time slice z
              ifelse(diff(as.numeric(y))[z]<0 && any(diff(as.numeric(y))[-c(1:z)]==1) &&
                       min(which(diff(as.numeric(y))[-c(1:z)]==1))==x, TRUE, FALSE)
            })
          })
        ),
        2, sum, na.rm = T)
    }))
    
  )
  
  obs <- apply(reapTab, 2, function(x) 1-(x[2]/x[1]))
  # opar <- par(mfrow = c(2,1))
  # plot(rev(m_slices[1:length(obs)]), obs, type= 'o', xlim = rev(range(m_slices[1:length(obs)])), xaxt = "n", xlab = "")
  # axis(1, at = rev(m_slices[1:length(obs)]))
  # abline(v = rev(m_slices[1:length(obs)])[c(TRUE, FALSE)], lty = 3)
  
  exp <- apply(sapply(2:ncol(reapTab), function(x)
    sapply(1:ncol(reapTab), function(y) sum(reapTab[y+2,1:x])/sum(reapTab[1,1:x]))
  ),
  2, function(z) cumsum(rev(z)))[ncol(reapTab):1,]
  
  # invisible(apply(exp[-1,], 2, function(x) lines(rev(m_slices[2:(length(obs))]), x, col = "grey70")))
  # lines(rev(m_slices[2:length(obs)]), exp[-1,4], lwd = 2)
  # plot(rev(m_slices[2:length(obs)]), exp[-1,4]- obs[-length(obs)], type= "o")
  # abline(h = c(0, 0.05))
  
  out <- apply(exp[-1,], 2, function(x) x - obs[-length(obs)])
  # plot(rev(m_slices[2:length(obs)]), out[,4], type= 'o', xlim = rev(range(m_slices[1:length(obs)])), xaxt = "n", xlab = "")  
  # axis(1, at = rev(m_slices[1:length(obs)]))
  # abline(h = 0)
  # abline(v = rev(m_slices[1:length(obs)])[c(TRUE, FALSE)], lty = 3)
  
  list(
    NcountTab2
    ,
    out)
}, mc.cores = 5)

#
save(smpl, file = "output/smpl9_corrected.rda")

###############################################################################
smpl <- parallel::mclapply(1:100, function(y) {
  
  resampleList <- df_long_for_rar[mean_age<=28000,
  ][,time_slice := as.integer(cut(mean_age, t_slices))
  ][,time_slice := abs(time_slice - max(.SD[,"time_slice"]))][order(time_slice)]
  
  resampleCores   <- resampleList[, .SD[core %in% unique(core)[
    sample(1:length(unique(core)), min(resampleList[, length(unique(core)), by = time_slice]$V1))
  ]] , by = time_slice]
  
  resampleSamples <- resampleCores[,.SD[sample_id %in% unique(sample_id)[
    sample(1:length(unique(sample_id)), min(resampleCores[, length(unique(sample_id)), by = time_slice]$V1))
  ]], by = time_slice]
  
  resampleReads   <- resampleSamples[,num:=1:.N
  ][, cbind(.SD, dup = 1:read_counts), by = "num"
  ][,read_counts := 1,][,c("uniq_label", "time_slice")][,
                                                        .SD[sample(1:.N, min(resampleSamples[, sum(read_counts), by = time_slice]$V1))],
                                                        by = time_slice][,.(.N), by = c("uniq_label", "time_slice")]
  
  NcountTab       <- dcast(resampleReads, uniq_label ~ time_slice, value.var = "N")  %>%
    mutate(across(-uniq_label, ~if_else(is.na(.), 0, 1))) %>%
    right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label)), by = "uniq_label")
  
  NcountTab2 <- as_tibble(resampleReads) %>% group_by(uniq_label, time_slice) %>% summarise(n_count = sum(N)) %>% arrange(desc(time_slice)) %>%
    pivot_wider(values_from = n_count, names_from = time_slice) %>%
    right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label))) %>%
    arrange(uniq_label) %>% ungroup()
  
  reapTab         <- rbind(
    
    apply(t(apply(NcountTab[,-1], 1, function(x) {
      diff(as.numeric(x))
    })), 2, function(y) sum(y<0, na.rm = T)),
    
    apply(t(apply(NcountTab[,-1], 1, function(x) {
      sapply(1:length(diff(as.numeric(x))), function(y) {
        diff(as.numeric(x))[y]<0 & all(diff(as.numeric(x))[-c(1:y)] == 0)
      })})), 2, function(z) sum(z, na.rm = T)),
    
    do.call("rbind", lapply((ncol(NcountTab[,-1])-1):1, function(x) {                     ## reappear after x time slices
      apply(
        t(
          apply(NcountTab[,-1], 1, function(y) {                                          ## specific sample y
            sapply(1:length(diff(as.numeric(y))), function(z) {                           ## specific time slice z
              ifelse(diff(as.numeric(y))[z]<0 && any(diff(as.numeric(y))[-c(1:z)]==1) &&
                       min(which(diff(as.numeric(y))[-c(1:z)]==1))==x, TRUE, FALSE)
            })
          })
        ),
        2, sum, na.rm = T)
    }))
    
  )
  
  obs <- apply(reapTab, 2, function(x) 1-(x[2]/x[1]))
  # opar <- par(mfrow = c(2,1))
  # plot(rev(m_slices[1:length(obs)]), obs, type= 'o', xlim = rev(range(m_slices[1:length(obs)])), xaxt = "n", xlab = "")
  # axis(1, at = rev(m_slices[1:length(obs)]))
  # abline(v = rev(m_slices[1:length(obs)])[c(TRUE, FALSE)], lty = 3)
  
  exp <- apply(sapply(2:ncol(reapTab), function(x)
    sapply(1:ncol(reapTab), function(y) sum(reapTab[y+2,1:x])/sum(reapTab[1,1:x]))
  ),
  2, function(z) cumsum(rev(z)))[ncol(reapTab):1,]
  
  # invisible(apply(exp[-1,], 2, function(x) lines(rev(m_slices[2:(length(obs))]), x, col = "grey70")))
  # lines(rev(m_slices[2:length(obs)]), exp[-1,4], lwd = 2)
  # plot(rev(m_slices[2:length(obs)]), exp[-1,4]- obs[-length(obs)], type= "o")
  # abline(h = c(0, 0.05))
  
  out <- apply(exp[-1,], 2, function(x) x - obs[-length(obs)])
  # plot(rev(m_slices[2:length(obs)]), out[,4], type= 'o', xlim = rev(range(m_slices[1:length(obs)])), xaxt = "n", xlab = "")  
  # axis(1, at = rev(m_slices[1:length(obs)]))
  # abline(h = 0)
  # abline(v = rev(m_slices[1:length(obs)])[c(TRUE, FALSE)], lty = 3)
  
  list(
    NcountTab2
    ,
    out)
}, mc.cores = 5)

#
save(smpl, file = "output/smpl10_corrected.rda")

#####################################################################
# 4.3 - Save resampled data!
#####################################################################
load("output/smpl1_corrected.rda")
smpl1 <- smpl
load("output/smpl2_corrected.rda")
smpl2 <- smpl
load("output/smpl3_corrected.rda")
smpl3 <- smpl
load("output/smpl4_corrected.rda")
smpl4 <- smpl
load("output/smpl5_corrected.rda")
smpl5 <- smpl
load("output/smpl6_corrected.rda")
smpl6 <- smpl
load("output/smpl7_corrected.rda")
smpl7 <- smpl
load("output/smpl8_corrected.rda")
smpl8 <- smpl
load("output/smpl9_corrected.rda")
smpl9 <- smpl
load("output/smpl10_corrected.rda")
smpl10 <- smpl

# merge all
smpl_raw <- do.call(c, list(smpl1, smpl2, smpl3, smpl4, smpl5, smpl6, smpl7, smpl8, smpl9, smpl10))
# check them
smpl_raw[[1]][[1]]
smpl_raw[[101]][[1]]
smpl_raw[[201]][[1]]
smpl_raw[[301]][[1]]
smpl_raw[[401]][[1]]
smpl_raw[[501]][[1]]
smpl_raw[[601]][[1]]
smpl_raw[[701]][[1]]
smpl_raw[[801]][[1]]
smpl_raw[[901]][[1]]

# Save and use in next script!
save(smpl_raw, file = "output/resampled_data_final_1000iterations.rda")
save(smpl_raw, file = "resampled_data_final_1000iterations.rda")
