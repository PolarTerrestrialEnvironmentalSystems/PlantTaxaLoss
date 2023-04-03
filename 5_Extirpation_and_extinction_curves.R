# TAXALOSS Project - script to plot Figure 2 and supplementary information
# Script from Jeremy Courtin and Simeon Lisovski
# Script last update - 03.04.2023

# Layout of the script: 
# 0 - GENERAL OPTION - SETTING UP SCRIPT
# 1 - Make the abundance table
# 2 - Plots
# 2.1 - Plot for the supplementary material
# 2.2 - Plot for the main manuscript
# 3 - What to report in the text
# 4 - MEGAFAUNA

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
library(data.table)

# Import the resampled data
load("output/resampled_data_final_1000iterations.rda")

resampleAll <- smpl_raw
look <- resampleAll[[1]][[1]]

# Set age timeslice: 
age <- tibble(timeslice = c(13:0), age = c(1000, 3000, 5000, 7000, 9000, 11000, 13000, 15000, 17000, 19000, 21000, 23000, 25000, 27000))
df_long_for_rar <- read_delim("output/df_long_for_rar.csv", delim = ",", col_names = T) %>%
  subset(read_counts > 0)

df_long_for_rar %>% #subset(mean_age < 28001) %>% 
  dplyr::select(mean_age, sample_id) %>%
  #summarise(sum = sum(read_counts)) %>%
  distinct() 

taxa_info <- read_delim("output/final_dataset_before_resampling_family_merge_new.csv", delim = ",", col_names = T) %>%
  dplyr::select(identity, family, genus, species, scientific_name, group, nb_in_group) %>% mutate(uniq = 1:nrow(.)) %>% 
  mutate(level = ifelse(is.na(species), ifelse(is.na(genus), "family", "genus"), "species")) %>%
  mutate(genus = ifelse(is.na(genus), family, genus)) %>% mutate(species = ifelse(is.na(species), paste(genus, uniq, sep = " "), paste(species, uniq, sep = " "))) %>% mutate(identity = ifelse(identity == 0.9, "cand", "1")) %>%
  mutate(key = paste(nb_in_group, identity, family, scientific_name, sep = "_"))

###############################################################################
# 1 - Make the abundance table
###############################################################################
abund_table <- lapply(resampleAll, function(x) x[[1]])
abund_table <- lapply(abund_table, function(x) replace(x, is.na(x), 0))

# Change the info on NcountTab adn save as different outputs
NcountTab <- lapply(abund_table, function(x) x %>% pivot_longer(`13`:`0`) %>% subset(value > 0) %>% 
                      mutate(name = as.numeric(name)) %>% setNames(c("uniq_label", "time_slice", "N")) %>% 
                      arrange(time_slice) %>% as.data.table() %>%
                      dcast(uniq_label ~ time_slice, value.var = "N")  %>%
                      mutate(across(-uniq_label, ~if_else(is.na(.), 0, 1))) %>%
                      subset(grepl("_cand_", uniq_label)) %>% # to set if we want to calculate the extinction rate
                      #subset(grepl("_1_", uniq_label)) %>% # to set if we want to calculate the extirpation rate
                      right_join(tibble(uniq_label = unique(df_long_for_rar$uniq_label)), by = "uniq_label"))

nSAMPLE <- 1000

reapTab = list()
for (i in 1:nSAMPLE) {
  reapTab[[i]] <- rbind(
    apply(t(apply(NcountTab[[i]][,-1], 1, function(x) {
      diff(as.numeric(x))
    })), 2, function(y) sum(y<0, na.rm = T)),
    
    apply(t(apply(NcountTab[[i]][,-1], 1, function(x) {
      sapply(1:length(diff(as.numeric(x))), function(y) {
        diff(as.numeric(x))[y]<0 & all(diff(as.numeric(x))[-c(1:y)] == 0)
      })})), 2, function(z) sum(z, na.rm = T)),
    
    do.call("rbind", lapply((ncol(NcountTab[[i]][,-1])-1):1, function(x) {                     ## reappear after x time slices
      apply(
        t(
          apply(NcountTab[[i]][,-1], 1, function(y) {                                          ## specific sample y
            sapply(1:length(diff(as.numeric(y))), function(z) {                           ## specific time slice z
              ifelse(diff(as.numeric(y))[z]<0 && any(diff(as.numeric(y))[-c(1:z)]==1) &&
                       min(which(diff(as.numeric(y))[-c(1:z)]==1))==x, TRUE, FALSE)
            })
          })
        ),
        2, sum, na.rm = T)
    }))
    
  )
}

obs = list()
for (i in 1:nSAMPLE) {
  obs[[i]] <- rbind(apply(reapTab[[i]], 2, function(x) 1-(x[2]/x[1])))}

exp = list()
for (i in 1:nSAMPLE) {
  exp[[i]] <- rbind(apply(sapply(2:ncol(reapTab[[i]]), function(x) 
  sapply(1:ncol(reapTab[[i]]), function(y) sum(reapTab[[i]][y+2,1:x])/sum(reapTab[[i]][1,1:x]))
),
2, function(z) cumsum(rev(z)))[ncol(reapTab[[i]]):1,])}

out = list() 
for (i in 1:nSAMPLE) {
  out[[i]] <- rbind(apply(exp[[i]][-1,], 2, function(x) x - obs[[i]][-length(obs[[i]])]))}


# Check the maximum number of timeslice a disappeared taxon takes to reappear
max_nb_timeslice <- lapply(reapTab, function(x) x %>% as_tibble() %>% mutate(col = 15:1) %>% subset(col < 14) %>% pivot_longer(cols =  V1:V13) %>% dplyr::select(-name) %>% subset(value > 0) %>% arrange(desc(col)) %>% summarise(max_nb_timeslice = max(col)))
max_nb_timeslice <- do.call(rbind, max_nb_timeslice)
median(max_nb_timeslice$max_nb_timeslice)
# For all (no subset): here median maximum number of timeslice to reappear is 9 timeslice. Which mean 19000 cal. yrs BP

# if no subset
out_all <- out
save(out_all, file = "data_curve_all_extirpation_extinction.RData")

# if subset cand
#out_cand <- out
#save(out_cand, file = "data_curve_cand_extinction.RData")

# if subset db
# out_db <- out
# save(out_db, file = "data_curve_db_extirpation.RData")

###############################################################################
# 2 - Plots
###############################################################################
load("data_curve_all_extirpation_extinction.RData")
load("data_curve_cand_extinction.RData")
load("data_curve_db_extirpation.RData")

#Set which smpl you want to plot
smpl <- out_all
smpl <- out_cand
smpl <- out_db

t_slices <- seq(0, 32000, 2000)
m_slices <- t_slices + median(diff(t_slices))/2

slices <- c(1,2,3,4,5)

###############################################################################
# 2.1 - Plot for the supplementary material
###############################################################################
opar <- par(mfrow = c(length(slices),1), mar = c(3,3,1,1), oma = c(1,1,0,0), las = 1)
for(i in slices) {
  dat <- do.call("rbind", lapply(smpl, function(x) x[,i]))
  matplot(x = rev(m_slices[2:(ncol(smpl[[1]])+1)]), t(dat), lwd = 1, lty = 1, col = "grey90", pch = 16,
          xlab = "", ylab = "", ylim = c(-0.37, 0.49), xlim = c(25000, 1000), xaxt = "n")
  #xlab = "", ylab = "", ylim = c(-0.2, 0.2), xlim = c(25000, 1000), xaxt = "n")
  axis(1, at = m_slices[2:(ncol(smpl[[1]])+1)])
  abline(v = m_slices[2:(ncol(smpl[[1]])+1)], lty = 3)
  segments(x0 = rev(m_slices[2:(ncol(smpl[[1]])+1)]), y0 = apply(dat, 2, function(x) quantile(x, probs = 0.05, na.rm = T)),
           x1 = rev(m_slices[2:(ncol(smpl[[1]])+1)]), y1 = apply(dat, 2, function(x) quantile(x, probs = 0.95, na.rm = T)), lwd = 2)
  points(rev(m_slices[2:(ncol(smpl[[1]])+1)]), apply(dat, 2, function(x) quantile(x, probs = 0.5, na.rm = T)), pch = 16, cex = 1.5)
  abline(h = c(0), lty = c(2))
  # abline(h = c(0.05), lty = c(2))
  mtext(glue::glue("Slice - {c(11:1)[i]}"), 1, line = -1.5)
  text(rev(m_slices[2:(ncol(smpl[[1]])+1)]), 0.45, round(apply(do.call("rbind", lapply(smpl, function(x) x[,i])), 2, function(x) sum(x>0)/length(smpl))*100, 1))
}
par(opar)

# For one slice chosen based on the median value calculated before
dat <- do.call("rbind", lapply(smpl, function(x) x[,3])) # x[,3] is taking the first expectation curve measured from timeslices 25000 to 19000.
dat1 <- do.call("rbind", lapply(smpl, function(x) x[,3])) %>% as_tibble() %>% # x[,3] is taking the first expectation curve measured from timeslices 25000 to 19000.
  setNames(c("25000", "23000", "21000", "19000", "17000", "15000", "13000", "11000", "9000", "7000", "5000", "3000")) %>%
  pivot_longer(cols = c("25000", "23000", "21000", "19000", "17000", "15000", "13000", "11000", "9000", "7000", "5000", "3000"))

# Save the data for the figure
write_delim(dat1, "Figure2_data_save_all.csv", delim = ",", col_names = T)
write_delim(dat1, "Figure2_data_save_cand_extinction.csv", delim = ",", col_names = T)
write_delim(dat1, "Figure2_data_save_db_extirpation.csv", delim = ",", col_names = T)

dat1 <- read_delim("Figure2_data_save_cand_extinction.csv", delim = ",", col_names = T)

###############################################################################
# 2.2 - Plot for the main manuscript
###############################################################################
fig2_1 <- 
  dat1 %>%
  setNames(c("age", "ext_rate")) %>%
  group_by(age) %>% 
  mutate(median_rate = median(ext_rate), q95 = quantile(ext_rate, 0.95), q5 = quantile(ext_rate, 0.05),  age = as.numeric(age), 
         color_point = ifelse(median_rate > 0, ifelse(median_rate > 0.05, "red", "orange"), "blue"), 
         perc_sup0 = sum(ext_rate>0)/1000*100, perc_sup0 = paste(perc_sup0, "%",sep = "")) %>%
  ggplot(aes(x= age)) +
  geom_ribbon(aes(ymin = 0, ymax = q95), fill = "red", alpha=0.1) +
  geom_ribbon(aes(ymin = q5, ymax = 0), fill = "grey90") +
  #  geom_point(aes(y=ext_rate), color = "gray80", size=1) +
  geom_point(aes(y=median_rate), size=3) +
  geom_line(aes(y=median_rate), size=0.1) +
  geom_hline(yintercept=0, linetype = "dashed") +
  geom_text(aes(y = q95+0.01, label  = perc_sup0)) +
  scale_x_reverse(breaks=seq(0, 25000, 2000)) +
  xlab("Age (cal. yrs BP)") +
  ylab("Proportion of loss species relative to the expectation") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
fig2_1

###############################################################################
# 3 - What to report in the text
###############################################################################
dat1 %>%
  setNames(c("age", "ext_rate")) %>%
  group_by(age) %>% 
  mutate(median_rate = median(ext_rate), q95 = quantile(ext_rate, 0.95), q5 = quantile(ext_rate, 0.05),  age = as.numeric(age), 
         color_point = ifelse(median_rate > 0, ifelse(median_rate > 0.05, "red", "orange"), "blue"), 
         perc_sup0 = sum(ext_rate>0)/1000*100, perc_sup0 = paste(perc_sup0, "%",sep = "")) %>% 
  dplyr::select(age, median_rate, perc_sup0) %>% distinct()

################################################################################
# 4 - MEGAFAUNA
################################################################################
megafauna <- read_delim("Megafauna_data_/Persotable.csv", delim = ",", col_names = T) %>% subset(!is.na(Area)) 

megafauna %>% dplyr::select(`Scientific name`, Area, `Approximate extinction time`) %>%
  setNames(c("name", "area", "age")) %>% arrange(desc(age)) %>% arrange(desc(name)) %>%
  mutate(uniq = 1:length(name), nameuniq = paste(name, uniq, sep = ".")) %>%
  mutate(nameuniq = factor(nameuniq, levels=c("Homotherium serum.7", "Equus spp..8", "Saiga tatarica.1", "Coelodonta antiquitatis.9", "Panthera spelaea.2", 
                                              "Panthera spelaea.3", "Mammuthus primigenius.5", "Mammuthus primigenius.6", "Bison priscus.10", "Megaloceros giganteus.4"))) %>%
  ggplot(aes(x = age, y = nameuniq)) +
  geom_segment( aes(x=25000, xend=age, y=nameuniq, yend=nameuniq, color = area)) +
  geom_point(aes(color= area), size=4, alpha=0.6) +
  # geom_segment(aes(x = 25000, y = nameuniq, xend = age, yend = nameuniq, color = area), size = 3) + 
  scale_x_reverse(breaks=seq(0, 25000, 2000)) +
  xlab("Age (cal. yrs BP)") +
  ylab("Megafauna estimated frame of extinction") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")
