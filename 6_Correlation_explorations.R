# TAXALOSS Project - script to explore which parameters to use for correlation in Figure 2 - mainly for supplementary material
# Script from Jeremy Courtin
# Script last update - 03.04.2023

# Layout of the script: 
# 0 - GENERAL OPTION - SETTING UP SCRIPT
# 1 - Import all data
# 2 - Make correlation tests
# 2.1 - Check for auto-correlations
# 2.2 - Compare extirpation and extinction
# 2.3 - Correlation to betadiversity estimate (dissimilarity - replacement rate)
# 2.4 - Compare to megafauna extinction
# 2.5 - Compare to climate change
# Simulated data selection
# Pollen data selection
# Best fitted with select simulated and pollen
# 2.6 - All best explanatory factors together

###############################################################################
# 0 - GENERAL OPTION - SETTING UP SCRIPT
###############################################################################
rm(list = ls()) # Remove all the objects we created so far.

# Set up options for the rest of the script
options(stringsAsFactors=FALSE) # state that everything is reads as character 

# Load all needed packages
library(tidyverse)
library(ggplot2)

###############################################################################
# 1 - Import all data
###############################################################################
# Climate reconstruction - Simulated data
climate <- read_delim("Temperature_estimates/2023_simulated_aggregated_timeslice.csv", delim = ",", col_names = T)

climate_dif_2000 <- climate %>% select(time_slice, mean_tjuly_anomaly1000, median_tjuly_anomaly1000, mean_pjuly_anomaly1000, median_pjuly_anomaly1000) %>% 
  mutate(time_slice = time_slice-2000, tjulyplus_mean = mean_tjuly_anomaly1000, tjulyplus_median = median_tjuly_anomaly1000, 
         pjulyplus_mean = mean_pjuly_anomaly1000, pjulyplus_median = median_pjuly_anomaly1000) %>%
  select(-mean_tjuly_anomaly1000, -median_tjuly_anomaly1000, -mean_pjuly_anomaly1000, -median_pjuly_anomaly1000) %>% 
  left_join(select(climate, time_slice, mean_tjuly_anomaly1000, median_tjuly_anomaly1000, mean_pjuly_anomaly1000, median_pjuly_anomaly1000)) %>% subset(time_slice > -1) %>%
  mutate(tjul_mean_diff = tjulyplus_mean-mean_tjuly_anomaly1000, tjul_median_diff = tjulyplus_median-median_tjuly_anomaly1000,
         pjul_mean_diff = pjulyplus_mean-mean_pjuly_anomaly1000, pjul_median_diff = pjulyplus_median-median_pjuly_anomaly1000) %>% 
  select(time_slice, tjul_mean_diff, tjul_median_diff, pjul_mean_diff, pjul_median_diff) %>%
  mutate(abs_tjul_mean_diff = abs(tjul_mean_diff), abs_tjul_median_diff = abs(tjul_median_diff), abs_pjul_mean_diff = abs(pjul_mean_diff), abs_pjul_median_diff = abs(pjul_median_diff)) %>% 
  distinct() %>% mutate(age = time_slice) %>% select(-time_slice)

climate_dif_2000 %>% pivot_longer(tjul_mean_diff:abs_pjul_median_diff) %>%
  subset(age > 1999) %>%
  subset(name == "tjul_median_diff") %>%
  ggplot(aes(x = age)) +
  geom_line(aes(y = value, color = name)) +
  scale_x_reverse(breaks=seq(0, 21000, 3000)) +
  xlab("Age (cal. yrs BP)") +
  ylab("") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

# Climate reconstruction - pollen
reconst_climate <- read_delim("Temperature_estimates/Temperature_differences_Alaska_Siberia_pollen_based_TO.csv", delim = ",", col_names = T)

reconst_prec <- read_delim("Temperature_estimates/precipitation_differences_alaska_sibeia_pollen_based.csv", delim = ",", col_names = T)

reconst_temp_ano <- read_delim("Temperature_estimates/Temperature_anomalies_Alaska_Siberia_pollen_based_JC.csv", delim = ",", col_names = T)
reconst_climate %>% pivot_longer(MAT_mean_diff:`WA-PLS_median_diff_clean`) %>%
  subset(age > 1999) %>%
  ggplot(aes(x = age)) +
  geom_line(aes(y = value, color = name)) +
  scale_x_reverse(breaks=seq(0, 25000, 2000)) +
  xlab("Age (cal. yrs BP)") +
  ylab("Pollen based - reconstructed July temperature (C)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

# Plot pollen based reconstruction
plot <- read_delim("Temperature_estimates/TJuly_MAT_anomalies_differences_pollen_siberia_alaska_28-1ka.csv", delim = ",", col_names = T) %>% 
  pivot_longer(cols = `1000`:`25000`) %>% mutate(name = as.numeric(name)) %>% setNames(c("site", "age", "diff_temp"))

plot %>% 
  subset(age > 1999) %>%
  ggplot(aes(x = age)) +
  geom_line(aes(y = diff_temp, color = site)) +
  scale_x_reverse(breaks=seq(0, 25000, 2000)) +
  xlab("Age (cal. yrs BP)") +
  ylab("Pollen based - reconstructed July temperature (C)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

names(reconst_prec)
reconst_prec <- reconst_prec %>% setNames(c("age", "MAT_mean_diffp", "MAT_median_diffp", "MAT_mean_diff_cleanp", "MAT_median_diff_cleanp", 
                                            "WAPLS_mean_diffp", "WAPLS_median_diffp", "WAPLS_mean_diff_cleanp", "WAPLS_median_diff_cleanp"))
reconst_prec %>% pivot_longer(MAT_mean_diffp:WAPLS_median_diff_cleanp) %>%
  subset(age > 1999) %>%
  ggplot(aes(x = age)) +
  geom_line(aes(y = value, color = name)) +
  scale_x_reverse(breaks=seq(0, 25000, 2000)) +
  xlab("Age (cal. yrs BP)") +
  ylab("Pollen based - reconstructed precipitation") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

# Diversity measurements
div_2000 <- read_delim("Median_Diversity_Estimates.csv", delim = ",", col_names = T)
div_all <- read_delim("beta_diversity_1000iterations.csv", delim = ",", col_names = T)

# Plant extinction rates
p_loss_2000 <- read_delim("Figure2_data_save.csv", delim = ",", col_names = T) %>% setNames(c("name", "p_loss"))
p_extirpation_2000 <- read_delim("Figure2_data_save_db_extirpation.csv", delim = ",", col_names = T) %>% setNames(c("name", "p_extir"))
p_extinction_2000 <- read_delim("Figure2_data_save_cand_extinction.csv", delim = ",", col_names = T) %>% setNames(c("name", "p_extin"))

p_ext_2000 <- cbind(p_loss_2000, select(p_extirpation_2000, p_extir), select(p_extinction_2000, p_extin)) %>% as_tibble()

# Megafauna extinction
megafauna <- read_delim("Megafauna_data/megafauna.csv", delim = ",", col_names = T) %>% subset(!is.na(Area)) %>% mutate(timeslice_disappear_1000 = ceiling(`Approximate extinction time`/1000)*1000)
megafauna_2000 <- megafauna %>% select(timeslice_disappear_1000) %>% group_by(timeslice_disappear_1000) %>% summarise(mega_fauna = sum(timeslice_disappear_1000>0)) %>% setNames(c("age", "megafauna"))

# Check 2000 timeslices
div_raw <- div_all %>% arrange(age) %>% subset(age > 1000)
div_time_lagminus <- div_all %>% mutate(age = age-2000) %>% arrange(age) %>% subset(age > 1000) %>% setNames(c("age1", "Beta_rich_minus", "Beta_D_minus", "Beta_repl_minus")) %>% add_row(age1 = 1:1000) %>% mutate(age1 = ifelse(age1 < 1001, 25000, age1))
div_time_lagplus <- div_all %>% mutate(age = age+2000) %>% arrange(age) %>% subset(age > 1000 & age < 26000) %>% setNames(c("age2", "Beta_rich_plus", "Beta_D_plus", "Beta_repl_plus"))
p_ext_raw <- p_ext_2000 %>% arrange(name) %>% subset(name > 1000) %>% mutate(p_loss_pos = ifelse(p_loss<0, 0, p_loss), p_extir_pos = ifelse(p_extir<0, 0, p_extir), p_extin_pos = ifelse(p_extin<0, 0, p_extin))

rate_2000 <- cbind(p_ext_raw, div_raw, div_time_lagminus, div_time_lagplus) %>% select(-name, -age1, -age2) %>% as_tibble()
rate_2000 <- left_join(rate_2000, megafauna_2000) %>% left_join(climate_dif_2000) %>% left_join(reconst_climate) %>% left_join(reconst_prec)

###############################################################################
# 2 - Make correlation tests
###############################################################################
# set megafauna extinction as non NA
modeDat <- rate_2000 %>% mutate(megafauna_pos = ifelse(is.na(megafauna), 0, megafauna)) %>% dplyr::select(-megafauna)
write_delim(modeDat, "taxaloss_parameters_for_correlation.csv", delim = ",")
library("ggpubr")
library("cowplot")
library(nlme)

##### General plots #####
names(modeDat)
M <- cor(modeDat[apply(modeDat, 1, function(x) all(!is.na(x))),-1])
corrplot::corrplot(M)

# Select only the final factors
modeDat1 <- modeDat %>% dplyr::select(age, p_extir_pos, p_extin_pos, Beta_repl_minus, megafauna_pos, MAT_mean_diff_clean) %>% filter(!is.na(Beta_repl_minus))
M <- cor(modeDat1[apply(modeDat1, 1, function(x) all(!is.na(x))),-1])
corrplot::corrplot(M)

# also quick overview:
# quick overview
pairs(modeDat1)

# take only the likely best data
names(modeDat)
modedat_c <- modeDat %>% select(-MAT_mean_diff, -MAT_median_diff, -'WA-PLS_mean_diff', -'WA-PLS_median_diff', -MAT_mean_diffp, -MAT_median_diffp, -WAPLS_mean_diffp, -WAPLS_median_diffp)
modedat_c <- modedat_c %>% select(age, p_extir_pos, p_extin_pos, megafauna_pos,
                                  Beta_D, Beta_repl, Beta_D_minus, Beta_repl_minus, Beta_D_plus, Beta_repl_plus,
                                  tjul_mean_diff, tjul_median_diff, abs_tjul_mean_diff, abs_tjul_median_diff,
                                  pjul_mean_diff, pjul_median_diff, abs_pjul_mean_diff, abs_pjul_median_diff,
                                  MAT_mean_diff_clean, MAT_median_diff_clean, MAT_mean_diff_cleanp, MAT_median_diff_cleanp,
                                  `WA-PLS_mean_diff_clean`, `WA-PLS_median_diff_clean`, WAPLS_mean_diff_cleanp, WAPLS_median_diff_cleanp)
M <- cor(modedat_c[apply(modedat_c, 1, function(x) all(!is.na(x))),])
corrplot::corrplot(M)

###############################################################################
# 2.1 - Check for auto-correlations
###############################################################################
##### Set age as autocorrelation factor #####
library(lme4) # http://www.flutterbys.com.au/stats/course.html 
##### extirpation #####
dat_extir <- modeDat %>% select(-p_extin_pos) %>%
  group_split(age) %>% lapply(., function(x) x %>% mutate(sim = 1:nrow(x))) %>% bind_rows()

names(dat_extir)
dat1 <- glmer(p_extir_pos ~  Beta_repl_minus + MAT_mean_diff_clean + megafauna_pos + (age | sim), data = dat_extir, family = binomial(), 
              control = glmerControl(optimizer="bobyqa"))
summary(dat1) 
car::Anova(dat1, type=3)

# without autocorrelation control
dat2 <- glm(p_extir_pos ~  Beta_repl_minus + MAT_mean_diff_clean + megafauna_pos, dat_extir, family = "binomial")
summary(dat2) 
car::Anova(dat2, type=3)

##### extinction #####
dat_extin <- modeDat %>% select(-p_extir_pos) %>%
  group_split(age) %>% lapply(., function(x) x %>% mutate(sim = 1:nrow(x))) %>% bind_rows()

###############################################################################
# 2.2 - Compare extirpation and extinction
###############################################################################
modeDat %>%
  group_by(age) %>% #summarise(p_extir_pos = median(p_extir_pos), p_extin_pos = median(p_extin_pos)) %>% # OPTIONAL
  ggscatter(x = "p_extir_pos", y = "p_extin_pos", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          title = "spearman")

modeDat %>% 
  group_by(age) %>% #summarise(p_extir_pos = median(p_extir_pos), p_extin_pos = median(p_extin_pos)) %>% # OPTIONAL
  ggscatter(x = "p_extir_pos", y = "p_extin_pos", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            title = "pearson")

corr <- glm(p_extir_pos~p_extin_pos, modeDat, family = "binomial")
summary(corr)
with(summary(corr), 1 - deviance/null.deviance)

###############################################################################
# 2.3 - Correlation to betadiversity estimate (dissimilarity - replacement rate)
###############################################################################
##### extinction correlation #####
names(dat_extin)
corr <- glmer(p_extin_pos~Beta_repl+Beta_D+Beta_D_minus+Beta_repl_minus+Beta_D_plus+Beta_repl_plus+(age | sim), dat_extin, family = binomial(), 
            control = glmerControl(optimizer="bobyqa"))
summary(corr)
car::Anova(corr, type=3)

corr1 <- glm(p_extin_pos~Beta_repl+Beta_D+Beta_D_minus+Beta_repl_minus+Beta_D_plus+Beta_repl_plus, dat_extin, family = "binomial")
summary(corr1)
car::Anova(corr1, type=3) # similar output between the corrected model for age and the standard model
step(corr1, test="LRT")

##### extirpation correlation #####
corr <- glmer(p_extir_pos~Beta_repl+Beta_D+Beta_D_minus+Beta_repl_minus+Beta_D_plus+Beta_repl_plus+(age | sim), dat_extir, family = binomial(), 
              control = glmerControl(optimizer="bobyqa"))
summary(corr)
car::Anova(corr, type=3)

corr1 <- glm(p_extir_pos~Beta_repl+Beta_D+Beta_D_minus+Beta_repl_minus+Beta_D_plus+Beta_repl_plus, dat_extir, family = "binomial")
summary(corr1)
car::Anova(corr1, type=3) # similar output between the corrected model for age and the standard model
step(corr1, test="LRT")

###############################################################################
# 2.4 - Compare to megafauna extinction
###############################################################################
##### extirpation correlation #####
names(dat_extir)
dat_extir %>% group_by(age) %>% #mutate(p_extir_pos = median(p_extir_pos)) %>%
  ggscatter(x = "p_extir_pos", y = "megafauna_pos", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "spearman",
            title = "spearman")

dat_extir %>% group_by(age) %>% #mutate(p_extir_pos = median(p_extir_pos)) %>%
  ggscatter(x = "p_extir_pos", y = "megafauna_pos", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            title = "pearson")

corr <- glmer(p_extir_pos~megafauna_pos+(age | sim), dat_extir, family = binomial(), 
              control = glmerControl(optimizer="bobyqa"))
summary(corr)
car::Anova(corr, type=3)

corr1 <- glm(p_extir_pos~megafauna_pos, dat_extir, family = "binomial")
summary(corr1)
with(summary(corr1), 1 - deviance/null.deviance)
step(corr1, test="LRT")

##### extinction correlation #####
names(dat_extin)
dat_extin %>% 
  group_by(age) %>% mutate(p_extin_pos = median(p_extin_pos)) %>% # OPTIONAL
  ggscatter(x = "p_extin_pos", y = "megafauna_pos", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "spearman",
            title = "spearman")

dat_extin %>% 
  #group_by(age) %>% mutate(p_extin_pos = median(p_extin_pos)) %>% # OPTIONAL
  ggscatter(x = "p_extin_pos", y = "megafauna_pos", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            title = "pearson")

corr <- glmer(p_extin_pos~megafauna_pos+(age | sim), dat_extin, family = binomial(), 
              control = glmerControl(optimizer="bobyqa"))
summary(corr)
car::Anova(corr, type=3)

corr1 <- glm(p_extin_pos~megafauna_pos, dat_extin, family = "binomial")
summary(corr1)
with(summary(corr1), 1 - deviance/null.deviance)
step(corr1, test="LRT")

###############################################################################
# 2.5 - Compare to climate change
###############################################################################
###############################################################################
# Simulated data selection
###############################################################################
##### extincion correlation #####
names(dat_extin)
corr <- glmer(p_extin_pos~tjul_mean_diff+tjul_median_diff+pjul_mean_diff+pjul_median_diff+abs_tjul_mean_diff+abs_tjul_median_diff+abs_pjul_mean_diff+abs_pjul_median_diff+ (age | sim), dat_extin, family = binomial(), 
              control = glmerControl(optimizer="bobyqa"))
summary(corr)
car::Anova(corr, type=3)

corr1 <- glm(p_extin_pos~tjul_mean_diff+tjul_median_diff+pjul_mean_diff+pjul_median_diff+abs_tjul_mean_diff+abs_tjul_median_diff+abs_pjul_mean_diff+abs_pjul_median_diff, dat_extin, family = "binomial")
summary(corr1) # similar output between the corrected model for age and the standard model
with(summary(corr1), 1 - deviance/null.deviance)
step(corr1, test="LRT") # best fitted are tjul_mean_diff / pjul_mean_diff / abs_pjul_mean_diff / pjul_median_diff / tjul_median_diff / abs_pjul_median_diff

##### extirpation correlation ##### 
corr <- glmer(p_extir_pos~tjul_mean_diff+tjul_median_diff+pjul_mean_diff+pjul_median_diff+abs_tjul_mean_diff+abs_tjul_median_diff+abs_pjul_mean_diff+abs_pjul_median_diff+ (age | sim), dat_extir, family = binomial(), 
              control = glmerControl(optimizer="bobyqa"))
summary(corr)
car::Anova(corr, type=3)

corr1 <- glm(p_extir_pos~tjul_mean_diff+tjul_median_diff+pjul_mean_diff+pjul_median_diff+abs_tjul_mean_diff+abs_tjul_median_diff+abs_pjul_mean_diff+abs_pjul_median_diff, dat_extir, family = "binomial")
summary(corr1) # similar output between the corrected model for age and the standard model
with(summary(corr1), 1 - deviance/null.deviance)
step(corr1, test="LRT") # best fitted are abs_pjul_median_diff /  pjul_mean_diff / tjul_median_diff / pjul_median_diff / abs_pjul_mean_diff / tjul_mean_diff / abs_tjul_median_diff / abs_tjul_mean_diff

# So we choose to keep for further tests the top 5 present in both extinction and extirpation:
# tjul_mean_diff / pjul_mean_diff / abs_pjul_median_diff / tjul_median_diff / pjul_median_diff

###############################################################################
# Pollen data selection
###############################################################################
##### extincion correlation #####
names(dat_extin)
corr <- glmer(p_extin_pos~MAT_mean_diff_clean+MAT_median_diff_clean+`WA-PLS_mean_diff_clean`+`WA-PLS_median_diff_clean`+
                MAT_mean_diff_cleanp+MAT_median_diff_cleanp+WAPLS_mean_diff_cleanp+WAPLS_median_diff_cleanp+(age | sim), dat_extin, family = binomial(), 
              control = glmerControl(optimizer="bobyqa"))
summary(corr)
car::Anova(corr, type=3)

corr1 <- glm(p_extin_pos~MAT_mean_diff_clean+MAT_median_diff_clean+`WA-PLS_mean_diff_clean`+`WA-PLS_median_diff_clean`+
               MAT_mean_diff_cleanp+MAT_median_diff_cleanp+WAPLS_mean_diff_cleanp+WAPLS_median_diff_cleanp, dat_extin, family = "binomial")
summary(corr1) # similar output between the corrected model for age and the standard model
with(summary(corr1), 1 - deviance/null.deviance)
step(corr1, test="LRT") # best fitted are WAPLS_mean_diff_cleanp / MAT_median_diff_cleanp / MAT_mean_diff_cleanp / `WA-PLS_median_diff_clean` / WAPLS_median_diff_cleanp / MAT_mean_diff_clean / MAT_median_diff_clean / `WA-PLS_mean_diff_clean`

##### extirpation correlation #####
names(dat_extir)
corr <- glmer(p_extir_pos~MAT_mean_diff_clean+MAT_median_diff_clean+`WA-PLS_mean_diff_clean`+`WA-PLS_median_diff_clean`+
                MAT_mean_diff_cleanp+MAT_median_diff_cleanp+WAPLS_mean_diff_cleanp+WAPLS_median_diff_cleanp+(age | sim), dat_extir, family = binomial(), 
              control = glmerControl(optimizer="bobyqa"))
summary(corr)
car::Anova(corr, type=3)

corr1 <- glm(p_extir_pos~MAT_mean_diff_clean+MAT_median_diff_clean+`WA-PLS_mean_diff_clean`+`WA-PLS_median_diff_clean`+
               MAT_mean_diff_cleanp+MAT_median_diff_cleanp+WAPLS_mean_diff_cleanp+WAPLS_median_diff_cleanp, dat_extir, family = "binomial")
summary(corr1) # similar output between the corrected model for age and the standard model
with(summary(corr1), 1 - deviance/null.deviance)
step(corr1, test="LRT") # best fitted are WAPLS_median_diff_cleanp / WAPLS_mean_diff_cleanp / MAT_median_diff_cleanp / MAT_mean_diff_cleanp / `WA-PLS_median_diff_clean`/ `WA-PLS_mean_diff_clean`/ MAT_median_diff_clean / MAT_mean_diff_clean

# So we choose to keep for further tests the top 5 present in both extinction and extirpation: we also prefer temperature over precipitation
# MAT_mean_diff_clean / MAT_median_diff_clean / WAPLS_median_diff_cleanp / WAPLS_mean_diff_cleanp / MAT_median_diff_cleanp / MAT_mean_diff_cleanp

###############################################################################
# Best fitted with select simulated and pollen
###############################################################################
dat_extin_s <- dat_extin %>% select(p_extin_pos, age, sim, 
                                    tjul_mean_diff, pjul_mean_diff, abs_pjul_median_diff, tjul_median_diff, pjul_median_diff,
                                    MAT_mean_diff_clean, MAT_median_diff_clean, WAPLS_median_diff_cleanp, WAPLS_mean_diff_cleanp, MAT_median_diff_cleanp, MAT_mean_diff_cleanp)
dat_extir_s <- dat_extir %>% select(p_extir_pos, age, sim, 
                                    tjul_mean_diff, pjul_mean_diff, abs_pjul_median_diff, tjul_median_diff, pjul_median_diff,
                                    MAT_mean_diff_clean, MAT_median_diff_clean, WAPLS_median_diff_cleanp, WAPLS_mean_diff_cleanp, MAT_median_diff_cleanp, MAT_mean_diff_cleanp)
##### extincion correlation #####
names(dat_extin_s)
corr <- glmer(p_extin_pos~tjul_mean_diff+pjul_mean_diff+abs_pjul_median_diff+tjul_median_diff+pjul_median_diff+
              MAT_mean_diff_clean+MAT_median_diff_clean+WAPLS_median_diff_cleanp+WAPLS_mean_diff_cleanp+MAT_median_diff_cleanp+MAT_mean_diff_cleanp+(age | sim), dat_extin_s, family = binomial(), 
              control = glmerControl(optimizer="bobyqa"))
summary(corr)
car::Anova(corr, type=3)

corr1 <- glm(p_extin_pos~tjul_mean_diff+pjul_mean_diff+abs_pjul_median_diff+tjul_median_diff+pjul_median_diff+
             MAT_mean_diff_clean+MAT_median_diff_clean+WAPLS_median_diff_cleanp+WAPLS_mean_diff_cleanp+MAT_median_diff_cleanp+MAT_mean_diff_cleanp, dat_extin_s, family = "binomial")
summary(corr1) # similar output between the corrected model for age and the standard model
with(summary(corr1), 1 - deviance/null.deviance)
step(corr1, test="LRT")

##### extirpation correlation #####
names(dat_extir_s)
corr <- glmer(p_extir_pos~tjul_mean_diff+pjul_mean_diff+abs_pjul_median_diff+tjul_median_diff+pjul_median_diff+
                MAT_mean_diff_clean+MAT_median_diff_clean+WAPLS_median_diff_cleanp+WAPLS_mean_diff_cleanp+MAT_median_diff_cleanp+MAT_mean_diff_cleanp+(age | sim), dat_extir_s, family = binomial(), 
              control = glmerControl(optimizer="bobyqa"))
summary(corr)
car::Anova(corr, type=3)

corr1 <- glm(p_extir_pos~tjul_mean_diff+pjul_mean_diff+abs_pjul_median_diff+tjul_median_diff+pjul_median_diff+
               MAT_mean_diff_clean+MAT_median_diff_clean+WAPLS_median_diff_cleanp+WAPLS_mean_diff_cleanp+MAT_median_diff_cleanp+MAT_mean_diff_cleanp, dat_extir_s, family = "binomial")
summary(corr1) # similar output between the corrected model for age and the standard model
with(summary(corr1), 1 - deviance/null.deviance)
step(corr1, test="LRT")

###############################################################################
# 2.6 - All best explanatory factors together
###############################################################################
##### extirpation correlation #####
dat_extir_f <- dat_extir %>% select(p_extir_pos, age, sim, MAT_median_diff_clean, Beta_repl_minus, megafauna_pos)
names(dat_extir_f)
corr <- glmer(p_extir_pos~MAT_median_diff_clean+Beta_repl_minus+megafauna_pos+(age | sim), dat_extir_f, family = binomial(), 
              control = glmerControl(optimizer="bobyqa"))
summary(corr)
car::Anova(corr, type=3)

#Single LGM:
# Temperature
corr <- glm(p_extir_pos~MAT_median_diff_clean, dat_extir_f, family = "binomial")
summary(corr)
car::Anova(corr, type=3)

# Replacement
corr <- glm(p_extir_pos~Beta_repl_minus, dat_extir_f, family = "binomial")
summary(corr)
car::Anova(corr, type=3)

# Megafauna
corr <- glm(p_extir_pos~megafauna_pos, dat_extir_f, family = "binomial")
summary(corr)
car::Anova(corr, type=3)

##### extincion correlation #####
dat_extin_f <- dat_extin %>% select(p_extin_pos, age, sim, MAT_median_diff_clean, Beta_repl_minus, megafauna_pos)
names(dat_extin_f)
corr <- glmer(p_extin_pos~MAT_median_diff_clean+Beta_repl_minus+megafauna_pos+(age | sim), dat_extin_f, family = binomial(), 
              control = glmerControl(optimizer="bobyqa"))
summary(corr)
car::Anova(corr, type=3)

#Single LGM:
# Temperature
corr <- glm(p_extin_pos~MAT_median_diff_clean, dat_extin_f, family = "binomial")
summary(corr)
car::Anova(corr, type=3)

# Replacement
corr <- glm(p_extin_pos~Beta_repl_minus, dat_extin_f, family = "binomial")
summary(corr)
car::Anova(corr, type=3)

# Megafauna
corr <- glm(p_extin_pos~megafauna_pos, dat_extin_f, family = "binomial")
summary(corr)
car::Anova(corr, type=3)






