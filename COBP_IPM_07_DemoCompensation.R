#//////////////////////////
# Integral Projection Models for Oenothera coloradensis: Analysis for demographic compensation
# Alice Stears
# 12 March 2022
#/////////////////////////

library(tidyverse)
library(ipmr)
library(rstatix)

# read in vital rate models
source("./analysis_scripts/01_VitalRateModels.R")
# read in data
dat_all <- read.csv(file = "../Processed_Data/allDat_plus_contSeedlings.csv")

for (i in 1:length(unique(dat_all$Site))) {
  dat_now <- dat_all[dat_all$Site == unique(dat_all$Site)[i],]
  ## vital rate models for each sites, including density dependence and all env. covariates 
  ## Survival ($s(z)$)
  survDat <- dat_now[dat_now$flowering==0,]
  # logistic glm with log-transformed size_t
  survMod_e_dd <- glm(survives_tplus1 ~ log_LL_t + tMean_grow_C  + tMean_winter_C  + precipWaterYr_cm + N_all_plot , data = survDat, family = "binomial")
  ## Growth ($G(z',z)$)
  sizeMod_e_dd <- lm(log_LL_tplus1 ~ log_LL_t  + tMean_grow_C  + tMean_winter_C  + precipWaterYr_cm + N_all_plot , data = dat_now)
  ## Number of seeds produced, according to plant size ($b(z)$)
  seedDat <- dat_now[dat_now$flowering == 1,]
  seedMod_e_dd <- MASS::glm.nb(Num_seeds ~ log_LL_t  + tMean_grow_C  + tMean_winter_C  + precipWaterYr_cm + N_all_plot , data = seedDat)
  ## Flowering probability ($p_b(z)$)
  flwrMod_e_dd <- glm(flowering ~ log_LL_t + I(log_LL_t^2)  + tMean_grow_C  + tMean_winter_C  + precipWaterYr_cm + N_all_plot , data = dat_now, family = "binomial")
  ## Distribution of recruit size ($c_o(z')$)
  recD = dat_now[dat_now$seedling==1,]
  recMod_env_dd <- lm(log_LL_t ~ tMean_grow_C  + tMean_winter_C  + precipWaterYr_cm + N_all_plot, data = recD)
  if (i == 1) {
    dc_mods <- data.frame("flwr" = coef(flwrMod_e_dd),
                          "surv" = c(coef(survMod_e_dd)[1:2], NA, coef(survMod_e_dd)[3:6]), 
                          "grow" = c(coef(sizeMod_e_dd)[1:2],NA,coef(sizeMod_e_dd)[3:6]), 
                          "seed" = c(coef(seedMod_e_dd)[1:2], NA, coef(seedMod_e_dd)[3:6]), 
                          "rec" = c(coef(recMod_env_dd)[1], NA, NA, coef(recMod_env_dd)[2:5]), 
                          "Site" = unique(dat_all$Site)[i])
  } else {
    dc_mods <- rbind(dc_mods, 
                     data.frame("flwr" = coef(flwrMod_e_dd),
                                "surv" = c(coef(survMod_e_dd)[1:2], NA, coef(survMod_e_dd)[3:6]), 
                                "grow" = c(coef(sizeMod_e_dd)[1:2],NA,coef(sizeMod_e_dd)[3:6]), 
                                "seed" = c(coef(seedMod_e_dd)[1:2], NA, coef(seedMod_e_dd)[3:6]), 
                                "rec" = c(coef(recMod_env_dd)[1], NA, NA, coef(recMod_env_dd)[2:5]), 
                                "Site" = unique(dat_all$Site)[i]))
  }
}

dc_mods$vitalRate <- rep(x = c("Intercept", "log_LL_t", "I(log_LL_t^2)", "tMean_grow_C", "tMean_winter_C", "precipWaterYr_cm", "N_all_plot"), times = 6)
# if the effect is 'NA', then replace with '0'
dc_mods[is.na(dc_mods)] <- 0
# use 'corr' to make a correlation between vital rate functions for each env. covariate
# correlation for mean growing season temp
tMean_grow_Cor <- cor(dc_mods[dc_mods$vitalRate == "tMean_grow_C",1:5], method = "pearson")
# get the p-values
tMean_grow_Cor_pMat <- rstatix::cor_pmat(dc_mods[dc_mods$vitalRate == "tMean_grow_C",1:5])

## use only tMean--there was not enough data to estimate correlation between each 
# correlation for mean winter temp 
tMean_winter_Cor <- cor(dc_mods[dc_mods$vitalRate == "tMean_winter_C",1:5])
## doesn't work--not significant in enough models

# correlation for mean growing season temp
precip_Cor <- cor(dc_mods[dc_mods$vitalRate == "precipWaterYr_cm",1:5])
## doesn't work--not significant in enough models

# do simulation of correlations of tMean coefficients to determine if there are significant differences
# make a normal dist. from which to pull the random coefficients
mean <- mean(as.matrix(dc_mods[dc_mods$vitalRate == "tMean_grow_C",1:5]), na.rm = TRUE)
sd <- sd(as.matrix(dc_mods[dc_mods$vitalRate == "tMean_grow_C",1:5]), na.rm = TRUE)
corList <- list()
for (i in 1:10000) {
  tempDF <- data.frame("flwr" = rnorm(6, mean = mean, sd = sd), 
             "surv" = rnorm(6, mean = mean, sd = sd),
             "grow" = rnorm(6, mean = mean, sd = sd), 
             "seed" = rnorm(6, mean = mean, sd = sd), 
             "rec" = rnorm(6, mean = mean, sd = sd))
  corList[[i]] <- cor(tempDF)
}
# count number of negative correlations in each matrix
numNegCorrs <- sapply(X = corList, FUN = function(x) sum(x < 0))/2
numNegMean <- mean(numNegCorrs)
numNegSd <- sd(numNegCorrs)

# the number of negative correlations in the actual data is 5 (for tMean_grow)
pnorm(q = 5, mean = numNegMean, sd = numNegSd)
