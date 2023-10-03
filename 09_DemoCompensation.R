#//////////////////////////
# Integral Projection Models for Oenothera coloradensis: Analysis for demographic compensation
# Alice Stears
# 12 March 2022
#/////////////////////////

library(tidyverse)
library(ipmr)
library(rstatix)

# load IPMs 
source("./analysis_scripts/06_IPMs_S_X.R")

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
