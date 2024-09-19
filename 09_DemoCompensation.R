#//////////////////////////
# Integral Projection Models for Oenothera coloradensis: Analysis for demographic compensation
# Alice Stears
# 12 March 2022
#/////////////////////////

library(tidyverse)
library(ipmr)
library(rstatix)

# load IPMs 
source("./COBP_manuscriptAnalysis/06_IPMs_S_X.R")

# use 'corr' to make a correlation between vital rate functions for each env. covariate
# correlation for mean growing season temp
tMean_grow_Cor <- cor(dc_mods[dc_mods$vitalRate == "tMean_grow_C",1:5], method = "spearman")
# get the p-values
tMean_grow_Cor_pMat <- rstatix::cor_pmat(dc_mods[dc_mods$vitalRate == "tMean_grow_C",1:5])
# are five negative correlations (0 significant negative correlations)
# are five positive correlations (1 significant positive correlation)

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
set.seed("01121993")
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

# count number of positive correlations in each matrix
numPosCorrs <- sapply(X = corList, FUN = function(x) sum(x > 0))/2
numPosMean <- mean(numNegCorrs)
numPosSd <- sd(numNegCorrs)
# the number of negative correlations in the actual data is 5, although none are significant (for tMean_grow)
pnorm(q = 5, mean = numNegMean, sd = numNegSd)
# the number of significant negative correlations in the actual data is 0
pnorm(q = 0, mean = numNegMean, sd = numNegSd)

# the number of significant positive correlations in the actual data is 1 (for tMean_grow)
pnorm(q = 1, mean = numPosMean, sd = numPosSd)
# the number of non-significant positive correlations in the actual data is 5 (for tMean_grow)
pnorm(q = 5, mean = numPosMean, sd = numPosSd)

