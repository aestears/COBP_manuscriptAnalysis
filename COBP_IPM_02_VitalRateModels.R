#///////////////////////////////////////////////////////
# Integral Projection Models for Oenothera coloradensis
# Part 2: Vital Rate Models for IPMs
# Alice Stears
# 3 December 2021
#///////////////////////////////////////////////////////
#### load packages ####
library(lme4) 
library(MASS)
#### load data from the previous script ####
# (COBP_IPM_01_dataPrep.R)
source("./analysis_scripts/COBP_IPM_01_dataPrep.R")

#### Vital Rate Models for Deterministic, non-density-dependent IPM with all data####
### Make vital rate models 
## Survival ($s(z)$)
# subset the data to exclude flowering individuals
survDat <- dat[dat$flowering==0 | is.na(dat$flowering),]
# logistic glm with log-transformed size_t
survMod <- glm(survives_tplus1 ~ log_LL_t , data = survDat, family = binomial)
summary(survMod)
# plot model results 
plot(survives_tplus1 ~ log_LL_t, data = survDat)
newdata <- data.frame("log_LL_t" = seq(from = min(survDat$log_LL_t, na.rm = TRUE), 
                                       to = max(survDat$log_LL_t, na.rm = TRUE),
                                       length.out = 100))
lines(x = newdata$log_LL_t, y = predict(object = survMod, newdata =  newdata, type = "response"), col = "red")

## Growth ($G(z',z)$)
# lm w/ log-transformed size_t and size_t+1
sizeMod <- lm(log_LL_tplus1 ~ log_LL_t , data = dat)
summary(sizeMod)
# plot model results
plot(log_LL_tplus1 ~ log_LL_t, data = dat)
newdata <- data.frame("log_LL_t" = seq(from = min(dat$log_LL_t, na.rm = TRUE), 
                                       to = max(dat$log_LL_t, na.rm = TRUE),
                                       length.out = 100))
lines(x = newdata$log_LL_t, y = predict(object = sizeMod, newdata =  newdata), col = "red")
lines(x = c(-1,4), y = c(-1,4), col = "darkgrey", lty = 2)

## Number of seeds produced, according to plant size ($b(z)$)
# using size in current year (no. of seeds/plant, for those that flowered ~ size_t)
seedDat <- dat[dat$flowering == 1,]
# fit a negative binomial glm (poisson was overdispersed)
seedMod_t <- MASS::glm.nb(Num_seeds ~ log_LL_t , data = seedDat)
summary(seedMod_t)
# plot model results
plot(Num_seeds ~ log_LL_t, data = seedDat)
newdata <- data.frame("log_LL_t" = seq(from = min(seedDat$log_LL_t, na.rm = TRUE), 
                                       to = max(seedDat$log_LL_t, na.rm = TRUE),
                                       length.out = 100))
lines(x = newdata$log_LL_t, y = predict(object = seedMod_t, newdata =  newdata, type = "response"), col = "red")

## Flowering probability ($p_b(z)$)
# using size in current year (w/ squared term)
# logistic glm with log-transformed size_t
flwrMod_t <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = dat, family = binomial)))
summary(flwrMod_t)
# plot model results 
plot(flowering ~ log_LL_t, data = dat)
newdata <- data.frame("log_LL_t" = seq(from = min(dat$log_LL_t, na.rm = TRUE), 
                                       to = max(dat$log_LL_t, na.rm = TRUE),
                                       length.out = 100))
lines(x = newdata$log_LL_t, y = predict(object = flwrMod_t, newdata =  newdata, type = "response"), col = "red")

## Distribution of recruit size ($c_o(z')$)
# subset the data
recD <- dat[dat$age == 0 & is.na(dat$age) == FALSE,]
# fit the model
recMod <- lm(log_LL_t ~ 1, data = recD)
summary(recMod)
#plot the results
hist(recD$log_LL_t)
abline(v = recMod$coefficients, col = "blue", lwd = 2)

## Probability of a seedling in year *t* establishing to a rosette in year *t+1* ($p_{estab}$)
# data.source is 'estabs' 
# calculate the probability of establishment value
p.estab.est <- sum(discDat[discDat$Seedling_t == 1, "Recruit_tplus1"], na.rm = TRUE)  / # number of new continuous stage recruits in t+1
  sum(discDat["Seedling_t"], na.rm = TRUE)  # number of seedlings in year t
#sum(estabs$P_estab, na.rm = TRUE)/sum(is.na(estabs$P_estab)==FALSE)


## Probability that a seed from the seedbank in year t will germinate to a seedling in year t+1 ($outSB$--is the 'germ.rt')
outSB.est <- germ.rt
  #sum(discDat[discDat$SeedBank_t == 1, "Seedling_tplus1"], na.rm = TRUE) / # number of seedlings in year t+1
  #sum(discDat[discDat$SeedBank_t == 1, "SeedBank_t"], na.rm = TRUE) # number of seedbank seeds in year t

## Probability that a seed from the seedbank in year t will stay in the seedbank in year t+1 ($staySB$)--Burgess, 2005 shows that rate of viability doesn't really decrease much with time
# (1 - germ.rt) * 0.9
staySB.est <- sum(discDat[discDat$SeedBank_t == 1, "SeedBank_tplus1"], na.rm = TRUE) / # number of seedbank seeds in year t+1
  sum(discDat[discDat$SeedBank_t == 1, "SeedBank_t"], na.rm = TRUE) # number of seedbank seeds in year t

## Probability that a seed produced by an adult plant in year t will enter the seedbank in year t+1 ($goSB$)
# viab.rt (1 - germ.rt) 
goSB.est <- viab.rt * (1 - germ.rt)
  #sum(discDat[discDat$NewSeeds_t == 1, "SeedBank_tplus1"], na.rm = TRUE) / # number of seedbank seeds in year t+1
  #sum(discDat[discDat$NewSeeds_t == 1, "NewSeeds_t"], na.rm = TRUE) # number of reproductive seeds in year t

## Probability that a seed from a plant in year t will go directly to the seedling stage ($goSdlng$)
# use the germination rate, since it doesn't seem to change much with age (Burgess, Hild & Shaw, 2005)
# viab.rt * germ.rt
 goSdlng.est <- viab.rt * germ.rt
# sum(discDat[discDat$NewSeeds_t == 1, "Seedling_tplus1"], na.rm = TRUE) / # number of seedbank seeds in year t+1
#   sum(discDat[discDat$NewSeeds_t == 1, "NewSeeds_t"], na.rm = TRUE) # number of seedbank seeds in year t

#### Vital Rate Models for deterministic, density-independent IPM with FIRST HALF OF DATA ####
## subset data for years 2018-2019
dat_first <- dat[dat$Year %in% c(2018),]

# subset the data to exclude flowering individuals
survDat_first <- dat_first[dat_first$flowering==0 | is.na(dat_first$flowering),]
# logistic glm with log-transformed size_t
survMod_first <- glm(survives_tplus1 ~ log_LL_t , data = survDat_first, family = binomial)

## Growth ($G(z',z)$)
# lm w/ log-transformed size_t and size_t+1
sizeMod_first <- lm(log_LL_tplus1 ~ log_LL_t , data = dat_first)

## Number of seeds produced, according to plant size ($b(z)$)
# using size in current year (no. of seeds/plant, for those that flowered ~ size_t)
seedDat_first <- dat_first[dat_first$flowering == 1,]
# fit a negative binomial glm (poisson was overdispersed)
seedMod_first <- MASS::glm.nb(Num_seeds ~ log_LL_t , data = seedDat_first)

## Flowering probability ($p_b(z)$)
# using size in current year (w/ squared term)
# logistic glm with log-transformed size_t
flwrMod_first <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = dat_first, family = binomial)))

## Distribution of recruit size ($c_o(z')$)
# subset the data
recD_first <- dat[dat$age == 0 & is.na(dat$age) == FALSE & dat$Year == 2019,]
# fit the model
recMod_first <- lm(log_LL_t ~ 1, data = recD_first)

#### Vital Rate Models for deterministic, density-independent IPM with SECOND HALF OF DATA ####
## subset data for years 2018-2019
dat_second <- dat[dat$Year %in% c(2018),]

# subset the data to exclude flowering individuals
survDat_second <- dat_second[dat_second$flowering==0 | is.na(dat_second$flowering),]
# logistic glm with log-transformed size_t
survMod_second <- glm(survives_tplus1 ~ log_LL_t , data = survDat_second, family = binomial)

## Growth ($G(z',z)$)
# lm w/ log-transformed size_t and size_t+1
sizeMod_second <- lm(log_LL_tplus1 ~ log_LL_t , data = dat_second)

## Number of seeds produced, according to plant size ($b(z)$)
# using size in current year (no. of seeds/plant, for those that flowered ~ size_t)
seedDat_second <- dat_second[dat_second$flowering == 1,]
# fit a negative binomial glm (poisson was overdispersed)
seedMod_second <- MASS::glm.nb(Num_seeds ~ log_LL_t , data = seedDat_second)

## Flowering probability ($p_b(z)$)
# using size in current year (w/ squared term)
# logistic glm with log-transformed size_t
flwrMod_second <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = dat_second, family = binomial)))

## Distribution of recruit size ($c_o(z')$)
# subset the data
recD_second <- dat[dat$age == 0 & is.na(dat$age) == FALSE & dat$Year == 2020,]
# fit the model
recMod_second <- lm(log_LL_t ~ 1, data = recD_second)
#### Vital Rate Models for Deterministic, density-dependent IPM with all data####
### Make vital rate models 
## Survival ($s(z)$)
# logistic glm with log-transformed size_t
survMod_dd <- glm(survives_tplus1 ~ log_LL_t + N_all_t , data = survDat, family = binomial)
summary(survMod_dd)
# plot model results 

## Growth ($G(z',z)$)
# lm w/ log-transformed size_t and size_t+1
sizeMod_dd <- lm(log_LL_tplus1 ~ log_LL_t + N_all_t , data = dat)
summary(sizeMod_dd)

## Number of seeds produced, according to plant size ($b(z)$)
# using size in current year (no. of seeds/plant, for those that flowered ~ size_t)
# fit poisson glm (for count data)
seedMod_dd <- MASS::glm.nb(Num_seeds ~ log_LL_t + N_all_t, data = seedDat)
summary(seedMod_dd)

## Flowering probability ($p_b(z)$)
# using size in current year (w/ squared term)
# logistic glm with log-transformed size_t
flwrMod_dd <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) + N_all_t, data = dat, family = binomial)))
summary(flwrMod_dd)

## Distribution of recruit size ($c_o(z')$)
# fit the model
recMod_dd <- lm(log_LL_t ~ 1 + N_all_t, data = recD)
summary(recMod_dd)

#### vital rate models for stochastic, density-independent IPM for all data with random effect for site and environmental covariates ####
### fit the vital rate models including environmental variables and random intercept for plot
## Survival ($s(z)$)
# data: survDat
# logistic glm with log-transformed size_t
survMod_e <- glmer(survives_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s  + SoilTemp_grow_C_s + tMean_grow_C_s + precipWaterYr_cm_s + (1|Plot_ID), data = survDat, family = binomial)
summary(survMod_e)
# try excluding the precip data (is not significant in the previous model)
survMod_e_1 <- glmer(survives_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s + SoilTemp_grow_C_s + tMean_grow_C_s + (1|Plot_ID),  data = survDat, family = binomial)
summary(survMod_e_1)
## survMod_e_1 is the better fit, use this model
survMod_env <- survMod_e_1

## Growth ($G(z',z)$)
# lm w/ log-transformed size_t and size_t+1
sizeMod_e <- lmer(log_LL_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s  + precipWaterYr_cm_s + tMean_grow_C_s +  (1|Plot_ID), data = dat)
summary(sizeMod_e)
# remove soilTemp in winter and growing season
sizeMod_e_1 <- lmer(log_LL_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s +  tMean_grow_C_s + precipWaterYr_cm_s + (1|Plot_ID), data = dat)
summary(sizeMod_e_1)
sizeMod_env <- sizeMod_e_1

## Number of seeds produced, according to plant size ($b(z)$)
# using size in current year (no. of seeds/plant, for those that flowered ~ size_t)
# use only plants that flowered 
# dat: seedDat
# fit poisson glm (for count data)
seedMod_e <- lme4::glmer.nb(Num_seeds ~ log_LL_t + SoilMoisture_m3m3_s + SoilTemp_winter_C_s + SoilTemp_grow_C_s + tMean_grow_C_s + precipWaterYr_cm_s + (1|Plot_ID), data = seedDat)
summary(seedMod_e)
# remove soil moisture, soil temp grow, soil temp winter (not significant)
seedMod_e_1 <- lme4::glmer.nb(Num_seeds ~ log_LL_t + tMean_grow_C_s + precipWaterYr_cm_s + (1|Plot_ID) , data = seedDat, family = poisson)
summary(seedMod_e_1)
seedMod_env <- seedMod_e_1

## Flowering probability ($p_b(z)$)
# using size in current year (w/ squared term)
# logistic glm with log-transformed size_t
flwrMod_e <- glmer(flowering ~ log_LL_t + I(log_LL_t^2) + SoilMoisture_m3m3_s + SoilTemp_grow_C_s + tMean_grow_C_s + precipWaterYr_cm_s + (1|Plot_ID), data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(flwrMod_e)
flwrMod_env <- flwrMod_e

## Distribution of recruit size ($c_o(z')$)
# subset the data
# data: recD 
# fit the model
recMod_env <- lmer(log_LL_t ~ SoilMoisture_m3m3_s + (1|Plot_ID), data = recD)
summary(recMod_env)
#plot the results

## Probability of a seedling in year *t* establishing to a rosette in year *t+1* ($p_{estab}$)
#p.estab.est

## Probability that a seed from the seedbank in year t will germinate to a seedling in year t+1 ($outSB$)
#outSB.est 

## Probability that a seed from the seedbank in year t will stay in the seedbank in year t+1 ($staySB$)
#staySB.est 

## Probability that a seed produced by an adult plant in year t will enter the seedbank in year t+1 ($goSB$)
#goSB.est

## Probability that a seed from a plant in year t will go directly to the seedling stage ($goSdlng$)
#goSdlng.est

#### Vital Rate Models for the Stochastic, density-dependent IPM for all data with random effect for site and environmental covariates ####
## Survival ($s(z)$)
# data: survDat 
dat$log_LL_t_s <- scale(dat$log_LL_t)
# logistic glm with log-transformed size_t
survMod_e_dd <- glmer(survives_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s  + SoilTemp_grow_C_s + tMean_grow_C_s + precipWaterYr_cm_s + N_all_t + (1|Plot_ID), data = survDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(survMod_e_dd)
# try excluding the precip data (is not significant in the previous model)
survMod_e_dd_1 <- glmer(survives_tplus1 ~ log_LL_t + tMean_grow_C_s + N_all_t +  (1|Plot_ID),  data = survDat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(survMod_e_dd_1)
## survMod_e_1 is the better fit, use this model
survMod_env_dd <- survMod_e_dd_1

## Growth ($G(z',z)$)
# lm w/ log-transformed size_t and size_t+1
sizeMod_e_dd <- lmer(log_LL_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s + SoilTemp_winter_C_s + SoilTemp_grow_C_s + tMean_grow_C_s + precipWaterYr_cm_s + N_all_t + (1|Plot_ID), data = dat)
summary(sizeMod_e_dd)
# remove soilTemp in winter and growing season; drop precip
sizeMod_e_dd_1 <- lmer(log_LL_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s +  tMean_grow_C_s + N_all_t + (1|Plot_ID), data = dat)
summary(sizeMod_e_dd_1)
sizeMod_env_dd <- sizeMod_e_dd_1

## Number of seeds produced, according to plant size ($b(z)$)
# data: seedDat 
# fit poisson glm (for count data)
seedMod_e_dd <- lme4::glmer.nb(Num_seeds ~ log_LL_t + SoilMoisture_m3m3_s + SoilTemp_winter_C_s + SoilTemp_grow_C_s + tMean_grow_C_s + precipWaterYr_cm_s +  N_all_t + (1|Plot_ID), data = seedDat, family = poisson)
summary(seedMod_e_dd)
# remove soil moisture, soil temp grow, soil temp winter (not significant)
seedMod_e_dd_1 <- lme4::glmer.nb(Num_seeds ~ log_LL_t + tMean_grow_C_s + precipWaterYr_cm_s + (1|Plot_ID) , data = seedDat, family = poisson)
summary(seedMod_e_dd_1)
## use the final model w/ no density dependence
seedMod_env_dd <- seedMod_e_dd_1

## Flowering probability ($p_b(z)$)
# using size in current year (w/ squared term)
# logistic glm with log-transformed size_t
flwrMod_e_dd <- glmer(flowering ~ log_LL_t + I(log_LL_t^2)  + tMean_grow_C_s  +  N_all_t + (1|Plot_ID), data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(flwrMod_e_dd)
# remove density dependence
flwrMod_e_dd_1 <- glmer(flowering ~ log_LL_t + I(log_LL_t^2)  + tMean_grow_C_s  + (1|Plot_ID), data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(flwrMod_e_dd_1)
# final model has no density dependence
flwrMod_env_dd <- flwrMod_e_dd

## Distribution of recruit size ($c_o(z')$)
# data = recD
# fit the model
recMod_env_dd <- lmer(log_LL_t ~ SoilMoisture_m3m3_s + N_all_t + (1|Plot_ID), data = recD)
summary(recMod_env_dd)
#plot the results

## Probability of a seedling in year *t* establishing to a rosette in year *t+1* ($p_{estab}$)
#p.estab.est

## Probability that a seed from the seedbank in year t will germinate to a seedling in year t+1 ($outSB$)
#outSB.est 

## Probability that a seed from the seedbank in year t will stay in the seedbank in year t+1 ($staySB$)
#staySB.est 

## Probability that a seed produced by an adult plant in year t will enter the seedbank in year t+1 ($goSB$)
#goSB.est

## Probability that a seed from a plant in year t will go directly to the seedling stage ($goSdlng$)
#goSdlng.est 

#### Models by site (deterministic, no env covariates, DI) ####
### Make vital rate models (put inside a for-loop; store models in a list)
# make an empty list to hold models
det_DI_mods <- list()
det_DI_mods[1:6] <- NA
# make a vector of site names
siteNames <- unique(dat$Site)
# start for-loop to fit models
for (i in 1:length(siteNames)) {
  # make a list to hold models for this site
  temp_mod_list <- list()
  temp_mod_list[1:5] <- NA
  # get the name of the site
  site_now <- siteNames[i]
  
  ## fit survival model
  # get data only for plants that didn't flower
  survDat_now <- dat[dat$flowering==0 | is.na(dat$flowering) & 
                       dat$Site == site_now,]
  # fit model and store in list
  temp_mod_list[[1]] <- glm(survives_tplus1 ~ log_LL_t , data = survDat_now, family = binomial)
  names(temp_mod_list)[1] <- paste0("surv")
  
  ## fit growth model
  # get all data, but just for this site
  allDat_now <- dat[dat$Site == site_now,]
  # fit model and store in a list
  temp_mod_list[[2]] <- lm(log_LL_tplus1 ~ log_LL_t , data = allDat_now)
  names(temp_mod_list)[2] <- paste0("growth")
  
  ## number of seeds produced, according to plant size
  seedDat_now <- dat[dat$flowering == 1 & dat$Site == site_now,]
  # fit the model and store it
  temp_mod_list[[3]] <- glm(Num_seeds ~ log_LL_t , data = seedDat_now, family = poisson)
  names(temp_mod_list)[3] <- paste0("seedProduction")
  
  ## probability of flowering
  temp_mod_list[[4]] <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = allDat_now, family = binomial)))
  names(temp_mod_list)[4] <- paste0("flowering")
  
  ## distribution of recruit size
  recD_now <- dat[dat$age == 0 & is.na(dat$age) == FALSE & 
                    dat$Site == site_now,]
  temp_mod_list[[5]] <- lm(log_LL_t ~ 1, data = recD_now)
  names(temp_mod_list)[5] <- paste0("recruitDist")
  
  ## store the temporary model list in the appropriate slot of the 'det_DI_mods' list
  det_DI_mods[[i]] <- temp_mod_list
  names(det_DI_mods)[i] <- site_now
  }
## use uniform p.estab, outSB, staySB, and goSB estimates

#### Models by site (deterministic, no env covariates, DD) ####
### Make vital rate models (put inside a for-loop; store models in a list)
# make an empty list to hold models
det_DD_mods <- list()
det_DD_mods[1:6] <- NA
# make a vector of site names
siteNames <- unique(dat$Site)
# start for-loop to fit models
for (i in 1:length(siteNames)) {
  # make a list to hold models for this site
  temp_mod_list <- list()
  temp_mod_list[1:5] <- NA
  # get the name of the site
  site_now <- siteNames[i]
  
  ## fit survival model
  # get data only for plants that didn't flower
  survDat_now <- dat[dat$flowering==0 | is.na(dat$flowering) & 
                       dat$Site == site_now,]
  # fit model and store in list
  temp_mod_list[[1]] <- glm(survives_tplus1 ~ log_LL_t + N_all_t, data = survDat_now, family = binomial)
  names(temp_mod_list)[1] <- paste0("surv")
  
  ## fit growth model
  # get all data, but just for this site
  allDat_now <- dat[dat$Site == site_now,]
  # fit model and store in a list
  temp_mod_list[[2]] <- lm(log_LL_tplus1 ~ log_LL_t + N_all_t, data = allDat_now)
  names(temp_mod_list)[2] <- paste0("growth")
  
  ## number of seeds produced, according to plant size
  seedDat_now <- dat[dat$flowering == 1 & dat$Site == site_now,]
  # fit the model and store it
  temp_mod_list[[3]] <- glm(Num_seeds ~ log_LL_t , data = seedDat_now, family = poisson)
  names(temp_mod_list)[3] <- paste0("seedProduction")
  
  ## probability of flowering
  temp_mod_list[[4]] <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = allDat_now, family = binomial)))
  names(temp_mod_list)[4] <- paste0("flowering")
  
  ## distribution of recruit size
  recD_now <- dat[dat$age == 0 & is.na(dat$age) == FALSE & 
                    dat$Site == site_now,]
  temp_mod_list[[5]] <- lm(log_LL_t ~ 1 + N_all_t, data = recD_now)
  names(temp_mod_list)[5] <- paste0("recruitDist")
  
  ## store the temporary model list in the appropriate slot of the 'det_DI_mods' list
  det_DD_mods[[i]] <- temp_mod_list
  names(det_DD_mods)[i] <- site_now
}
## use uniform p.estab, outSB, staySB, and goSB estimates

#### Models by site (stochastic, yes env covariates, DI) ####
### Make vital rate models (put inside a for-loop; store models in a list)
## for now, using the same model structure that we did for the whole-dataset models (i.e. using the covariates that best fit for larger models in each plot-level model)
# make an empty list to hold models
stoch_DI_mods <- list()
stoch_DI_mods[1:6] <- NA
# make a vector of site names
siteNames <- unique(dat$Site)
# start for-loop to fit models
for (i in 1:length(siteNames)) {
  # make a list to hold models for this site
  temp_mod_list <- list()
  temp_mod_list[1:5] <- NA
  # get the name of the site
  site_now <- siteNames[i]
  
  ## fit survival model
  # get data only for plants that didn't flower
  survDat_now <- dat[dat$flowering==0 | is.na(dat$flowering) & 
                       dat$Site == site_now,]
  # fit model and store in list
  temp_mod_list[[1]] <- glm(survives_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s + tMean_grow_C_s +  
                              SoilTemp_grow_C_s , 
                            data = survDat_now, 
                            family = binomial)
  names(temp_mod_list)[1] <- paste0("surv")
  
  ## fit growth model
  # get all data, but just for this site
  allDat_now <- dat[dat$Site == site_now,]
  # fit model and store in a list
  temp_mod_list[[2]] <- lm(log_LL_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s + tMean_grow_C_s , 
                           data = allDat_now)
  names(temp_mod_list)[2] <- paste0("growth")
  
  ## number of seeds produced, according to plant size
  seedDat_now <- dat[dat$flowering == 1 & dat$Site == site_now,]
  # fit the model and store it
  temp_mod_list[[3]] <- glm(Num_seeds ~ log_LL_t + tMean_grow_C_s , 
                            data = seedDat_now, 
                            family = poisson)
  names(temp_mod_list)[3] <- paste0("seedProduction")
  
  ## probability of flowering
  temp_mod_list[[4]] <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) + SoilMoisture_m3m3_s + tMean_grow_C_s + precipWaterYr_cm_s, 
                                              data = allDat_now, 
                                              family = binomial)))
  names(temp_mod_list)[4] <- paste0("flowering")
  
  ## distribution of recruit size
  recD_now <- dat[dat$age == 0 & is.na(dat$age) == FALSE & 
                    dat$Site == site_now,]
  temp_mod_list[[5]] <- lm(log_LL_t ~ 1 + SoilMoisture_m3m3_s, data = recD_now)
  names(temp_mod_list)[5] <- paste0("recruitDist")
  
  ## store the temporary model list in the appropriate slot of the 'det_DI_mods' list
  stoch_DI_mods[[i]] <- temp_mod_list
  names(stoch_DI_mods)[i] <- site_now
}
## use uniform p.estab, outSB, staySB, and goSB estimates 

#### Models by site (stochastic, yes env covariates, DD) ####
### Make vital rate models (put inside a for-loop; store models in a list)
## for now, using the same model structure that we did for the whole-dataset models (i.e. using the covariates that best fit for larger models in each plot-level model)
# make an empty list to hold models
stoch_DD_mods <- list()
stoch_DD_mods[1:6] <- NA
# make a vector of site names
siteNames <- unique(dat$Site)
# start for-loop to fit models
for (i in 1:length(siteNames)) {
  # make a list to hold models for this site
  temp_mod_list <- list()
  temp_mod_list[1:5] <- NA
  # get the name of the site
  site_now <- siteNames[i]
  
  ## fit survival model
  # get data only for plants that didn't flower
  survDat_now <- dat[dat$flowering==0 | is.na(dat$flowering) & 
                       dat$Site == site_now,]
  # fit model and store in list
  temp_mod_list[[1]] <- glm(survives_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s + tMean_grow_C_s +  
                              SoilTemp_grow_C_s + N_all_t, 
                            data = survDat_now, 
                            family = binomial)
  names(temp_mod_list)[1] <- paste0("surv")
  
  ## fit growth model
  # get all data, but just for this site
  allDat_now <- dat[dat$Site == site_now,]
  # fit model and store in a list
  temp_mod_list[[2]] <- lm(log_LL_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s + tMean_grow_C_s +  
                              N_all_t, 
                           data = allDat_now)
  names(temp_mod_list)[2] <- paste0("growth")
  
  ## number of seeds produced, according to plant size
  seedDat_now <- dat[dat$flowering == 1 & dat$Site == site_now,]
  # fit the model and store it
  temp_mod_list[[3]] <- glm(Num_seeds ~ log_LL_t + tMean_grow_C_s + precipWaterYr_cm_s , 
                            data = seedDat_now, 
                            family = poisson)
  names(temp_mod_list)[3] <- paste0("seedProduction")
  
  ## probability of flowering
  temp_mod_list[[4]] <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) + SoilMoisture_m3m3_s  + tMean_grow_C_s + precipWaterYr_cm_s, 
                                              data = allDat_now, 
                                              family = binomial)))
  names(temp_mod_list)[4] <- paste0("flowering")
  
  ## distribution of recruit size
  recD_now <- dat[dat$age == 0 & is.na(dat$age) == FALSE & 
                    dat$Site == site_now,]
  temp_mod_list[[5]] <- lm(log_LL_t ~ 1 + SoilMoisture_m3m3_s + N_all_t, data = recD_now)
  names(temp_mod_list)[5] <- paste0("recruitDist")
  
  ## store the temporary model list in the appropriate slot of the 'det_DI_mods' list
  stoch_DD_mods[[i]] <- temp_mod_list
  names(stoch_DD_mods)[i] <- site_now
}
## use uniform p.estab, outSB, staySB, and goSB estimates 
