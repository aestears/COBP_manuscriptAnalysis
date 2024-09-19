#///////////////////////////////////////////////////////
# Integral Projection Models for Oenothera coloradensis
# Part 2: Vital Rate Models for IPMs
# Alice Stears
# 3 December 2021
#///////////////////////////////////////////////////////
#### load packages ####
library(tidyverse)
library(lme4) 
library(MASS)
#### load data from the previous script ####
# (COBP_IPM_01_dataPrep.R)
#source("./analysis_scripts/COBP_IPM_01_dataPrep.R")
dat_all <- read.csv(file = "../Processed_Data/ZenodoSubmission/allDat_plus_contSeedlings.csv")
discDat <- read.csv(file = "../Processed_Data/ZenodoSubmission/discreteStageData.csv")

#### calculate seedbank vital rates ####
## data from Burgess, Hild & Shaw, 2005--NOT using data from MT place, only from FEWAFB
# Seeds per Capsule (total no., viable and inviable)
seed_per_cap <- mean(c(2.4, 1.0))
# Capsule Viability(%) (percentage of capsules that are viable--contain >=1 seed)
capsule_viab.rt <- mean(c(81, 61, 54))/100
# Seed Viability (%) (percentage of seeds in a viable capsule that are viable)
seed_viab.rt <- mean(c(1.9/2.4,  1/1))

## calculate the rate at which a seed produced in a capsule in year t is viable (probability of a viable capsule * probability that a seed inside a viable capsule is viable) 
total_seed_viab.rt <- capsule_viab.rt * seed_viab.rt

# (data in SeedBagGreenhouseSeedlings.csv)--seeds from previous year
germ.rt.ours <- .03
#data from (Burgess, Hild & Shaw, 2005)--seedbank seed viability/germination rate doesn't seem to change much over time-- only from WAFB sites
germ.rt.Burgess <- mean(c(13.0, 12, 8.3, 7.0, 5.3)/(45 * seed_per_cap))

# germination rate from the Burgess paper incorporates both viability and germination. To isolate just the germination rate, divide the Burgess germination rate by the viability rate
germ.rt_temp <-germ.rt.Burgess/total_seed_viab.rt  
# because the Burgess germ.rate was estimated in a greenhouse, and we are confident that field rates are lower, multiply teh Burgess germ.rate by .80
germ.rt <- germ.rt_temp * .8

# the viability rate (proportion of seeds produced by an adult plant that are viable) is the 'total_seed_viab.rt' derived from the Burgess paper results
viab.rt <- total_seed_viab.rt

#### Vital Rate Models for IPM A and IPM B####
### Deterministic, density-independent IPM with all data + continuous seedlings ###
## the dataset is called 'dat_all'
 ### Make vital rate models 
 ## Survival ($s(z)$)
 # subset the data to exclude flowering individuals
 survDat_all <- dat_all[dat_all$flowering==0 | is.na(dat_all$flowering),]
 # logistic glm with log-transformed size_t
 survMod_all <- glm(survives_tplus1 ~ log_LL_t , data = survDat_all, family = binomial)
 summary(survMod_all)
 # plot model results 
 plot(survives_tplus1 ~ log_LL_t, data = survDat_all)
 newdata <- data.frame("log_LL_t" = seq(from = min(survDat_all$log_LL_t, na.rm = TRUE), 
                                        to = max(survDat_all$log_LL_t, na.rm = TRUE),
                                        length.out = 100))
 lines(x = newdata$log_LL_t, y = predict(object = survMod_all, newdata =  newdata, type = "response"), col = "red")
 
 ## Growth ($G(z',z)$)
 # lm w/ log-transformed size_t and size_t+1
 sizeMod_all <- lm(log_LL_tplus1 ~ log_LL_t , data = dat_all)
 summary(sizeMod_all)
 # plot model results
 plot(log_LL_tplus1 ~ log_LL_t, data = dat_all)
 newdata <- data.frame("log_LL_t" = seq(from = min(dat_all$log_LL_t, na.rm = TRUE), 
                                        to = max(dat_all$log_LL_t, na.rm = TRUE),
                                        length.out = 100))
 lines(x = newdata$log_LL_t, y = predict(object = sizeMod_all, newdata =  newdata), col = "red")
 lines(x = c(-1,4), y = c(-1,4), col = "darkgrey", lty = 2)
 
 ## Number of seeds produced, according to plant size ($b(z)$)
 # using size in current year (no. of seeds/plant, for those that flowered ~ size_t)
 seedDat_all <- dat_all[dat_all$flowering == 1,]
 # fit a negative binomial glm (poisson was overdispersed)
 seedMod_all <- MASS::glm.nb(Num_seeds ~ log_LL_t , data = seedDat_all)
 summary(seedMod_all)
 # plot model results
 plot(Num_seeds ~ log_LL_t, data = seedDat_all)
 newdata <- data.frame("log_LL_t" = seq(from = min(seedDat_all$log_LL_t, na.rm = TRUE), 
                                        to = max(seedDat_all$log_LL_t, na.rm = TRUE),
                                        length.out = 100))
 lines(x = newdata$log_LL_t, y = predict(object = seedMod_all, newdata =  newdata, type = "response"), col = "red")
 
 ## Flowering probability ($p_b(z)$)
 # using size in current year (w/ squared term)
 # logistic glm with log-transformed size_t
 flwrMod_all <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = dat_all, family = binomial)))
 summary(flwrMod_all)
 # plot model results 
 plot(flowering ~ log_LL_t, data = dat_all)
 newdata <- data.frame("log_LL_t" = seq(from = min(dat_all$log_LL_t, na.rm = TRUE), 
                                        to = max(dat_all$log_LL_t, na.rm = TRUE),
                                        length.out = 100))
 lines(x = newdata$log_LL_t, y = predict(object = flwrMod_all, newdata =  newdata, type = "response"), col = "red")
 
 ## Distribution of recruit size ($c_o(z')$)
 # subset the data
 recD_all <- dat_all[dat_all$seedling == 1,]
 # plot the data
 hist((recD_all$log_LL_t))
 plot(density(recD_all$log_LL_t), ylim = c(0,5))
 
 recMod_all <- lm(log_LL_t ~ 1, data = recD_all)
 summary(recMod_all)
 #plot the results
 abline(v = recMod_all$coefficients, col = "blue", lwd = 2)
 
## Probability that a seed from the seedbank in year t will germinate to a seedling in year t+1 ($outSB$--is the 'germ.rt')
 outSB_all <- germ.rt * .9 # (???)
 
## Probability that a seed from the seedbank in year t will stay in the seedbank in year t+1 ($staySB$)--Burgess, 2005 shows that rate of viability doesn't really decrease much with time
 # (1 - germ.rt) * 0.9
 staySB_all<- (1-germ.rt) * .9
 
 ## Probability that a seed produced by an adult plant in year t will enter the seedbank in year t+1 ($goSB$)
 # viab.rt (1 - germ.rt) 
 goSB_all <- viab.rt * (1 - germ.rt)
 
 ## Probability that a seed from a plant in year t will go directly to the continous stage ($goCont$)
 # use the germination rate, since it doesn't seem to change much with age (Burgess, Hild & Shaw, 2005)
 # viab.rt * germ.rt
 goCont_all <- viab.rt * germ.rt
 
 #### Vital Rate Models for Deterministic, density-dependent IPM with all data + continuous seedlings####
 ## the dataset is called 'dat_all'
 ### Make vital rate models 
 ## Survival ($s(z)$)
 # subset the data to exclude flowering individuals
 survDat_N <- dat_all[dat_all$flowering==0 | is.na(dat_all$flowering),]
 # logistic glm with log-transformed size_t
 survMod_N <- glm(survives_tplus1 ~ log_LL_t + N_all_plot_t , data = survDat_N, family = binomial)
 summary(survMod_N)
 # plot model results 
 plot(survives_tplus1 ~ log_LL_t, data = survDat_N)
 newdata <- data.frame("log_LL_t" = seq(from = min(survDat_N$log_LL_t, na.rm = TRUE), 
                                        to = max(survDat_N$log_LL_t, na.rm = TRUE),
                                        length.out = 100),  
                       "N_all_plot_t" = seq(from = min(survDat_N$N_all_plot_t, na.rm = TRUE),
                                     to = max(survDat_N$N_all_plot_t, na.rm = TRUE), 
                                     length.out = 100))
 lines(x = newdata$log_LL_t, y = predict(object = survMod_N, newdata =  newdata, type = "response"), col = "red")
 
 ## Growth ($G(z',z)$)
 # lm w/ log-transformed size_t and size_t+1
 sizeMod_N <- lm(log_LL_tplus1 ~ log_LL_t + N_all_plot_t, data = dat_all)
 summary(sizeMod_N)
 # plot model results
 plot(log_LL_tplus1 ~ log_LL_t, data = dat_all)
 newdata <- data.frame("log_LL_t" = seq(from = min(dat_all$log_LL_t, na.rm = TRUE), 
                                        to = max(dat_all$log_LL_t, na.rm = TRUE),
                                        length.out = 100),
                       "N_all_plot_t" = seq(from = min(dat_all$N_all_plot_t, na.rm = TRUE),
                                     to = max(dat_all$N_all_plot_t, na.rm = TRUE),
                                     length.out = 100))
 lines(x = newdata$log_LL_t, y = predict(object = sizeMod_N, newdata =  newdata), col = "red")
 lines(x = c(-1,4), y = c(-1,4), col = "darkgrey", lty = 2)
 
 ## Number of seeds produced, according to plant size ($b(z)$)
 # using size in current year (no. of seeds/plant, for those that flowered ~ size_t)
 seedDat_N <- dat_all[dat_all$flowering == 1,]
 # fit a negative binomial glm (poisson was overdispersed)
 seedMod_N <- MASS::glm.nb(Num_seeds ~ log_LL_t, data = seedDat_N)
 summary(seedMod_N)
 # plot model results
 plot(Num_seeds ~ log_LL_t, data = seedDat_all)
 newdata <- data.frame("log_LL_t" = seq(from = min(seedDat_N$log_LL_t, na.rm = TRUE), 
                                        to = max(seedDat_N$log_LL_t, na.rm = TRUE),
                                        length.out = 100))
 lines(x = newdata$log_LL_t, y = predict(object = seedMod_N, newdata =  newdata, type = "response"), col = "red")
 
 ## Flowering probability ($p_b(z)$)
 # using size in current year (w/ squared term)
 # logistic glm with log-transformed size_t
 flwrMod_N <- suppressWarnings(glm(flowering ~ log_LL_t + I(log_LL_t^2) + N_all_plot_t,
                                    data = dat_all, family = binomial))
 summary(flwrMod_N)
 # plot model results 
 plot(flowering ~ log_LL_t, data = dat_all)
 newdata <- data.frame("log_LL_t" = seq(from = min(dat_all$log_LL_t, na.rm = TRUE), 
                                        to = max(dat_all$log_LL_t, na.rm = TRUE),
                                        length.out = 100),
                       "N_all_plot_t" = seq(from = min(dat_all$N_all_plot_t, na.rm = TRUE),
                                     to = max(dat_all$N_all_plot_t, na.rm = TRUE), length.out = 100))

 lines(x = newdata$log_LL_t, y = predict(object = flwrMod_N, newdata =  newdata, type = "response"), col = "red")
 
 ## Distribution of recruit size ($c_o(z')$)
 # subset the data
 recD_all <- dat_all[dat_all$seedling == 1,]
 # fit the model
 recMod_N <- lm(log_LL_t ~ 1, data = recD_all)
 summary(recMod_N)
 #plot the results
 hist(recD_all$log_LL_t)
 abline(v = recMod_all$coefficients, col = "blue", lwd = 2)
 
 ## Probability that a seed from the seedbank in year t will germinate to a seedling in year t+1 ($outSB$--is the 'germ.rt')
 outSB_all <- germ.rt * .9
 
 ## Probability that a seed from the seedbank in year t will stay in the seedbank in year t+1 ($staySB$)--Burgess, 2005 shows that rate of viability doesn't really decrease much with time
 # (1 - germ.rt) * 0.9
 staySB_all<- (1-germ.rt) * .9
 
 ## Probability that a seed produced by an adult plant in year t will enter the seedbank in year t+1 ($goSB$)
 # viab.rt (1 - germ.rt) 
 goSB_all <- viab.rt * (1 - germ.rt)
 
 ## Probability that a seed from a plant in year t will go directly to the continous stage ($goCont$)
 # use the germination rate, since it doesn't seem to change much with age (Burgess, Hild & Shaw, 2005)
 # viab.rt * germ.rt
 goCont_all <- viab.rt * germ.rt

 #### Models by site (no env. co-variates, density independent, all data + continuous seedlings) ####
 ### Make vital rate models (put inside a for-loop; store models in a list)
 # make an empty list to hold models
 det_DI_mods <- list()
 det_DI_mods[1:6] <- NA
 # make a vector of site names
 siteNames <- unique(dat_all$Site)
 # start for-loop to fit models
 for (i in 1:length(siteNames)) {
   # make a list to hold models for this site
   temp_mod_list <- list()
   temp_mod_list[1:5] <- NA
   # get the name of the site
   site_now <- siteNames[i]
   
   ## fit survival model
   # get data only for plants that didn't flower
   survDat_now <- dat_all[(dat_all$flowering==0 | is.na(dat_all$flowering)) & 
                        dat_all$Site == site_now,]
   # fit model and store in list
   temp_mod_list[[1]] <- glm(survives_tplus1 ~ log_LL_t , data = survDat_now, family = binomial)
   names(temp_mod_list)[1] <- paste0("surv")
   
   ## fit growth model
   # get all data, but just for this site
   allDat_now <- dat_all[dat_all$Site == site_now,]
   # fit model and store in a list
   temp_mod_list[[2]] <- lm(log_LL_tplus1 ~ log_LL_t , data = allDat_now)
   names(temp_mod_list)[2] <- paste0("growth")
   
   ## number of seeds produced, according to plant size
   seedDat_now <- dat_all[dat_all$flowering == 1 & dat_all$Site == site_now,]
   # fit the model and store it
   temp_mod_list[[3]] <- glm(Num_seeds ~ log_LL_t , data = seedDat_now, family = poisson)
   names(temp_mod_list)[3] <- paste0("seedProduction")
   
   ## probability of flowering
   temp_mod_list[[4]] <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = allDat_now, family = binomial)))
   names(temp_mod_list)[4] <- paste0("flowering")
   
   ## distribution of recruit size
   recD_now <- dat_all[dat_all$age == 0 & is.na(dat_all$age) == FALSE & 
                     dat_all$Site == site_now,]
   temp_mod_list[[5]] <- lm(log_LL_t ~ 1, data = recD_now)
   names(temp_mod_list)[5] <- paste0("recruitDist")
   
   ## store the temporary model list in the appropriate slot of the 'det_DI_mods' list
   det_DI_mods[[i]] <- temp_mod_list
   names(det_DI_mods)[i] <- site_now
 }
 ## use uniform p.estab, outSB, staySB, and goSB estimates
 
 #### Models by site (no env covariates, density dependent, with continuous seedlings) ####
 ### Make vital rate models (put inside a for-loop; store models in a list)
 # make an empty list to hold models
 det_DD_mods <- list()
 det_DD_mods[1:6] <- NA
 # make a vector of site names
 siteNames <- unique(dat_all$Site)
 # start for-loop to fit models
 for (i in 1:length(siteNames)) {
   # make a list to hold models for this site
   temp_mod_list <- vector(mode = "list", length = 5)
   # get the name of the site
   site_now <- siteNames[i]
   
   ## fit survival model
   # get data only for plants that didn't flower
   survDat_now <- dat_all[(dat_all$flowering==0 | is.na(dat_all$flowering)) & 
                        dat_all$Site == site_now,]
   # fit model and store in list
   temp_mod_list[[1]] <- glm(survives_tplus1 ~ log_LL_t + N_Site_t, data = survDat_now, family = binomial)
   names(temp_mod_list)[1] <- paste0("surv")
   
   ## fit growth model
   # get all data, but just for this site
   allDat_now <- dat_all[dat_all$Site == site_now,]
   
   # fit model and store in a list
   temp_mod_list[[2]] <- lm(log_LL_tplus1 ~ log_LL_t + N_Site_t, data = allDat_now)
   names(temp_mod_list)[2] <- paste0("growth")
   
   ## number of seeds produced, according to plant size
   seedDat_now <- dat_all[dat_all$flowering == 1 & dat_all$Site == site_now,]
   # fit the model and store it
   temp_mod_list[[3]] <- glm(Num_seeds ~ log_LL_t , data = seedDat_now, family = poisson)
   names(temp_mod_list)[3] <- paste0("seedProduction")
   
   ## probability of flowering
   temp_mod_list[[4]] <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = allDat_now, family = binomial)))
   names(temp_mod_list)[4] <- paste0("flowering")
   
   ## distribution of recruit size
   recD_now <- dat_all[dat_all$age == 0 & is.na(dat_all$age) == FALSE & 
                     dat_all$Site == site_now,]
   temp_mod_list[[5]] <- lm(log_LL_t ~ 1 + N_Site_t, data = recD_now)
   names(temp_mod_list)[5] <- paste0("recruitDist")
   
   ## store the temporary model list in the appropriate slot of the 'det_DI_mods' list
   det_DD_mods[[i]] <- temp_mod_list
   names(det_DD_mods)[i] <- site_now
 }
 
 # 
#### Vital Rate Models for deterministic, density-independent IPM with FIRST HALF OF DATA (including cont. seedlings) #### 
## subset data for years 2018-2019
dat_first <- dat_all[dat_all$Year %in% c(2018),]

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
recD_first <- dat_all[dat_all$seedling == 1  & dat_all$Year == 2019,]
# fit the model
recMod_first <- lm(log_LL_t ~ 1, data = recD_first)

#### Vital Rate Models for deterministic, density-independent IPM with SECOND HALF OF DATA (including cont. seedlings) ####
## subset data for years 2018-2019
dat_second <- dat_all[dat_all$Year %in% c(2019),]

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
recD_second <- dat_all[dat_all$seedling == 1 & dat_all$Year == 2020,]
# fit the model
recMod_second <- lm(log_LL_t ~ 1, data = recD_second)

#------------------ vital rate models for IPMR functions 

#### vital rate models for stochastic, density-independent IPM for all data with random effect for site and environmental covariates ####
### fit the vital rate models including environmental variables and random intercept for plot
## Survival ($s(z)$)
# data: survDat
# logistic glm with log-transformed size_t
survMod_e <- glmer(survives_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s  + SoilTemp_grow_C_s + tMean_grow_C_s + precipWaterYr_cm_s + (1|Plot_ID), data = survDat_all, family = binomial)
summary(survMod_e)
# try excluding the precip data (is not significant in the previous model)
survMod_e_1 <- glmer(survives_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s + SoilTemp_grow_C_s + tMean_grow_C_s + (1|Plot_ID),  data = survDat_all, family = binomial)
summary(survMod_e_1)
## survMod_e_1 is the better fit, use this model
survMod_env <- survMod_e_1

## Growth ($G(z',z)$)
# lm w/ log-transformed size_t and size_t+1
sizeMod_e <- lmer(log_LL_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s  + precipWaterYr_cm_s + tMean_grow_C_s +  (1|Plot_ID), data = dat_all)
summary(sizeMod_e)
# remove soilTemp in winter and growing season
sizeMod_e_1 <- lmer(log_LL_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s +  tMean_grow_C_s + precipWaterYr_cm_s + (1|Plot_ID), data = dat_all)
summary(sizeMod_e_1)
sizeMod_env <- sizeMod_e_1

## Number of seeds produced, according to plant size ($b(z)$)
# using size in current year (no. of seeds/plant, for those that flowered ~ size_t)
# use only plants that flowered 
# dat: seedDat
# fit poisson glm (for count data)
seedMod_e <- lme4::glmer.nb(Num_seeds ~ log_LL_t + SoilMoisture_m3m3_s + SoilTemp_winter_C_s + SoilTemp_grow_C_s + tMean_grow_C_s + precipWaterYr_cm_s + (1|Plot_ID), data = seedDat_all)
summary(seedMod_e)
# remove soil moisture, soil temp grow, soil temp winter (not significant)
seedMod_e_1 <- lme4::glmer.nb(Num_seeds ~ log_LL_t + tMean_grow_C_s + precipWaterYr_cm_s + (1|Plot_ID) , data = seedDat_all, family = poisson)
summary(seedMod_e_1)
seedMod_env <- seedMod_e_1

## Flowering probability ($p_b(z)$)
# using size in current year (w/ squared term)
# logistic glm with log-transformed size_t
flwrMod_e <- glmer(flowering ~ log_LL_t + I(log_LL_t^2) + SoilMoisture_m3m3_s + SoilTemp_grow_C_s + tMean_grow_C_s + precipWaterYr_cm_s + (1|Plot_ID), data = dat_all, family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(flwrMod_e)
flwrMod_env <- flwrMod_e

## Distribution of recruit size ($c_o(z')$)
# subset the data
# data: recD 
# fit the model
recMod_env <- lmer(log_LL_t ~ SoilMoisture_m3m3_s + (1|Plot_ID), data = recD_all)
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

#### Models by site (deterministic, no env covariates, DI) ####
### Make vital rate models (put inside a for-loop; store models in a list)
# make an empty list to hold models
det_DI_mods <- list()
det_DI_mods[1:6] <- NA
# make a vector of site names
siteNames <- unique(dat_all$Site)
# start for-loop to fit models
for (i in 1:length(siteNames)) {
  # make a list to hold models for this site
  temp_mod_list <- list()
  temp_mod_list[1:5] <- NA
  # get the name of the site
  site_now <- siteNames[i]
  
  ## fit survival model
  # get data only for plants that didn't flower
  survDat_now <- dat_all[dat_all$flowering==0 | is.na(dat_all$flowering) & 
                       dat_all$Site == site_now,]
  # fit model and store in list
  temp_mod_list[[1]] <- glm(survives_tplus1 ~ log_LL_t , data = survDat_now, family = binomial)
  names(temp_mod_list)[1] <- paste0("surv")
  
  ## fit growth model
  # get all data, but just for this site
  allDat_now <- dat_all[dat_all$Site == site_now,]
  # fit model and store in a list
  temp_mod_list[[2]] <- lm(log_LL_tplus1 ~ log_LL_t , data = allDat_now)
  names(temp_mod_list)[2] <- paste0("growth")
  
  ## number of seeds produced, according to plant size
  seedDat_now <- dat_all[dat_all$flowering == 1 & dat_all$Site == site_now,]
  # fit the model and store it
  temp_mod_list[[3]] <- glm(Num_seeds ~ log_LL_t , data = seedDat_now, family = poisson)
  names(temp_mod_list)[3] <- paste0("seedProduction")
  
  ## probability of flowering
  temp_mod_list[[4]] <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = allDat_now, family = binomial)))
  names(temp_mod_list)[4] <- paste0("flowering")
  
  ## distribution of recruit size
  recD_now <- dat_all[dat_all$age == 0 & is.na(dat_all$age) == FALSE & 
                    dat_all$Site == site_now,]
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
siteNames <- unique(dat_all$Site)
# start for-loop to fit models
for (i in 1:length(siteNames)) {
  # make a list to hold models for this site
  temp_mod_list <- list()
  temp_mod_list[1:5] <- NA
  # get the name of the site
  site_now <- siteNames[i]
  
  ## fit survival model
  # get data only for plants that didn't flower
  survDat_now <- dat_all[dat_all$flowering==0 | is.na(dat_all$flowering) & 
                       dat_all$Site == site_now,]
  # fit model and store in list
  temp_mod_list[[1]] <- glm(survives_tplus1 ~ log_LL_t + N_all_plot_t, data = survDat_now, family = binomial)
  names(temp_mod_list)[1] <- paste0("surv")
  
  ## fit growth model
  # get all data, but just for this site
  allDat_now <- dat_all[dat_all$Site == site_now,]
  # fit model and store in a list
  temp_mod_list[[2]] <- lm(log_LL_tplus1 ~ log_LL_t + N_all_plot_t, data = allDat_now)
  names(temp_mod_list)[2] <- paste0("growth")
  
  ## number of seeds produced, according to plant size
  seedDat_now <- dat_all[dat_all$flowering == 1 & dat_all$Site == site_now,]
  # fit the model and store it
  temp_mod_list[[3]] <- glm(Num_seeds ~ log_LL_t , data = seedDat_now, family = poisson)
  names(temp_mod_list)[3] <- paste0("seedProduction")
  
  ## probability of flowering
  temp_mod_list[[4]] <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = allDat_now, family = binomial)))
  names(temp_mod_list)[4] <- paste0("flowering")
  
  ## distribution of recruit size
  recD_now <- dat_all[dat_all$age == 0 & is.na(dat_all$age) == FALSE & 
                    dat_all$Site == site_now,]
  temp_mod_list[[5]] <- lm(log_LL_t ~ 1 + N_all_plot_t, data = recD_now)
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
siteNames <- unique(dat_all$Site)
# start for-loop to fit models
for (i in 1:length(siteNames)) {
  # make a list to hold models for this site
  temp_mod_list <- list()
  temp_mod_list[1:5] <- NA
  # get the name of the site
  site_now <- siteNames[i]
  
  ## fit survival model
  # get data only for plants that didn't flower
  survDat_now <- dat_all[dat_all$flowering==0 | is.na(dat_all$flowering) & 
                       dat_all$Site == site_now,]
  # fit model and store in list
  temp_mod_list[[1]] <- glm(survives_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s + tMean_grow_C_s +  
                              SoilTemp_grow_C_s , 
                            data = survDat_now, 
                            family = binomial)
  names(temp_mod_list)[1] <- paste0("surv")
  
  ## fit growth model
  # get all data, but just for this site
  allDat_now <- dat_all[dat_all$Site == site_now,]
  # fit model and store in a list
  temp_mod_list[[2]] <- lm(log_LL_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s + tMean_grow_C_s , 
                           data = allDat_now)
  names(temp_mod_list)[2] <- paste0("growth")
  
  ## number of seeds produced, according to plant size
  seedDat_now <- dat_all[dat_all$flowering == 1 & dat_all$Site == site_now,]
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
  recD_now <- dat_all[dat_all$age == 0 & is.na(dat_all$age) == FALSE & 
                    dat_all$Site == site_now,]
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
siteNames <- unique(dat_all$Site)
# start for-loop to fit models
for (i in 1:length(siteNames)) {
  # make a list to hold models for this site
  temp_mod_list <- list()
  temp_mod_list[1:5] <- NA
  # get the name of the site
  site_now <- siteNames[i]
  
  ## fit survival model
  # get data only for plants that didn't flower
  survDat_now <- dat_all[dat_all$flowering==0 | is.na(dat_all$flowering) & 
                       dat_all$Site == site_now,]
  # fit model and store in list
  temp_mod_list[[1]] <- glm(survives_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s + tMean_grow_C_s +  
                              SoilTemp_grow_C_s + N_all_plot_t, 
                            data = survDat_now, 
                            family = binomial)
  names(temp_mod_list)[1] <- paste0("surv")
  
  ## fit growth model
  # get all data, but just for this site
  allDat_now <- dat_all[dat_all$Site == site_now,]
  # fit model and store in a list
  temp_mod_list[[2]] <- lm(log_LL_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s +  
                              N_all_plot_t, 
                           data = allDat_now)
  names(temp_mod_list)[2] <- paste0("growth")
  
  ## number of seeds produced, according to plant size
  seedDat_now <- dat_all[dat_all$flowering == 1 & dat_all$Site == site_now,]
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
  recD_now <- dat_all[dat_all$age == 0 & is.na(dat_all$age) == FALSE & 
                    dat_all$Site == site_now,]
  temp_mod_list[[5]] <- lm(log_LL_t ~ 1 + SoilMoisture_m3m3_s + N_all_plot_t, data = recD_now)
  names(temp_mod_list)[5] <- paste0("recruitDist")
  
  ## store the temporary model list in the appropriate slot of the 'det_DI_mods' list
  stoch_DD_mods[[i]] <- temp_mod_list
  names(stoch_DD_mods)[i] <- site_now
}
## use uniform p.estab, outSB, staySB, and goSB estimates 


## find out how many individuals are biennial 
# calculate age at death for each individual
test <- dat_all %>% 
  group_by(Location, Site, Plot_ID, Quadrant, ID) %>% 
  summarize(ageAtDeath = max(age),
            diedDuringStudy = (survives_tplus1 == 0)) %>% 
  group_by(Location, Site, Plot_ID, Quadrant, ID) %>% 
  summarize(diedDuringStudy = sum(diedDuringStudy), 
            ageAtDeath = max(ageAtDeath)) %>% 
  filter(ageAtDeath == 1 & 
           diedDuringStudy == 1) 
nrow(test)

test2 <- dat_all %>% 
  filter(!is.na(age)) %>% 
  dplyr::select(Location, Site, Plot_ID, Quadrant, ID) %>% 
  unique()
nrow(test2)

dat_all$uniqueID <- paste0(dat_all$Plot_ID, "_", dat_all$Quadrant, "_", dat_all$ID)
# look at cohort recruited in 2018 
RecruitIDs_2018  <- dat_all %>% 
  filter(Year == 2018 & 
           seedling == 1)

nrow(RecruitIDs_2018)  ## 1802 individuals that were recruited in 2018

survIDs <- RecruitIDs_2018 %>% 
  filter(survives_tplus1 == 1)
test <- dat_all %>% 
  filter(uniqueID %in% survIDs$uniqueID) %>% # get those individuals that were recruited in 2018 and survived past 2018
  filter(survives_tplus1 == 0) %>% # filter to the year when these individuals died
  filter(Year == 2019)# filter for those that died in 2020 (had a 0 for survives_tplus1 in 2019)

length(test$uniqueID)
## 351 individuals recruited in 2018 died after 2 years (20% were biennial)
  
## how many of the 2018 individuals flowered? 
test2 <- dat_all %>% 
  filter(uniqueID %in% survIDs$uniqueID) %>% 
  filter(flowering == 1)
  ## 24/1802 individuals from 2018 flowered in the time we observed them (1.3%)

## how many of the 2018 *biennials flowered?* 
test3 <- dat_all %>% 
  filter(uniqueID %in% test$uniqueID)%>% 
  filter(flowering == 1)
# only 5 of the 351 biennials flowered (1.4% of the biennial population)


## what is the average age of flowering individuals 
flwr <- dat_all %>% 
  filter(flowering == 1) 

mean(flwr$age, na.rm = TRUE)

  
