#//////////////////////////
# Integral Projection Models for Oenothera coloradensis: Analysis for negative density dependence
# Alice Stears
# 11 March 2022
#/////////////////////////

# load required packages
library(tidyverse)
library(MASS)
library(ipmr)
library(ggpubr)

# load data from script 1
dat_all <- read.csv(file = "../Processed_Data/allDat_plus_contSeedlings.csv")

### IPMs for this analysis are in the "03_IPMs_C_N.R" file 
#### read in data saved to file ####
fileLoc <- "./intermediate_analysis_Data/site_level_IPMs_allYears/"

## site-level DI IPM matrices
IPMs_C_H <- readRDS(file = paste0(fileLoc,"/IPMs_C_H.RDS"))
## site-level DI bootstrap CI data
IPMs_C_H_bootCI_lambdas <- readRDS(file = paste0(fileLoc,"/IPMs_C_H_bootCI_lambdas.RDS"))
IPMs_C_H_bootCI_params <- readRDS(file = paste0(fileLoc,"/IPMs_C_H_bootCI_params.RDS"))

## site-level DD IPM matrices
IPMs_I_N <- readRDS(file = paste0(fileLoc,"/IPMs_I_N.RDS"))
## site-level DD bootstrap CI data
IPMs_I_N_bootCI_lambdas <- readRDS(file = paste0(fileLoc, "/IPMs_I_N_bootCI_lambdas.RDS"))
IPMs_I_N_bootCI_params <- readRDS(file = paste0(fileLoc, "/IPMs_I_N_bootCI_params.RDS"))

#### calculate 95% confidence intervals for lambdas ####
## DI lambda CIs
## calculate 95% CI for lambdas
DI_lams_CIs <- as.data.frame(lapply(X = IPMs_C_H_bootCI_lambdas, FUN = function(x) 
  c(log(mean(x) - 1.96*(sd(x)/sqrt(1000))),log(mean(x) + 1.96*(sd(x)/sqrt(1000))))
  ))
names(DI_lams_CIs) <- unique(dat_all$Site)

## DD lambda CIs
DD_lams_CIs <- as.data.frame(lapply(X = IPMs_I_N_bootCI_lambdas, FUN = function(x) 
  c(log(mean(x) - 1.96*(sd(x)/sqrt(1000))),log(mean(x) + 1.96*(sd(x)/sqrt(1000))))
))
names(DD_lams_CIs) <- unique(dat_all$Site)

#### compare the vital rate models w/ and w/out DD--use AIC ####
## Density Independent
# list to hold models
subPop_modList_DI <- list()
subPop_AIClist_DI <- list()
for (i in 1:length(unique(dat_all$Site))) {
  ## get data for this 'current' site
  dat_now <- dat_all[dat_all$Site == unique(dat_all$Site)[i],]
  
  ## fit vital rate models
  ## Survival ($s(z)$)
  survDat_now <- dat_now[dat_now$flowering==0 | is.na(dat_now$flowering),]
  survMod_now <- glm(survives_tplus1 ~ log_LL_t , data = survDat_now, family = binomial)
  ## Growth ($G(z',z)$)
  sizeMod_now <- lm(log_LL_tplus1 ~ log_LL_t , data = dat_now)
  ## Number of seeds produced, according to plant size ($b(z)$)
  seedDat_now <- dat_now[dat_now$flowering==1,]
  # fit poisson glm (for count data)
  seedMod_now <- MASS::glm.nb(Num_seeds ~ log_LL_t , data = seedDat_now)
  ## Flowering probability ($p_b(z)$)
  flwrMod_now <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = dat_now, family = binomial)))
  ## Distribution of recruit size ($c_o(z')$)
  # subset the data
  recD_now <- dat_now[dat_now$seedling == 1,]
  # fit the model
  recMod_now <- lm(log_LL_t ~ 1, data = recD_now)
  
  ## put in the parameter list (paramCont)
  modelList <- list(
    surv     = survMod_now, 
    growth   = sizeMod_now,
    seed     = seedMod_now,
    flwr     = flwrMod_now, 
    recs     = recMod_now              
  )
  
  ## put AICs for each model in a list
  AIClist <- list(
    surv     = AIC(survMod_now),  
    growth   = AIC(sizeMod_now),
    seed     = AIC(seedMod_now),
    flwr     = AIC(flwrMod_now), 
    recs     = AIC(recMod_now)
  )
  subPop_modList_DI[[i]] <- modelList
  subPop_AIClist_DI[[i]] <- AIClist
}
names(subPop_modList_DI) <- unique(dat_all$Site)
names(subPop_AIClist_DI) <- unique(dat_all$Site)

## Density Dependent
# list to hold models
subPop_modList_DD <- list()
subPop_AIClist_DD <- list()
for (i in 1:length(unique(dat_all$Site))) {
  ## get data for this 'current' site
  dat_now <- dat_all[dat_all$Site == unique(dat_all$Site)[i],]
  
  ## fit vital rate models
  ## Survival ($s(z)$)
  survDat_now <- dat_now[dat_now$flowering==0 | is.na(dat_now$flowering),]
  survMod_now <- glm(survives_tplus1 ~ log_LL_t + N_Site_t, data = survDat_now, family = binomial)
  ## Growth ($G(z',z)$)
  sizeMod_now <- lm(log_LL_tplus1 ~ log_LL_t + N_Site_t , data = dat_now)
  ## Number of seeds produced, according to plant size ($b(z)$)
  seedDat_now <- dat_now[dat_now$flowering==1,]
  # fit poisson glm (for count data)
  seedMod_now <- MASS::glm.nb(Num_seeds ~ log_LL_t + N_Site_t, data = seedDat_now)
  ## Flowering probability ($p_b(z)$)
  flwrMod_now <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) + N_Site_t, data = dat_now, family = binomial)))
  ## Distribution of recruit size ($c_o(z')$)
  # subset the data
  recD_now <- dat_now[dat_now$seedling == 1,]
  # fit the model
  recMod_now <- lm(log_LL_t ~ 1 + N_Site_t, data = recD_now)
  
  ## put in the parameter list (paramCont)
  modelList <- list(
    surv     = survMod_now, # growth 
    growth   = sizeMod_now,
    seed     = seedMod_now,
    flwr     = flwrMod_now, 
    recs     = recMod_now              
  )
  
  ## put AICs for each model in a list
  AIClist <- list(
    surv     = AIC(survMod_now),  
    growth   = AIC(sizeMod_now),
    seed     = AIC(seedMod_now),
    flwr     = AIC(flwrMod_now), 
    recs     = AIC(recMod_now)
  )
  
  subPop_modList_DD[[i]] <- modelList
  subPop_AIClist_DD[[i]] <- AIClist
}
names(subPop_modList_DD) <- unique(dat_all$Site)
names(subPop_AIClist_DD) <- unique(dat_all$Site)

## compare AICs
subPop_AIC_compare <- rbind(as.data.frame(subPop_AIClist_DI), 
                            as.data.frame(subPop_AIClist_DD))
subPop_AIC_compare[3,] <- 
  subPop_AIC_compare[1,] - subPop_AIC_compare[2,]

subPop_AIC_compare$type <- c("DI", "DD", "DI-DD")
## these data go into Table 4

#### compare lambdas and subpop N ####
# read in site-level IPMs
site_IPMs_DI <- readRDS("./intermediate_analysis_Data/site_level_IPMs_allYears/site_level_DI_IPMs.RDS")
site_IPMs_DD <- readRDS("./intermediate_analysis_Data/site_level_IPMs_allYears/site_level_DD_IPMs.RDS")

# calculate lambdas
siteDI_lams <- sapply(site_IPMs_DI, function(x) as.numeric(eigen(x)$values[1]))
siteDI_lams <- data.frame("Site" = names(siteDI_lams), "log_lambda" = log(siteDI_lams), "type" = "DI")
siteDD_lams <- sapply(site_IPMs_DD, function(x) as.numeric(eigen(x)$values[1]))
siteDD_lams <- data.frame("Site" = names(siteDD_lams), "log_lambda" = log(siteDD_lams), "type" = "DD")

site_lams <- rbind(siteDI_lams, siteDD_lams)

## get lambdas from each year and each site
byYear_DI_lams <- data.frame(
  "Site" = c("Unnamed_Creek", "Diamond_Creek", "Crow_Creek", "Meadow", "HQ3", "HQ5"), 
  "Year" = as.factor(c(rep(2018, length.out = 6), rep(2019, length.out = 6))),
  "lambda" = c(0.16, 1.21, 0.62, 0.38, 0.51, 0.82, 0.26, 0.48, 0.46, 0.20, -0.10, -0.19)
                             )

N_site <- unique(dat[,c("Site", "Year", "N_all_t")]) %>% 
  group_by(Site, Year) %>% 
  summarize(N_all_t = sum(N_all_t))

# calculate N_all_tplus1
N_site$N_all_tplus1 <- NA
for (i in 1:length(unique(N_site$Site))) {
  N_site[N_site$Site == unique(N_site$Site)[i], "N_all_tplus1"] <- 
    c(N_site[N_site$Site == unique(N_site$Site)[i], ]$N_all_t[2:3], NA)
}

lam_N <- left_join(site_lams, N_site)

byYear_DI_lams <- byYear_DI_lams %>% 
  left_join(N_site)
# add data from Floyd thesis 
byYear_DI_lams <- rbind( byYear_DI_lams, data.frame(Site = c("Crow_Creek", "Unnamed_Creek", "Diamond_Creek"), 
           Year = as.factor(1993),
           lambda = log(c(1.36, 1.02, 1.24)),
           N_all_t = c(435, 550, 590),
           N_all_tplus1 = c(368, 581, 658)),
data.frame(Site = c("Crow_Creek", "Unnamed_Creek", "Diamond_Creek"), 
           Year = as.factor(1992),
           lambda = log(c(1.23, 2.04, 1.67)),
           N_all_t = c(233, 355, 375),
           N_all_tplus1 = c(435, 550, 590)))

# model comparing log(lambda) and N_all_t, w/ a fixed effect of Site to account for different pop sizes between sites
mod <- (lm(lambda ~ N_all_t + Site, data = byYear_DI_lams))
newdata <- data.frame(N_all_t = rep(seq(min(byYear_DI_lams$N_all_t), max(byYear_DI_lams$N_all_t), length.out = 20), 6), Site = c(rep(unique(byYear_DI_lams$Site)[1],length.out = 20), 
                rep(unique(byYear_DI_lams$Site)[2],length.out = 20), 
                rep(unique(byYear_DI_lams$Site)[3],length.out = 20), 
                rep(unique(byYear_DI_lams$Site)[4],length.out = 20), 
                rep(unique(byYear_DI_lams$Site)[5],length.out = 20), 
                rep(unique(byYear_DI_lams$Site)[6],length.out = 20))
                      )
modPreds <- predict(mod, newdata = newdata, interval = "confidence")
modPreds <- as.data.frame(modPreds)
modPreds$N_all_t <- newdata$N_all_t
modPreds$Site <- newdata$Site
modPreds <- modPreds %>% 
  group_by(N_all_t) %>% 
  summarize(fit = mean(fit), 
            lwr = mean(lwr), 
            upr = mean(upr))
# plot DI subpop-level lambdas for each transition according to size in year_t
logLambda_nt_figure <- ggplot(data = byYear_DI_lams) + 
  #geom_ribbon(aes(x = log(N_all_t), ymin = lwr, ymax = upr), data = modPreds, col = "lightgrey", fill = "lightgrey", alpha = .5) +
  #geom_line(aes(x = log(N_all_t), y = fit), data = modPreds, lty = 2, lwd = 1) +
  geom_point(aes(x = log(N_all_t), y = lambda, col = Site)) + 
  geom_smooth(aes(x = log(N_all_t), y = lambda, col = Site), se = FALSE, method = "lm", lty = 1, lwd = .75) +
  xlab(expression(paste("log(N ", italic(t),")"))) +
  scale_color_brewer(type = "qual", palette = "Dark2") +
  ylab(expression(paste("log(", lambda,"), year ", italic(t)))) +
  theme_classic()

# get slopes of line for each site
slopeDist <- c(coefficients(lm(lambda ~ log(N_all_t), data = byYear_DI_lams[byYear_DI_lams$Site == "Crow_Creek",]))[2],
coefficients(lm(lambda ~ log(N_all_t), data = byYear_DI_lams[byYear_DI_lams$Site == "Diamond_Creek",]))[2],
coefficients(lm(lambda ~ log(N_all_t), data = byYear_DI_lams[byYear_DI_lams$Site == "HQ3",]))[2],
coefficients(lm(lambda ~ log(N_all_t), data = byYear_DI_lams[byYear_DI_lams$Site == "HQ5",]))[2],
coefficients(lm(lambda ~ log(N_all_t), data = byYear_DI_lams[byYear_DI_lams$Site == "Meadow",]))[2],
coefficients(lm(lambda ~ log(N_all_t), data = byYear_DI_lams[byYear_DI_lams$Site == "Unnamed_Creek",]))[2])
names(slopeDist) <- c("Crow_Creek", "Diamond_Creek", "HQ3", "HQ5", "Meadow", "Unnamed_Creek")

logLambda_nt_distFig <- ggplot() +
  geom_ribbon(aes(x = seq(-3.25,1.5,length.out = 60)[41:60], 
                  ymin = 0, 
                  ymax = dnorm(x = c(seq(-3.25,1.5,length.out = 60)), mean = mean(slopeDist), sd = sd(slopeDist))[41:60]), colour = "grey80", fill = "grey80") +
  geom_vline(aes(xintercept = 0), lty = 2, col = "darkgrey") +
  geom_line(aes(x = seq(-3.25,1.5,length.out = 60), y = dnorm(x = c(seq(-3.25,1.5,length.out = 60)), mean = mean(slopeDist), sd = sd(slopeDist)))) + 
  geom_rug((aes(x = slopeDist, col = names(slopeDist))), lwd = 1) + 
  xlab(expression(paste("slope of log(", lambda, ") ~ log(N ",italic(t), ")"))) + 
  ylab("Probability") + 
  theme_classic() +
  scale_colour_brewer(palette = "Dark2", guide = "none") 
# probability of a value higher than 0
pnorm(0, mean = mean(slopeDist), sd = sd(slopeDist), lower.tail = FALSE)

## plot subpop size in year t+1 as a function of subpop size in year t
ntplus1_nt_figure <- ggplot(data = byYear_DI_lams) +
  geom_point(aes(x = log(N_all_t), y = log(N_all_tplus1/N_all_t), col = Site)) +
  geom_smooth(aes(x =  log(N_all_t), y = log(N_all_tplus1/N_all_t), col = Site), se = FALSE, method = "lm", lwd = .75) +
  theme_classic()  +
  scale_color_brewer(type = "qual", palette = "Dark2") +
  xlab(expression(paste("log(N ", italic(t), ")"))) + 
  ylab(expression(paste("log(N ", italic(t+1), ") / log(N ", italic(t), ")")))

slopeDist_2 <- c(coefficients(lm(log(N_all_tplus1/N_all_t) ~ log(N_all_t), data = byYear_DI_lams[byYear_DI_lams$Site == "Crow_Creek",]))[2],
               coefficients(lm(log(N_all_tplus1/N_all_t) ~ log(N_all_t), data = byYear_DI_lams[byYear_DI_lams$Site == "Diamond_Creek",]))[2],
               coefficients(lm(log(N_all_tplus1/N_all_t) ~ log(N_all_t), data = byYear_DI_lams[byYear_DI_lams$Site == "HQ3",]))[2],
               coefficients(lm(log(N_all_tplus1/N_all_t) ~ log(N_all_t), data = byYear_DI_lams[byYear_DI_lams$Site == "HQ5",]))[2],
               coefficients(lm(log(N_all_tplus1/N_all_t) ~ log(N_all_t), data = byYear_DI_lams[byYear_DI_lams$Site == "Meadow",]))[2],
               coefficients(lm(log(N_all_tplus1/N_all_t) ~ log(N_all_t), data = byYear_DI_lams[byYear_DI_lams$Site == "Unnamed_Creek",]))[2])
names(slopeDist_2) <- c("Crow_Creek", "Diamond_Creek", "HQ3", "HQ5", "Meadow", "Unnamed_Creek")

ntplus1_nt_distFig <- ggplot() +
  geom_ribbon(aes(x = seq(-3,.5,length.out = 50)[43:50], 
                  ymin = 0, 
                  ymax = dnorm(x = c(seq(-3,.5,length.out = 50)), mean = mean(slopeDist_2), sd = sd(slopeDist_2))[43:50]), colour = "grey80", fill = "grey80") + 
  geom_vline(aes(xintercept = 0), lty = 2, col = "darkgrey") +
  geom_line(aes(x = seq(-3,.5,length.out = 50), y = dnorm(x = c(seq(-3,.5,length.out = 50)), mean = mean(slopeDist_2), sd = sd(slopeDist_2)))) + 
  geom_rug(aes(x = slopeDist_2, color = names(slopeDist_2)), lwd = 1) + 
  xlab(expression(paste("slope of log(N ",italic(t+1)," / N ", italic(t) ,") ~ log(N", italic(t), ")"))) + 
  ylab("Probability") + 
  theme_classic() +
  scale_colour_brewer(palette = "Dark2", guide = "none") 
# probability of a value higher than 0
pnorm(0, mean = mean(slopeDist_2), sd = sd(slopeDist_2), lower.tail = FALSE)

ggarrange(logLambda_nt_figure, ntplus1_nt_figure, logLambda_nt_distFig, ntplus1_nt_distFig, ncol = 2, nrow = 2, common.legend = TRUE, heights = c(1,.75), legend = "right", align = "v", labels = c("A", "B", "C", "D"))
