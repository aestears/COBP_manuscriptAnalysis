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
fileLoc <- "./intermediate_analysis_Data/"

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
## site-level IPM matrices are read in previously in this script as "IPMs_C_H" (DI models) and 
# "IPMs_I_N" (DD models)

# calculate lambdas
siteDI_lams <- sapply(IPMs_C_H, function(x) as.numeric(eigen(x)$values[1]))
siteDI_lams <- data.frame("Site" = names(siteDI_lams), "log_lambda" = log(siteDI_lams), "type" = "DI")
siteDD_lams <- sapply(IPMs_I_N, function(x) as.numeric(eigen(x)$values[1]))
siteDD_lams <- data.frame("Site" = names(siteDD_lams), "log_lambda" = log(siteDD_lams), "type" = "DD")

site_lams <- rbind(siteDI_lams, siteDD_lams)

## get lambdas from each year and each site
# have to read in IPMs from file ("IPMs_CC_HH" are DI models for 2018-2019; "IPMs_II_NN" are DI models for 2019-2020)
fileLoc <- "./intermediate_analysis_Data/"
## site-level DI IPM matrices for 2018-2019
IPMs_CC_HH <- readRDS(file = paste0(fileLoc,"/IPMs_CC_HH.RDS"))
## site-level DI IPM matrices for 2019-2020
IPMs_II_NN <- readRDS(file = paste0(fileLoc,"/IPMs_II_NN.RDS"))

## calculate lambdas
# 2018-2019
siteDI_firstYr_lams <- sapply(IPMs_CC_HH, function(x) as.numeric(eigen(x$KMatrix)$values[1]))
siteDI_firstYr_lams <- data.frame("Site" = names(siteDI_firstYr_lams), "log_lambda" = unname(log(siteDI_firstYr_lams)), "type" = "DI", "Year" = 2018)
# 2019-2020
siteDI_secondYr_lams <- sapply(IPMs_II_NN, function(x) as.numeric(eigen(x$KMatrix)$values[1]))
siteDI_secondYr_lams <- data.frame("Site" = names(siteDI_secondYr_lams), "log_lambda" = unname(log(siteDI_secondYr_lams)), "type" = "DI", "Year" = 2019)

byYear_DI_lams <- rbind(siteDI_firstYr_lams, siteDI_secondYr_lams)
# update names for IPMs (don't need year identifer, since it's in a separate column)
byYear_DI_lams[byYear_DI_lams$Site %in% c("Crow_Creek_18_19", "Crow_Creek_19_20"), "Site"] <- "Crow_Creek"
byYear_DI_lams[byYear_DI_lams$Site %in% c("Diamond_Creek_18_19", "Diamond_Creek_19_20"), "Site"] <- "Diamond_Creek"
byYear_DI_lams[byYear_DI_lams$Site %in% c("Unnamed_Creek_18_19", "Unnamed_Creek_19_20"), "Site"] <- "Unnamed_Creek"
byYear_DI_lams[byYear_DI_lams$Site %in% c("HQ5_18_19", "HQ5_19_20"), "Site"] <- "HQ5"
byYear_DI_lams[byYear_DI_lams$Site %in% c("HQ3_18_19", "HQ3_19_20"), "Site"] <- "HQ3"
byYear_DI_lams[byYear_DI_lams$Site %in% c("Meadow_18_19", "Meadow_19_20"), "Site"] <- "Meadow"

# get values for population size for each site in each year
N_site <- unique(dat_all[,c("Site", "Year", "N_Site_t")]) 

# calculate N_all_tplus1
N_site$N_Site_tplus1 <- NA
for (i in 1:length(unique(N_site$Site))) {
  N_site[N_site$Site == unique(N_site$Site)[i], "N_Site_tplus1"] <- 
    c(N_site[N_site$Site == unique(N_site$Site)[i], ]$N_Site_t[2:3], NA)
}

lam_N <- left_join(site_lams, N_site)

byYear_DI_lams <- byYear_DI_lams %>% 
  left_join(N_site)

# add data from Floyd thesis 
byYear_DI_lams <- rbind( byYear_DI_lams, data.frame(Site = c("Crow_Creek", "Unnamed_Creek", "Diamond_Creek"), 
           log_lambda = log(c(1.36, 1.02, 1.24)),
           type = "DI",
           Year = as.factor(1993),
           N_Site_t = c(435, 550, 590),
           N_Site_tplus1 = c(368, 581, 658)),
data.frame(Site = c("Crow_Creek", "Unnamed_Creek", "Diamond_Creek"), 
           log_lambda = log(c(1.23, 2.04, 1.67)),
           type = "DI",
           Year = as.factor(1992),
           N_Site_t = c(233, 355, 375),
           N_Site_tplus1 = c(435, 550, 590)))

# model comparing log(lambda) and N_Site_t, w/ a fixed effect of Site to account for different pop sizes between sites
mod <- (lm(log_lambda ~ N_Site_t + Site, data = byYear_DI_lams))
newdata <- data.frame(N_Site_t = rep(seq(min(byYear_DI_lams$N_Site_t), max(byYear_DI_lams$N_Site_t), length.out = 20), 6), Site = c(rep(unique(byYear_DI_lams$Site)[1],length.out = 20), 
                rep(unique(byYear_DI_lams$Site)[2],length.out = 20), 
                rep(unique(byYear_DI_lams$Site)[3],length.out = 20), 
                rep(unique(byYear_DI_lams$Site)[4],length.out = 20), 
                rep(unique(byYear_DI_lams$Site)[5],length.out = 20), 
                rep(unique(byYear_DI_lams$Site)[6],length.out = 20))
                      )
modPreds <- predict(mod, newdata = newdata, interval = "confidence")
modPreds <- as.data.frame(modPreds)
modPreds$N_Site_t <- newdata$N_Site_t
modPreds$Site <- newdata$Site
modPreds <- modPreds %>% 
  group_by(N_Site_t) %>% 
  summarize(fit = mean(fit), 
            lwr = mean(lwr), 
            upr = mean(upr))
# plot DI subpop-level lambdas for each transition according to size in year_t


# get slopes of line for each site
slopeDist <- c(coefficients(lm(log_lambda ~ log(N_Site_t), data = byYear_DI_lams[byYear_DI_lams$Site == "Crow_Creek",]))[2],
               coefficients(lm(log_lambda ~ log(N_Site_t), data = byYear_DI_lams[byYear_DI_lams$Site == "Diamond_Creek",]))[2],
               coefficients(lm(log_lambda ~ log(N_Site_t), data = byYear_DI_lams[byYear_DI_lams$Site == "HQ3",]))[2],
               coefficients(lm(log_lambda ~ log(N_Site_t), data = byYear_DI_lams[byYear_DI_lams$Site == "HQ5",]))[2],
               coefficients(lm(log_lambda ~ log(N_Site_t), data = byYear_DI_lams[byYear_DI_lams$Site == "Meadow",]))[2],
               coefficients(lm(log_lambda ~ log(N_Site_t), data = byYear_DI_lams[byYear_DI_lams$Site == "Unnamed_Creek",]))[2])
names(slopeDist) <- c("Crow_Creek", "Diamond_Creek", "HQ3", "HQ5", "Meadow", "Unnamed_Creek")
slopeDist <- data.frame("Site" = names(slopeDist), "slope" = round(slopeDist,1),
                        "x" = c(5.8, 6.4, 5.4, 6.9, 4.7, 6.5),
                        "y" = c(.645, .8, .4, .6, .5, .33))

(logLambda_nt_figure <- ggplot(data = byYear_DI_lams) + 
  #geom_ribbon(aes(x = log(N_Site_t), ymin = lwr, ymax = upr), data = modPreds, col = "lightgrey", fill = "lightgrey", alpha = .5) +
  #geom_line(aes(x = log(N_Site_t), y = fit), data = modPreds, lty = 2, lwd = 1) +
  geom_smooth(aes(x = log(N_Site_t), y = log_lambda), se = FALSE, method = "lm", lty = 2, linewidth = .75, col = "grey") +
  geom_point(aes(x = log(N_Site_t), y = log_lambda, col = Site)) + 
  geom_smooth(aes(x = log(N_Site_t), y = log_lambda, col = Site), se = FALSE, method = "lm", lty = 1, linewidth = .75) +
  xlab(expression(paste("log(N ", italic(t),")"))) +
  scale_color_brewer(type = "qual", palette = "Dark2") +
  ylab(expression(paste("log(", lambda,"), year ", italic(t)))) +
  theme_classic() + 
  geom_text(data = slopeDist, aes(x = x, y = y, col = Site, label = slope), show.legend = FALSE, fontface = "bold")) 

## plot subpop size in year t+1 as a function of subpop size in year t
slopeDist_2 <- c(coefficients(lm((log(N_Site_tplus1)/log(N_Site_t)) ~ log(N_Site_t), data = byYear_DI_lams[byYear_DI_lams$Site == "Crow_Creek",]))[2],
                 coefficients(lm((log(N_Site_tplus1)/log(N_Site_t)) ~ log(N_Site_t), data = byYear_DI_lams[byYear_DI_lams$Site == "Diamond_Creek",]))[2],
                 coefficients(lm((log(N_Site_tplus1)/log(N_Site_t)) ~ log(N_Site_t), data = byYear_DI_lams[byYear_DI_lams$Site == "HQ3",]))[2],
                 coefficients(lm((log(N_Site_tplus1)/log(N_Site_t)) ~ log(N_Site_t), data = byYear_DI_lams[byYear_DI_lams$Site == "HQ5",]))[2],
                 coefficients(lm((log(N_Site_tplus1)/log(N_Site_t)) ~ log(N_Site_t), data = byYear_DI_lams[byYear_DI_lams$Site == "Meadow",]))[2],
                 coefficients(lm((log(N_Site_tplus1)/log(N_Site_t)) ~ log(N_Site_t), data = byYear_DI_lams[byYear_DI_lams$Site == "Unnamed_Creek",]))[2])
names(slopeDist_2) <- c("Crow_Creek", "Diamond_Creek", "HQ3", "HQ5", "Meadow", "Unnamed_Creek")

slopeDist_2 <- data.frame("Site" = names(slopeDist_2), "slope" = round(slopeDist_2,1),
                        "x" = c(5.75, 6.57, 5.4, 7.3, 4.55, 6.65),
                        "y" = c(1.09, 0.975, .98, 1.04, 0.97, 1.03))


(ntplus1_nt_figure <- ggplot(data = byYear_DI_lams) +
    geom_smooth(aes(x = log(N_Site_t), y =  (log(N_Site_tplus1)/log(N_Site_t))), se = FALSE, method = "lm", lty = 2, linewidth = .75, col = "grey") +
  geom_point(aes(x = log(N_Site_t), y = (log(N_Site_tplus1)/log(N_Site_t)), col = Site)) +
  geom_smooth(aes(x =  log(N_Site_t), y = (log(N_Site_tplus1)/log(N_Site_t)), col = Site), se = FALSE, method = "lm", lwd = .75) +
  theme_classic()  +
  scale_color_brewer(type = "qual", palette = "Dark2") +
  xlab(expression(paste("log(N ", italic(t), ")"))) + 
  ylab(expression(paste("log(N ", italic(t+1), ") / log(N ", italic(t), ")")))+ 
  geom_text(data = slopeDist_2, aes(x = x, y = y, col = Site, label = slope), 
            show.legend = FALSE, fontface = "bold")
) 



DDfigure <- ggarrange(logLambda_nt_figure, ntplus1_nt_figure, ncol = 2, nrow = 1, common.legend = TRUE,  legend = "right", align = "v", labels = c("A", "B"))

ggsave("./figures/densityDependenceFigure.pdf", plot = DDfigure)
