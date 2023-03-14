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

#### save the data to file ####
## site-level DI IPM matrices
saveRDS(site_IPMs_DI, file = '/Users/astears/COBP_project/site_level_IPMs/site_level_DI_IPMs')
## site-level DI bootstrap CI data
saveRDS(siteDI_bootCI_lambdas, file = "/Users/astears/COBP_project/site_level_IPMs/site_level_DI_bootCI_lambdas")
saveRDS(siteDI_bootCI_params, file = "/Users/astears/COBP_project/site_level_IPMs/site_level_DI_bootCI_params")


## site-level DD IPM matrices
saveRDS(site_IPMs_DD, file = '/Users/astears/COBP_project/site_level_IPMs/site_level_DD_IPMs')
## site-level DD bootstrap CI data
saveRDS(siteDD_bootCI_lambdas, file = "/Users/astears/COBP_project/site_level_IPMs/site_level_DD_bootCI_lambdas")
saveRDS(siteDD_bootCI_params, file = "/Users/astears/COBP_project/site_level_IPMs/site_level_DD_bootCI_params")

#### calculate 95% confidence intervals for lambdas ####
## DI lambda CIs
## calculate 95% CI for lambdas
DI_lams_CIs <- as.data.frame(lapply(X = siteDI_bootCI_lambdas, FUN = function(x) 
  c(log(mean(x) - 1.96*(sd(x)/sqrt(1000))),log(mean(x) + 1.96*(sd(x)/sqrt(1000))))
  ))
names(DI_lams_CIs) <- unique(dat_all$Site)

## DD lambda CIs
DD_lams_CIs <- as.data.frame(lapply(X = siteDD_bootCI_lambdas, FUN = function(x) 
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

#### project each Population forward in time  ####
## forward 100 years, w/ and w/out DD, including demographic stochasticity (random effect of plot?), starting pop.dist is the stable size dist from the discrete IPM

### Density-Independent simulations
## for Soapstone
#### Deterministic, density-independent IPM with CONTINOUS SEEDLINGS ####
## fit vital rate models
datSoap <- dat_all[dat_all$Location == "Soapstone",]
## Survival ($s(z)$)
survDat <- datSoap[datSoap$flowering == 0,]
survMod_proj <- glmer(survives_tplus1 ~ log_LL_t + (1|Plot_ID), data = survDat, family = binomial)
## Growth ($G(z',z)$)
sizeMod_proj <- lmer(log_LL_tplus1 ~ log_LL_t +  (1|Plot_ID), data = datSoap)
## Number of seeds produced, according to plant size ($b(z)$)
seedDat <- datSoap[datSoap$flowering==1,]
seedMod_proj <- lme4::glmer.nb(Num_seeds ~ log_LL_t +  (1|Plot_ID), data = seedDat)
## Flowering probability ($p_b(z)$)
# using size in current year (w/ squared term)
# logistic glm with log-transformed size_t
flwrMod_proj <- glmer(flowering ~ log_LL_t + I(log_LL_t^2) +  (1|Plot_ID), data = datSoap, family = binomial, control = glmerControl(optimizer = "bobyqa"))
## Distribution of recruit size ($c_o(z')$)
recD <- datSoap[datSoap$seedling==1,]
recMod_proj <- lmer(log_LL_t ~ 1 +  (1|Plot_ID), data = recD)

data_list <- list(
  g_int     = fixef(sizeMod_proj)[1], # growth 
  g_slope   = fixef(sizeMod_proj)[2],
  g_sd      = sd(residuals(sizeMod_proj)),
  s_int     = fixef(survMod_proj)[1], # survival
  s_slope   = fixef(survMod_proj)[2],
  p_b_int   = fixef(flwrMod_proj)[1], #probability of flowering
  p_b_slope = fixef(flwrMod_proj)[2],
  p_b_slope_2 = fixef(flwrMod_proj)[3],
  b_int   = fixef(seedMod_proj)[1], #seed production
  b_slope = fixef(seedMod_proj)[2],
  c_o_int    = fixef(recMod_proj), #recruit size distribution
  c_o_sd    = sd(residuals(recMod_proj)), 
  outSB  = outSB_all,
  staySB = staySB_all,
  goSB   = goSB_all, 
  goCont = goCont_all                  
)

## now get a list of values for the random effects
# random values for the growth model
g_r_int <- unlist(ranef(sizeMod_proj))
# random values for the survival model
s_r_int <- unlist(ranef(survMod_proj))
# random values for the flowering model
p_b_r_int <- unlist(ranef(flwrMod_proj))
# random values for the seed production model
b_r_int <- unlist(ranef(seedMod_proj))
# random values for the recruit size model
c_o_r_int <- unlist(ranef(recMod_proj))

# name the list of random effects
nms <- paste("r_", 1:9, sep = "")

names(g_r_int)   <- paste('g_', nms, sep = "")
names(s_r_int)   <- paste('s_', nms, sep = "")
names(p_b_r_int) <- paste('p_b_', nms, sep = "")
names(b_r_int) <- paste('b_', nms, sep = "")
names(c_o_r_int) <- paste('c_o_', nms, sep = "")

# Each set of parameters is converted to a named list. The names should match
# the variables referenced in each define_kernel() call.
g_params   <- as.list(g_r_int)
s_params   <- as.list(s_r_int)
p_b_params <- as.list(p_b_r_int)
b_params <- as.list(b_r_int)
c_o_params <- as.list(c_o_r_int)

# add them all together using c()
all_params_list <- c(data_list, g_params, s_params, p_b_params, b_params, c_o_params)

# initial population state -- (IPM matrix is called "mat_all_Soap")
## calculate the stable size distribution
tempStableDist<-  Re(eigen(mat_all_Soap)$vectors[,1])/sum( Re(eigen(mat_all_Soap)$vectors[,1])) 
init_size_state <- tempStableDist[2:501]

## iterate through this IPM 1000 times
# make a list to hold the output
soapProjIPMs_DI <- list()
for (i in 1:1000) {temp_IPM <- init_ipm(sim_gen   = "general", 
                                        di_dd     = "di", 
                                        det_stoch = "stoch",
                                        kern_param = "kern") %>% 
  define_kernel(
    name          = "P_plot",
    formula       =(1-p_b_plot) * s_plot * g_plot * d_size,
    
    s_plot            = 1/(1 + exp(-(s_int + s_slope * size_1 + s_r_plot))),
    g_plot            = dnorm(size_2, g_mu., g_sd), 
    g_mu.         = g_int + g_slope * size_1 + g_r_plot, 
    p_b_plot          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2) + p_b_r_plot))),
    
    family        = "CC",
    data_list     = all_params_list,
    states        = list(c('size')),
    uses_par_sets = TRUE,
    par_set_indices = list(plot = 1:9),
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", target = "g_plot")
  )  %>% 
  define_kernel(
    name          = "F_plot", 
    formula       = goCont. * (p_b_plot * b_plot * c_o_plot * d_size),
    
    p_b_plot          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2) + p_b_r_plot))),
    b_plot            = exp(b_int + b_slope * size_1 + b_r_plot),
    c_o_plot          = dnorm(size_2, mean =c_o_mu., sd = c_o_sd ),
    c_o_mu.        = c_o_int + c_o_r_plot,
    goCont.       = goCont,
    
    family        = "CC",
    data_list     = all_params_list,
    states        = list(c('size')),
    uses_par_sets = TRUE,
    par_set_indices = list(plot = 1:9),
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", target = "c_o_plot")
  ) %>% define_kernel(
    name          = "seedbank_to_continuous_plot", 
    formula       = outSB. * c_o_plot ,
    
    c_o_plot      = dnorm(size_2, mean = c_o_mu., sd = c_o_sd ),
    c_o_mu.        = c_o_int + c_o_r_plot,
    outSB.       = outSB,
    
    family        = "DC",
    data_list     = all_params_list,
    states        = list(c('size')),
    uses_par_sets = TRUE,
    par_set_indices = list(plot = 1:9),
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", target = "c_o_plot")
  ) %>% define_kernel(
    name          = "seedbank_to_seedbank", 
    formula       = staySB.,
    
    staySB.       = staySB,
    
    family        = "DD",
    data_list     = all_params_list,
    states        = list(c('b')),
    uses_par_sets = FALSE,
    evict_cor     = FALSE
  )   %>% define_kernel(
    name          = "continuous_to_seedbank_plot", 
    formula       = goSB.  * (p_b_plot * b_plot * d_size),
    
    p_b_plot          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2) + p_b_r_plot))),
    b_plot            = exp(b_int + b_slope * size_1 + b_r_plot),
    goSB.       = goSB,
    
    family        = "CD",
    data_list     = all_params_list,
    states        = list(c('size', 'b')),
    uses_par_sets = TRUE,
    par_set_indices = list(plot = 1:9),
    evict_cor     = FALSE
  )  %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P_plot", "F_plot", "seedbank_to_continuous_plot", 
                       "seedbank_to_seedbank", "continuous_to_seedbank_plot"), 
      int_rule = rep("midpoint", 5),
      state_start = c("size", "size", "b", "b", "size"), 
      state_end = c("size", "size", "size", "b", "b")
    )
  ) %>% 
  define_domains(
    size = c(
      min(dat_all$log_LL_t, na.rm = TRUE) * 1.2, # lower bound (L)
      max(dat_all$log_LL_t, na.rm = TRUE) * 1.2, # upper bound (U)
      500 # number of mesh points
    )
  ) %>% 
  define_pop_state(
    n_size = init_size_state,
    n_b = tempStableDist[1], 
  ) %>% 
  make_ipm( iterate = TRUE, iterations = 100,
            normalize_pop_size = TRUE,
            kernel_seq = sample(1:9, 100, replace = TRUE)
  )

## get the vector of lambdas
IPM_out_list <- list(
  stochLambda = suppressMessages(lambda(temp_IPM, type_lambda = "stochastic", burn_in = .15)),
  iterationLambdas = temp_IPM$pop_state$lambda,
  sizeDists <- temp_IPM$pop_state$n_size,
  SBsize <- temp_IPM$pop_state$n_b
)

soapProjIPMs_DI[[i]] <- IPM_out_list
}


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
