#//////////////////////////
# Integral Projection Models for Oenothera coloradensis: Projecting populations forward in time
# Alice Stears
# 08 December 2021
#/////////////////////////
library(lme4)
library(tidyverse)
library(ipmr)
library(MASS)
library(ggpubr)

# load data 
dat_all <- read.csv(file = "../Processed_Data/allDat_plus_contSeedlings.csv")

#### fit basic IPMs using ipmr ####
# models for each population
# Soapstone data
datSoap <- dat_all[dat_all$Location=="Soapstone",]
survDat_soap <- datSoap[datSoap$flowering==0 | is.na(datSoap$flowering),]
survMod_soap <- glm(survives_tplus1 ~ log_LL_t , data = survDat_soap, family = binomial)
## Growth ($G(z',z)$)
# lm w/ log-transformed size_t and size_t+1
sizeMod_soap <- lm(log_LL_tplus1 ~ log_LL_t , data = datSoap)
## Number of seeds produced, according to plant size ($b(z)$)
seedDat_soap <- datSoap[datSoap$flowering == 1,]
# fit a negative binomial glm (poisson was overdispersed)
seedMod_soap <- MASS::glm.nb(Num_seeds ~ log_LL_t , data = seedDat_soap)
## Flowering probability ($p_b(z)$)
flwrMod_soap <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = datSoap, family = binomial)))
## Distribution of recruit size ($c_o(z')$)
recD_soap <- datSoap[datSoap$seedling == 1,]
recMod_soap <- lm(log_LL_t ~ 1, data = recD_soap)

# FEWAFB data
datBase <- dat_all[dat_all$Location=="FEWAFB",]
survDat_Base <- datBase[datBase$flowering==0 | is.na(datBase$flowering),]
survMod_Base <- glm(survives_tplus1 ~ log_LL_t , data = survDat_Base, family = binomial)
## Growth ($G(z',z)$)
# lm w/ log-transformed size_t and size_t+1
sizeMod_Base <- lm(log_LL_tplus1 ~ log_LL_t , data = datBase)
## Number of seeds produced, according to plant size ($b(z)$)
seedDat_Base <- datBase[datBase$flowering == 1,]
# fit a negative binomial glm (poisson was overdispersed)
seedMod_Base <- MASS::glm.nb(Num_seeds ~ log_LL_t , data = seedDat_Base)
## Flowering probability ($p_b(z)$)
flwrMod_Base <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = datBase, family = binomial)))
## Distribution of recruit size ($c_o(z')$)
recD_Base <- datBase[datBase$seedling == 1,]
recMod_Base <- lm(log_LL_t ~ 1, data = recD_Base)

#seedbank params
germ.rt <- 0.16
viab.rt <- 0.58
surv.seeds <- 0.9
outSB <- germ.rt * surv.seeds 
staySB <- (1-germ.rt) * surv.seeds
goCont <- viab.rt * germ.rt 
goSB <- viab.rt * (1 - germ.rt)

### Deterministic, DI IPM for Base data  
data_list <- list(
  g_int     = coef(sizeMod_Base)[1], # growth 
  g_slope   = coef(sizeMod_Base)[2],
  g_sd      = sd(residuals(sizeMod_Base)),
  s_int     = coef(survMod_Base)[1], # survival
  s_slope   = coef(survMod_Base)[2],
  p_b_int   = coef(flwrMod_Base)[1], #probability of flowering
  p_b_slope = coef(flwrMod_Base)[2],
  p_b_slope_2 = coef(flwrMod_Base)[3],
  b_int   = coef(seedMod_Base)[1], #seed production
  b_slope = coef(seedMod_Base)[2],
  c_o_mu    = coef(recMod_Base), #recruit size distribution
  c_o_sd    = sd(residuals(recMod_Base)), 
  outSB  = outSB_all,
  staySB = staySB_all,
  goSB   = goSB_all, 
  goCont = goCont_all                  
)

# inital population state
init_size_state <- runif(500)

base_IPM <- init_ipm(sim_gen   = "general", 
                     di_dd     = "di", 
                     det_stoch = "det") %>% 
  define_kernel(
    name          = "P",
    formula       =(1-p_b.) * s. * g. * d_size,
    
    s.            = 1/(1 + exp(-(s_int + s_slope * size_1))),
    g.            = dnorm(size_2, g_mu., g_sd), 
    g_mu.         = g_int + g_slope * size_1, 
    p_b.          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2)))),
    
    family        = "CC",
    data_list     = data_list,
    states        = list(c('size')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "g.")
  ) %>% 
  define_kernel(
    name          = "F", 
    formula       = goCont. * (p_b. * b. * c_o. * d_size),
    
    p_b.          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2)))),
    b.            = exp(b_int + b_slope * size_1),
    c_o.          = dnorm(size_2, mean =c_o_mu, sd = c_o_sd ),
    goCont.       = goCont,
    
    family        = "CC",
    data_list     = data_list,
    states        = list(c('size')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "c_o.")
  ) %>% define_kernel(
    name          = "seedbank_to_continuous", 
    formula       = outSB. * c_o. ,
    
    c_o.          = dnorm(size_2, mean =c_o_mu, sd = c_o_sd ),
    outSB.       = outSB,
    
    family        = "DC",
    data_list     = data_list,
    states        = list(c('size')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "c_o.")
  ) %>% define_kernel(
    name          = "seedbank_to_seedbank", 
    formula       = staySB.,
    
    staySB.       = staySB,
    
    family        = "DD",
    data_list     = data_list,
    states        = list(c('b')),
    uses_par_sets = FALSE,
    evict_cor     = FALSE
  )   %>% define_kernel(
    name          = "continuous_to_seedbank", 
    formula       = goSB.  * (p_b. * b. * d_size),
    
    p_b.          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2)))),
    b.            = exp(b_int + b_slope * size_1),
    goSB.       = goSB,
    
    family        = "CD",
    data_list     = data_list,
    states        = list(c('size', 'b')),
    uses_par_sets = FALSE,
    evict_cor     = FALSE
  )  %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "F", "seedbank_to_continuous", "seedbank_to_seedbank", "continuous_to_seedbank"), 
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
    n_size = runif(500),
    n_b = 400, 
    
  ) %>% 
  make_ipm(
    normalize_pop_size = FALSE
    #iterations = 100
  )

ipmr::lambda(base_IPM)

### Deterministic, DI IPM for Soapstone data 
data_list <- list(
  g_int     = coef(sizeMod_soap)[1], # growth 
  g_slope   = coef(sizeMod_soap)[2],
  g_sd      = sd(residuals(sizeMod_soap)),
  s_int     = coef(survMod_soap)[1], # survival
  s_slope   = coef(survMod_soap)[2],
  p_b_int   = coef(flwrMod_soap)[1], #probability of flowering
  p_b_slope = coef(flwrMod_soap)[2],
  p_b_slope_2 = coef(flwrMod_soap)[3],
  b_int   = coef(seedMod_soap)[1], #seed production
  b_slope = coef(seedMod_soap)[2],
  c_o_mu    = coef(recMod_soap), #recruit size distribution
  c_o_sd    = sd(residuals(recMod_soap)), 
  outSB  = outSB_all,
  staySB = staySB_all,
  goSB   = goSB_all, 
  goCont = goCont_all                  
)

# inital population state
init_size_state <- runif(500)

soapstone_IPM <- init_ipm(sim_gen   = "general", 
                          di_dd     = "di", 
                          det_stoch = "det") %>% 
  define_kernel(
    name          = "P",
    formula       =(1-p_b.) * s. * g. * d_size,
    
    s.            = 1/(1 + exp(-(s_int + s_slope * size_1))),
    g.            = dnorm(size_2, g_mu., g_sd), 
    g_mu.         = g_int + g_slope * size_1, 
    p_b.          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2)))),
    
    family        = "CC",
    data_list     = data_list,
    states        = list(c('size')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "g.")
  ) %>% 
  define_kernel(
    name          = "F", 
    formula       = goCont. * (p_b. * b. * c_o. * d_size),
    
    p_b.          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2)))),
    b.            = exp(b_int + b_slope * size_1),
    c_o.          = dnorm(size_2, mean =c_o_mu, sd = c_o_sd ),
    goCont.       = goCont,
    
    family        = "CC",
    data_list     = data_list,
    states        = list(c('size')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "c_o.")
  ) %>% define_kernel(
    name          = "seedbank_to_continuous", 
    formula       = outSB. * c_o. ,
    
    c_o.          = dnorm(size_2, mean =c_o_mu, sd = c_o_sd ),
    outSB.       = outSB,
    
    family        = "DC",
    data_list     = data_list,
    states        = list(c('size')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "c_o.")
  ) %>% define_kernel(
    name          = "seedbank_to_seedbank", 
    formula       = staySB.,
    
    staySB.       = staySB,
    
    family        = "DD",
    data_list     = data_list,
    states        = list(c('b')),
    uses_par_sets = FALSE,
    evict_cor     = FALSE
  )   %>% define_kernel(
    name          = "continuous_to_seedbank", 
    formula       = goSB.  * (p_b. * b. * d_size),
    
    p_b.          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2)))),
    b.            = exp(b_int + b_slope * size_1),
    goSB.       = goSB,
    
    family        = "CD",
    data_list     = data_list,
    states        = list(c('size', 'b')),
    uses_par_sets = FALSE,
    evict_cor     = FALSE
  )  %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "F", "seedbank_to_continuous", "seedbank_to_seedbank", "continuous_to_seedbank"), 
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
    n_size = runif(500),
    n_b = 400, 
    
  ) %>% 
  make_ipm(
    normalize_pop_size = FALSE
    #iterations = 100
  )

ipmr::lambda(soapstone_IPM)

#### projections ####
# each pop forward 100 years, 100 times, including density dependence, random effect of Site (demographic stochasticity) and environmental variation

## soapstone w/ site ranef
## fit vital rate models
datSoap <- dat_all[dat_all$Location == "Soapstone",]
## Survival ($s(z)$)
survDat <- datSoap[datSoap$flowering == 0,]
survMod_proj <- glmer(survives_tplus1 ~ log_LL_t + N_all_t +tMean_grow_C_s  + (1|Site), data = survDat, family = binomial, glmerControl(optimizer = "bobyqa"))
## Growth ($G(z',z)$)
sizeMod_proj <- lmer(log_LL_tplus1 ~ log_LL_t +  N_all_t + tMean_grow_C_s  + (1|Site), data = datSoap)
## Number of seeds produced, according to plant size ($b(z)$)
seedDat <- datSoap[datSoap$flowering==1,]
seedMod_proj <- MASS::glm.nb(Num_seeds ~ log_LL_t +  tMean_grow_C_s + precipWaterYr_cm_s + N_all_t , data = seedDat)
## Flowering probability ($p_b(z)$)
# using size in current year (w/ squared term)
# logistic glm with log-transformed size_t
flwrMod_proj <- glmer(flowering ~ log_LL_t + I(log_LL_t^2) +  tMean_grow_C_s + precipWaterYr_cm_s + N_all_t + (1|Site), data = datSoap, family = binomial, control = glmerControl(optimizer = "bobyqa"))

## Distribution of recruit size ($c_o(z')$)
recD <- datSoap[datSoap$seedling==1,]
recMod_proj <- lm(log_LL_t ~ 1 +  N_all_t + tMean_grow_C_s + precipWaterYr_cm_s, data = recD)

data_list <- list(
  g_int     = fixef(sizeMod_proj)[1], # growth 
  g_slope   = fixef(sizeMod_proj)[2],
  g_dd      = fixef(sizeMod_proj)[3],
  g_tMean   = fixef(sizeMod_proj)[4],
  g_sd      = sd(residuals(sizeMod_proj)),
  s_int     = fixef(survMod_proj)[1], # survival
  s_slope   = fixef(survMod_proj)[2],
  s_dd      = fixef(survMod_proj)[3],
  s_tMean   = fixef(survMod_proj)[4],
  p_b_int   = fixef(flwrMod_proj)[1], #probability of flowering
  p_b_slope = fixef(flwrMod_proj)[2],
  p_b_slope_2 = fixef(flwrMod_proj)[3],
  p_b_tMean = fixef(flwrMod_proj)[4],
  p_b_precip = fixef(flwrMod_proj)[5],
  p_b_dd    = fixef(flwrMod_proj)[6],
  b_int   = coef(seedMod_proj)[1], #seed production
  b_slope = coef(seedMod_proj)[2],
  b_tMean = coef(seedMod_proj)[3],
  b_precip = coef(seedMod_proj)[4],
  b_dd    = coef(seedMod_proj)[5],
  c_o_int    = coef(recMod_proj)[1], #recruit size distribution
  c_o_dd   = coef(recMod_proj)[2],
  c_o_tMean   = coef(recMod_proj)[3],
  c_o_precip   = coef(recMod_proj)[4],
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

# name the list of random effects
nms <- paste("r_", 1:3, sep = "")

names(g_r_int)   <- paste('g_', nms, sep = "")
names(s_r_int)   <- paste('s_', nms, sep = "")
names(p_b_r_int) <- paste('p_b_', nms, sep = "")

# Each set of parameters is converted to a named list. The names should match
# the variables referenced in each define_kernel() call.
g_params   <- as.list(g_r_int)
s_params   <- as.list(s_r_int)
p_b_params <- as.list(p_b_r_int)

# add them all together using c()
all_params_list <- c(data_list, g_params, s_params, p_b_params)

### make environmental parameters
env_params <- list(
  tMean_mu = mean(unique(dat_all$tMean_grow_C_s), na.rm = TRUE),
  tMean_sd = sd(unique(dat_all$tMean_grow_C_s), na.rm = TRUE),
  precip_mu  = mean(unique(dat_all$precipWaterYr_cm_s), na.rm = TRUE),
  precip_sd  = sd(unique(dat_all$precipWaterYr_cm_s), na.rm = TRUE)
)
# define a wrapper function that samples from these distributions
sample_env <- function(env_params) {
  tMean_now  <- rnorm(1, mean = env_params$tMean_mu, sd = env_params$tMean_sd)
  precip_now <- rnorm(1, mean = env_params$precip_mu, sd  = env_params$precip_sd)
  
  out        <- list(tMean = tMean_now, precip = precip_now)
  
  return(out)
  
}

# initial population state -- (IPM matrix is called "mat_all_Soap")
## calculate the stable size distribution
#tempStableDist<-  Re(eigen(mat_all_Soap)$vectors[,1])/sum( Re(eigen(mat_all_Soap)$vectors[,1])) 
init_size_state <- (c(right_ev(soapstone_IPM)$b_w, right_ev(soapstone_IPM)$size_w) / sum(c(right_ev(soapstone_IPM)$b_w, right_ev(soapstone_IPM)$size_w)))[2:501]
init_sb_state <- (c(right_ev(soapstone_IPM)$b_w, right_ev(soapstone_IPM)$size_w) / sum(c(right_ev(soapstone_IPM)$b_w, right_ev(soapstone_IPM)$size_w)))[1]

## make elas and sens matrices
v.dot.w_test <- sum(stable.dist_test * repro.val_test)*h
# calculate the sensitivity function (whole kernel)
sens_test <- outer(repro.val_test,stable.dist_test, '*')/(v.dot.w_test)
# calculate the elasticity function (whole kernel)
elas_test <- matrix(as.vector(sens_test)*as.vector(mat_all_DI)/lambda(contSeedlings_IPM),nrow=501)
## use these instead...


## iterate through this IPM 1000 times
# make a list to hold the output
soapProjIPMs <- list()
for (i in 1:10) {
  temp_IPM <- init_ipm(sim_gen   = "general", 
                       di_dd     = "dd", 
                       det_stoch = "stoch",
                       kern_param = "param") %>% 
    define_kernel(
      name          = "P_site",
      formula       =(1-p_b_site) * s_site * g_site * d_size,
      
      s_site            = 1/(1 + exp(-(s_int + s_slope * size_1 + s_r_site + s_dd * sum(n_size_t) + s_tMean * tMean))),
      g_site            = dnorm(size_2, g_mu., g_sd), 
      g_mu.         = g_int + g_slope * size_1 + g_r_site + g_dd * sum(n_size_t) + g_tMean * tMean, 
      p_b_site          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2) + p_b_r_site + p_b_dd * sum(n_size_t) + p_b_tMean * tMean + p_b_precip * precip))),
      
      family        = "CC",
      data_list     = all_params_list,
      states        = list(c('size')),
      uses_par_sets = TRUE,
      par_set_indices = list(site = 1:3),
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions("norm", target = "g_site")
    )  %>% 
    define_kernel(
      name          = "F_site", 
      formula       = goCont. * (p_b_site * b. * c_o. * d_size),
      
      p_b_site          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2) + p_b_r_site + p_b_dd * sum(n_size_t) + p_b_tMean * tMean + p_b_precip * precip))),
      b.            = exp(b_int + b_slope * size_1 + b_dd * sum(n_size_t) + b_tMean * tMean + b_precip * precip),
      c_o.         = dnorm(size_2, mean =c_o_mu., sd = c_o_sd ),
      c_o_mu.        = c_o_int + c_o_dd * sum(n_size_t) + c_o_tMean * tMean + c_o_precip * precip,
      goCont.       = goCont,
      
      family        = "CC",
      data_list     = all_params_list,
      states        = list(c('size')),
      uses_par_sets = TRUE,
      par_set_indices = list(site = 1:3),
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions("norm", target = "c_o.")
    ) %>% define_kernel(
      name          = "seedbank_to_continuous", 
      formula       = outSB. * c_o. ,
      
      c_o. =        dnorm(size_2, mean =c_o_mu., sd = c_o_sd ),
      c_o_mu.        = c_o_int + c_o_dd * sum(n_size_t) + c_o_tMean * tMean + c_o_precip * precip,
      outSB.       = outSB,
      
      family        = "DC",
      data_list     = all_params_list,
      states        = list(c('size')),
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions("norm", target = "c_o.")
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
      name          = "continuous_to_seedbank_site", 
      formula       = goSB.  * (p_b_site * b. * d_size),
      
      p_b_site          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2) + p_b_r_site + p_b_dd * sum(n_size_t) + p_b_tMean * tMean + p_b_precip * precip))),
      b.            = exp(b_int + b_slope * size_1 + b_dd * sum(n_size_t) + b_tMean * tMean + b_precip * precip),
      goSB.       = goSB,
      
      family        = "CD",
      data_list     = all_params_list,
      states        = list(c('size', 'b')),
      uses_par_sets = TRUE,
      par_set_indices = list(site = 1:3),
      evict_cor     = FALSE
    )  %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("P_site", "F_site", "seedbank_to_continuous", 
                         "seedbank_to_seedbank", "continuous_to_seedbank_site"), 
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
    define_env_state( env_covs   = sample_env(env_params),
      data_list  = list(env_params = env_params,
                        sample_env = sample_env)
    ) %>% 
    define_pop_state(
      n_size = init_size_state,
      n_b = init_sb_state, 
    ) %>% 
    make_ipm( iterate = TRUE, iterations = 100,
              normalize_pop_size = FALSE,
              kernel_seq = sample(1:3, 100, replace = TRUE)
    )
  
  ## get the vector of lambdas
  IPM_out_list <- list(
    stochLambda = suppressMessages(lambda(temp_IPM, type_lambda = "stochastic", burn_in = .15)),
    iterationLambdas = temp_IPM$pop_state$lambda,
    sizeDists <- temp_IPM$pop_state$n_size,
    SBsize <- temp_IPM$pop_state$n_b
  )
  
  soapProjIPMs[[i]] <- IPM_out_list
}

## base w/ site ranef
## fit vital rate models
## fit vital rate models
datBase <- dat_all[dat_all$Location == "FEWAFB",]
## Survival ($s(z)$)
survDat <- datBase[datBase$flowering == 0,]
survMod_proj <- glmer(survives_tplus1 ~ log_LL_t + N_all_t +tMean_grow_C_s  + (1|Site), data = survDat, family = binomial, glmerControl(optimizer = "bobyqa"))
## Growth ($G(z',z)$)
sizeMod_proj <- lmer(log_LL_tplus1 ~ log_LL_t +  N_all_t + tMean_grow_C_s  + (1|Site), data = datBase)
## Number of seeds produced, according to plant size ($b(z)$)
seedDat <- datBase[datBase$flowering==1,]
seedMod_proj <- MASS::glm.nb(Num_seeds ~ log_LL_t +  tMean_grow_C_s + precipWaterYr_cm_s + N_all_t , data = seedDat)
## Flowering probability ($p_b(z)$)
# using size in current year (w/ squared term)
# logistic glm with log-transformed size_t
flwrMod_proj <- glmer(flowering ~ log_LL_t + I(log_LL_t^2) +  tMean_grow_C_s + precipWaterYr_cm_s + N_all_t + (1|Site), data = datBase, family = binomial, control = glmerControl(optimizer = "bobyqa"))

## Distribution of recruit size ($c_o(z')$)
recD <- datBase[datBase$seedling==1,]
recMod_proj <- lm(log_LL_t ~ 1 +  N_all_t + tMean_grow_C_s + precipWaterYr_cm_s, data = recD)

data_list <- list(
  g_int     = fixef(sizeMod_proj)[1], # growth 
  g_slope   = fixef(sizeMod_proj)[2],
  g_dd      = fixef(sizeMod_proj)[3],
  g_tMean   = fixef(sizeMod_proj)[4],
  g_sd      = sd(residuals(sizeMod_proj)),
  s_int     = fixef(survMod_proj)[1], # survival
  s_slope   = fixef(survMod_proj)[2],
  s_dd      = fixef(survMod_proj)[3],
  s_tMean   = fixef(survMod_proj)[4],
  p_b_int   = fixef(flwrMod_proj)[1], #probability of flowering
  p_b_slope = fixef(flwrMod_proj)[2],
  p_b_slope_2 = fixef(flwrMod_proj)[3],
  p_b_tMean = fixef(flwrMod_proj)[4],
  p_b_precip = fixef(flwrMod_proj)[5],
  p_b_dd    = fixef(flwrMod_proj)[6],
  b_int   = coef(seedMod_proj)[1], #seed production
  b_slope = coef(seedMod_proj)[2],
  b_tMean = coef(seedMod_proj)[3],
  b_precip = coef(seedMod_proj)[4],
  b_dd    = coef(seedMod_proj)[5],
  c_o_int    = coef(recMod_proj)[1], #recruit size distribution
  c_o_dd   = coef(recMod_proj)[2],
  c_o_tMean   = coef(recMod_proj)[3],
  c_o_precip   = coef(recMod_proj)[4],
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

# name the list of random effects
nms <- paste("r_", 1:3, sep = "")

names(g_r_int)   <- paste('g_', nms, sep = "")
names(s_r_int)   <- paste('s_', nms, sep = "")
names(p_b_r_int) <- paste('p_b_', nms, sep = "")

# Each set of parameters is converted to a named list. The names should match
# the variables referenced in each define_kernel() call.
g_params   <- as.list(g_r_int)
s_params   <- as.list(s_r_int)
p_b_params <- as.list(p_b_r_int)

# add them all together using c()
all_params_list <- c(data_list, g_params, s_params, p_b_params)

### make environmental parameters
env_params <- list(
  tMean_mu = mean(unique(dat_all$tMean_grow_C_s), na.rm = TRUE),
  tMean_sd = sd(unique(dat_all$tMean_grow_C_s), na.rm = TRUE),
  precip_mu  = mean(unique(dat_all$precipWaterYr_cm_s), na.rm = TRUE),
  precip_sd  = sd(unique(dat_all$precipWaterYr_cm_s), na.rm = TRUE)
)
# define a wrapper function that samples from these distributions
sample_env <- function(env_params) {
  tMean_now  <- rnorm(1, mean = env_params$tMean_mu, sd = env_params$tMean_sd)
  precip_now <- rnorm(1, mean = env_params$precip_mu, sd  = env_params$precip_sd)
  
  out        <- list(tMean = tMean_now, precip = precip_now)
  
  return(out)
  
}

# initial population state -- (IPM matrix is called "mat_all_Soap")
## calculate the stable size distribution
#tempStableDist<-  Re(eigen(mat_all_Soap)$vectors[,1])/sum( Re(eigen(mat_all_Soap)$vectors[,1])) 
init_size_state <- (c(right_ev(base_IPM)$b_w, right_ev(base_IPM)$size_w) / sum(c(right_ev(base_IPM)$b_w, right_ev(base_IPM)$size_w)))[2:501]
init_sb_state <- (c(right_ev(base_IPM)$b_w, right_ev(base_IPM)$size_w) / sum(c(right_ev(base_IPM)$b_w, right_ev(base_IPM)$size_w)))[1]

## iterate through this IPM 1000 times
# make a list to hold the output
baseProjIPMs <- list()
for (i in 1:10) {
  temp_IPM <- init_ipm(sim_gen   = "general", 
                       di_dd     = "dd", 
                       det_stoch = "stoch",
                       kern_param = "param") %>% 
    define_kernel(
      name          = "P_site",
      formula       =(1-p_b_site) * s_site * g_site * d_size,
      
      s_site            = 1/(1 + exp(-(s_int + s_slope * size_1 + s_r_site + s_dd * sum(n_size_t) + s_tMean * tMean))),
      g_site            = dnorm(size_2, g_mu., g_sd), 
      g_mu.         = g_int + g_slope * size_1 + g_r_site + g_dd * sum(n_size_t) + g_tMean * tMean, 
      p_b_site          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2) + p_b_r_site + p_b_dd * sum(n_size_t) + p_b_tMean * tMean + p_b_precip * precip))),
      
      family        = "CC",
      data_list     = all_params_list,
      states        = list(c('size')),
      uses_par_sets = TRUE,
      par_set_indices = list(site = 1:3),
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions("norm", target = "g_site")
    )  %>% 
    define_kernel(
      name          = "F_site", 
      formula       = goCont. * (p_b_site * b. * c_o. * d_size),
      
      p_b_site          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2) + p_b_r_site + p_b_dd * sum(n_size_t) + p_b_tMean * tMean + p_b_precip * precip))),
      b.            = exp(b_int + b_slope * size_1 + b_dd * sum(n_size_t) + b_tMean * tMean + b_precip * precip),
      c_o.         = dnorm(size_2, mean =c_o_mu., sd = c_o_sd ),
      c_o_mu.        = c_o_int + c_o_dd * sum(n_size_t) + c_o_tMean * tMean + c_o_precip * precip,
      goCont.       = goCont,
      
      family        = "CC",
      data_list     = all_params_list,
      states        = list(c('size')),
      uses_par_sets = TRUE,
      par_set_indices = list(site = 1:3),
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions("norm", target = "c_o.")
    ) %>% define_kernel(
      name          = "seedbank_to_continuous", 
      formula       = outSB. * c_o. ,
      
      c_o. =        dnorm(size_2, mean =c_o_mu., sd = c_o_sd ),
      c_o_mu.        = c_o_int + c_o_dd * sum(n_size_t) + c_o_tMean * tMean + c_o_precip * precip,
      outSB.       = outSB,
      
      family        = "DC",
      data_list     = all_params_list,
      states        = list(c('size')),
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions("norm", target = "c_o.")
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
      name          = "continuous_to_seedbank_site", 
      formula       = goSB.  * (p_b_site * b. * d_size),
      
      p_b_site          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2) + p_b_r_site + p_b_dd * sum(n_size_t) + p_b_tMean * tMean + p_b_precip * precip))),
      b.            = exp(b_int + b_slope * size_1 + b_dd * sum(n_size_t) + b_tMean * tMean + b_precip * precip),
      goSB.       = goSB,
      
      family        = "CD",
      data_list     = all_params_list,
      states        = list(c('size', 'b')),
      uses_par_sets = TRUE,
      par_set_indices = list(site = 1:3),
      evict_cor     = FALSE
    )  %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("P_site", "F_site", "seedbank_to_continuous", 
                         "seedbank_to_seedbank", "continuous_to_seedbank_site"), 
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
    define_env_state( env_covs   = sample_env(env_params),
                      data_list  = list(env_params = env_params,
                                        sample_env = sample_env)
    ) %>% 
    define_pop_state(
      n_size = init_size_state,
      n_b = init_sb_state, 
    ) %>% 
    make_ipm( iterate = TRUE, iterations = 100,
              normalize_pop_size = FALSE,
              kernel_seq = sample(1:3, 100, replace = TRUE)
    )
  
  ## get the vector of lambdas
  IPM_out_list <- list(
    stochLambda = suppressMessages(lambda(temp_IPM, type_lambda = "stochastic", burn_in = .15)),
    iterationLambdas = temp_IPM$pop_state$lambda,
    sizeDists <- temp_IPM$pop_state$n_size,
    SBsize <- temp_IPM$pop_state$n_b
  )
  
  baseProjIPMs[[i]] <- IPM_out_list
}

#### projections with drier/hotter climate ####
# Soapstone
## fit vital rate models
datSoap <- dat_all[dat_all$Location == "Soapstone",]
## Survival ($s(z)$)
survDat <- datSoap[datSoap$flowering == 0,]
survMod_proj <- glmer(survives_tplus1 ~ log_LL_t + N_all_t +tMean_grow_C_s  + (1|Site), data = survDat, family = binomial, glmerControl(optimizer = "bobyqa"))
## Growth ($G(z',z)$)
sizeMod_proj <- lmer(log_LL_tplus1 ~ log_LL_t +  N_all_t + tMean_grow_C_s  + (1|Site), data = datSoap)
## Number of seeds produced, according to plant size ($b(z)$)
seedDat <- datSoap[datSoap$flowering==1,]
seedMod_proj <- MASS::glm.nb(Num_seeds ~ log_LL_t +  tMean_grow_C_s + precipWaterYr_cm_s + N_all_t , data = seedDat)
## Flowering probability ($p_b(z)$)
# using size in current year (w/ squared term)
# logistic glm with log-transformed size_t
flwrMod_proj <- glmer(flowering ~ log_LL_t + I(log_LL_t^2) +  tMean_grow_C_s + precipWaterYr_cm_s + N_all_t + (1|Site), data = datSoap, family = binomial, control = glmerControl(optimizer = "bobyqa"))
## Distribution of recruit size ($c_o(z')$)
recD <- datSoap[datSoap$seedling==1,]
recMod_proj <- lm(log_LL_t ~ 1 +  N_all_t + tMean_grow_C_s + precipWaterYr_cm_s, data = recD)

data_list <- list(
  g_int     = fixef(sizeMod_proj)[1], # growth 
  g_slope   = fixef(sizeMod_proj)[2],
  g_dd      = fixef(sizeMod_proj)[3],
  g_tMean   = fixef(sizeMod_proj)[4],
  g_sd      = sd(residuals(sizeMod_proj)),
  s_int     = fixef(survMod_proj)[1], # survival
  s_slope   = fixef(survMod_proj)[2],
  s_dd      = fixef(survMod_proj)[3],
  s_tMean   = fixef(survMod_proj)[4],
  p_b_int   = fixef(flwrMod_proj)[1], #probability of flowering
  p_b_slope = fixef(flwrMod_proj)[2],
  p_b_slope_2 = fixef(flwrMod_proj)[3],
  p_b_tMean = fixef(flwrMod_proj)[4],
  p_b_precip = fixef(flwrMod_proj)[5],
  p_b_dd    = fixef(flwrMod_proj)[6],
  b_int   = coef(seedMod_proj)[1], #seed production
  b_slope = coef(seedMod_proj)[2],
  b_tMean = coef(seedMod_proj)[3],
  b_precip = coef(seedMod_proj)[4],
  b_dd    = coef(seedMod_proj)[5],
  c_o_int    = coef(recMod_proj)[1], #recruit size distribution
  c_o_dd   = coef(recMod_proj)[2],
  c_o_tMean   = coef(recMod_proj)[3],
  c_o_precip   = coef(recMod_proj)[4],
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

# name the list of random effects
nms <- paste("r_", 1:3, sep = "")

names(g_r_int)   <- paste('g_', nms, sep = "")
names(s_r_int)   <- paste('s_', nms, sep = "")
names(p_b_r_int) <- paste('p_b_', nms, sep = "")

# Each set of parameters is converted to a named list. The names should match
# the variables referenced in each define_kernel() call.
g_params   <- as.list(g_r_int)
s_params   <- as.list(s_r_int)
p_b_params <- as.list(p_b_r_int)

# add them all together using c()
all_params_list <- c(data_list, g_params, s_params, p_b_params)

### make environmental parameters
env_params <- list(
  tMean_mu = mean((unique(dat_all$tMean_grow_C_s)*1.1), na.rm = TRUE),
  tMean_sd = sd((unique(dat_all$tMean_grow_C_s)), na.rm = TRUE),
  precip_mu  = mean((unique(dat_all$precipWaterYr_cm_s)*1.1), na.rm = TRUE),
  precip_sd  = sd(unique(dat_all$precipWaterYr_cm_s), na.rm = TRUE)
)
# define a wrapper function that samples from these distributions
sample_env <- function(env_params) {
  tMean_now  <- rnorm(1, mean = env_params$tMean_mu, sd = env_params$tMean_sd)
  precip_now <- rnorm(1, mean = env_params$precip_mu, sd  = env_params$precip_sd)
  
  out        <- list(tMean = tMean_now, precip = precip_now)
  
  return(out)
  
}

# initial population state -- (IPM matrix is called "mat_all_Soap")
## calculate the stable size distribution
#tempStableDist<-  Re(eigen(mat_all_Soap)$vectors[,1])/sum( Re(eigen(mat_all_Soap)$vectors[,1])) 
init_size_state <- (c(right_ev(soapstone_IPM)$b_w, right_ev(soapstone_IPM)$size_w) / sum(c(right_ev(soapstone_IPM)$b_w, right_ev(soapstone_IPM)$size_w)))[2:501]
init_sb_state <- (c(right_ev(soapstone_IPM)$b_w, right_ev(soapstone_IPM)$size_w) / sum(c(right_ev(soapstone_IPM)$b_w, right_ev(soapstone_IPM)$size_w)))[1]

## iterate through this IPM 1000 times
# make a list to hold the output
soapProjIPMs_hot <- list()
for (i in 1:10) {
  temp_IPM <- init_ipm(sim_gen   = "general", 
                       di_dd     = "dd", 
                       det_stoch = "stoch",
                       kern_param = "param") %>% 
    define_kernel(
      name          = "P_site",
      formula       =(1-p_b_site) * s_site * g_site * d_size,
      
      s_site            = 1/(1 + exp(-(s_int + s_slope * size_1 + s_r_site + s_dd * sum(n_size_t) + s_tMean * tMean))),
      g_site            = dnorm(size_2, g_mu., g_sd), 
      g_mu.         = g_int + g_slope * size_1 + g_r_site + g_dd * sum(n_size_t) + g_tMean * tMean, 
      p_b_site          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2) + p_b_r_site + p_b_dd * sum(n_size_t) + p_b_tMean * tMean + p_b_precip * precip))),
      
      family        = "CC",
      data_list     = all_params_list,
      states        = list(c('size')),
      uses_par_sets = TRUE,
      par_set_indices = list(site = 1:3),
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions("norm", target = "g_site")
    )  %>% 
    define_kernel(
      name          = "F_site", 
      formula       = goCont. * (p_b_site * b. * c_o. * d_size),
      
      p_b_site          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2) + p_b_r_site + p_b_dd * sum(n_size_t) + p_b_tMean * tMean + p_b_precip * precip))),
      b.            = exp(b_int + b_slope * size_1 + b_dd * sum(n_size_t) + b_tMean * tMean + b_precip * precip),
      c_o.         = dnorm(size_2, mean =c_o_mu., sd = c_o_sd ),
      c_o_mu.        = c_o_int + c_o_dd * sum(n_size_t) + c_o_tMean * tMean + c_o_precip * precip,
      goCont.       = goCont,
      
      family        = "CC",
      data_list     = all_params_list,
      states        = list(c('size')),
      uses_par_sets = TRUE,
      par_set_indices = list(site = 1:3),
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions("norm", target = "c_o.")
    ) %>% define_kernel(
      name          = "seedbank_to_continuous", 
      formula       = outSB. * c_o. ,
      
      c_o. =        dnorm(size_2, mean =c_o_mu., sd = c_o_sd ),
      c_o_mu.        = c_o_int + c_o_dd * sum(n_size_t) + c_o_tMean * tMean + c_o_precip * precip,
      outSB.       = outSB,
      
      family        = "DC",
      data_list     = all_params_list,
      states        = list(c('size')),
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions("norm", target = "c_o.")
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
      name          = "continuous_to_seedbank_site", 
      formula       = goSB.  * (p_b_site * b. * d_size),
      
      p_b_site          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2) + p_b_r_site + p_b_dd * sum(n_size_t) + p_b_tMean * tMean + p_b_precip * precip))),
      b.            = exp(b_int + b_slope * size_1 + b_dd * sum(n_size_t) + b_tMean * tMean + b_precip * precip),
      goSB.       = goSB,
      
      family        = "CD",
      data_list     = all_params_list,
      states        = list(c('size', 'b')),
      uses_par_sets = TRUE,
      par_set_indices = list(site = 1:3),
      evict_cor     = FALSE
    )  %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("P_site", "F_site", "seedbank_to_continuous", 
                         "seedbank_to_seedbank", "continuous_to_seedbank_site"), 
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
    define_env_state( env_covs   = sample_env(env_params),
                      data_list  = list(env_params = env_params,
                                        sample_env = sample_env)
    ) %>% 
    define_pop_state(
      n_size = init_size_state,
      n_b = init_sb_state, 
    ) %>% 
    make_ipm( iterate = TRUE, iterations = 100,
              normalize_pop_size = FALSE,
              kernel_seq = sample(1:3, 100, replace = TRUE)
    )
  
  ## get the vector of lambdas
  IPM_out_list <- list(
    stochLambda = suppressMessages(lambda(temp_IPM, type_lambda = "stochastic", burn_in = .15)),
    iterationLambdas = temp_IPM$pop_state$lambda,
    sizeDists <- temp_IPM$pop_state$n_size,
    SBsize <- temp_IPM$pop_state$n_b
  )
  
  soapProjIPMs_hot[[i]] <- IPM_out_list
}


## Base
## fit vital rate models
datBase <- dat_all[dat_all$Location == "FEWAFB",]
## Survival ($s(z)$)
survDat <- datBase[datBase$flowering == 0,]
survMod_proj <- glmer(survives_tplus1 ~ log_LL_t + N_all_t +tMean_grow_C_s  + (1|Site), data = survDat, family = binomial, glmerControl(optimizer = "bobyqa"))
## Growth ($G(z',z)$)
sizeMod_proj <- lmer(log_LL_tplus1 ~ log_LL_t +  N_all_t + tMean_grow_C_s  + (1|Site), data = datBase)
## Number of seeds produced, according to plant size ($b(z)$)
seedDat <- datBase[datBase$flowering==1,]
seedMod_proj <- MASS::glm.nb(Num_seeds ~ log_LL_t +  tMean_grow_C_s + precipWaterYr_cm_s + N_all_t , data = seedDat)
## Flowering probability ($p_b(z)$)
# using size in current year (w/ squared term)
# logistic glm with log-transformed size_t
flwrMod_proj <- glmer(flowering ~ log_LL_t + I(log_LL_t^2) +  tMean_grow_C_s + precipWaterYr_cm_s + N_all_t + (1|Site), data = datBase, family = binomial, control = glmerControl(optimizer = "bobyqa"))
## Distribution of recruit size ($c_o(z')$)
recD <- datBase[datBase$seedling==1,]
recMod_proj <- lm(log_LL_t ~ 1 +  N_all_t + tMean_grow_C_s + precipWaterYr_cm_s, data = recD)

data_list <- list(
  g_int     = fixef(sizeMod_proj)[1], # growth 
  g_slope   = fixef(sizeMod_proj)[2],
  g_dd      = fixef(sizeMod_proj)[3],
  g_tMean   = fixef(sizeMod_proj)[4],
  g_sd      = sd(residuals(sizeMod_proj)),
  s_int     = fixef(survMod_proj)[1], # survival
  s_slope   = fixef(survMod_proj)[2],
  s_dd      = fixef(survMod_proj)[3],
  s_tMean   = fixef(survMod_proj)[4],
  p_b_int   = fixef(flwrMod_proj)[1], #probability of flowering
  p_b_slope = fixef(flwrMod_proj)[2],
  p_b_slope_2 = fixef(flwrMod_proj)[3],
  p_b_tMean = fixef(flwrMod_proj)[4],
  p_b_precip = fixef(flwrMod_proj)[5],
  p_b_dd    = fixef(flwrMod_proj)[6],
  b_int   = coef(seedMod_proj)[1], #seed production
  b_slope = coef(seedMod_proj)[2],
  b_tMean = coef(seedMod_proj)[3],
  b_precip = coef(seedMod_proj)[4],
  b_dd    = coef(seedMod_proj)[5],
  c_o_int    = coef(recMod_proj)[1], #recruit size distribution
  c_o_dd   = coef(recMod_proj)[2],
  c_o_tMean   = coef(recMod_proj)[3],
  c_o_precip   = coef(recMod_proj)[4],
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

# name the list of random effects
nms <- paste("r_", 1:3, sep = "")

names(g_r_int)   <- paste('g_', nms, sep = "")
names(s_r_int)   <- paste('s_', nms, sep = "")
names(p_b_r_int) <- paste('p_b_', nms, sep = "")

# Each set of parameters is converted to a named list. The names should match
# the variables referenced in each define_kernel() call.
g_params   <- as.list(g_r_int)
s_params   <- as.list(s_r_int)
p_b_params <- as.list(p_b_r_int)

# add them all together using c()
all_params_list <- c(data_list, g_params, s_params, p_b_params)

### make environmental parameters
env_params <- list(
  tMean_mu = mean((unique(dat_all$tMean_grow_C_s)*1.1), na.rm = TRUE),
  tMean_sd = sd((unique(dat_all$tMean_grow_C_s)), na.rm = TRUE),
  precip_mu  = mean((unique(dat_all$precipWaterYr_cm_s)*1.1), na.rm = TRUE),
  precip_sd  = sd(unique(dat_all$precipWaterYr_cm_s), na.rm = TRUE)
)
# define a wrapper function that samples from these distributions
sample_env <- function(env_params) {
  tMean_now  <- rnorm(1, mean = env_params$tMean_mu, sd = env_params$tMean_sd)
  precip_now <- rnorm(1, mean = env_params$precip_mu, sd  = env_params$precip_sd)
  
  out        <- list(tMean = tMean_now, precip = precip_now)
  
  return(out)
  
}

# initial population state -- (IPM matrix is called "mat_all_Soap")
## calculate the stable size distribution
#tempStableDist<-  Re(eigen(mat_all_Soap)$vectors[,1])/sum( Re(eigen(mat_all_Soap)$vectors[,1])) 
init_size_state <- init_size_state <- (c(right_ev(base_IPM)$b_w, right_ev(base_IPM)$size_w) / sum(c(right_ev(base_IPM)$b_w, right_ev(base_IPM)$size_w)))[2:501]
init_sb_state <- (c(right_ev(base_IPM)$b_w, right_ev(base_IPM)$size_w) / sum(c(right_ev(base_IPM)$b_w, right_ev(base_IPM)$size_w)))[1]

## iterate through this IPM 1000 times
# make a list to hold the output
baseProjIPMs_hot <- list()
for (i in 1:10) {
  temp_IPM <- init_ipm(sim_gen   = "general", 
                       di_dd     = "dd", 
                       det_stoch = "stoch",
                       kern_param = "param") %>% 
    define_kernel(
      name          = "P_site",
      formula       =(1-p_b_site) * s_site * g_site * d_size,
      
      s_site            = 1/(1 + exp(-(s_int + s_slope * size_1 + s_r_site + s_dd * sum(n_size_t) + s_tMean * tMean))),
      g_site            = dnorm(size_2, g_mu., g_sd), 
      g_mu.         = g_int + g_slope * size_1 + g_r_site + g_dd * sum(n_size_t) + g_tMean * tMean, 
      p_b_site          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2) + p_b_r_site + p_b_dd * sum(n_size_t) + p_b_tMean * tMean + p_b_precip * precip))),
      
      family        = "CC",
      data_list     = all_params_list,
      states        = list(c('size')),
      uses_par_sets = TRUE,
      par_set_indices = list(site = 1:3),
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions("norm", target = "g_site")
    )  %>% 
    define_kernel(
      name          = "F_site", 
      formula       = goCont. * (p_b_site * b. * c_o. * d_size),
      
      p_b_site          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2) + p_b_r_site + p_b_dd * sum(n_size_t) + p_b_tMean * tMean + p_b_precip * precip))),
      b.            = exp(b_int + b_slope * size_1 + b_dd * sum(n_size_t) + b_tMean * tMean + b_precip * precip),
      c_o.         = dnorm(size_2, mean =c_o_mu., sd = c_o_sd ),
      c_o_mu.        = c_o_int + c_o_dd * sum(n_size_t) + c_o_tMean * tMean + c_o_precip * precip,
      goCont.       = goCont,
      
      family        = "CC",
      data_list     = all_params_list,
      states        = list(c('size')),
      uses_par_sets = TRUE,
      par_set_indices = list(site = 1:3),
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions("norm", target = "c_o.")
    ) %>% define_kernel(
      name          = "seedbank_to_continuous", 
      formula       = outSB. * c_o. ,
      
      c_o. =        dnorm(size_2, mean =c_o_mu., sd = c_o_sd ),
      c_o_mu.        = c_o_int + c_o_dd * sum(n_size_t) + c_o_tMean * tMean + c_o_precip * precip,
      outSB.       = outSB,
      
      family        = "DC",
      data_list     = all_params_list,
      states        = list(c('size')),
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions("norm", target = "c_o.")
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
      name          = "continuous_to_seedbank_site", 
      formula       = goSB.  * (p_b_site * b. * d_size),
      
      p_b_site          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2) + p_b_r_site + p_b_dd * sum(n_size_t) + p_b_tMean * tMean + p_b_precip * precip))),
      b.            = exp(b_int + b_slope * size_1 + b_dd * sum(n_size_t) + b_tMean * tMean + b_precip * precip),
      goSB.       = goSB,
      
      family        = "CD",
      data_list     = all_params_list,
      states        = list(c('size', 'b')),
      uses_par_sets = TRUE,
      par_set_indices = list(site = 1:3),
      evict_cor     = FALSE
    )  %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("P_site", "F_site", "seedbank_to_continuous", 
                         "seedbank_to_seedbank", "continuous_to_seedbank_site"), 
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
    define_env_state( env_covs   = sample_env(env_params),
                      data_list  = list(env_params = env_params,
                                        sample_env = sample_env)
    ) %>% 
    define_pop_state(
      n_size = init_size_state,
      n_b = init_sb_state, 
    ) %>% 
    make_ipm( iterate = TRUE, iterations = 100,
              normalize_pop_size = FALSE,
              kernel_seq = sample(1:3, 100, replace = TRUE)
    )
  
  ## get the vector of lambdas
  IPM_out_list <- list(
    stochLambda = suppressMessages(lambda(temp_IPM, type_lambda = "stochastic", burn_in = .15)),
    iterationLambdas = temp_IPM$pop_state$lambda,
    sizeDists <- temp_IPM$pop_state$n_size,
    SBsize <- temp_IPM$pop_state$n_b
  )
  
  baseProjIPMs_hot[[i]] <- IPM_out_list
}


#### plot the projected lambdas ####
# soapProjIPMs; soapProj_IPMs_hot; baseProjIPMs; baseProjIPMs_hot

# stoch lambdas for soapstone
# for normal climate scenario
soapStochLams <- sapply(X = soapProjIPMs, FUN = function(x) x$stochLambda)
soapIteratelams <- as.data.frame(sapply(X = soapProjIPMs, FUN = function(x) x$iterationLambdas))
names(soapIteratelams) <- paste0("Iteration_",1:10)
soapIteratelams$Year <- 1:100
soapIteratelams <- pivot_longer(soapIteratelams, cols = 1:10, names_to = "Iteration", values_to = "iterateLam")
soapIteratelams$iterateLam <- log(soapIteratelams$iterateLam)
# for hotter climate scenario
soapStochLams_hot <- sapply(X = soapProjIPMs_hot, FUN = function(x) x$stochLambda)
soapIteratelams_hot <- as.data.frame(sapply(X = soapProjIPMs_hot, 
                                            FUN = function(x) x$iterationLambdas))
names(soapIteratelams_hot) <- paste0("Iteration_",1:10)
soapIteratelams_hot$Year <- 1:100
soapIteratelams_hot <- pivot_longer(soapIteratelams_hot, cols = 1:10, names_to = "Iteration", values_to = "iterateLam")
soapIteratelams_hot$iterateLam <- log(soapIteratelams_hot$iterateLam)

# stoch lambdas for base
# for normal climate scenario
baseStochLams <- sapply(X = baseProjIPMs, FUN = function(x) x$stochLambda)
baseIteratelams <- as.data.frame(sapply(X = baseProjIPMs, FUN = function(x) x$iterationLambdas))
names(baseIteratelams) <- paste0("Iteration_",1:10)
baseIteratelams$Year <- 1:100
baseIteratelams <- pivot_longer(baseIteratelams, cols = 1:10, names_to = "Iteration", values_to = "iterateLam")
baseIteratelams$iterateLam <- log(baseIteratelams$iterateLam)
# for hotter climate scenario
baseStochLams_hot <- sapply(X = baseProjIPMs_hot, FUN = function(x) x$stochLambda)
baseIteratelams_hot <- as.data.frame(sapply(X = baseProjIPMs_hot, 
                                            FUN = function(x) x$iterationLambdas))
names(baseIteratelams_hot) <- paste0("Iteration_",1:10)
baseIteratelams_hot$Year <- 1:100
baseIteratelams_hot <- pivot_longer(baseIteratelams_hot, cols = 1:10, names_to = "Iteration", values_to = "iterateLam")
baseIteratelams_hot$iterateLam <- log(baseIteratelams_hot$iterateLam)

soapIteratelams$Climate <- "Reference"
soapIteratelams$Population <- "Soapstone"
soapIteratelams_hot$Climate <- "hot/dry"
soapIteratelams_hot$Population <- "Soapstone"
baseIteratelams$Climate <- "Reference"
baseIteratelams$Population <- "FEWAFB"
baseIteratelams_hot$Climate <- "hot/dry"
baseIteratelams_hot$Population <- "FEWAFB"

# put all of the iteration lambdas in the same d.f
iterateLams <- rbind(baseIteratelams, soapIteratelams, baseIteratelams_hot, soapIteratelams_hot)

# add mean stochastic lambdas to this d.f 
stochLams_long <- rbind(data.frame("Population" = "Soapstone", "Climate" = "Reference", "logLambda_s" = soapStochLams), 
                        data.frame("Population" = "Soapstone", "Climate" = "hot/dry", "logLambda_s" = soapStochLams_hot), 
                        data.frame("Population" = "FEWAFB", "Climate" = "Reference", "logLambda_s" = baseStochLams), 
                        data.frame("Population" = "FEWAFB", "Climate" = "hot/dry", "logLambda_s" = baseStochLams_hot))

stochLams <- data.frame("Population" = c("Soapstone", "Soapstone", "FEWAFB", "FEWAFB"), 
                        "Climate" = c("Reference", "hot/dry","Reference", "hot/dry"), 
                        "meanlogLambda_s" = c(mean(soapStochLams), mean(soapStochLams_hot), mean(baseStochLams), mean(baseStochLams_hot)))

stochLams_long <- left_join(stochLams_long, stochLams)
iterateLams <- left_join(iterateLams, stochLams)
# calculate mean iteration lambda for each pop/climate/year
meanItLams <- iterateLams %>% 
  group_by(Population, Climate, Year) %>% 
  summarize(meanIterateLam = mean(iterateLam))
iterateLams <- left_join(iterateLams, meanItLams)

# plot of iteration lambdas
iterationLamFig <- ggplot(data = iterateLams) +
  geom_line(aes(x = Year, y = iterateLam, col = Climate, lty = Iteration), alpha = 0.1) + 
  geom_line(aes(x = Year, y = meanIterateLam, col = Climate)) +
  facet_wrap(.~Population) + 
  scale_linetype_manual(values = c(1,1,1,1,1,1,1,1,1,1), guide = "none") +
  geom_hline(aes(yintercept = 0), col = "grey50", lty = 2) +
  #geom_hline(aes(yintercept = meanlogLambda_s, col = Climate)) +
  ylim(c(-.2,1.25)) +
  theme_classic() + 
  ylab(expression(paste("log(", lambda, ")")))

# density of stochastic lambdas 
stochLamDensitFig <- ggplot(data = stochLams_long) +
  geom_density(aes(logLambda_s, col = Climate)) + 
  geom_vline(aes(xintercept = meanlogLambda_s, col = Climate), lty = 2) + 
  facet_wrap(.~Population, ncol = 1) +
  scale_color_discrete(guide = "none") +
  xlab(expression(paste("log(", lambda[s], ")"))) +
  theme_classic() + 
  theme(plot.margin = margin(50, 15, 50, 15))

## plot of population sizes (
# soapstone 
soapN <- as.data.frame(sapply(soapProjIPMs, function(x) colSums(x[[3]])))
# get mean final population size
soapN_tot <- mean(as.numeric(soapN[100,]))
soapN_hot <- as.data.frame(sapply(soapProjIPMs_hot, function(x) colSums(x[[3]]) ))
# get mean final population size
soapN_hot_tot <- mean(as.numeric(soapN_hot[100,]))
# base
baseN <- as.data.frame(sapply(baseProjIPMs, function(x) colSums(x[[3]])))
baseN_tot <- mean(as.numeric(baseN[100,]))
baseN_hot <- as.data.frame(sapply(baseProjIPMs_hot, function(x) colSums(x[[3]])))
baseN_hot_tot <- mean(as.numeric(baseN_hot[100,]))

names(soapN) <- paste0("Iteration_", 1:10)
soapN$meanN <- apply(soapN, MARGIN = 1, FUN = mean )
soapN$Population <- "Soapstone"
soapN$Climate <- "Reference"
soapN$Year <- 1:101
soapN <- pivot_longer(soapN, cols = paste0("Iteration_", 1:10), names_to = "Iteration", values_to = "N")
names(soapN_hot) <- paste0("Iteration_", 1:10)
soapN_hot$meanN <- apply(soapN_hot, MARGIN = 1, FUN = mean )
soapN_hot$Population <- "Soapstone"
soapN_hot$Climate <- "hot/dry"
soapN_hot$Year <- 1:101
soapN_hot <- pivot_longer(soapN_hot, cols = paste0("Iteration_", 1:10), names_to = "Iteration", values_to = "N")
names(baseN) <- paste0("Iteration_", 1:10)
baseN$meanN <- apply(baseN, MARGIN = 1, FUN = mean )
baseN$Population <- "FEWAFB"
baseN$Climate <- "Reference"
baseN$Year <- 1:101
baseN <- pivot_longer(baseN, cols = paste0("Iteration_", 1:10), names_to = "Iteration", values_to = "N")
names(baseN_hot) <- paste0("Iteration_", 1:10)
baseN_hot$meanN <- apply(baseN_hot, MARGIN = 1, FUN = mean )
baseN_hot$Population <- "FEWAFB"
baseN_hot$Climate <- "hot/dry"
baseN_hot$Year <- 1:101
baseN_hot <- pivot_longer(baseN_hot, cols = paste0("Iteration_", 1:10), names_to = "Iteration", values_to = "N")

datN <- rbind(soapN, soapN_hot, baseN, baseN_hot)
datN$meanFinalN <- c(rep(soapN_tot, length.out = nrow(soapN)), 
                     rep(soapN_hot_tot, length.out = nrow(soapN_hot)),
                     rep(baseN_tot, length.out = nrow(baseN)),
                     rep(baseN_hot_tot, length.out = nrow(baseN_hot)))

popSizeFig <- ggplot(datN) +
  geom_line(aes(x = Year, y = N, col = Climate, lty = Iteration), alpha = 0.2) + 
  geom_line(aes(x = Year, y = meanN, col = Climate)) +
  facet_wrap(.~Population) +
  scale_linetype_manual(guide = "none", values = rep(1, length.out = 10)) +
  ylab("Population Size") +
  theme_classic()

#### make a multipanel figure ####
vertPlots <- ggarrange(iterationLamFig, popSizeFig, ncol = 1, nrow = 2, align = "v", legend = "bottom", common.legend = TRUE, labels = c("A", "B"), heights= c(1.5,1))
ggarrange(vertPlots,  stochLamDensitFig, ncol = 2, nrow = 1, widths = c(1.5, .75), labels = c(NA, "C"), label.y = 0.9, label.x = .115)


## save current projections
saveRDS(soapProjIPMs, file = "/Users/Alice/Dropbox/Grad School/Research/Oenothera coloradensis project/COBP_analysis/intermediate_analysis_Data/Projections/SoapstoneNormalClim_projections.RDS")
saveRDS(soapProjIPMs_hot, file = "/Users/Alice/Dropbox/Grad School/Research/Oenothera coloradensis project/COBP_analysis/intermediate_analysis_Data/Projections/SoapstoneHotterClim_projections.RDS")

saveRDS(baseProjIPMs, file = "/Users/Alice/Dropbox/Grad School/Research/Oenothera coloradensis project/COBP_analysis/intermediate_analysis_Data/Projections/BaseNormalClim_projections.RDS")
saveRDS(baseProjIPMs_hot, file = "/Users/Alice/Dropbox/Grad School/Research/Oenothera coloradensis project/COBP_analysis/intermediate_analysis_Data/Projections/BaseHotterClim_projections.RDS")
