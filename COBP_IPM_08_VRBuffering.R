#//////////////////////////
# Integral Projection Models for Oenothera coloradensis: Analysis for Vital Rate Buffering
# Alice Stears
# 11 March 2022
#/////////////////////////

library(tidyverse)
library(popbio) 

# load data from script 01
dat_all <- read.csv(file = "../Processed_Data/allDat_plus_contSeedlings.csv")
source("./analysis_scripts/01_VitalRateModels.R")

#### calculate IPMs (for each site and each transition) ####
## code is also in script "05_IPMs_CC_NN.R"
IPMs_CC_HH <- readRDS("./intermediate_analysis_Data/site_level_IPMs_eachYear/IPMs_CC_HH.RDS")
IPMs_II_NN <- readRDS("./intermediate_analysis_Data/site_level_IPMs_eachYear/IPMs_II_NN.RDS")


#### calculate w/ ipmr ####
### IPM CC-HH ###
### DI IPM for each site--first half of data, discrete ###
# Define the lower and upper integration limit
L <-  1.2 * min(dat_all$log_LL_t, na.rm = TRUE) # minimum size
U <-  1.2 * max(dat_all$log_LL_t, na.rm = TRUE) # maximum size

n <-500 # bins

# These are the parameters for the discrete stages
outSB <- outSB_all #SB to continuous stage
staySB <- staySB_all # staying in SB
goCont <- goCont_all # seeds become continuous right away (without going to the seed bank) 
goSB <- goSB_all # seeds go to the seedbank
surv.seeds <-  0.9 # survival of seeds

## make an empty list to hold the IPM kernels 
IPMs_CC_HH_ipmr <- list()
for (i in 1:length(unique(dat_all$Site))) {
  ## get data for this 'current' site
  dat_now <- dat_all[dat_all$Site == unique(dat_all$Site)[i] # data for this site
                     & dat_all$Year == 2018# data for this year
                     ,]
  
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
  recD_now <- dat_all[dat_all$seedling == 1 & dat_all$Year == 2019,]
  # fit the model
  recMod_now <- lm(log_LL_t ~ 1, data = recD_now)
  
  ## put in the parameter list (paramCont)
  param_cont <- list(
    g_int     = coef(sizeMod_now)[1], # growth 
    g_slope   = coef(sizeMod_now)[2],
    g_sd      = summary(sizeMod_now)$sigma,
    s_int     = coef(survMod_now)[1], # survival
    s_slope   = coef(survMod_now)[2],
    p_b_int   = coef(flwrMod_now)[1], #probability of flowering
    p_b_slope = coef(flwrMod_now)[2],
    p_b_slope_2 = coef(flwrMod_now)[3],
    b_int   = coef(seedMod_now)[1], #seed production
    b_slope = coef(seedMod_now)[2],
    c_o_mu    = coef(recMod_now), #recruit size distribution
    c_o_sd    = summary(recMod_now)$sigma,
    outSB  = outSB_all,
    staySB = staySB_all,
    goSB   = goSB_all, 
    goCont = goCont_all                  
  )
  # inital population state
  init_size_state <- runif(500)
  
  ipm_temp <- init_ipm(sim_gen   = "general", 
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
      data_list     = param_cont,
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
      data_list     = param_cont,
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
      data_list     = param_cont,
      states        = list(c('size')),
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions("norm", "c_o.")
    ) %>% define_kernel(
      name          = "seedbank_to_seedbank", 
      formula       = staySB.,
      
      staySB.       = staySB,
      
      family        = "DD",
      data_list     = param_cont,
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
      data_list     = param_cont,
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
      n_b = 200, 
      
    ) %>% 
    make_ipm(
      normalize_pop_size = FALSE,
      iterations = 200,
      return_all_envs = TRUE
    )
  
  IPMs_CC_HH_ipmr[[i]] <- ipm_temp
}
names(IPMs_CC_HH_ipmr) <- paste0(unique(dat_all$Site),"_18_19")
lambdas_IPMs_CC_HH_ipmr <- sapply(IPMs_CC_HH_ipmr, FUN = function(x) ipmr::lambda(x))

### IPMs II-NN ###
### DI IPM for each site--second half of data, discrete ###
# Define the lower and upper integration limit
L <-  1.2 * min(dat_all$log_LL_t, na.rm = TRUE) # minimum size
U <-  1.2 * max(dat_all$log_LL_t, na.rm = TRUE) # maximum size

n <-500 # bins

# These are the parameters for the discrete stages
outSB <- outSB_all #SB to continuous stage
staySB <- staySB_all # staying in SB
goCont <- goCont_all # seeds become continuous right away (without going to the seed bank) 
goSB <- goSB_all # seeds go to the seedbank
surv.seeds <-  0.9 # survival of seeds

## make an empty list to hold the IPM kernels 
IPMs_II_NN_ipmr <- list()
for (i in 1:length(unique(dat_all$Site))) {
  ## get data for this 'current' site
  dat_now <- dat_all[dat_all$Site == unique(dat_all$Site)[i] # data for this site
                     & dat_all$Year == 2019# data for this year
                     ,]
  
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
  recD_now <- dat_all[dat_all$seedling == 1 & dat_all$Year == 2020,]
  # fit the model
  recMod_now <- lm(log_LL_t ~ 1, data = recD_now)
  
  ## put in the parameter list (paramCont)
  param_cont <- list(
    g_int     = coef(sizeMod_now)[1], # growth 
    g_slope   = coef(sizeMod_now)[2],
    g_sd      = summary(sizeMod_now)$sigma,
    s_int     = coef(survMod_now)[1], # survival
    s_slope   = coef(survMod_now)[2],
    p_b_int   = coef(flwrMod_now)[1], #probability of flowering
    p_b_slope = coef(flwrMod_now)[2],
    p_b_slope_2 = coef(flwrMod_now)[3],
    b_int   = coef(seedMod_now)[1], #seed production
    b_slope = coef(seedMod_now)[2],
    c_o_mu    = coef(recMod_now), #recruit size distribution
    c_o_sd    = summary(recMod_now)$sigma,
    outSB  = outSB_all,
    staySB = staySB_all,
    goSB   = goSB_all, 
    goCont = goCont_all                  
  )  # initial population state
  init_size_state <- runif(500)
  
  ipm_temp <- init_ipm(sim_gen   = "general", 
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
      data_list     = param_cont,
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
      data_list     = param_cont,
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
      data_list     = param_cont,
      states        = list(c('size')),
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions("norm", "c_o.")
    ) %>% define_kernel(
      name          = "seedbank_to_seedbank", 
      formula       = staySB.,
      
      staySB.       = staySB,
      
      family        = "DD",
      data_list     = param_cont,
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
      data_list     = param_cont,
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
      n_b = 200, 
      
    ) %>% 
    make_ipm(
      normalize_pop_size = FALSE,
      iterations = 200
    )
  
  IPMs_II_NN_ipmr[[i]] <- ipm_temp
}
names(IPMs_II_NN_ipmr) <- paste0(unique(dat_all$Site),"_19_20")
lambdas_IPMs_II_NN_ipmr <- sapply(IPMs_II_NN_ipmr, FUN = function(x) ipmr::lambda(x))
megaMat_IPMs_II_NN_ipmr <- lapply(IPMs_II_NN_ipmr, FUN = function(x) ipmr::format_mega_kernel(ipm = x, 
           mega_mat = c(seedbank_to_seedbank, continuous_to_seedbank, seedbank_to_continuous, P + F))$mega_matrix)

#### Calculate corrected standard deviation of Vital Rates #### 

# put all of the IPMs into one big list
allIPMs <- c(IPMs_II_NN, IPMs_CC_HH)

# get all of the respective vital rate matrices in their own lists
allKmat <- lapply(allIPMs, FUN = function(x) x$KMatrix)
allGmat <- lapply(allIPMs, FUN = function(x) x$GMatrix)
allSmat <- lapply(allIPMs, FUN = function(x) x$SMatrix)
allFmat <- lapply(allIPMs, FUN = function(x) x$FMatrix)
allStaySB <- lapply(allIPMs, FUN = function(x) x$staySB_vec)
allLeaveSB <- lapply(allIPMs, FUN = function(x) x$leaveSB_vec)
allGoSb <- lapply(allIPMs, FUN = function(x) x$goSB_vec)

# In this section we first extract the underlying vital rates, to then calculate their corrected standard deviation (according to McDonald et al. 2017). The corrections are required because we want to be able to include in the same analysis vital rates such as survival and growth (bounded between 0 and 1) and fecundity (bounded only by 0). The variance of 0-1 vital rates is constrained by a lower and upper limit, therefore these vital rates have to be transformed to free variance from this constraint. 

### Survival
# use the survival function w/ parameters for each IPM to calculate the survival vector for each IPM (want just survival, not survival - probability of flowering)
S.fun <- function(z, paramCont) {
  mu.surv=paramCont$s_int + paramCont$s_slope *z
  return(1/(1 + exp(-(mu.surv))))
}
surv <- sapply(allIPMs, FUN = function(x) 
  S.fun(z = meshp, paramCont = x$params)
)

plot(0,0, xlim = c(0,500), ylim = c(0,1), type = "n")
for (i in 1:length(allIPMs)) {
  lines(surv[,i], col = i)
}

# get the mean value of the vital rate over all of the IPMs (surv.mu is a vector with the mean survival rate for each cell of the IPM)
surv.mu <- apply(surv, 1, mean)
# get the standard deviation of survival over all of the IPMs 
surv.sd <- apply(surv, 1, sd)
# get the corrected sd (use a logit transformation, since it is a probability)
corr.surv.sd <- apply(car::logit(surv,adjust = 0.001), 1 ,sd) # McDonald et al. (2017) used logit transformation on 0-1 vital rates

### Probability of Flowering (Pb(z))
# use the flowering function w/ parameters for each IPM to calculate the P(flowering) vector for each IPM 
# PROBABILITY OF FLOWERING 
FL.fun <- function(z, paramCont) {
  mu.fl = paramCont$p_b_int + paramCont$p_b_slope*z +  paramCont$p_b_slope_2 * (z^2)
  return(1/(1+ exp(-(mu.fl))))
}

flowering <- lapply(allIPMs, FUN = function(x) 
  FL.fun(z = meshp, paramCont = x$params)
)
# make into a full matrix

# continuous part
flowering.cont <- lapply(flowering, FUN = function(x)
  as.matrix(diag(x)))
# add discrete part (contribution of flowering probability to the seedbank)
flowering.disc <- lapply(flowering, FUN = function(x)
  matrix(c(0, x), nrow = 1))

# combine into one matrix
for (i in 1:length(flowering.cont)) {
  flowering[[i]] <- rbind(flowering.disc[[i]], cbind(rep(0, length.out = 500), flowering.cont[[i]]))
}

# don't actually have to do this, since this is a probability, which we can do a logit-transformation on
# # are there any values of 0 in the F matrices? 
# lapply(flowering, FUN = function(x) sum(x==0))
# # Since the correction that we must apply to fecundities is the
# # log-transformation, we encounter problems when the the vital rate has a value
# # of 0. The following if-else loops are needed to deal with fecundity matrices
# # with "0" entries where they should have a positive value (0s are a problem
# # because log transformation can't be done on the value "0"). Therefore, we add
# # a small value (0.01) to the fecundity value which is 0, over all the years of
# # the study. This way, we keep a biological meaning, and we stay consistent.
# flowering <- lapply(flowering, FUN = function(x)
#   replace(x, list = which(x==0, arr.ind = TRUE), values = .0001)
#   )
# Corrected standard deviation for flowering probability (again according to McDonald et al. 2017)
# mean flowering probability
flowering.mu <- apply(simplify2array(flowering), 1:2, mean) # get a matrix that is an average of all of the different Fmats from each IPM
# sd of fecundity
flowering.sd <- apply(simplify2array(flowering), 1:2, sd)
# get the corrected sd (use a logit transformation, since it is a probability)
corr.flowering.sd <- apply(car::logit(simplify2array(flowering),adjust = 0.001), 1:2 ,sd) # McDonald et al. (2017) used logit transformation on 0-1 vital rates

### Seed production (b(z))
# use the flowering function w/ parameters for each IPM to calculate the seed production vector for each IPM 
# SEED PRODUCTION
SDP.fun <- function(z, paramCont) {
  mu.fps=exp(paramCont$b_int + paramCont$b_slope *z)
  return(mu.fps)
}

seedProd <- lapply(allIPMs, FUN = function(x) 
  SDP.fun(z = meshp, paramCont = x$params)
)

# make into a full matrix
# continuous part
seedProd.cont <- lapply(seedProd, FUN = function(x)
  as.matrix(diag(x)))
# add discrete part (contribution of flowering probability to the seedbank)
seedProd.disc <- lapply(seedProd, FUN = function(x)
  matrix(c(0, x), nrow = 1))

# combine into one matrix
for (i in 1:length(seedProd.cont)) {
  seedProd[[i]] <- rbind(seedProd.disc[[i]], cbind(rep(0, length.out = 500), seedProd.cont[[i]]))
}

# are there any values of 0 in the seed production matrices? 
lapply(seedProd, FUN = function(x) sum(x==0))
# Since the correction that we must apply to fecundities is the
# log-transformation, we encounter problems when the the vital rate has a value
# of 0. The following if-else loops are needed to deal with fecundity matrices
# with "0" entries where they should have a positive value (0s are a problem
# because log transformation can't be done on the value "0"). Therefore, we add
# a small value (0.01) to the fecundity value which is 0, over all the years of
# the study. This way, we keep a biological meaning, and we stay consistent.
seedProd <- lapply(seedProd, FUN = function(x)
  replace(x, list = which(x==0, arr.ind = TRUE), values = .0001)
)
# Corrected standard deviation for flowering probability (again according to McDonald et al. 2017)
# mean seedProd probability
seedProd.mu <- apply(simplify2array(seedProd), 1:2, mean) # get a matrix that is an average of all of the different Fmats from each IPM
# sd of fecundity
seedProd.sd <- apply(simplify2array(seedProd), 1:2, sd)
# get the corrected sd (use a logit transformation, since it is a probability)
corr.seedProd.sd <- apply(log(simplify2array(seedProd)), 1:2 ,sd) # McDonald et al. (2017) used logit transformation on 0-1 vital rates

### Growth. %%% I think this is correct? treat it the same way as fecundity? (i.e. values for the entire matrix, don't sum across columns like for survival)
# use the flowering function w/ parameters for each IPM to calculate the seed production vector for each IPM 
# GROWTH (we assume a constant variance)
GR.fun <- function(z,zz, paramCont){
  growth.mu = paramCont$g_int + paramCont$g_slope*z
  return(dnorm(zz, mean = growth.mu, sd = paramCont$g_sd))
}
growth <- lapply(allIPMs, FUN = function(x) 
  h * t(outer(meshp,meshp,GR.fun, paramCont = x$params)) # Growth
)
# mean of growth
growth.mean <- apply(simplify2array(growth), 1:2, mean)
# sd of growth
growth.sd <- apply(simplify2array(growth), 1:2, sd)
# corrected sd of growth (use a logit transformation, since it is a probability)
corr.growth.sd <- apply(car::logit(simplify2array(allGmat), adjust=0.001), 1:2, sd)

# make the growth vector
growthvector <- rep(0, length(MatMeanG[1,]))

for (i in 1:(length(MatMeanG[1,])-1)) {
      growthvector[i] <- MatMeanG[i+1,i] # a vector with only the values of the mean growth matrix that we need
    }
  
# MatMeanG and growthvector will be needed later, for the sensitivity calculations
### staying in the seedbank # values are identical, so don't bother

### leaving the seedbank
# mean of leaving SB
leaveSB.mean <- apply(simplify2array(allLeaveSB), 1:2, mean)
# sd of leaving SB
leaveSB.sd <- apply(simplify2array(allLeaveSB), 1:2, sd)
# corrected sd of leaving SB (use a logit transformation, since it is a probability)
corr.leaveSB.sd <- apply(car::logit(simplify2array(allLeaveSB), adjust=0.001), 1:2, sd) 

### going to the seedbank
# mean of going to SB
goSB.mean <- apply(simplify2array(allGoSb), 1:2, mean)
# sd of going to SB
goSB.sd <- apply(simplify2array(allGoSb), 1:2, sd)
# corrected sd of going to SB (use a log transformation, since it is NOT a probability)
corr.goSB.sd <- apply(log(simplify2array(allGoSb)), 1:2, sd)

##### 3. SENSITIVITY ###################################################

# Calculation of the corrected sensitivity. The first step is the calculation 
# according to Silvertown and Franco (2004), and then we apply a correction according to McDonald et al. (2017)
# calculate the sensitivity based on the matrix from IPM B (called "mat_al) (IPM for all sites and both transitions) (%%% I think this is correct? not too sure?)
# the matrix is called "mat_all_DI" from the "02_IPMs_A_B.R" script

## calculate sensitivity for the entire matrix

# calculate the eigen vectors for the mean matrix (according to Ellner, Childs & Rees 2016--code in their book, pg. 96)
# stable stage distribution
w.z <- Re(eigen(mat_all_DI)$vectors[,1])

# reproductive value
v.z1 <- Re(eigen(t(mat_all_DI))$vectors[,1])

# lambda
lambda <- Re(eigen(mat_all_DI)$values[1])
# calculate meshpoint size "h"
h <- diff(meshp)[1]
# hand-calculate kernel sensitivity 
S <- outer(v.z1, w.z, "*")/sum(v.z1*  w.z * h)

# calculate elasticity as a gut check (should sum to one)
E <- S * (mat_all_DI/h) / lambda
sum(E) * h^2 # (sums to 1!)

## now, calculate vital rate sensitivity following code from Ellner, Childs and Rees 2016 book (pg. 98-100)
## below: functions from IPM calculation
paramCont=list(NULL)
# survival model is called 'survMod_all'
paramCont[[1]]=as.matrix(coef(survMod_all)) # save coefficients 
# growth model is called 'sizeMod_all'
paramCont[[2]]=cbind(as.matrix(coef(sizeMod_all)),sd(residuals(sizeMod_all))) # the third column is for the standard deviation of growth 
# seedling size distribution is a uniform distribution (of exp(size_2)) with a min of 0.1 and a max 0f 3
paramCont[[3]]= cbind(as.matrix(coef(recMod_all)), sd(residuals(recMod_all)))
# model for probability of flowering is flwrMod_all
paramCont[[4]]=as.matrix(coef(flwrMod_all))
# model for seed production per plant (if reproductive) is seedMod_all
paramCont[[5]]=as.matrix(coef(seedMod_all))
# name the paramCont list to keep track of coefficients
names(paramCont) <- c("survival", "growth", "recruitDist", "flowering", "seedProduction")

# Construct an IPM kernel K using the parameters we obtained from the models
# define the vital rate functions
# SURVIVAL:
S.fun <- function(z, paramCont) {
  mu.surv=paramCont[["survival"]]["(Intercept)",] + paramCont[["survival"]]["log_LL_t",]*z
  return(1/(1 + exp(-(mu.surv))))
}
# GROWTH (we assume a constant variance)
GR.fun <- function(z,z1, paramCont){
  growth.mu = paramCont[["growth"]]["(Intercept)",1] + paramCont[["growth"]]["log_LL_t",1]*z
  return(dnorm(z1, mean = growth.mu, sd = paramCont[["growth"]][1,2]))
}
## SEEDLING SIZES (same approach as in growth function)
SDS.fun <- function(z1, paramCont){
  rec_mu <- paramCont[["recruitDist"]][1]
  rec_sd <- paramCont[["recruitDist"]][2]
  return(dnorm(z1, mean = rec_mu, sd = rec_sd))
}
# PROBABILITY OF FLOWERING 
FL.fun <- function(z, paramCont) {
  mu.fl = paramCont[["flowering"]][1,] + paramCont[["flowering"]][2,]*z +  paramCont[["flowering"]][3,]* (z^2)
  return(1/(1+ exp(-(mu.fl))))
}
# SEED PRODUCTION
SDP.fun <- function(z, paramCont) {
  mu.fps=exp(paramCont[["seedProduction"]][1,1] + paramCont[["seedProduction"]][2,1]*z)
  return(mu.fps)
}

# Define the lower and upper integration limit
L <-  1.2 * min(dat_all$log_LL_t, na.rm = TRUE) # minimum size
U <-  1.2 * max(dat_all$log_LL_t, na.rm = TRUE) # maximum size

n <-500 # bins

# These are the parameters for the discrete stages
outSB <- outSB_all #SB to continuous stage
staySB <- staySB_all # staying in SB
goCont <- goCont_all # seeds become continuous right away (without going to the seed bank) 
goSB <- goSB_all # seeds go to the seedbank
surv.seeds <-  0.9 # survival of seeds

#### Survival
# Ellner book code
# function: $\frac{\delta\lambda}{\delta s(z_0)} = \frac{\int v(z')(1-Pb(z_o))G(z',z_0)w(z_0)dz'}{\int v(z)w(z)dz}$, where $v(z)$ is the reproductive value distribution and $w(z_0)$ is the stable size distribution
dK_by_ds_z1z <- outer (meshp, meshp,
                      function (z1, z, paramCont) {
                        GR.fun(z, z1, paramCont) * (1 - FL.fun(z, paramCont))
                      }, paramCont
                        )
## %%% then, add a row and column of zeros to represent the seedbank?? otherwise the sensitivity matrix is one row and one column larger?
dK_by_ds_z1z <- cbind(0,rbind(0,dK_by_ds_z1z))
# "then, multiply this element-wise by the sensitivity function from earlier and
# sum over the rows to do the integration, remembering to multiply by the
# meshwidth"
s.sens.z <- apply(S * dK_by_ds_z1z, 2, sum) * h


#code from Maria:  sensitivity of population growth rate (i.e. lambda) to changes in survival rates
# VSS on survival (correction suggested by McDonald et al., to account for 0-1 boundaries in vital rates such as survival and growth)
VSS.surv <- rep(0, length(mat_all_DI[,1]))

# θ*(1-θ)/λ * (dλ/dθ)
# add a '0' to the beginning of the surv.mu and surv.sd vectors (for the seedbank)
surv.mu <- c(0,surv.mu)
surv.sd <- c(0, surv.sd)

VSS.surv <- s.sens.z * ((surv.mu*(1-surv.mu))/lambda)

#### 3b) Fecundity
# For fecundity rates the VSS transformation corresponds to the elasticities
elast     <- elasticity(MatMean)
elast.fec <- elast[which(fec.mu != 0)] # keep only the elasticities of fecundity values
sens.fec= sensitivity(MatMean)
sens.fec <- sens.fec[which(fec.mu != 0)]

#### 3c) Growth
# sensitivity of lambda to growth according to Silvertown and Franco 2004

if(Species_name[k] == "Suricata_suricatta")
{
  sens.growth <- rep(0,length(MatMean[,1]))
  
  for(i in 1:length(MatMean[,1])-1)
  {
    sens.growth[i]<- abs(S[i,i]*(-surv.mu[i])+S[i+1,i]*surv.mu[i]) 
  }
  sens.growth[3] <- sens.growth[2]
  sens.growth[2] <- abs(S[1,1]*(-surv.mu[1])+S[3,1]*surv.mu[1]) # insert here growth from stage 1 to stage 3
}else{
  sens.growth <- rep(0,length(MatMean[,1])-1)
  for (i in 1:length(MatMean[,1])-1)
  {
    sens.growth[i] <- abs(S[i,i]*(-surv.mu[i])+S[i+1,i]*surv.mu[i]) 
  }
  
}
# VSS on growth (following McDonald 2017)
VSS.growth <- rep(0,length(MatMean[,1])-1)
if(sum(MatMeanG) != (dim(MatMeanG)[1]*dim(MatMeanG)[2]))
{
  
  if(Species_name[k] == "Suricata_suricatta")
  {
    VSS.growth    <- rep(0,3)
    VSS.growth[1] <- sens.growth[1]*growthvector[1]*(1-growthvector[1])/lambda(MatMean)
    VSS.growth[2] <- sens.growth[2]*growthvector[2]*(1-growthvector[2])/lambda(MatMean)
    VSS.growth[3] <- sens.growth[3]*growthvector[3]*(1-growthvector[3])/lambda(MatMean)
  }else{
    for(i in 1:length(MatMean[,1])-1)
    {
      VSS.growth[i] <- sens.growth[i]*growthvector[i]*(1-growthvector[i])/lambda(MatMean)
    }
  }
}




temp=data.frame(mean=surv.mu,sd=surv.sd,sd.corr=corr.surv.sd,vss=VSS.surv,vr="surv",sens=sens.surv)


temp=rbind(temp,data.frame(mean=fec.mu[fec.mu>0],sd=fec.sd,sd.corr=corr.fec.sd,vss=elast.fec,vr="fec",sens=sens.fec))



if(sum(diag(MatMeanU)[-dim(MatMeanU)[1]]) > 0){
  
  temp=rbind(temp,data.frame(mean=growth.mean,sd=growth.sd,sd.corr=corr.growth.sd,vss=VSS.growth,vr="gr",sens=sens.growth)) 
  
} 
# calculate lambda (pop.growth rate) and store it
lambdas[k] <- lambda(MatMean)

temp$species=as.character(Species_name[k])
mean.var=rbind(mean.var,temp)


# IPMR sensitivity calculation method  ------------------------------------

#from code from the ipmr_esa workshop materials from Sam Levin's GitHub
#https://github.com/aestears/ipmr_esa ## Perturbation analyses Simple analyses
#of asymptotic dynamics are rarely the final piece in using IPMs. `ipmr` does
#not contain additional higher level functions to compute, for example,
#sensitivity or life expectancy. However, it does contain a variety of helpers
#to extract and format the quantities you need to compute those from an IPM.
#This next section will introduce some of these helpers and show you how to use
#them to compute a few frequently used demographic quantities. The first
#analysis we will run is computing sensitivity and elasticity. We will compute
#it at both the kernel level and the vital rate level. Lower level perturbations
#analyses (_i.e._ parameter level) are possible using slight modifications to
#the latter analysis.

### Sensitivity First, we will do kernel level perturbations. This is a partial
#derivative of $\lambda$ with respect to localized perturbation at $z_0, z'_0$
#in the iteration kernel, $K(z',z)$.

#$S(z'_0, z_0) = \frac{\delta\lambda}{\delta K(z'_0, z_0)} = \frac{v(z'_0)w(z_0)}{\langle v, w \rangle}$

# The only piece of information we have not already extracted from our IPM is
# the $dz$ term, which is required to compute the denominator. We can extract
# that using the `int_mesh()` function like so:

mesh_info <- ipmr::int_mesh(IPMs_CC_HH_ipmr$Crow_Creek_18_19)
d_size       <- mesh_info$d_size



# We can now right a function that takes $v(z')$, $w(z)$, and $dz$ as arguments, and compute the kernel sensitivity surface:

  sens <- function(v_z, w_z, d_z) {
    
    outer(v_z, w_z) / sum(v_z * w_z * d_z)
    
  }
  # calculate eigenvectors
  v_size <- left_ev(IPMs_CC_HH_ipmr$Crow_Creek_18_19)
  w_size <- right_ev(IPMs_CC_HH_ipmr$Crow_Creek_18_19)
  
 ipmr_sens <- sens(v_size$size_v, w_size$size_w, d_size)
 
# compare to popbio package
 # get undiscretized iteration matrix
 ipmr_Crow_MegaMat <- ipmr::format_mega_kernel(ipm = (IPMs_CC_HH_ipmr$Crow_Creek_18_19), 
                          mega_mat = c(seedbank_to_seedbank, continuous_to_seedbank, seedbank_to_continuous, P + F))$mega_matrix
 popbio_sens <- popbio::sensitivity(ipmr_Crow_MegaMat[2:201, 2:201])

 w_popbio <- popbio::eigen.analysis(ipmr_Crow_MegaMat)$stable.stage
 v_popbio <- popbio::eigen.analysis(ipmr_Crow_MegaMat)$repro.value
 
 w_test <- eigen(ipmr_Crow_MegaMat)$vectors
 v_test <- Conj(solve(eigen(ipmr_Crow_MegaMat)$vectors))
 
 sens_test <- Re(v_test[1,] %*% t(w_test[,1]))
 
 elas_test <- (1/(Re(eigen(ipmr_Crow_MegaMat)$values[1]))) * sens_test * ipmr_Crow_MegaMat
 elas_popbio <- popbio::eigen.analysis(ipmr_Crow_MegaMat)$elasticities
 
 sum(elas_popbio)
 sum(elas_test) 
 ## difference is negligible (rounding errors, probably)

### Lower level perturbations
                                                                                                                                                                                                              
#Function value perturbations allow us to ask questions like *what happens if we
#tweak survival of size $z_0$ individuals?* or *what happens if we increase
#propagule number for the largest individuals?* Function value perturbations are
#expensive to compute by brute force, but fortunately, analytical formulae exist
#to help us out. The general formula for sensitivity from above can also be
#written:

  #1. $\frac{\delta \lambda(\epsilon)}{\delta \epsilon}\Bigg\rvert_{\epsilon =
  #0}  = \frac{\int \int v(z')C(z',z)w(z) dz' dz}{\int v(z)w(z) dz}$ where
  #$C(z',z)$ is a perturbation kernel, and $v$ and $w$ are the left and right
  #eigenvectors, respectively. The perturbation function $\delta_{z_0}(z)$
  #introduces a localized perturbation at size $z_0$. We will start with the
  #example of perturbing survival for $z_0 \in [L, U]$. Recalling that

# 2. $K(z',z) = P(z',z) + F(z',z) = s(z)G(z',z) + r_p(z)r_n(z)r_d(z')r_g$,
# we want to know what happens to $\lambda$ when $s(z)$ is changed. We can write this as 

#3. $K(z',z) + \epsilon \delta_{z_0}(z)G(z',z)$. 

#This gives us the following perturbation kernel:

# 4. $C(z',z) = \delta_{z_0}(z)G(z',z)$.

# Substituting Eq 4 into Eq 1, we get:

# 5. $\frac{\delta \lambda}{\delta s(z_0)} = \frac{\int v(z')G(z',z_0)w(z_0)dz'}{\int v(z)w(z)dz}$

# This looks pretty nasty, but we can re-write is using operator notation, which
# excludes the $z$ and $z'$s from the variables. It looks like this:

# 6. $\frac{(vG) \: \circ \: w}{\langle v,w \rangle}$.

# This still looks nasty, but we can parse it into more manageable language as
# follows: the change in $\lambda$ induced by a small change in $s(z)$ near
# $z_0$ depends on the fraction of individuals of size $z_0$ that would be
# affected (given by $w$), their size after they grow (given by $G$), and their
# reproductive value after growth (given by $v$). This is scaled by the total
# reproductive value of the population (given by $\langle v,w \rangle$).

# Once we have the sensitivity written out, we can also write out the elasticity. The equation for it from above becomes

# $\frac{\delta \text{ log } \lambda}{\delta \text{ log } s(z)} = \frac{s\: \circ \: vG \: \circ \: w}{\lambda \langle v,w \rangle}$

# Next, we will implement this using some helpers from `ipmr`. 

### Sensitivity and elasticity code

# The next helper function that we will introduce is called `vital_rate_funs()`.
# It extracts the function values for each vital rate from an IPM object. We
# have already extracted the left and right eigenvectors for the sensitivity
# analysis above, so we are ready to proceed.

# NB: for simple IPMs `vital_rate_funs()` always returns an $m \times m$
# function representing every value of the function for $z,z'$. These functions
# have not been integrated yet, so they can be used directly in the calculations
# of sensitivity and elasticity, but care must be taken when using them for
# other applications.
 
 ### Sensitivity of lambda to survival 
 
 vr_funs <- vital_rate_funs(IPMs_CC_HH_ipmr$Crow_Creek_18_19)
 
 G   <- vr_funs$P$g.
 s   <- vr_funs$P$s.
 Pb <- vr_funs$P$p_b.
 
 # add rows and columns for the seedbank (are all 0s in this case)
 #G <- cbind(0, rbind(0, G))
 #s <- cbind(0, rbind(0, s))
 #Pb <- cbind(0, rbind(0, Pb))
 
 # Every row of the "s" object will be identical, because s is only a function
 # of z. Therefore, we take the first row to get its univariate form. (same for Pb)
 s   <- s[1, ]
 Pb <- Pb[1,]
 lambda <- ipmr::lambda(IPMs_CC_HH_ipmr$Crow_Creek_18_19)
 
 # this v isn't standardized... is that a problem? also doesn't have values for seedbank...? 
 v_size <- #c(left_ev(IPMs_CC_HH_ipmr$Crow_Creek_18_19)$b_v , 
   left_ev(IPMs_CC_HH_ipmr$Crow_Creek_18_19)$size_v
 #)
 w_size <- #c(right_ev(IPMs_CC_HH_ipmr$Crow_Creek_18_19)$b_w, 
   right_ev(IPMs_CC_HH_ipmr$Crow_Creek_18_19)$size_w
 #)
 d_size <- int_mesh(IPMs_CC_HH_ipmr$Crow_Creek_18_19)$d_size
 
 sens_s <- (left_mult(G * (1-Pb), v_size) * w_size) / (sum(v_size * w_size* d_size))
 elas_s <- (s * left_mult(G* (1-Pb), v_size) * w_size) / (lambda * sum(v_size * w_size * d_size))
 
 plot(sens_s, type = 'l')
 lines(1:200, elas_s, lty = 2)
 
 ### check with a brute force method
 ### Manually compute an intercept sensitivity
 
 # We will start by extracting the vital rate expressions for the model and
 # having a peak at them so that we can get a feel for how `ipmr` handles these.
 # Since we are going to rebuild the model, we also want to pull out the
 # `proto_ipm` object so that we modify that each time. It is stored in the
 # `ipm$proto_ipm` slot of the IPM object. Once we are familiar with them, we
 # can start our perturbing.
 
 
 base_proto <-IPMs_CC_HH_ipmr$Crow_Creek_18_19$proto_ipm
 
 vr_exprs <- vital_rate_exprs(IPMs_CC_HH_ipmr$Crow_Creek_18_19)
 
 print(vr_exprs)
 
 
 # These are raw expressions, not text strings. If we wanted to, for example,
 # perturb the survival function, we could use the following code:

 # new_fun_form() is necessary when assigning new values to a given vital rate expression. 
 # The "why" isn't important for this exercise, we just have to remember to use it. 
 
 vital_rate_exprs(base_proto, kernel = "P", "s") <- new_fun_form( 1/(1 + exp(-(s_int + s_slope * size_1))) + pert)
 
 # This is great for interactive use, but it would be better to have something
 # we can program with. That way, we do not have to re-type this expression for
 # each vital rate. Below is a function that will programatically append the
 # perturbation term. It uses `parse_expr` from `rlang`  instead of
 # `base::parse`, as the latter does not quite produce the format we need.

 
 library(rlang)
 
 append_pert <- function(vr_expr, type = c("sens", "elas")) {
   
   fun <- switch(type, sens = "+", elas = "*")
   
   parse_expr(paste(deparse(vr_expr), fun, "pert"))
   
 }
 
 # The next piece of the puzzle is to select the perturbation magnitude, and add
 # it as a parameter to our data list so that `make_ipm()` can find it when goes
 # to build the model. We will use a very small magnitude for sensitivity.
 
 parameters(base_proto) <- list(pert = 0.0001)
 
 # Now, for the final piece - a function that wraps up the whole process of perturbing a single vital rate and returns the value.
 
 
 sens_vr <- function(proto_ipm, vital_rate, vr_exprs, pert_magnitude, init_lambda) {
   
   # Create the new vital rate expression
   
   new_vr_expr <- append_pert(vr_exprs[[vital_rate]], "sens")
   
   # The vital_rate_exprs setting function requires the name of the kernel
   # that the vital rate appears in. This step scans the proto_ipm object
   # and pulls out the kernel_id that contains the vital rate. Some vital
   # rates may appear in multiple kernels, so this generalizes over that case.
   
   kern_ind    <- vapply(proto_ipm$params, 
                         function(x, vr_nm) {
                           any(vr_nm %in% names(x$vr_text))
                         },
                         logical(1L),
                         vr_nm = vital_rate)
   
   kern_nms    <- proto_ipm$kernel_id[kern_ind]
   
   # Now, we loop over all the kernels that have the vital rate of interest
   # and set the new expression. NB: The "!!" in new_fun_form ensures that the
   # the expression we created is passed properly. The details of why this matters
   # aren't very important, we just have to remember to use this syntax whenever
   # we are trying to program with the function. You DO NOT need to use this
   # when using new_fun_form interactively, like we did above.
   
   for(i in seq_along(kern_nms)) {
     
     vital_rate_exprs(proto_ipm, 
                      kern_nms[i], 
                      vital_rate) <- new_fun_form(!! new_vr_expr)
     
   }
   
   parameters(proto_ipm) <- list(pert = pert_magnitude)
   
   # Rebuild the IPM, then calculate the change of lambda by the perturbation
   # magnitude.
   
   lambda <- make_ipm(proto_ipm, iterations = 100) %>%
     ipmr::lambda()
   
   out <- (lambda - init_lambda) / pert_magnitude
   
   return(out)
   
 }
 
 current_lambda <- ipmr::lambda(IPMs_CC_HH_ipmr$Crow_Creek_18_19)
 
 s_sens <- sens_vr(base_proto, "s", vr_exprs, 0.0001, current_lambda)
 
 
 # We will now create a `data.frame` to store the results for all of the vital
 # rates, loop over them, and plot the results. We will exclude $G(z', z)$ and
 # $r_d(z')$ for this example.

 
 vr_exprs$G <- vr_exprs$r_d <- NULL
 
 all_sens <- data.frame(vr = names(vr_exprs),
                        value = NA)
 
 for(i in seq_along(vr_exprs)) {
   
   all_sens$value[i] <- sens_vr(base_proto, 
                                all_sens$vr[i], 
                                vr_exprs,
                                0.0001, 
                                current_lambda)
   
 }
 
 
 barplot(all_sens$value, names.arg = all_sens$vr)
 
 ```
 
 
# #### read in elasticity data ####
# contParamElas <- readRDS("./intermediate_analysis_Data/allSiteAllYears_noDDnoEnv/continuousParamElasticity.RDS")
# discParamElas <- readRDS("./intermediate_analysis_Data/allSiteAllYears_noDDnoEnv/discreteParamElasticity.RDS")
# # reorder correctly
# contParamElas <- contParamElas[match(c("growth_(Intercept)", "growth_log_LL_t", "growth_stndDev", "survival_(Intercept)", "survival_log_LL_t","flowering_(Intercept)", "flowering_log_LL_t", "flowering_I(log_LL_t^2)","seedProduction_(Intercept)", "seedProduction_log_LL_t", "recruitDist_(Intercept)", "recruitDist_stndDev"), contParamElas$param_name),] 
# ## put elasticity in the vrDF
# vrDF$Elas <- contParamElas$elas_mean
# vrDF$Sens <- contParamElas$sens_mean
# 
# ## add data for sb params? 
# # vrDF[13,"mean"] <- mean(c(13.0, 12, 8.3, 7.0, 5.3)/(45 * seed_per_cap))
# # vrDF[13, "sd"] <- sd(c(13.0, 12, 8.3, 7.0, 5.3)/(45 * seed_per_cap))
# # 
# # vrDF[14,"mean"] <- mean(c(81,61,54)/100 * 0.89)
# # vrDF[14, "sd"] <- sd(c(81,61,54)/100 * 0.89)
# # 
# # vrDF[15, "mean"] <- 0.5
# # vrDF[15, "sd"] <- 0.175
# # simulate the sb params using a uniform distribution 
# rownames(vrDF[13:15,]) <- c("germ.rt", "surv.rt", "viab.rt")
# vrDF[13:15,"quarDisp"] <- c(.75-.25)/(.75+.25)
# 
# vrDF[13:15,"Elas"] <- discParamElas$elas_mean
# vrDF[13:15, "Sens"] <- discParamElas$sens_mean
# 
# # vrDF$CV <- vrDF$sd/vrDF$mean * 100
 
# #### calculate the correlation between CV and elasticity ####
# cor.test((vrDF$quarDisp[c(1:3,5:12)]), abs(vrDF$Elas[c(1:3,5:12)]), method = "pearson")
# ggplot(data = vrDF) +
#   geom_smooth(aes(x = (quarDisp), y = abs(Elas)),data = vrDF, method = "lm", lty = 2, col = "red", alpha = .3, se = FALSE) + 
#   geom_smooth(aes(x = (quarDisp), y = abs(Elas)),data = vrDF[c(1:3,5:12),], method = "lm", lty = 2, col = "blue", alpha = .3, se = FALSE) +
#   geom_point(aes(x = (quarDisp), y = abs(Elas))) +
#   xlab(c("Quartile coefficient of dispersion")) +
#   ylab(c("|Vital Rate Param. Elasticity|")) + 
#   theme_classic() 
# 
# # plot(x = (vrDF$quarDisp), y = vrDF$Elas, 
# #      xlab = "Vital Rate Param. Coefficient of Variation (CV)", 
# #      ylab = "Vital Rate Param. Elasticity", 
# #      pch = 16, 
# #      col = "grey20")
# # abline(lm(vrDF$Elas ~ vrDF$quarDisp), col = "red", lty = 2)
# # abline(lm(Elas ~ quarDisp, data = vrDF[c(1:3,5:12),]), col = "blue", lty = 2)
# # text(x =650, y = -1.5, labels = "r = 0.116 \n P-value = 0.68 \n(t = 0.42, df = 13)")
# 
# #### fit IPMs to get lambdas ####
# outSB <- outSB_all #SB to continuous stage
# staySB <- staySB_all # staying in SB
# goCont <- goCont_all # seeds become continuous right away (without going to the seed bank) 
# goSB <- goSB_all # seeds go to the seedbank
# surv.seeds <-  0.9 # survival of seeds
# # SURVIVAL:
# S.fun <- function(z, paramCont) {
#   mu.surv=paramCont$s_int + paramCont$s_slope *z
#   return(1/(1 + exp(-(mu.surv))))
# }
# # GROWTH (we assume a constant variance)
# GR.fun <- function(z,zz, paramCont){
#   growth.mu = paramCont$g_int + paramCont$g_slope*z
#   return(dnorm(zz, mean = growth.mu, sd = paramCont$g_sd))
# }
# ## SEEDLING SIZES (same approach as in growth function)
# SDS.fun <- function(zz, paramCont){
#   rec_mu <- paramCont$c_o_mu
#   rec_sd <- paramCont$c_o_sd
#   return(dnorm(zz, mean = rec_mu, sd = rec_sd))
# }
# # PROBABILITY OF FLOWERING 
# FL.fun <- function(z, paramCont) {
#   mu.fl = paramCont$p_b_int + paramCont$p_b_slope*z +  paramCont$p_b_slope_2 * (z^2)
#   return(1/(1+ exp(-(mu.fl))))
# }
# # SEED PRODUCTION
# SDP.fun <- function(z, paramCont) {
#   mu.fps=exp(paramCont$b_int + paramCont$b_slope *z)
#   return(mu.fps)
# }
# # for first half of data
# subPop_first_mats <- list()
# for (i in 1:length(subPop_first_VRs)) {
#   paramsNow <- subPop_first_VRs[[i]]
#   ## fit the IPM
#   K <- array(0,c(n+1,n+1))
#   # Setting up the kernels
#   b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
#   meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
#   h=(U-L)/n # bin width 
#   # Survival and growth 
#   S <- diag(S.fun(meshp, paramCont = paramsNow)) # Survival # put survival probabilities in the diagonal of the matrix
#   G <- h * t(outer(meshp,meshp,GR.fun, paramCont = paramsNow)) # Growth
#   #Recruits distribution (seeds recruited from the seedbank into the continuous stage)
#   c_o <- h * matrix(rep(SDS.fun(meshp, paramCont = paramsNow),n),n,n,byrow=F)
#   #Probability of flowering
#   Pb = (FL.fun(meshp, paramCont = paramsNow))
#   #Number of seeds produced according to adult size
#   b_seed = (SDP.fun(meshp, paramCont = paramsNow))
#   FecALL= Pb * b_seed
#   # update the 'S' matrix by multiplying it by (1-Pb), since this is a monocarpic perennial
#   S_new <- S * (1-Pb)
#   # Control for eviction:
#   G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
#   c_o <- c_o/matrix(as.vector(apply(c_o,2,sum)),nrow=n,ncol=n,byrow=TRUE)
#   # make the continuous part of the P matrix
#   Pkernel.cont <- as.matrix(G %*% S_new)
#   # seedbank (first column of your K)
#   Pkernel.seedbank = c(staySB, outSB*c_o[,1]) # seeds survive and go to continuous
#   # Make the full P kernel
#   Pkernel <- cbind(Pkernel.seedbank,rbind(rep(0,length(meshp)),Pkernel.cont)) # discrete component
#   ## make the F kernel
#   Fkernel.cont <-  as.matrix(goCont * ((c_o) %*% diag(FecALL))) 
#   Fkernel.discr  <- matrix(c(0, goSB * (FecALL)), nrow = 1)
#   # multiply the cont_to_disc distribution by the binwidth (h)
#   Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))
#   
#   mat <-Pkernel+Fkernel
#   eigenMat <- eigen(mat)
#   
#   subPop_first_mats[[i]] <- mat
# }
# names(subPop_first_mats) <- names(subPop_first_VRs)
# 
# # for second half of data
# subPop_second_mats <- list()
# for (i in 1:length(subPop_second_VRs)) {
#   paramsNow <- subPop_second_VRs[[i]]
#   ## fit the IPM
#   K <- array(0,c(n+1,n+1))
#   # Setting up the kernels
#   b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
#   meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
#   h=(U-L)/n # bin width 
#   # Survival and growth 
#   S <- diag(S.fun(meshp, paramCont = paramsNow)) # Survival # put survival probabilities in the diagonal of the matrix
#   G <- h * t(outer(meshp,meshp,GR.fun, paramCont = paramsNow)) # Growth
#   #Recruits distribution (seeds recruited from the seedbank into the continuous stage)
#   c_o <- h * matrix(rep(SDS.fun(meshp, paramCont = paramsNow),n),n,n,byrow=F)
#   #Probability of flowering
#   Pb = (FL.fun(meshp, paramCont = paramsNow))
#   #Number of seeds produced according to adult size
#   b_seed = (SDP.fun(meshp, paramCont = paramsNow))
#   FecALL= Pb * b_seed
#   # update the 'S' matrix by multiplying it by (1-Pb), since this is a monocarpic perennial
#   S_new <- S * (1-Pb)
#   # Control for eviction:
#   G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
#   c_o <- c_o/matrix(as.vector(apply(c_o,2,sum)),nrow=n,ncol=n,byrow=TRUE)
#   # make the continuous part of the P matrix
#   Pkernel.cont <- as.matrix(G %*% S_new)
#   # seedbank (first column of your K)
#   Pkernel.seedbank = c(staySB, outSB*c_o[,1]) # seeds survive and go to continuous
#   # Make the full P kernel
#   Pkernel <- cbind(Pkernel.seedbank,rbind(rep(0,length(meshp)),Pkernel.cont)) # discrete component
#   ## make the F kernel
#   Fkernel.cont <-  as.matrix(goCont * ((c_o) %*% diag(FecALL))) 
#   Fkernel.discr  <- matrix(c(0, goSB * (FecALL)), nrow = 1)
#   # multiply the cont_to_disc distribution by the binwidth (h)
#   Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))
#   
#   mat <-Pkernel+Fkernel
#   eigenMat <- eigen(mat)
#   
#   subPop_second_mats[[i]] <- mat
# }
# names(subPop_second_mats) <- names(subPop_second_VRs)
# 
# # calculate lambdas
# log(sapply(subPop_first_mats, function(x) as.numeric(eigen(x)$values[1])))
# 
# log(sapply(subPop_second_mats, function(x) as.numeric(eigen(x)$values[1])))
