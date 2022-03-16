#//////////////////////////
# Integral Projection Models for Oenothera coloradensis: Analysis for demographic compensation
# Alice Stears
# 12 March 2022
#/////////////////////////

library(tidyverse)
library(ipmr)
library(rstatix)

# read in data
dat_all <- read.csv(file = "/Users/astears/COBP_project/allDat_plus_contSeedlings.csv")

for (i in 1:length(unique(dat_all$Site))) {
  dat_now <- dat_all[dat_all$Site == unique(dat_all$Site)[i],]
  ## vital rate models for each sites, including density dependence and all env. covariates 
  ## Survival ($s(z)$)
  survDat <- dat_now[dat_now$flowering==0,]
  # logistic glm with log-transformed size_t
  survMod_e_dd <- glm(survives_tplus1 ~ log_LL_t + tMean_grow_C  + tMean_winter_C  + precipWaterYr_cm + N_all_t , data = survDat, family = binomial)
  ## Growth ($G(z',z)$)
  sizeMod_e_dd <- lm(log_LL_tplus1 ~ log_LL_t  + tMean_grow_C  + tMean_winter_C  + precipWaterYr_cm + N_all_t , data = dat_now)
  ## Number of seeds produced, according to plant size ($b(z)$)
  seedDat <- dat_now[dat_now$flowering == 1,]
  seedMod_e_dd <- MASS::glm.nb(Num_seeds ~ log_LL_t  + tMean_grow_C  + tMean_winter_C  + precipWaterYr_cm + N_all_t , data = seedDat)
  ## Flowering probability ($p_b(z)$)
  flwrMod_e_dd <- glm(flowering ~ log_LL_t + I(log_LL_t^2)  + tMean_grow_C  + tMean_winter_C  + precipWaterYr_cm + N_all_t , data = dat_now, family = binomial)
  ## Distribution of recruit size ($c_o(z')$)
  recD = dat_now[dat_now$seedling==1,]
  recMod_env_dd <- lm(log_LL_t ~ tMean_grow_C  + tMean_winter_C  + precipWaterYr_cm + N_all_t, data = recD)
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

dc_mods$vitalRate <- rep(x = c("Intercept", "log_LL_t", "I(log_LL_t^2)", "tMean_grow_C", "tMean_winter_C", "precipWaterYr_cm", "N_all_t"), times = 6)
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


#### fit IMPs #### 
# use parameters in the 'dc_mods' d.f
for (i in unique(dc_mods$Site)) {
  # get the data for this site
  dat_now <- dc_mods[dc_mods$Site == i,]
  # get the vital rates into the right format
  paramNow <- list(
    g_int = dat_now[dat_now$vitalRate == "Intercept","grow"], 
    g_slope = dat_now[dat_now$vitalRate == "log_LL_t","grow"], 
    g_dd  = dat_now[dat_now$vitalRate == "N_all_t","grow"], 
    g_tGrow = dat_now[dat_now$vitalRate == "tMean_grow","grow"], 
    g_sd = 0.4837,
    s_int = dat_now[dat_now$vitalRate == "Intercept","surv"], 
    s_slope = dat_now[dat_now$vitalRate == "log_LL_t","surv"], 
    s_dd = dat_now[dat_now$vitalRate == "N_all_t","surv"], 
    s_tGrow = dat_now[dat_now$vitalRate == "tMean_grow","surv"], 
    p_b_int = dat_now[dat_now$vitalRate == "Intercept","flwr"], 
    p_b_slope = dat_now[dat_now$vitalRate == "log_LL_t","flwr"], 
    p_b_slope_2 = dat_now[dat_now$vitalRate == "I(log_LL_t^2)","flwr"], 
    p_b_dd = dat_now[dat_now$vitalRate == "N_all_t","flwr"], 
    p_b_tGrow = dat_now[dat_now$vitalRate == "tMean_grow_C","flwr"], 
    p_b_tWinter = dat_now[dat_now$vitalRate == "tMean_winter_C","flwr"], 
    b_int = dat_now[dat_now$vitalRate == "Intercept","seed"], 
    b_slope = dat_now[dat_now$vitalRate == "log_LL_t","seed"], 
    b_dd = dat_now[dat_now$vitalRate == "N_all_t","seed"], 
    b_tGrow = dat_now[dat_now$vitalRate == "tMean_grow_C","seed"], 
    b_tWinter = dat_now[dat_now$vitalRate == "tMean_winter_C","seed"], 
    c_o_int = dat_now[dat_now$vitalRate == "Intercept","rec"],
    c_o_dd = dat_now[dat_now$vitalRate == "N_all_t","rec"],
    c_o_tGrow = dat_now[dat_now$vitalRate == "tMean_grow_C","rec"],
    c_o_tWinter = dat_now[dat_now$vitalRate == "tMean_winter_C","rec"],
    c_o_sd = 0.7684952
  )
}
param_list_N <- list(
  g_int     = coef(sizeMod_N)[1], # growth 
  g_slope   = coef(sizeMod_N)[2],
  g_dd      = coef(sizeMod_N)[3],
  g_sd      = summary(sizeMod_N)$sigma,
  s_int     = coef(survMod_N)[1], # survival
  s_slope   = coef(survMod_N)[2],
  s_dd      = coef(survMod_N)[3],
  p_b_int   = coef(flwrMod_N)[1], #probability of flowering
  p_b_slope = coef(flwrMod_N)[2],
  p_b_slope_2 = coef(flwrMod_N)[3],
  p_b_dd    = coef(flwrMod_N)[4],
  b_int   = coef(seedMod_N)[1], #seed production
  b_slope = coef(seedMod_N)[2],
  c_o_mu    = coef(recMod_N), #recruit size distribution
  c_o_sd    = summary(recMod_N)$sigma,
  outSB  = outSB_all,
  staySB = staySB_all,
  goSB   = goSB_all, 
  goCont = goCont_all                  
)

# Construct an IPM kernel K using the parameters we obtained from the models

# First define the elements that make up the IPM (vital rates):
# SURVIVAL:
S.fun <- function(z, N_all) {
  mu.surv=param_list_N$s_int + param_list_N$s_slope *z + param_list_N$s_dd * N_all
  return(1/(1 + exp(-(mu.surv))))
}

# GROWTH (we assume a constant variance)
GR.fun <- function(z,zz, N_all){
  growth.mu = param_list_N$g_int + param_list_N$g_slope *z + param_list_N$g_dd * N_all
  return(dnorm(zz, mean = growth.mu, sd = param_list_N$g_sd))
}

## SEEDLING SIZES (same approach as in growth function)
SDS.fun <- function(zz){
  rec_mu <- param_list_N$c_o_mu
  rec_sd <- param_list_N$c_o_sd
  return(dnorm(zz, mean = rec_mu, sd = rec_sd))
}

# PROBABILITY OF FLOWERING 
FL.fun <- function(z, N_all) {
  mu.fl = param_list_N$p_b_int + param_list_N$p_b_slope*z +  param_list_N$p_b_slope_2 * (z^2) + param_list_N$p_b_dd * N_all
  return(1/(1+ exp(-(mu.fl))))
}

# SEED PRODUCTION
SDP.fun <- function(z) {
  mu.fps=exp(param_list_N$b_int + param_list_N$b_slope *z)
  return(mu.fps)
}

# Define the lower and upper integration limit
L <-  1.2 * min(dat_all$log_LL_t, na.rm = TRUE) # minimum size
U <-  1.2 * max(dat_all$log_LL_t, na.rm = TRUE) # maximum size

n <-500 # bins

# These are the parameters for the discrete stages
# I usually only have seed banks (SB), but now I added a seedling stage

outSB <- outSB_all #SB to continuous stage
staySB <- staySB_all # staying in SB
goCont <- goCont_all # seeds become continuous right away (without going to the seed bank) 
goSB <- goSB_all # seeds go to the seedbank
surv.seeds <-  0.9 # survival of seeds

K <- array(0,c(n+1,n+1))

# I recommend you set i = 1, set n low to say 10 

# Setting up the kernels
b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint

h=(U-L)/n # bin width 

# Survival and growth 
S <- diag(S.fun(meshp, N_all = 500)) # Survival # put survival probabilities in the diagonal of the matrix
G <- h * t(outer(meshp,meshp,GR.fun, N_all = 500)) # Growth
# G <- t(outer(meshp,meshp,GR.fun)) # Growth

#Recruits distribution (seeds recruited from the seedbank into the continuous stage)
c_o <- h * matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
# c_o <- matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)

#Probability of flowering
Pb = (FL.fun(meshp, N_all = 500))

#Number of seeds produced according to adult size
b_seed = (SDP.fun(meshp))

FecALL= Pb * b_seed

# update the 'S' matrix by multiplying it by (1-Pb), since this is a monocarpic perennial
S_new <- S * (1-Pb)

# Control for eviction:
# this is equivalent to redistributing evicted sizes evenly among existing size classes 
G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
c_o <- c_o/matrix(as.vector(apply(c_o,2,sum)),nrow=n,ncol=n,byrow=TRUE)

# make the continuous part of the P matrix
Pkernel.cont <- as.matrix(G %*% S_new)
# multiply the continuous kernel by the binwidth (h)
# Pkernel.cont <- h * Pkernel.cont

# seedbank (first column of your K)
Pkernel.seedbank = c(staySB, outSB*c_o[,1]) # seeds survive and go to continuous

# Make the full P kernel
Pkernel <- cbind(Pkernel.seedbank,rbind(rep(0,length(meshp)),Pkernel.cont)) # discrete component

## make the F kernel
Fkernel.cont <-  as.matrix(goCont * ((c_o) %*% diag(FecALL))) # the size of seedlings that go into the seed bank from each continuous size class
# multiply the continuous kernel by the binwidth (h)
#Fkernel.cont <- as.matrix(h * Fkernel.cont)

Fkernel.discr  <- matrix(c(0, goSB * (FecALL)), nrow = 1)
# multiply the cont_to_disc distribution by the binwidth (h)
#Fkernel.discr <- as.matrix(h * Fkernel.discr)
Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))

mat_all_DD <-Pkernel+Fkernel

eigenMat <- base::eigen(mat_all_DD)
# get the lambda
lam_all_DD <-  eigenMat$values[1]

eigenMat$vectors