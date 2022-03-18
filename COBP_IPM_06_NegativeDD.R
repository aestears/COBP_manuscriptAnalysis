#//////////////////////////
# Integral Projection Models for Oenothera coloradensis: Analysis for negative density dependence
# Alice Stears
# 11 March 2022
#/////////////////////////

# load required packages
library(tidyverse)
library(MASS)
library(ipmr)

# load data from script 1
dat_all <- read.csv(file = "/Users/astears/COBP_project/allDat_plus_contSeedlings.csv")

#### DI IPM for each site ####
## write vital rate functions
# Define the lower and upper integration limit
L <-  1.2 * min(dat_all$log_LL_t, na.rm = TRUE) # minimum size
U <-  1.2 * max(dat_all$log_LL_t, na.rm = TRUE) # maximum size

n <-500 # bins

germ.rt <-  0.1629526
viab.rt <- 0.5852778
surv.seeds <-  0.9 

# These are the parameters for the discrete stages
outSB <- germ.rt * surv.seeds 
staySB <- (1-germ.rt) * surv.seeds
goSB <- viab.rt * (1 - germ.rt)
goCont <- viab.rt * germ.rt

## make an empty list to hold the IPM kernels 
site_IPMs_DI <- list()
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
  paramCont <- list(
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
    outSB  = outSB,
    staySB = staySB,
    goSB   = goSB, 
    goCont = goCont                 
  )
  # SURVIVAL:
  S.fun <- function(z) {
    mu.surv=paramCont$s_int + paramCont$s_slope *z
    return(1/(1 + exp(-(mu.surv))))
  }
  # GROWTH (we assume a constant variance)
  GR.fun <- function(z,zz){
    growth.mu = paramCont$g_int + paramCont$g_slope*z
    return(dnorm(zz, mean = growth.mu, sd = paramCont$g_sd))
  }
  ## SEEDLING SIZES (same approach as in growth function)
  SDS.fun <- function(zz){
    rec_mu <- paramCont$c_o_mu
    rec_sd <- paramCont$c_o_sd
    return(dnorm(zz, mean = rec_mu, sd = rec_sd))
  }
  # PROBABILITY OF FLOWERING 
  FL.fun <- function(z) {
    mu.fl = paramCont$p_b_int + paramCont$p_b_slope*z +  paramCont$p_b_slope_2 * (z^2)
    return(1/(1+ exp(-(mu.fl))))
  }
  # SEED PRODUCTION
  SDP.fun <- function(z) {
    mu.fps=exp(paramCont$b_int + paramCont$b_slope *z)
    return(mu.fps)
  }
  
  ## fit the IPM
  K <- array(0,c(n+1,n+1))
  # Setting up the kernels
  b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
  meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
  h=(U-L)/n # bin width 
  # Survival and growth 
  S <- diag(S.fun(meshp)) # Survival # put survival probabilities in the diagonal of the matrix
  G <- h * t(outer(meshp,meshp,GR.fun)) # Growth
  #Recruits distribution (seeds recruited from the seedbank into the continuous stage)
  c_o <- h * matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
  #Probability of flowering
  Pb = (FL.fun(meshp))
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
  # seedbank (first column of your K)
  Pkernel.seedbank = c(staySB, outSB*c_o[,1]) # seeds survive and go to continuous
  # Make the full P kernel
  Pkernel <- cbind(Pkernel.seedbank,rbind(rep(0,length(meshp)),Pkernel.cont)) # discrete component
  ## make the F kernel
  Fkernel.cont <-  as.matrix(goCont * ((c_o) %*% diag(FecALL))) # the size of seedlings that go into the seed bank from each continuous size class
  Fkernel.discr  <- matrix(c(0, goSB * (FecALL)), nrow = 1)
  # multiply the cont_to_disc distribution by the binwidth (h)
  Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))
  
  mat <-Pkernel+Fkernel
  
  eigenMat <- eigen(mat)
  
  site_IPMs_DI[[i]] <- mat
}
names(site_IPMs_DI) <- paste0(unique(dat_all$Site))
lambdas_site_DI <- sapply(site_IPMs_DI, FUN = function(x) eigen(x)$values[1])

## estimate bootstrap confidence intervals
# make an empty list to hold the estimate results (parameters, labmda)
siteDI_bootCI_lambdas <- list()

for (i in 1:length(unique(dat_all$Site))) {
  # get the data just for this site
  dat_now <- dat_all[dat_all$Site == unique(dat_all$Site)[i],]
  # make a vector to hold all of the lambdas from the resampling runs
  siteDI_lambdas <- as.vector(rep(0, length.out = 1000))
  # Now, we refit the vital rate models, and refit the IPM
  for(j in 1:1000) {
    ## sample continuous data
    sample_ind <- seq(1, nrow(dat_now), by = 1)
    
    boot_data_ind   <- sample(sample_ind, size = length(sample_ind), replace = TRUE)
    
    boot_data <- dat_now[boot_data_ind,]
    
    ## fit vital rate models
    ## Survival ($s(z)$)
    survDat_now <- boot_data[boot_data$flowering==0 | is.na(boot_data$flowering),]
    survMod_now <- glm(survives_tplus1 ~ log_LL_t, data = survDat_now, family = binomial)
    ## Growth ($G(z',z)$)
    sizeMod_now <- lm(log_LL_tplus1 ~ log_LL_t, data = boot_data)
    ## Number of seeds produced, according to plant size ($b(z)$)
    seedDat_now <- boot_data[boot_data$flowering==1,]
    # fit poisson glm (for count data)
    seedMod_now <- MASS::glm.nb(Num_seeds ~ log_LL_t , data = seedDat_now)
    ## Flowering probability ($p_b(z)$)
    flwrMod_now <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2), data = boot_data, family = binomial)))
    ## Distribution of recruit size ($c_o(z')$)
    # subset the data
    recD_now <- boot_data[boot_data$seedling == 1,]
    # fit the model
    recMod_now <- lm(log_LL_t ~ 1, data = recD_now)
    
    # update the parameter list with parameters from this iteration
    paramCont <- list(
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
      outSB  = outSB,
      staySB = staySB,
      goSB   = goSB, 
      goCont = goCont                  
    )
    S.fun <- function(z) {
      mu.surv=paramCont$s_int + paramCont$s_slope *z 
      return(1/(1 + exp(-(mu.surv))))
    }
    # GROWTH (we assume a constant variance)
    GR.fun <- function(z,zz){
      growth.mu = paramCont$g_int + paramCont$g_slope *z 
      return(dnorm(zz, mean = growth.mu, sd = paramCont$g_sd))
    }
    ## SEEDLING SIZES (same approach as in growth function)
    SDS.fun <- function(zz){
      rec_mu <- paramCont$c_o_mu
      rec_sd <- paramCont$c_o_sd
      return(dnorm(zz, mean = rec_mu, sd = rec_sd))
    }
    # PROBABILITY OF FLOWERING 
    FL.fun <- function(z) {
      mu.fl = paramCont$p_b_int + paramCont$p_b_slope*z +  paramCont$p_b_slope_2 * (z^2) 
      return(1/(1+ exp(-(mu.fl))))
    }
    # SEED PRODUCTION
    SDP.fun <- function(z) {
      mu.fps=exp(paramCont$b_int + paramCont$b_slope *z)
      return(mu.fps)
    }
    
    ## fit the IPM
    K <- array(0,c(n+1,n+1))
    # Setting up the kernels
    b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
    meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
    h=(U-L)/n # bin width 
    # Survival and growth 
    S <- diag(S.fun(meshp)) # Survival # put survival probabilities in the diagonal of the matrix
    G <- h * t(outer(meshp,meshp,GR.fun)) # Growth
    #Recruits distribution (seeds recruited from the seedbank into the continuous stage)
    c_o <- h * matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
    #Probability of flowering
    Pb = (FL.fun(meshp))
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
    # seedbank (first column of your K)
    Pkernel.seedbank = c(staySB, outSB*c_o[,1]) # seeds survive and go to continuous
    # Make the full P kernel
    Pkernel <- cbind(Pkernel.seedbank,rbind(rep(0,length(meshp)),Pkernel.cont)) # discrete component
    ## make the F kernel
    Fkernel.cont <-  as.matrix(goCont * ((c_o) %*% diag(FecALL))) # the size of seedlings that go into the seed bank from each continuous size class
    Fkernel.discr  <- matrix(c(0, goSB * (FecALL)), nrow = 1)
    # multiply the cont_to_disc distribution by the binwidth (h)
    Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))
    mat <-Pkernel+Fkernel
    
    # save the lambda
    siteDI_lambdas[j] <-base::eigen(mat)$values[1]
    if (j == 1000) {
      siteDI_bootCI_lambdas[[i]] <- siteDD_lambdas
    }
    # save the parameters
    if (j == 1 & i == 1){
      siteDI_bootCI_params <- as.data.frame(c(paramCont, "subPop" =unique(dat_all$Site)[i]))
    } else {
      siteDI_bootCI_params <- rbind(siteDI_bootCI_params, as.data.frame(c(paramCont, "subPop" =unique(dat_all$Site)[i])))
    }
  }
}

#### DD IPM for each site ####

# Define the lower and upper integration limit
L <-  1.2 * min(dat_all$log_LL_t, na.rm = TRUE) # minimum size
U <-  1.2 * max(dat_all$log_LL_t, na.rm = TRUE) # maximum size

n <-500 # bins

# These are the parameters for the discrete stages
# I usually only have seed banks (SB), but now I added a seedling stage
# 
# outSB <- outSB_all #SB to continuous stage
# staySB <- staySB_all # staying in SB
# goCont <- goCont_all # seeds become continuous right away (without going to the seed bank) 
# goSB <- goSB_all # seeds go to the seedbank
# surv.seeds <-  0.9 # survival of seeds

## make an empty list to hold the IPM kernels 
site_IPMs_DD <- list()
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
  seedMod_now <- MASS::glm.nb(Num_seeds ~ log_LL_t + N_Site_t , data = seedDat_now)
  ## Flowering probability ($p_b(z)$)
  flwrMod_now <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) + N_Site_t, data = dat_now, family = binomial)))
  ## Distribution of recruit size ($c_o(z')$)
  # subset the data
  recD_now <- dat_now[dat_now$seedling == 1,]
  # fit the model
  recMod_now <- lm(log_LL_t ~ 1 + N_Site_t, data = recD_now)
  
  ## put in the parameter list (paramCont)
  paramCont <- list(
    g_int     = coef(sizeMod_now)[1], # growth 
    g_slope   = coef(sizeMod_now)[2],
    g_sd      = summary(sizeMod_now)$sigma,
    g_dd      = coef(sizeMod_now)[3],
    s_int     = coef(survMod_now)[1], # survival
    s_slope   = coef(survMod_now)[2],
    s_dd      = coef(survMod_now)[3],
    p_b_int   = coef(flwrMod_now)[1], #probability of flowering
    p_b_slope = coef(flwrMod_now)[2],
    p_b_slope_2 = coef(flwrMod_now)[3],
    p_b_dd    = coef(flwrMod_now)[4],
    b_int   = coef(seedMod_now)[1], #seed production
    b_slope = coef(seedMod_now)[2],
    b_dd   = coef(seedMod_now)[3],
    c_o_int    = coef(recMod_now)[1], #recruit size distribution
    c_o_dd  = coef(recMod_now)[2],
    c_o_sd    = summary(recMod_now)$sigma,
    outSB  = outSB,
    staySB = staySB,
    goSB   = goSB, 
    goCont = goCont                  
  )
  S.fun <- function(z, N_all) {
    mu.surv=paramCont$s_int + paramCont$s_slope *z + paramCont$s_dd * N_all
    return(1/(1 + exp(-(mu.surv))))
  }
  # GROWTH (we assume a constant variance)
  GR.fun <- function(z,zz, N_all){
    growth.mu = paramCont$g_int + paramCont$g_slope *z + paramCont$g_dd * N_all
    return(dnorm(zz, mean = growth.mu, sd = paramCont$g_sd))
  }
  ## SEEDLING SIZES (same approach as in growth function)
  SDS.fun <- function(zz, N_all){
    rec_mu <- paramCont$c_o_int + paramCont$c_o_dd * N_all
    rec_sd <- paramCont$c_o_sd
    return(dnorm(zz, mean = rec_mu, sd = rec_sd))
  }
  # PROBABILITY OF FLOWERING 
  FL.fun <- function(z, N_all) {
    mu.fl = paramCont$p_b_int + paramCont$p_b_slope*z +  paramCont$p_b_slope_2 * (z^2) + paramCont$p_b_dd * N_all
    return(1/(1+ exp(-(mu.fl))))
  }
  # SEED PRODUCTION
  SDP.fun <- function(z, N_all) {
    mu.fps=exp(paramCont$b_int + paramCont$b_slope *z + paramCont$b_dd * N_all)
    return(mu.fps)
  }
  
  ## fit the IPM
  K <- array(0,c(n+1,n+1))
  # Setting up the kernels
  b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
  meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
  h=(U-L)/n # bin width 
  # Survival and growth 
  S <- diag(S.fun(meshp, N_all = 500)) # Survival # put survival probabilities in the diagonal of the matrix
  G <- h * t(outer(meshp,meshp,GR.fun, N_all = 500)) # Growth
  #Recruits distribution (seeds recruited from the seedbank into the continuous stage)
  c_o <- h * matrix(rep(SDS.fun(meshp, N_all = 500),n),n,n,byrow=F)
  #Probability of flowering
  Pb = (FL.fun(meshp, N_all = 500))
  #Number of seeds produced according to adult size
  b_seed = (SDP.fun(meshp, N_all = 500))
  FecALL= Pb * b_seed
  # update the 'S' matrix by multiplying it by (1-Pb), since this is a monocarpic perennial
  S_new <- S * (1-Pb)
  # Control for eviction:
  # this is equivalent to redistributing evicted sizes evenly among existing size classes 
  G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  c_o <- c_o/matrix(as.vector(apply(c_o,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  # make the continuous part of the P matrix
  Pkernel.cont <- as.matrix(G %*% S_new)
  # seedbank (first column of your K)
  Pkernel.seedbank = c(staySB, outSB*c_o[,1]) # seeds survive and go to continuous
  # Make the full P kernel
  Pkernel <- cbind(Pkernel.seedbank,rbind(rep(0,length(meshp)),Pkernel.cont)) # discrete component
  ## make the F kernel
  Fkernel.cont <-  as.matrix(goCont * ((c_o) %*% diag(FecALL))) # the size of seedlings that go into the seed bank from each continuous size class
  Fkernel.discr  <- matrix(c(0, goSB * (FecALL)), nrow = 1)
  # multiply the cont_to_disc distribution by the binwidth (h)
  Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))
  mat <-Pkernel+Fkernel
  
  site_IPMs_DD[[i]] <- mat
}

names(site_IPMs_DD) <- paste0(unique(dat_all$Site))

lambdas_site_DD <- sapply(site_IPMs_DD, FUN = function(x) eigen(x)$values[1])

## estimate bootstrap confidence intervals
## estimate bootstrap confidence intervals for each model parameter for each site
# make an empty list to hold the estimate results (parameters, labmda)
siteDD_bootCI_lambdas <- list()
for (i in 1:length(unique(dat_all$Site))) {
  # get the data just for this site
  dat_now <- dat_all[dat_all$Site == unique(dat_all$Site)[i],]
  # make a vector to hold all of the lambdas from the resampling runs
  siteDD_lambdas <- as.vector(rep(0, length.out = 1000))
  # Now, we refit the vital rate models, and refit the IPM
  for(j in 1:1000) {
    ## sample continuous data
    sample_ind <- seq(1, nrow(dat_now), by = 1)
    
    boot_data_ind   <- sample(sample_ind, size = length(sample_ind), replace = TRUE)
    
    boot_data <- dat_now[boot_data_ind,]
    
    ## fit vital rate models
    ## Survival ($s(z)$)
    survDat_now <- boot_data[boot_data$flowering==0 | is.na(boot_data$flowering),]
    survMod_now <- glm(survives_tplus1 ~ log_LL_t + N_Site_t, data = survDat_now, family = binomial)
    ## Growth ($G(z',z)$)
    sizeMod_now <- lm(log_LL_tplus1 ~ log_LL_t + N_Site_t , data = boot_data)
    ## Number of seeds produced, according to plant size ($b(z)$)
    seedDat_now <- boot_data[boot_data$flowering==1,]
    # fit poisson glm (for count data)
    seedMod_now <- MASS::glm.nb(Num_seeds ~ log_LL_t + N_Site_t, data = seedDat_now)
    ## Flowering probability ($p_b(z)$)
    flwrMod_now <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) + N_Site_t, data = boot_data, family = binomial)))
    ## Distribution of recruit size ($c_o(z')$)
    # subset the data
    recD_now <- boot_data[boot_data$seedling == 1,]
    # fit the model
    recMod_now <- lm(log_LL_t ~ 1 + N_Site_t, data = recD_now)
    
    # update the parameter list with parameters from this iteration
    paramCont <- list(
      g_int     = coef(sizeMod_now)[1], # growth 
      g_slope   = coef(sizeMod_now)[2],
      g_sd      = summary(sizeMod_now)$sigma,
      g_dd      = coef(sizeMod_now)[3],
      s_int     = coef(survMod_now)[1], # survival
      s_slope   = coef(survMod_now)[2],
      s_dd      = coef(survMod_now)[3],
      p_b_int   = coef(flwrMod_now)[1], #probability of flowering
      p_b_slope = coef(flwrMod_now)[2],
      p_b_slope_2 = coef(flwrMod_now)[3],
      p_b_dd    = coef(flwrMod_now)[4],
      b_int   = coef(seedMod_now)[1], #seed production
      b_slope = coef(seedMod_now)[2],
      b_dd   = coef(seedMod_now)[3],
      c_o_int    = coef(recMod_now)[1], #recruit size distribution
      c_o_dd  = coef(recMod_now)[2],
      c_o_sd    = summary(recMod_now)$sigma,
      outSB  = outSB,
      staySB = staySB,
      goSB   = goSB, 
      goCont = goCont                  
    )
    S.fun <- function(z, N_all) {
      mu.surv=paramCont$s_int + paramCont$s_slope *z + paramCont$s_dd * N_all
      return(1/(1 + exp(-(mu.surv))))
    }
    # GROWTH (we assume a constant variance)
    GR.fun <- function(z,zz, N_all){
      growth.mu = paramCont$g_int + paramCont$g_slope *z + paramCont$g_dd * N_all
      return(dnorm(zz, mean = growth.mu, sd = paramCont$g_sd))
    }
    ## SEEDLING SIZES (same approach as in growth function)
    SDS.fun <- function(zz, N_all){
      rec_mu <- paramCont$c_o_int + paramCont$c_o_dd * N_all
      rec_sd <- paramCont$c_o_sd
      return(dnorm(zz, mean = rec_mu, sd = rec_sd))
    }
    # PROBABILITY OF FLOWERING 
    FL.fun <- function(z, N_all) {
      mu.fl = paramCont$p_b_int + paramCont$p_b_slope*z +  paramCont$p_b_slope_2 * (z^2) + paramCont$p_b_dd * N_all
      return(1/(1+ exp(-(mu.fl))))
    }
    # SEED PRODUCTION
    SDP.fun <- function(z, N_all) {
      mu.fps=exp(paramCont$b_int + paramCont$b_slope *z + paramCont$b_dd * N_all)
      return(mu.fps)
    }
    
    ## fit the IPM
    K <- array(0,c(n+1,n+1))
    # Setting up the kernels
    b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
    meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
    h=(U-L)/n # bin width 
    # Survival and growth 
    S <- diag(S.fun(meshp, N_all = 500)) # Survival # put survival probabilities in the diagonal of the matrix
    G <- h * t(outer(meshp,meshp,GR.fun, N_all = 500)) # Growth
    #Recruits distribution (seeds recruited from the seedbank into the continuous stage)
    c_o <- h * matrix(rep(SDS.fun(meshp, N_all = 500),n),n,n,byrow=F)
    #Probability of flowering
    Pb = (FL.fun(meshp, N_all = 500))
    #Number of seeds produced according to adult size
    b_seed = (SDP.fun(meshp, N_all = 500))
    FecALL= Pb * b_seed
    # update the 'S' matrix by multiplying it by (1-Pb), since this is a monocarpic perennial
    S_new <- S * (1-Pb)
    # Control for eviction:
    # this is equivalent to redistributing evicted sizes evenly among existing size classes 
    G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
    c_o <- c_o/matrix(as.vector(apply(c_o,2,sum)),nrow=n,ncol=n,byrow=TRUE)
    # make the continuous part of the P matrix
    Pkernel.cont <- as.matrix(G %*% S_new)
    # seedbank (first column of your K)
    Pkernel.seedbank = c(staySB, outSB*c_o[,1]) # seeds survive and go to continuous
    # Make the full P kernel
    Pkernel <- cbind(Pkernel.seedbank,rbind(rep(0,length(meshp)),Pkernel.cont)) # discrete component
    ## make the F kernel
    Fkernel.cont <-  as.matrix(goCont * ((c_o) %*% diag(FecALL))) # the size of seedlings that go into the seed bank from each continuous size class
    Fkernel.discr  <- matrix(c(0, goSB * (FecALL)), nrow = 1)
    # multiply the cont_to_disc distribution by the binwidth (h)
    Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))
    mat <-Pkernel+Fkernel
    
    # save the lambda
    siteDD_lambdas[j] <-base::eigen(mat)$values[1]
    if (j == 1000) {
      siteDD_bootCI_lambdas[[i]] <- siteDD_lambdas
    }
    # save the parameters
    if (j == 1 & i == 1){
      siteDD_bootCI_params <- as.data.frame(c(paramCont, "subPop" =unique(dat_all$Site)[i]))
    } else {
      siteDD_bootCI_params <- rbind(siteDD_bootCI_params, as.data.frame(c(paramCont, "subPop" =unique(dat_all$Site)[i])))
    }
  }
}
siteDD_bootCI_lamdas[[6]] <- tempDD_6

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
siteDI_lams <- sapply(site_IPMs_DI, function(x) as.numeric(eigen(x)$values[1]))
siteDI_lams <- data.frame("Site" = names(siteDI_lams), "log_lambda" = log(siteDI_lams), "type" = "DI")
siteDD_lams <- sapply(site_IPMs_DD, function(x) as.numeric(eigen(x)$values[1]))
siteDD_lams <- data.frame("Site" = names(siteDD_lams), "log_lambda" = log(siteDD_lams), "type" = "DD")

site_lams <- rbind(siteDI_lams, siteDD_lams)

meanN_subPop <- N_site %>% 
  group_by(Site) %>% 
  summarize("meanN_t" = mean(N_Site_t))

lam_N <- left_join(site_lams, meanN_subPop)

ggplot(data = lam_N) + 
  geom_point(aes(x = meanN_t, y = log_lambda, col = type)) + 
  geom_smooth(aes(x = meanN_t, y = log_lambda, col = type), se = FALSE, method = "lm") +
  theme_classic(
    
  )
