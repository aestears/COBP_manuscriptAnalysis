#//////////////////////////
# Integral Projection Models for Oenothera coloradensis: Analysis for negative density dependence
# Alice Stears
# 11 March 2022
#/////////////////////////

# load required packages
library(tidyverse)
library(MASS)

# load data from script 2
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
  seedMod_now <- MASS::glm.nb(Num_seeds ~ log_LL_t , data = seedDat_now)
  ## Flowering probability ($p_b(z)$)
  flwrMod_now <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) + N_Site_t, data = dat_now, family = binomial)))
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
    c_o_mu    = coef(recMod_now), #recruit size distribution
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
  SDS.fun <- function(zz){
    rec_mu <- paramCont$c_o_mu
    rec_sd <- paramCont$c_o_sd
    return(dnorm(zz, mean = rec_mu, sd = rec_sd))
  }
  # PROBABILITY OF FLOWERING 
  FL.fun <- function(z, N_all) {
    mu.fl = paramCont$p_b_int + paramCont$p_b_slope*z +  paramCont$p_b_slope_2 * (z^2) + paramCont$p_b_dd * N_all
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
  S <- diag(S.fun(meshp, N_all = 500)) # Survival # put survival probabilities in the diagonal of the matrix
  G <- h * t(outer(meshp,meshp,GR.fun, N_all = 500)) # Growth
  
  #Recruits distribution (seeds recruited from the seedbank into the continuous stage)
  c_o <- h * matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
  
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

for (i in 6:6){#length(unique(dat_all$Site))) {
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
    seedMod_now <- MASS::glm.nb(Num_seeds ~ log_LL_t , data = seedDat_now)
    ## Flowering probability ($p_b(z)$)
    flwrMod_now <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) + N_Site_t, data = boot_data, family = binomial)))
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
      c_o_mu    = coef(recMod_now), #recruit size distribution
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
    SDS.fun <- function(zz){
      rec_mu <- paramCont$c_o_mu
      rec_sd <- paramCont$c_o_sd
      return(dnorm(zz, mean = rec_mu, sd = rec_sd))
    }
    # PROBABILITY OF FLOWERING 
    FL.fun <- function(z, N_all) {
      mu.fl = paramCont$p_b_int + paramCont$p_b_slope*z +  paramCont$p_b_slope_2 * (z^2) + paramCont$p_b_dd * N_all
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
    S <- diag(S.fun(meshp, N_all = 500)) # Survival # put survival probabilities in the diagonal of the matrix
    G <- h * t(outer(meshp,meshp,GR.fun, N_all = 500)) # Growth
    #Recruits distribution (seeds recruited from the seedbank into the continuous stage)
    c_o <- h * matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
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

## save the data to file
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
