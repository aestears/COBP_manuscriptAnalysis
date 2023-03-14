#//////////////////////////
# Integral Projection Models for Oenothera coloradensis: Code for IPMs CC - NN
# Alice Stears
# 08 December 2021
#/////////////////////////

#### load vital rate models from previous script ####
source("./analysis_scripts/COBP_IPM_02_VitalRateModels.R")

#### IPM CC-HH ####
### DI IPM for each site--first half of data, discrete ###
# Define the lower and upper integration limit
L <-  1.2 * min(dat_all$log_LL_t, na.rm = TRUE) # minimum size
U <-  1.2 * max(dat_all$log_LL_t, na.rm = TRUE) # maximum size

n <-200 # bins

# These are the parameters for the discrete stages
outSB <- outSB_all #SB to continuous stage
staySB <- staySB_all # staying in SB
goCont <- goCont_all # seeds become continuous right away (without going to the seed bank) 
goSB <- goSB_all # seeds go to the seedbank
surv.seeds <-  0.9 # survival of seeds

## make an empty list to hold the IPM kernels 
IPMs_CC_HH <- list()
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
    outSB  = outSB_all,
    staySB = staySB_all,
    goSB   = goSB_all, 
    goCont = goCont_all                  
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
  # G <- t(outer(meshp,meshp,GR.fun)) # Growth
  #Recruits distribution (seeds recruited from the seedbank into the continuous stage)
  c_o <- h * matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
  # c_o <- matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
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
  
  IPMs_CC_HH[[i]] <- mat
}
names(IPMs_CC_HH) <- paste0(unique(dat_all$Site))
lambdas_IPMs_CC_HH <- sapply(IPMs_CC_HH, FUN = function(x) eigen(x)$values[1])

## estimate 95% Bootstrap CIs 
IPMs_CC_HH_bootCI_lambdas <- list()

for (i in 1:length(unique(dat_all$Site))) {
  # get the data just for this site
  dat_now <- dat_all[dat_all$Site == unique(dat_all$Site)[i]
                     & dat_all$Year == 2018# data for this year
                     ,]
  # make a vector to hold all of the lambdas from the resampling runs
  siteDI_lambdas <- as.vector(rep(0, length.out = 1000))
  # Now, we refit the vital rate models, and refit the IPM
  for(j in 1:100) {
    ## sample continuous data
    sample_ind <- seq(1, nrow(dat_now), by = 1)
    
    boot_data_ind   <- sample(sample_ind, size = length(sample_ind), replace = TRUE)
    
    boot_data <- dat_now[boot_data_ind,]
    
    ## fit vital rate models
    ## Survival ($s(z)$)
    survDat_now <- boot_data[boot_data$flowering==0 | is.na(boot_data$flowering),]
    survMod_now <- glm(survives_tplus1 ~ log_LL_t , data = survDat_now, family = binomial)
    ## Growth ($G(z',z)$)
    sizeMod_now <- lm(log_LL_tplus1 ~ log_LL_t , data = boot_data)
    ## Number of seeds produced, according to plant size ($b(z)$)
    seedDat_now <- boot_data[boot_data$flowering==1,]
    # fit poisson glm (for count data)
    seedMod_now <- MASS::glm.nb(Num_seeds ~ log_LL_t , data = seedDat_now)
    ## Flowering probability ($p_b(z)$)
    flwrMod_now <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = boot_data, family = binomial)))
    ## Distribution of recruit size ($c_o(z')$)
    # subset the data
    recD_now <- dat_all[dat_all$seedling == 1 & dat_all$Year == 2019,]
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
      outSB  = outSB_all,
      staySB = staySB_all,
      goSB   = goSB_all, 
      goCont = goCont_all                  
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
    # G <- t(outer(meshp,meshp,GR.fun)) # Growth
    #Recruits distribution (seeds recruited from the seedbank into the continuous stage)
    c_o <- h * matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
    # c_o <- matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
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
    if (j == 100) {
      IPMs_CC_HH_bootCI_lambdas[[i]] <- siteDI_lambdas
    }
    # save the parameters
    if (j == 1 & i == 1){
      IPMs_CC_HH_bootCI_params <- as.data.frame(c(paramCont, "subPop" =unique(dat_all$Site)[i]))
    } else {
      IPMs_CC_HH_bootCI_params <- rbind(IPMs_CC_HH_bootCI_params, as.data.frame(c(paramCont, "subPop" =unique(dat_all$Site)[i])))
    }
  }
}
IPMs_CC_HH_bootCI <- lapply(IPMs_CC_HH_bootCI_lambdas,
                          FUN = function(x)
                            log(c((mean(x[1:100]) - 1.96*sd(x[1:100])/sqrt(100)),
                                  (mean(x[1:100]) + 1.96*sd(x[1:100])/sqrt(100))))
)
names(IPMs_CC_HH_bootCI) <- unique(dat_all$Site)
means <- lapply(IPMs_CC_HH_bootCI_lambdas, FUN = function(x) log(mean(x[1:100])))
names(means) <- unique(dat_all$Site)

#### IPMs II-NN #### 
### DI IPM for each site--second half of data, discrete ###
# Define the lower and upper integration limit
L <-  1.2 * min(dat_all$log_LL_t, na.rm = TRUE) # minimum size
U <-  1.2 * max(dat_all$log_LL_t, na.rm = TRUE) # maximum size

n <-200 # bins

# These are the parameters for the discrete stages
outSB <- outSB_all #SB to continuous stage
staySB <- staySB_all # staying in SB
goCont <- goCont_all # seeds become continuous right away (without going to the seed bank) 
goSB <- goSB_all # seeds go to the seedbank
surv.seeds <-  0.9 # survival of seeds

## make an empty list to hold the IPM kernels 
IPMs_II_NN <- list()
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
    outSB  = outSB_all,
    staySB = staySB_all,
    goSB   = goSB_all, 
    goCont = goCont_all                  
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
  # G <- t(outer(meshp,meshp,GR.fun)) # Growth
  #Recruits distribution (seeds recruited from the seedbank into the continuous stage)
  c_o <- h * matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
  # c_o <- matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
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
  
  IPMs_II_NN[[i]] <- mat
}
names(IPMs_II_NN) <- paste0(unique(dat_all$Site))
lambdas_IPMs_II_NN <- sapply(IPMs_II_NN, FUN = function(x) eigen(x)$values[1])

## estimate 95% Bootstrap CIs 
IPMs_II_NN_bootCI_lambdas <- list()

for (i in c(1:4,6)) {
  # get the data just for this site
  dat_now <- dat_all[dat_all$Site == unique(dat_all$Site)[i]
                     & dat_all$Year == 2019# data for this year
                     ,]
  # make a vector to hold all of the lambdas from the resampling runs
  siteDI_lambdas <- as.vector(rep(0, length.out = 100))
  # Now, we refit the vital rate models, and refit the IPM
  for(j in 1:100) {
    ## sample continuous data
    sample_ind <- seq(1, nrow(dat_now), by = 1)
    
    boot_data_ind   <- sample(sample_ind, size = length(sample_ind), replace = TRUE)
    
    boot_data <- dat_now[boot_data_ind,]
    
    ## fit vital rate models
    ## Survival ($s(z)$)
    survDat_now <- boot_data[boot_data$flowering==0 | is.na(boot_data$flowering),]
    survMod_now <- glm(survives_tplus1 ~ log_LL_t , data = survDat_now, family = binomial)
    ## Growth ($G(z',z)$)
    sizeMod_now <- lm(log_LL_tplus1 ~ log_LL_t , data = boot_data)
    ## Number of seeds produced, according to plant size ($b(z)$)
    seedDat_now <- boot_data[boot_data$flowering==1,]
    # fit poisson glm (for count data)
    seedMod_now <- MASS::glm.nb(Num_seeds ~ log_LL_t , data = seedDat_now)
    ## Flowering probability ($p_b(z)$)
    flwrMod_now <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = boot_data, family = binomial)))
    ## Distribution of recruit size ($c_o(z')$)
    # subset the data
    recD_now <- dat_all[dat_all$seedling == 1 & dat_all$Year == 2019,]
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
      outSB  = outSB_all,
      staySB = staySB_all,
      goSB   = goSB_all, 
      goCont = goCont_all                  
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
    # G <- t(outer(meshp,meshp,GR.fun)) # Growth
    #Recruits distribution (seeds recruited from the seedbank into the continuous stage)
    c_o <- h * matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
    # c_o <- matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
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
    if (j == 100) {
      IPMs_II_NN_bootCI_lambdas[[i]] <- siteDI_lambdas
    }
    # save the parameters
    if (j == 1 & i == 1){
      IPMs_II_NN_bootCI_params <- as.data.frame(c(paramCont, "subPop" =unique(dat_all$Site)[i]))
    } else {
      IPMs_II_NN_bootCI_params <- rbind(IPMs_II_NN_bootCI_params, as.data.frame(c(paramCont, "subPop" =unique(dat_all$Site)[i])))
    }
  }
}
site_secondDI_CI <- lapply(IPMs_II_NN_bootCI_lambdas,
                           FUN = function(x)
                             log(c((mean(x[1:100]) - 1.96*sd(x[1:100])/sqrt(100)),
                                   (mean(x[1:100]) + 1.96*sd(x[1:100])/sqrt(100))))
)
names(site_secondDI_CI) <- unique(dat_all$Site)
lapply(IPMs_II_NN_bootCI_lambdas, 
       FUN = function(x) log(mean(x)))


#### save the data to file ####
fileLoc <- "./intermediate_analysis_Data/site_level_IPMs_eachYear/"
## site-level DI IPM matrices
saveRDS(IPMs_CC_HH, file = paste0(fileLoc,"/IPMs_CC_HH.RDS"))
## site-level DI bootstrap CI data
saveRDS(IPMs_CC_HH_bootCI_lambdas, file = paste0(fileLoc,"/IPMs_CC_HH_bootCI_lambdas.RDS"))
saveRDS(IPMs_CC_HH_bootCI_params, file = paste0(fileLoc,"/IPMs_CC_HH_bootCI_params.RDS"))

## site-level DD IPM matrices
saveRDS(IPMs_II_NN, file = paste0(fileLoc,"/IPMs_II_NN.RDS"))
## site-level DD bootstrap CI data
saveRDS(IPMs_II_NN_bootCI_lambdas, file = paste0(fileLoc, "/IPMs_II_NN_bootCI_lambdas.RDS"))
saveRDS(IPMs_II_NN_bootCI_params, file = paste0(fileLoc, "/IPMs_II_NN_bootCI_params.RDS"))