#//////////////////////////
# Integral Projection Models for Oenothera coloradensis: Analysis for Vital Rate Buffering
# Alice Stears
# 11 March 2022
#/////////////////////////

library(tidyverse)

# load data from script 1
dat_all <- read.csv(file = "../Processed_Data/allDat_plus_contSeedlings.csv")

#### fit vital rate models for 2018-2019, each subpopulation ####
subPop_first_VRs <- list()
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
  seedMod_now <- MASS::glm.nb(Num_seeds ~ log_LL_t , data = seedDat_now)
  ## Flowering probability ($p_b(z)$)
  flwrMod_now <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = dat_now, family = binomial)))
  ## Distribution of recruit size ($c_o(z')$)
  recD_now <- dat_all[dat_all$seedling == 1 & dat_all$Year == 2019,]
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
  subPop_first_VRs[[i]] <- paramCont
}
names(subPop_first_VRs) <- paste0(unique(dat_all$Site))

#### fit vital rate models for 2019-2020, each subpopulation ####
subPop_second_VRs <- list()
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
  seedMod_now <- MASS::glm.nb(Num_seeds ~ log_LL_t , data = seedDat_now)
  ## Flowering probability ($p_b(z)$)
  flwrMod_now <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = dat_now, family = binomial)))
  ## Distribution of recruit size ($c_o(z')$)
  recD_now <- dat_all[dat_all$seedling == 1 & dat_all$Year == 2020,]
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
  subPop_second_VRs[[i]] <- paramCont
}
names(subPop_second_VRs) <- paste0(unique(dat_all$Site))

#### combine the vital rate params for comparison ####
temp <- as.data.frame(sapply(subPop_first_VRs, function(x) data.frame(x)))
names(temp)  <- paste0(names(subPop_second_VRs), "_first")
temp2 <- as.data.frame(sapply(subPop_second_VRs, function(x) data.frame(x)))
names(temp2)  <- paste0(names(subPop_second_VRs),"_second")
vrDF <- cbind(temp, temp2)
vrDF <- as.data.frame(apply(vrDF, 2, unlist))

vrDF$mean <- apply(X = vrDF, MARGIN = 1, FUN = mean)
vrDF$sd <- apply(X = vrDF[,1:12], MARGIN = 1, FUN = sd)
vrDF$CV <- vrDF$sd/vrDF$mean * 100
# calculate the "quartile coefficient of dispersion" https://en.wikipedia.org/wiki/Quartile_coefficient_of_dispersion (Q3-Q1)/(Q3+Q1)
quantiles <- as.data.frame(apply(vrDF[,1:12], MARGIN = 1, FUN = function(x) quantile(x, probs = c(.25, .75))))
quartile_coef_disp <- apply(quantiles, MARGIN = 2, FUN = function(x) (x[2] - x[1])/(x[2] + x[1]))
vrDF$quarDisp <- quartile_coef_disp

# simulate variance of sb parameters as normal dists w/ mean of .5 and sd of .175??
#vrDF[13:16,"sd"]<- .175
#vrDF[14:16,"mean"] <- .5
## drop seedbank params, since we don't have variability data for them
vrDF <- vrDF[1:12,] 

#### read in elasticity data ####
contParamElas <- readRDS("./intermediate_analysis_Data/allSiteAllYears_noDDnoEnv/continuousParamElasticity.RDS")
discParamElas <- readRDS("./intermediate_analysis_Data/allSiteAllYears_noDDnoEnv/discreteParamElasticity.RDS")
# reorder correctly
contParamElas <- contParamElas[match(c("growth_(Intercept)", "growth_log_LL_t", "growth_stndDev", "survival_(Intercept)", "survival_log_LL_t","flowering_(Intercept)", "flowering_log_LL_t", "flowering_I(log_LL_t^2)","seedProduction_(Intercept)", "seedProduction_log_LL_t", "recruitDist_(Intercept)", "recruitDist_stndDev"), contParamElas$param_name),] 
## put elasticity in the vrDF
vrDF$Elas <- contParamElas$elas_mean
vrDF$Sens <- contParamElas$sens_mean

## add data for sb params? 
# vrDF[13,"mean"] <- mean(c(13.0, 12, 8.3, 7.0, 5.3)/(45 * seed_per_cap))
# vrDF[13, "sd"] <- sd(c(13.0, 12, 8.3, 7.0, 5.3)/(45 * seed_per_cap))
# 
# vrDF[14,"mean"] <- mean(c(81,61,54)/100 * 0.89)
# vrDF[14, "sd"] <- sd(c(81,61,54)/100 * 0.89)
# 
# vrDF[15, "mean"] <- 0.5
# vrDF[15, "sd"] <- 0.175
# simulate the sb params using a uniform distribution 
rownames(vrDF[13:15,]) <- c("germ.rt", "surv.rt", "viab.rt")
vrDF[13:15,"quarDisp"] <- c(.75-.25)/(.75+.25)

vrDF[13:15,"Elas"] <- discParamElas$elas_mean
vrDF[13:15, "Sens"] <- discParamElas$sens_mean

# vrDF$CV <- vrDF$sd/vrDF$mean * 100
#### calculate the correlation between CV and elasticity ####
cor.test((vrDF$quarDisp[c(1:3,5:12)]), abs(vrDF$Elas[c(1:3,5:12)]), method = "pearson")
ggplot(data = vrDF) +
  geom_smooth(aes(x = (quarDisp), y = abs(Elas)),data = vrDF, method = "lm", lty = 2, col = "red", alpha = .3, se = FALSE) + 
  geom_smooth(aes(x = (quarDisp), y = abs(Elas)),data = vrDF[c(1:3,5:12),], method = "lm", lty = 2, col = "blue", alpha = .3, se = FALSE) +
  geom_point(aes(x = (quarDisp), y = abs(Elas))) +
  xlab(c("Quartile coefficient of dispersion")) +
  ylab(c("|Vital Rate Param. Elasticity|")) + 
  theme_classic() 

# plot(x = (vrDF$quarDisp), y = vrDF$Elas, 
#      xlab = "Vital Rate Param. Coefficient of Variation (CV)", 
#      ylab = "Vital Rate Param. Elasticity", 
#      pch = 16, 
#      col = "grey20")
# abline(lm(vrDF$Elas ~ vrDF$quarDisp), col = "red", lty = 2)
# abline(lm(Elas ~ quarDisp, data = vrDF[c(1:3,5:12),]), col = "blue", lty = 2)
# text(x =650, y = -1.5, labels = "r = 0.116 \n P-value = 0.68 \n(t = 0.42, df = 13)")

#### fit IPMs to get lambdas ####
outSB <- outSB_all #SB to continuous stage
staySB <- staySB_all # staying in SB
goCont <- goCont_all # seeds become continuous right away (without going to the seed bank) 
goSB <- goSB_all # seeds go to the seedbank
surv.seeds <-  0.9 # survival of seeds
# SURVIVAL:
S.fun <- function(z, paramCont) {
  mu.surv=paramCont$s_int + paramCont$s_slope *z
  return(1/(1 + exp(-(mu.surv))))
}
# GROWTH (we assume a constant variance)
GR.fun <- function(z,zz, paramCont){
  growth.mu = paramCont$g_int + paramCont$g_slope*z
  return(dnorm(zz, mean = growth.mu, sd = paramCont$g_sd))
}
## SEEDLING SIZES (same approach as in growth function)
SDS.fun <- function(zz, paramCont){
  rec_mu <- paramCont$c_o_mu
  rec_sd <- paramCont$c_o_sd
  return(dnorm(zz, mean = rec_mu, sd = rec_sd))
}
# PROBABILITY OF FLOWERING 
FL.fun <- function(z, paramCont) {
  mu.fl = paramCont$p_b_int + paramCont$p_b_slope*z +  paramCont$p_b_slope_2 * (z^2)
  return(1/(1+ exp(-(mu.fl))))
}
# SEED PRODUCTION
SDP.fun <- function(z, paramCont) {
  mu.fps=exp(paramCont$b_int + paramCont$b_slope *z)
  return(mu.fps)
}
# for first half of data
subPop_first_mats <- list()
for (i in 1:length(subPop_first_VRs)) {
  paramsNow <- subPop_first_VRs[[i]]
  ## fit the IPM
  K <- array(0,c(n+1,n+1))
  # Setting up the kernels
  b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
  meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
  h=(U-L)/n # bin width 
  # Survival and growth 
  S <- diag(S.fun(meshp, paramCont = paramsNow)) # Survival # put survival probabilities in the diagonal of the matrix
  G <- h * t(outer(meshp,meshp,GR.fun, paramCont = paramsNow)) # Growth
  #Recruits distribution (seeds recruited from the seedbank into the continuous stage)
  c_o <- h * matrix(rep(SDS.fun(meshp, paramCont = paramsNow),n),n,n,byrow=F)
  #Probability of flowering
  Pb = (FL.fun(meshp, paramCont = paramsNow))
  #Number of seeds produced according to adult size
  b_seed = (SDP.fun(meshp, paramCont = paramsNow))
  FecALL= Pb * b_seed
  # update the 'S' matrix by multiplying it by (1-Pb), since this is a monocarpic perennial
  S_new <- S * (1-Pb)
  # Control for eviction:
  G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  c_o <- c_o/matrix(as.vector(apply(c_o,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  # make the continuous part of the P matrix
  Pkernel.cont <- as.matrix(G %*% S_new)
  # seedbank (first column of your K)
  Pkernel.seedbank = c(staySB, outSB*c_o[,1]) # seeds survive and go to continuous
  # Make the full P kernel
  Pkernel <- cbind(Pkernel.seedbank,rbind(rep(0,length(meshp)),Pkernel.cont)) # discrete component
  ## make the F kernel
  Fkernel.cont <-  as.matrix(goCont * ((c_o) %*% diag(FecALL))) 
  Fkernel.discr  <- matrix(c(0, goSB * (FecALL)), nrow = 1)
  # multiply the cont_to_disc distribution by the binwidth (h)
  Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))
  
  mat <-Pkernel+Fkernel
  eigenMat <- eigen(mat)
  
  subPop_first_mats[[i]] <- mat
}
names(subPop_first_mats) <- names(subPop_first_VRs)

# for second half of data
subPop_second_mats <- list()
for (i in 1:length(subPop_second_VRs)) {
  paramsNow <- subPop_second_VRs[[i]]
  ## fit the IPM
  K <- array(0,c(n+1,n+1))
  # Setting up the kernels
  b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
  meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
  h=(U-L)/n # bin width 
  # Survival and growth 
  S <- diag(S.fun(meshp, paramCont = paramsNow)) # Survival # put survival probabilities in the diagonal of the matrix
  G <- h * t(outer(meshp,meshp,GR.fun, paramCont = paramsNow)) # Growth
  #Recruits distribution (seeds recruited from the seedbank into the continuous stage)
  c_o <- h * matrix(rep(SDS.fun(meshp, paramCont = paramsNow),n),n,n,byrow=F)
  #Probability of flowering
  Pb = (FL.fun(meshp, paramCont = paramsNow))
  #Number of seeds produced according to adult size
  b_seed = (SDP.fun(meshp, paramCont = paramsNow))
  FecALL= Pb * b_seed
  # update the 'S' matrix by multiplying it by (1-Pb), since this is a monocarpic perennial
  S_new <- S * (1-Pb)
  # Control for eviction:
  G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  c_o <- c_o/matrix(as.vector(apply(c_o,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  # make the continuous part of the P matrix
  Pkernel.cont <- as.matrix(G %*% S_new)
  # seedbank (first column of your K)
  Pkernel.seedbank = c(staySB, outSB*c_o[,1]) # seeds survive and go to continuous
  # Make the full P kernel
  Pkernel <- cbind(Pkernel.seedbank,rbind(rep(0,length(meshp)),Pkernel.cont)) # discrete component
  ## make the F kernel
  Fkernel.cont <-  as.matrix(goCont * ((c_o) %*% diag(FecALL))) 
  Fkernel.discr  <- matrix(c(0, goSB * (FecALL)), nrow = 1)
  # multiply the cont_to_disc distribution by the binwidth (h)
  Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))
  
  mat <-Pkernel+Fkernel
  eigenMat <- eigen(mat)
  
  subPop_second_mats[[i]] <- mat
}
names(subPop_second_mats) <- names(subPop_second_VRs)

# calculate lambdas
log(sapply(subPop_first_mats, function(x) as.numeric(eigen(x)$values[1])))

log(sapply(subPop_second_mats, function(x) as.numeric(eigen(x)$values[1])))
