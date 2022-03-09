#//////////////////////////
# Integral Projection Models for Oenothera coloradensis: simple model using code from Ellner and Rees
# Alice Stears
# 08 December 2021
#/////////////////////////

#### load vital rate models from previous script ####
source("./analysis_scripts/COBP_IPM_02_VitalRateModels.R")

#### hand-calculated model for continuous-ized seedlings and discrete seedbank ####
# Density Independent 

### Define the demographic functions and parameters ###
## (Code below is modified from Maria Paniw)

# Empty list to save model coefficients 
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

# First define the elements that make up the IPM (vital rates):

# SURVIVAL:
S.fun <- function(z) {
  mu.surv=paramCont[["survival"]]["(Intercept)",] + paramCont[["survival"]]["log_LL_t",]*z
  return(1/(1 + exp(-(mu.surv))))
}

# GROWTH (we assume a constant variance)
GR.fun <- function(z,zz){
  growth.mu = paramCont[["growth"]]["(Intercept)",1] + paramCont[["growth"]]["log_LL_t",1]*z
  return(dnorm(zz, mean = growth.mu, sd = paramCont[["growth"]][1,2]))
}

## SEEDLING SIZES (same approach as in growth function)
SDS.fun <- function(zz){
  rec_mu <- paramCont[["recruitDist"]][1]
  rec_sd <- paramCont[["recruitDist"]][2]
  return(dnorm(zz, mean = rec_mu, sd = rec_sd))
}

# PROBABILITY OF FLOWERING 
FL.fun <- function(z) {
  mu.fl = paramCont[["flowering"]][1,] + paramCont[["flowering"]][2,]*z +  paramCont[["flowering"]][3,]* (z^2)
  return(1/(1+ exp(-(mu.fl))))
}

# SEED PRODUCTION
SDP.fun <- function(z) {
  mu.fps=exp(paramCont[["seedProduction"]][1,1] + paramCont[["seedProduction"]][2,1]*z)
  return(mu.fps)
}

# Second, put together the kernels - four kernels for four years:
years <- 1 # for right now, treating everything as one year

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

mat_all_DI <-Pkernel+Fkernel

eigenMat <- eigen(mat_all_DI)
# get the lambda
eigenMat$values[1]

eigenMat$vectors


## calculate the uncertainty using bootstrapping
# make a vector to hold all of the lambdas from the resampling runs
allDI_lambdas <- numeric(200L)
# make a list to hold the model parameters
paramCont_now <- list()
allDI_params <- list()
# set seedbank parameters, which don't change
outSB <- outSB_all #SB to continuous stage
staySB <- staySB_all # staying in SB
goCont <- goCont_all # seeds become continuous right away (without going to the seed bank) 
goSB <- goSB_all # seeds go to the seedbank
surv.seeds <-  0.9 # survival of seeds
# Now, we refit the vital rate models, and refit the IPM
for(i in 1:200) {
  ## sample continuous data
  sample_ind <- seq(1, nrow(dat_all), by = 1)
  
  boot_data_ind   <- sample(sample_ind, size = 6000, replace = TRUE)
  
  boot_data <- dat_all[boot_data_ind,]
  
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
  recD_now <- boot_data[boot_data$seedling == 1,]
  # fit the model
  recMod_now <- lm(log_LL_t ~ 1, data = recD_now)
 
  # update the parameter list with parameters from this iteration
  # survival model is called 'survMod_all'
  paramCont_now[[1]]=as.matrix(coef(survMod_now)) # save coefficients 
  # growth model is called 'sizeMod_all'
  paramCont_now[[2]]=cbind(as.matrix(coef(sizeMod_now)),sd(residuals(sizeMod_now))) # the third column is for the standard deviation of growth 
  # seedling size distribution is a uniform distribution (of exp(size_2)) with a min of 0.1 and a max 0f 3
  paramCont_now[[3]]= cbind(as.matrix(coef(recMod_now)), sd(residuals(recMod_now)))
  # model for probability of flowering is flwrMod_all
  paramCont_now[[4]]=as.matrix(coef(flwrMod_now))
  # model for seed production per plant (if reproductive) is seedMod_all
  paramCont_now[[5]]=as.matrix(coef(seedMod_now))
  
  # name the paramCont list to keep track of coefficients
  names(paramCont) <- c("survival", "growth", "recruitDist", "flowering", "seedProduction")
  
  # Insert the new vital rate models into the proto_ipm, then rebuild the IPM.
  K_now <- array(0,c(n+1,n+1))
  
  # I recommend you set i = 1, set n low to say 10 
  
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
  
  mat <-Pkernel+Fkernel
  
  eigenMat <- eigen(mat)
  # get the lambda and store it
  allDI_lambdas[i] <- eigenMat$values[1]
  
  allDI_params[[i]] <- paramCont_now
}

### plot the values
# get the parameters out
for (i in 1:200) {
  if (i == 1) {
    param_fig_DI <- data.frame("param" = "s_int", "value" = allDI_params[[i]][[1]][1,])
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "s_slope", "value" = allDI_params[[i]][[1]][2,]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "g_int", "value" = allDI_params[[i]][[2]][1,1]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "g_slope", "value" = allDI_params[[i]][[2]][2,1]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "g_sd", "value" = allDI_params[[i]][[2]][1,2]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "c_o_mu", "value" = allDI_params[[i]][[3]][1,1]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "c_o_sd", "value" = allDI_params[[i]][[3]][1,2]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "p_b_int", "value" = allDI_params[[i]][[4]][1,]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "p_b_slope", "value" = allDI_params[[i]][[4]][2,]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "p_b_slope_2", "value" = allDI_params[[i]][[4]][3,]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "b_int", "value" = allDI_params[[i]][[5]][1,]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "b_slope", "value" = allDI_params[[i]][[5]][2,]))
  } else {
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "s_int", "value" = allDI_params[[i]][[1]][1,]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "s_slope", "value" = allDI_params[[i]][[1]][2,]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "g_int", "value" = allDI_params[[i]][[2]][1,1]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "g_slope", "value" = allDI_params[[i]][[2]][2,1]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "g_sd", "value" = allDI_params[[i]][[2]][1,2]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "c_o_mu", "value" = allDI_params[[i]][[3]][1,1]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "c_o_sd", "value" = allDI_params[[i]][[3]][1,2]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "p_b_int", "value" = allDI_params[[i]][[4]][1,]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "p_b_slope", "value" = allDI_params[[i]][[4]][2,]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "p_b_slope_2", "value" = allDI_params[[i]][[4]][3,]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "b_int", "value" = allDI_params[[i]][[5]][1,]))
    param_fig_DI <- rbind(param_fig_DI, data.frame("param" = "b_slope", "value" = allDI_params[[i]][[5]][2,]))
  }
}

# calculate the mean values
test <- param_fig_DI %>% 
  group_by(param) %>% 
  summarise("mean_value" = mean(value))
param_fig <- left_join(param_fig_DI, test)
# plot the distribution of values
ggplot(data = param_fig) + 
  geom_density(mapping = aes(value)) +
  geom_vline(aes(xintercept = mean_value, col = param)) + 
  geom_vline(aes(xintercept = 0), col = "grey", lty = 2) +
  facet_wrap(~ param, scales = "free") +
  scale_color_discrete(guide = "none") + 
  theme_classic()

# compare lambdas
plot(density(as.numeric(allDI_lambdas)))
abline(v = eigenMat$values[1], col = 'red', lwd = 2, lty = 2)
abline(v = mean(allDI_lambdas), col = "blue", lwd = 2, lty = 2)

#### DI, all-site data for first half of data ####
### Define the demographic functions and parameters ###
# Empty list to save model coefficients 
paramCont_first=list(NULL)
# survival model is called 'survMod_all'
paramCont_first[[1]]=as.matrix(coef(survMod_first)) # save coefficients 
# growth model is called 'sizeMod_all'
paramCont_first[[2]]=cbind(as.matrix(coef(sizeMod_first)),sd(residuals(sizeMod_first))) # the third column is for the standard deviation of growth 
# seedling size distribution is a uniform distribution (of exp(size_2)) with a min of 0.1 and a max 0f 3
paramCont_first[[3]]= cbind(as.matrix(coef(recMod_first)), sd(residuals(recMod_first)))
# model for probability of flowering is flwrMod_first
paramCont_first[[4]]=as.matrix(coef(flwrMod_first))
# model for seed production per plant (if reproductive) is seedMod_first
paramCont_first[[5]]=as.matrix(coef(seedMod_first))
# name the paramCont_first list to keep track of coefficients
names(paramCont_first) <- c("survival", "growth", "recruitDist", "flowering", "seedProduction")

# Construct an IPM kernel K using the parameters we obtained from the models

# First define the elements that make up the IPM (vital rates):

# SURVIVAL:
S.fun <- function(z, paramCont) {
  mu.surv=paramCont[["survival"]]["(Intercept)",] + paramCont[["survival"]]["log_LL_t",]*z
  return(1/(1 + exp(-(mu.surv))))
}

# GROWTH (we assume a constant variance)
GR.fun <- function(z,zz, paramCont){
  growth.mu = paramCont[["growth"]]["(Intercept)",1] + paramCont[["growth"]]["log_LL_t",1]*z
  return(dnorm(zz, mean = growth.mu, sd = paramCont[["growth"]][1,2]))
}

## SEEDLING SIZES (same approach as in growth function)
SDS.fun <- function(zz, paramCont){
  rec_mu <- paramCont[["recruitDist"]][1]
  rec_sd <- paramCont[["recruitDist"]][2]
  return(dnorm(zz, mean = rec_mu, sd = rec_sd))
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

K <- array(0,c(n+1,n+1))

# Setting up the kernels
b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint

h=(U-L)/n # bin width 

# Survival and growth 
S <- diag(S.fun(meshp, paramCont = paramCont_first)) # Survival # put survival probabilities in the diagonal of the matrix
G <- h * t(outer(meshp,meshp,GR.fun, paramCont = paramCont_first)) # Growth
# G <- t(outer(meshp,meshp,GR.fun)) # Growth

#Recruits distribution (seeds recruited from the seedbank into the continuous stage)
c_o <- h * matrix(rep(SDS.fun(meshp, paramCont = paramCont_first),n),n,n,byrow=F)
# c_o <- matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)

#Probability of flowering
Pb = (FL.fun(meshp, paramCont = paramCont_first))

#Number of seeds produced according to adult size
b_seed = (SDP.fun(meshp, paramCont = paramCont_first))

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
#Fkernel.discr <- as.matrix(h * Fkernel.discr)
Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))

mat_firstYear_DI <-Pkernel+Fkernel

eigenMat <- eigen(mat_firstYear_DI)
# get the lambda
lam_firstYear_DI <-  eigenMat$values[1]

eigenMat$vectors

#### DI, all-site data for second half of data ####
# Empty list to save model coefficients 
paramCont_second=list(NULL)
# survival model is called 'survMod_all'
paramCont_second[[1]]=as.matrix(coef(survMod_second)) # save coefficients 
# growth model is called 'sizeMod_all'
paramCont_second[[2]]=cbind(as.matrix(coef(sizeMod_second)),sd(residuals(sizeMod_second))) # the third column is for the standard deviation of growth 
# seedling size distribution is a uniform distribution (of exp(size_2)) with a min of 0.1 and a max 0f 3
paramCont_second[[3]]= cbind(as.matrix(coef(recMod_second)), sd(residuals(recMod_second)))
# model for probability of flowering is flwrMod_second
paramCont_second[[4]]=as.matrix(coef(flwrMod_second))
# model for seed production per plant (if reproductive) is seedMod_second
paramCont_second[[5]]=as.matrix(coef(seedMod_second))
# name the paramCont_second list to keep track of coefficients
names(paramCont_second) <- c("survival", "growth", "recruitDist", "flowering", "seedProduction")

# Construct an IPM kernel K using the parameters we obtained from the models
# First define the elements that make up the IPM (vital rates):

# SURVIVAL:
S.fun <- function(z, paramCont) {
  mu.surv=paramCont[["survival"]]["(Intercept)",] + paramCont[["survival"]]["log_LL_t",]*z
  return(1/(1 + exp(-(mu.surv))))
}

# GROWTH (we assume a constant variance)
GR.fun <- function(z,zz, paramCont){
  growth.mu = paramCont[["growth"]]["(Intercept)",1] + paramCont[["growth"]]["log_LL_t",1]*z
  return(dnorm(zz, mean = growth.mu, sd = paramCont[["growth"]][1,2]))
}

## SEEDLING SIZES (same approach as in growth function)
SDS.fun <- function(zz, paramCont){
  rec_mu <- paramCont[["recruitDist"]][1]
  rec_sd <- paramCont[["recruitDist"]][2]
  return(dnorm(zz, mean = rec_mu, sd = rec_sd))
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

K <- array(0,c(n+1,n+1))

# Setting up the kernels
b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint

h=(U-L)/n # bin width 

# Survival and growth 
S <- diag(S.fun(meshp, paramCont = paramCont_second)) # Survival # put survival probabilities in the diagonal of the matrix
G <- h * t(outer(meshp,meshp,GR.fun, paramCont = paramCont_second)) # Growth

#Recruits distribution (seeds recruited from the seedbank into the continuous stage)
c_o <- h * matrix(rep(SDS.fun(meshp, paramCont = paramCont_second),n),n,n,byrow=F)

#Probability of flowering
Pb = (FL.fun(meshp, paramCont = paramCont_second))

#Number of seeds produced according to adult size
b_seed = (SDP.fun(meshp, paramCont = paramCont_second))

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
#Fkernel.discr <- as.matrix(h * Fkernel.discr)
Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))

mat_secondYear_DI <-Pkernel+Fkernel

eigenMat <- eigen(mat_secondYear_DI)
# get the lambda
lam_secondYear_DI <- eigenMat$values[1]

eigenMat$vectors
#### hand-calculated model for continuous-ized seedlings and discrete seedbank --Density Dependent ####

### Define the demographic functions and parameters ###
## (Code below is modified from Maria Paniw)

# Empty list to save model coefficients 
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
eigenMat$values[1]

eigenMat$vectors

## calculate the uncertainty using bootstrapping
# make a vector to hold all of the lambdas from the resampling runs
allDD_lambdas <- numeric(200L)
# make a list to hold the model parameters
paramCont_now <- list()
allDD_params <- list()
# set seedbank parameters, which don't change
outSB <- outSB_all #SB to continuous stage
staySB <- staySB_all # staying in SB
goCont <- goCont_all # seeds become continuous right away (without going to the seed bank) 
goSB <- goSB_all # seeds go to the seedbank
surv.seeds <-  0.9 # survival of seeds
# Now, we refit the vital rate models, and refit the IPM
for(i in 1:200) {
  ## sample continuous data
  sample_ind <- seq(1, nrow(dat_all), by = 1)
  
  boot_data_ind   <- sample(sample_ind, size = 6000, replace = TRUE)
  
  boot_data <- dat_all[boot_data_ind,]
  
  ## fit vital rate models
  ## Survival ($s(z)$)
  survDat_now <- boot_data[boot_data$flowering==0 | is.na(boot_data$flowering),]
  survMod_now <- glm(survives_tplus1 ~ log_LL_t + N_Site_t, data = survDat_now, family = binomial)
  ## Growth ($G(z',z)$)
  sizeMod_now <- lm(log_LL_tplus1 ~ log_LL_t + N_Site_t, data = boot_data)
  ## Number of seeds produced, according to plant size ($b(z)$)
  seedDat_now <- boot_data[boot_data$flowering==1,]
  # fit  negative binomial model (for count data)
  seedMod_now <- MASS::glm.nb(Num_seeds ~ log_LL_t , data = seedDat_now)
  ## Flowering probability ($p_b(z)$)
  flwrMod_now <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) + N_Site_t , data = boot_data, family = binomial)))
  ## Distribution of recruit size ($c_o(z')$)
  # subset the data
  recD_now <- boot_data[boot_data$seedling == 1,]
  # fit the model
  recMod_now <- lm(log_LL_t ~ 1, data = recD_now)
  
  # update the parameter list with parameters from this iteration
  paramCont_now<- list(
    g_int     = coef(sizeMod_now)[1], # growth 
    g_slope   = coef(sizeMod_now)[2],
    g_dd      = coef(sizeMod_now)[3],
    g_sd      = summary(sizeMod_now)$sigma,
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
    outSB  = outSB_all,
    staySB = staySB_all,
    goSB   = goSB_all, 
    goCont = goCont_all                  
  )
  
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
  
  mat <-Pkernel+Fkernel
  
  eigenMat <- base::eigen(mat)
  # save the lambda
  allDD_lambdas[i] <- eigenMat$values[1]
  # save the parameters
  allDD_params[[i]] <- paramCont_now
}

### plot the values
# get the parameters out
# Plot the results
param_fig <- data.frame("param" = "g_int", 
                        "value" = sapply(X = allDD_params, FUN = function (x) x[[1]]))
param_fig <- rbind(param_fig,  data.frame("param" = "g_slope", 
                                          "value" = sapply(X = allDD_params,
                                                           FUN = function (x) x[[2]])))
param_fig <- rbind(param_fig,  data.frame("param" = "g_dd", 
                                          "value" = sapply(X = allDD_params,
                                                           FUN = function (x) x[[3]])))
param_fig <- rbind(param_fig, data.frame("param" = "g_sd", 
                                         "value" = sapply(X = allDD_params,
                                                          FUN = function (x) x[[4]])))
param_fig <- rbind(param_fig, data.frame("param" = "s_int", 
                                         "value" = sapply(X = allDD_params,
                                                          FUN = function (x) x[[5]])))
param_fig <- rbind(param_fig, data.frame("param" = "s_slope",
                                         "value" = sapply(X = allDD_params,
                                                          FUN = function (x) x[[6]])))
param_fig <- rbind(param_fig,  data.frame("param" = "s_dd", 
                                          "value" = sapply(X = allDD_params,
                                                           FUN = function (x) x[[7]])))
param_fig <- rbind(param_fig, data.frame("param" = "p_b_int", 
                                         "value" = sapply(X = allDD_params,
                                                          FUN = function (x) x[[8]])))
param_fig <- rbind(param_fig, data.frame("param" = "p_b_slope", 
                                         "value" = sapply(X = allDD_params,
                                                          FUN = function (x) x[[9]])))
param_fig <- rbind(param_fig, data.frame("param" = "p_b_slope_2", 
                                         "value" = sapply(X = allDD_params,
                                                          FUN = function (x) x[[10]])))
param_fig <- rbind(param_fig,  data.frame("param" = "p_b_dd", 
                                          "value" = sapply(X = allDD_params,
                                                           FUN = function (x) x[[11]])))
param_fig <- rbind(param_fig, data.frame("param" = "b_int", 
                                         "value" = sapply(X = allDD_params,
                                                          FUN = function (x) x[[12]])))
param_fig <- rbind(param_fig, data.frame("param" = "b_slope", 
                                         "value" = sapply(X = allDD_params,
                                                          FUN = function (x) x[[13]])))
param_fig <- rbind(param_fig, data.frame("param" = "c_o_mu", 
                                         "value" = sapply(X = allDD_params,
                                                          FUN = function (x) x[[14]])))
param_fig <- rbind(param_fig, data.frame("param" = "c_o_sd", 
                                         "value" = sapply(X = allDD_params,
                                                          FUN = function (x) x[[15]])))

# calculate the mean values
test <- param_fig%>% 
  group_by(param) %>% 
  summarise("mean_value" = mean(value))
param_fig <- left_join(param_fig, test)
# plot the distribution of values
ggplot(data = param_fig) + 
  geom_density(mapping = aes(value)) +
  geom_vline(aes(xintercept = mean_value, col = param)) + 
  geom_vline(aes(xintercept = 0), col = "grey", lty = 2) +
  facet_wrap(~ param, scales = "free") +
  scale_color_discrete(guide = "none") + 
  theme_classic()

# compare lambdas
plot(density(as.numeric(allDD_lambdas)))
abline(v = eigenMat$values[1], col = 'red', lwd = 2, lty = 2)
abline(v = mean(allDD_lambdas), col = "blue", lwd = 2, lty = 2)

#### DI IPM for each site ####
## write vital rate functions
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
  
  site_IPMs_DI[[i]] <- mat
}
names(site_IPMs_DI) <- paste0(unique(dat_all$Site))
lambdas_site_DI <- sapply(site_IPMs_DI, FUN = function(x) eigen(x)$values[1])

## estimate bootstrap confidence intervals for each model parameter for each site
# make an empty list to hold the estimate results (parameters, labmda)
siteDI_bootCIs <- list()
# set seedbank parameters, which don't change
outSB <- outSB_all #SB to continuous stage
staySB <- staySB_all # staying in SB
goCont <- goCont_all # seeds become continuous right away (without going to the seed bank) 
goSB <- goSB_all # seeds go to the seedbank
surv.seeds <-  0.9 # survival of seeds
for (i in 1:length(unique(dat_all$Site))) {
  # get the data just for this site
  dat_now <- dat_all[dat_all$Site == unique(dat_all$Site)[i],]
  # make a vector to hold all of the lambdas from the resampling runs
  siteDI_lambdas <- numeric(200L)
  # make a list to hold the model parameters (output)
  siteDI_params <- list()
  
  # Now, we refit the vital rate models, and refit the IPM
  for(j in 1:200) {
    ## sample continuous data
    sample_ind <- seq(1, nrow(dat_now), by = 1)
    
    boot_data_ind   <- sample(sample_ind, size = 300, replace = TRUE)
    
    boot_data <- dat_now[boot_data_ind,]
    
    ## fit vital rate models
    ## Survival ($s(z)$)
    survDat_now <- boot_data[boot_data$flowering==0 | is.na(boot_data$flowering),]
    survMod_now <- glm(survives_tplus1 ~ log_LL_t , data = survDat_now, family = binomial)
    ## Growth ($G(z',z)$)
    sizeMod_now <- lm(log_LL_tplus1 ~ log_LL_t , data = boot_data)
    ## Number of seeds produced, according to plant size ($b(z)$)
    seedDat_now <- boot_data[boot_data$flowering==1,]
    # fit  negative binomial model (for count data)
    seedMod_now <- MASS::glm.nb(Num_seeds ~ log_LL_t , data = seedDat_now)
    ## Flowering probability ($p_b(z)$)
    flwrMod_now <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2)  , data = boot_data, family = binomial)))
    ## Distribution of recruit size ($c_o(z')$)
    # subset the data
    recD_now <- boot_data[boot_data$seedling == 1,]
    # fit the model
    recMod_now <- lm(log_LL_t ~ 1, data = recD_now)
    
    # update the parameter list with parameters from this iteration
    paramCont<- list(
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
    # multiply the continuous kernel by the binwidth (h)
    # Pkernel.cont <- h * Pkernel.cont
    
    # seedbank (first column of your K)
    Pkernel.seedbank = c(staySB, outSB*c_o[,1]) # seeds survive and go to continuous
    
    # Make the full P kernel
    Pkernel <- cbind(Pkernel.seedbank,rbind(rep(0,length(meshp)),Pkernel.cont)) # discrete component
    
    ## make the F kernel
    Fkernel.cont <-  as.matrix(goCont * ((c_o) %*% diag(FecALL))) # the size of seedlings that go into the seed bank from each continuous size class
    Fkernel.discr  <- matrix(c(0, goSB * (FecALL)), nrow = 1)
    Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))
    
    mat <-Pkernel+Fkernel
    
    # save the lambda
    siteDI_lambdas[j] <-base::eigen(mat)$values[1]
    # save the parameters
    siteDI_params[[j]] <- paramCont
  }
  
  ### plot the values
  # get the parameters out
  # Plot the results
  param_fig <- data.frame("param" = "g_int", 
                          "value" = sapply(X = siteDI_params, FUN = function (x) x[[1]]))
  param_fig <- rbind(param_fig,  data.frame("param" = "g_slope", 
                                            "value" = sapply(X = siteDI_params,
                                                             FUN = function (x) x[[2]])))
  param_fig <- rbind(param_fig, data.frame("param" = "g_sd", 
                                           "value" = sapply(X = siteDI_params,
                                                            FUN = function (x) x[[3]])))
  param_fig <- rbind(param_fig, data.frame("param" = "s_int", 
                                           "value" = sapply(X = siteDI_params,
                                                            FUN = function (x) x[[4]])))
  param_fig <- rbind(param_fig, data.frame("param" = "s_slope",
                                           "value" = sapply(X = siteDI_params,
                                                            FUN = function (x) x[[5]])))
  param_fig <- rbind(param_fig, data.frame("param" = "p_b_int", 
                                           "value" = sapply(X = siteDI_params,
                                                            FUN = function (x) x[[6]])))
  param_fig <- rbind(param_fig, data.frame("param" = "p_b_slope", 
                                           "value" = sapply(X = siteDI_params,
                                                            FUN = function (x) x[[7]])))
  param_fig <- rbind(param_fig, data.frame("param" = "p_b_slope_2", 
                                           "value" = sapply(X = siteDI_params,
                                                            FUN = function (x) x[[8]])))
  param_fig <- rbind(param_fig, data.frame("param" = "b_int", 
                                           "value" = sapply(X = siteDI_params,
                                                            FUN = function (x) x[[9]])))
  param_fig <- rbind(param_fig, data.frame("param" = "b_slope", 
                                           "value" = sapply(X = siteDI_params,
                                                            FUN = function (x) x[[10]])))
  param_fig <- rbind(param_fig, data.frame("param" = "c_o_mu", 
                                           "value" = sapply(X = siteDI_params,
                                                            FUN = function (x) x[[11]])))
  param_fig <- rbind(param_fig, data.frame("param" = "c_o_sd", 
                                           "value" = sapply(X = siteDI_params,
                                                            FUN = function (x) x[[12]])))
  
  # save the parameters in 'param_fig' format, as well as lambdas
siteDI_bootCIs[[i]] <- list("params" = param_fig, "lambdas" = siteDI_lambdas)
}

#### DI IPM for each site--first half of data ####
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
site_IPMs_first_DI <- list()
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
  
  site_IPMs_first_DI[[i]] <- mat
}
names(site_IPMs_first_DI) <- paste0(unique(dat_all$Site))
lambdas_site_first_DI <- sapply(site_IPMs_first_DI, FUN = function(x) eigen(x)$values[1])

#### DI IPM for each site--second half of data ####
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
site_IPMs_second_DI <- list()
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
  
  site_IPMs_second_DI[[i]] <- mat
}
names(site_IPMs_second_DI) <- paste0(unique(dat_all$Site))
lambdas_site_second_DI <- sapply(site_IPMs_second_DI, FUN = function(x) eigen(x)$values[1])

#### DD IPM for each site ####
## define the vital rate models
# SURVIVAL:
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
    outSB  = outSB_all,
    staySB = staySB_all,
    goSB   = goSB_all, 
    goCont = goCont_all                  
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
