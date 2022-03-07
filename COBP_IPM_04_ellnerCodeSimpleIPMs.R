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
  
  # how Maria did it
  # # Get residual variance 
  # var.res= (paramCont[["growth"]][1,2])^2
  # # Density distribution function of the normal distribution
  # gr1 = sqrt(2*pi*var.res)
  # gr2 = ((zz-growth.mu)^2)/(2*var.res)
  # 
  # return(exp(-gr2)/gr1)
  
  # makes more sense to me... (got the same answer)
  return(dnorm(zz, mean = growth.mu, sd = paramCont[["growth"]][1,2]))
  
}

## SEEDLING SIZES (same approach as in growth function)

SDS.fun <- function(zz){
  # how Maria did it: 
  # sds.mu=(paramCont[[3]][year,"(Intercept)"])
  # 
  # # Get residual variance 
  # var.res= (paramCont[[3]][1,2])^2
  # # Density distribution function of the normal distribution
  # sds1 = sqrt(2*pi*var.res)
  # sds2 = ((zz-sds.mu)^2)/(2*var.res)
  # 
  # return(exp(-sds2)/sds1)
  
  # I'll do it differently, since our distribution of seedling size is different
  rec_mu <- paramCont[["recruitDist"]][1]
  rec_sd <- paramCont[["recruitDist"]][2]
  
  return(dnorm(zz, mean = rec_mu, sd = rec_sd))
  # try the uniform dist. 
  #return(dunif(exp(zz), min = 0.1, max = 3))
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

mat <-Pkernel+Fkernel

eigenMat <- eigen(mat)
# get the lambda
eigenMat$values[1]

eigenMat$vectors

## 
## 
# compare the P and F continuous kernels
# the continuous part of the P kernel is the same between the ipmr and by-hand models! 
Pkernel.cont - contSeedlings_IPM$sub_kernels$P
# the contiuous F kernel and the ipmr F kernela are the same! 
Fkernel.cont - contSeedlings_IPM$sub_kernels$F

# compare the distributions for disc transitions
plot(x = meshp, y = contSeedlings_IPM$sub_kernels$seedbank_to_continuous, type = 'l')
lines(x = meshp, y = mat[2:501,1])
lines(x = meshp, y = c( outSB*c_o[,1])/h)
# the same! 

plot(x = meshp, y = contSeedlings_IPM$sub_kernels$continuous_to_seedbank, type = 'l')
lines(x = meshp, y = Fkernel.discr[2:501])
lines(x = meshp, y = matrix(c( goSB * h * (FecALL)), nrow = 1))
# not the same :-( 

library(popbio)

popbio::lambda(mat)
eigenList <- eigen.analysis(mat)
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
  # try the uniform dist. 
  #return(dunif(exp(zz), min = 0.1, max = 3))
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

mat <-Pkernel+Fkernel

eigenMat <- base::eigen(mat)
# get the lambda
eigenMat$values[1]

eigenMat$vectors

## 
## 
# compare the P and F continuous kernels
# the continuous part of the P kernel is the same between the ipmr and by-hand models! 
Pkernel.cont - contSeedlings_IPM$sub_kernels$P
# the contiuous F kernel and the ipmr F kernela are the same! 
Fkernel.cont - contSeedlings_IPM$sub_kernels$F

# compare the distributions for disc transitions
plot(x = meshp, y = contSeedlings_IPM$sub_kernels$seedbank_to_continuous, type = 'l')
lines(x = meshp, y = mat[2:501,1])
lines(x = meshp, y = c( outSB*c_o[,1])/h)
# the same! 

plot(x = meshp, y = contSeedlings_IPM$sub_kernels$continuous_to_seedbank, type = 'l')
lines(x = meshp, y = Fkernel.discr[2:501])
lines(x = meshp, y = matrix(c( goSB * h * (FecALL)), nrow = 1))
# not the same :-( 

library(popbio)

popbio::lambda(mat)
eigenList <- eigen.analysis(mat)

