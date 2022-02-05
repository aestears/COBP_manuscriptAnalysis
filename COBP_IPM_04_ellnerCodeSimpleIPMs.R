#//////////////////////////
# Integral Projection Models for Oenothera coloradensis: simple model using code from Ellner and Rees
# Alice Stears
# 08 December 2021
#/////////////////////////

#### load vital rate models from previous script ####
source("./analysis_scripts/COBP_IPM_02_VitalRateModels.R")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Section 1 - Define the demographic functions and parameters ####
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the true parameter vector, each parameter is given a name so the formulae, below, are easier to
## read. We'll use 'model.term' scheme to name elements of the vector

m.par.true <- c(## survival
  surv.int  =  coef(survMod)[1],
  surv.z    =   coef(survMod)[2],
  ## flowering
  flow.int  = coef(flwrMod_t)[1],
  flow.z    = coef(flwrMod_t)[2],
  flow.z.2  = coef(flwrMod_t)[3],
  ## growth
  grow.int  =   coef(sizeMod)[1],
  grow.z    =   coef(sizeMod)[2],
  grow.sd   =   summary(sizeMod)$sigma,
  ## recruit size
  rcsz.int  =   coef(recMod), 
  rcsz.sd   =   summary(recMod)$sigma,
  ## seed production by size
  seed.int  =   coef(seedMod_t)[1],
  seed.z    =   coef(seedMod_t)[2],
  ## recruitment probability
  p.r       =   p.estab.simple)  
names(m.par.true) <- c("surv.int", "surv.z", "flow.int", "flow.z", "flow.z.2", "grow.int", "grow.z", "grow.sd", "rcsz.int", "rcsz.sd", "seed.int", "seed.z", "p.r")


## Growth function, given you are size z now returns the pdf of size z1 next time

G_z1z <- function(z1, z, m.par)
{
  mu <- m.par["grow.int"] + m.par["grow.z"] * z           # mean size next year
  sig <- m.par["grow.sd"]                                    # sd about mean
  p.den.grow <- dnorm(z1, mean = mu, sd = sig)             # pdf that you are size z1 given you were size z
  return(p.den.grow)
}

## Survival function, logistic regression
s_z <- function(z, m.par)
{
  #linear.p <- m.par["surv.int"] + m.par["surv.z"] * z  # linear predictor
  p <- 1/(1+exp(-( m.par["surv.int"] + m.par["surv.z"] * z))) # logistic transformation to probability
  return(p)
}

## Probability of flowering function, logistic regression

p_bz <- function(z, m.par)
{
  # linear.p <- m.par["flow.int"] + m.par["flow.z"] * z + m.par["flow.z.2"] * (z^2)     # linear predictor
  p <- 1/(1+exp(-(m.par["flow.int"] + m.par["flow.z"] * z + m.par["flow.z.2"] * (z^2))))  # logistic transformation to probability
  return(p)
}

## Seed production function

b_z <- function(z, m.par)
{
  N <- exp(m.par["seed.int"] + m.par["seed.z"] * z)    # seed production of a size z plant
  return(N)
}

## Recruit size pdf

c_0z1 <- function(z1, m.par)
{
  mu <- m.par["rcsz.int"]
  sig <- m.par["rcsz.sd"]
  p.deRecr <- dnorm(z1, mean = mu, sd = sig)             # pdf of a size z1 recruit
  return(p.deRecr)
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Section 2 - Functions to build IPM kernels P, F, and K ####
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the survival kernel
P_z1z <- function (z1, z, m.par) {
  
  return((1 - p_bz(z, m.par)) * s_z(z, m.par) * G_z1z(z1, z, m.par))
  
}

## Define the fecundity kernel
F_z1z <- function (z1, z, m.par) {
  
  return( p_bz(z, m.par) * b_z(z, m.par) * m.par["p.r"] * c_0z1(z1, m.par))
  
}

mk_K <- function(m, m.par, L, U) {
  
  # mesh points 
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h
  P <- h * (outer(meshpts, meshpts, P_z1z, m.par = m.par))
  F <- h * (outer(meshpts, meshpts, F_z1z, m.par = m.par))
  K <- P + F
  return(list(K = K, meshpts = meshpts, P = P, F = F))
}


# flexible size limits, defaults set for Oenothera model
mk_K_ceiling <- function(m, m.par, L, U, U1 = U) {
  # mesh points 
  h <- (U - L)/m;
  meshpts <- L + ((1:m) - 1/2) * h;
  P <- h * (outer(meshpts, pmin(meshpts,U1), P_z1z, m.par = m.par));
  F <- h * (outer(meshpts, pmin(meshpts,U1), F_z1z, m.par = m.par));
  K <- P + F;
  return(list(K = K, meshpts = meshpts, P = P, F = F))
}

#Function to calculate mean and variance in reproductive output for a size z

get_mean_var_Repr <- function(init.z,n.samp) {
  
  # initial population sizes and ages
  z   <- rep(init.z,1000)
  Repr.out <- NULL
  
  repeat {
    
    ## calculate population size
    pop.size <- length(z)
    
    ## generate binomial random number for the probability of flowering, where the probability of flowering
    ## depends on your size z, this is a vector of 0's and 1's, you get a 1 if you flower
    Repr <- rbinom(n=pop.size, prob=p_bz(z, m.par.true), size=1)
    
    ## number of plants that flowered
    num.Repr <- sum(Repr)
    
    if(num.Repr>0) {
      Seeds <- rpois(num.Repr, m.par.true["p.r"] * b_z(z[Repr==1],m.par.true))
      Repr.out <- c(Repr.out,Seeds)
    }
    
    ## generate new recruit sizes
    ## rnorm generated normally distributed random numbers
    Rcsz <- rep(init.z,100)
    
    ## for the non-reproductive plants generate random number for survival
    Surv <- rep(NA, pop.size)
    Surv[Repr==0] <- rbinom(n = pop.size - num.Repr, prob = s_z(z[Repr==0], m.par.true), size = 1)
    num.die <- sum(Surv==0, na.rm=TRUE)
    
    if(num.die>0) Repr.out <- c(Repr.out,rep(0,num.die))
    
    ## index for individuals that did not flower and survived
    i.subset <- which(Repr==0 & Surv==1)
    
    ## let them grow
    E.z1 <- m.par.true["grow.int"]+m.par.true["grow.z"]*z[i.subset]
    z1 <- rnorm(n = pop.size - num.Repr - num.die, mean = E.z1, sd = m.par.true["grow.sd"])
    
    z <- c(Rcsz, z1)
    
    if(length(Repr.out)>n.samp) break
    
  }
  
  return(c(mean(Repr.out),var(Repr.out),mean(Repr.out>0),var(Repr.out>0)))
  
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
#### Section 3 - Construct Kernels and projection population size ####
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nBigMatrix <- 250
min.size <- min(dat$log_LL_t, na.rm = TRUE) * 1.2 # make lower limit 20% smaller than actual data (multiply by 1.2 since it is a negative number)
max.size <- max(dat$log_LL_t, na.rm = TRUE) * 1.2 # make upper limit 20% larger than actual data


# so let's set the lower and upper limits for the size range at -2.65
# and 4.5 so slightly smaller/bigger than observed

IPM.true <- mk_K(m = 500, m.par = m.par.true, L = min.size, U = max.size)

lambda_simple_math <- Re(eigen(IPM.true$K)$values[1])
## lambda from this mathematical method is the same as the lambda from the ipmr method! 

meshpts <- IPM.true$meshpts

# find the stable size distribution and mean size
w.true <- Re(eigen(IPM.true$K)$vectors[, 1])
stable.z.dist.true <- w.true/sum(w.true)
mean.z.true <- sum(stable.z.dist.true * meshpts)

# find the stable flowering size distribution and mean flowering size
wb.true <- p_bz(meshpts, m.par.true) * w.true
stable.flowering.dist.true <- wb.true/sum(wb.true)
mean.flowering.z.true <- sum(stable.flowering.dist.true * meshpts)

dev.new(); 
par(mfrow = c(2,2))
## 1 - plot population density versus time...
pop_size <- dat %>% 
  group_by(Year) %>% 
  summarize(N_all = length(log_LL_t)) %>% 
  mutate(Year = as.numeric(as.character(Year)))
plot(pop_size$Year, pop_size$N_all, type = "l", xlab = "Time", ylab = "Population size", ylim = c(1600,4000))
lines(x = c(2018, 2019, 2020), y = c(pop_size[1,2], pop_size[1,2]*lambda_simple_math, (pop_size[1,2]*lambda_simple_math)*lambda_simple_math), col = "blue")

## 2 - plot mean size versus time...
mean_size <- dat %>% 
  group_by(Year) %>% 
  summarize(mean_LL = (sum(log_LL_t, na.rm = TRUE)/length(log_LL_t))) %>% 
  mutate(Year = as.numeric(as.character(Year)))
plot(x = mean_size$Year, y = mean_size$mean_LL, type = "l", xlab = "Time", ylab = "ln(Mean plant size)", ylim = c(0,3))
abline(h = mean.z.true, col = "red")

## 3 - plot mean flowering size versus time...
mean_flwr <- dat %>% 
  filter(flowering == 1) %>% 
  group_by(Year) %>% 
  summarize(mean_LL = (sum(log_LL_t, na.rm = TRUE)/length(log_LL_t))) %>% 
  mutate(Year = as.numeric(as.character(Year)))
plot(x = mean_flwr$Year, y = mean_flwr$mean_LL, type = "l", xlab = "Time", ylab = "ln(Mean flowering plant size)", ylim = c(0,3))
abline(h = mean.flowering.z.true, col = "red")

## 4 - plot of density estimates at time 50 and the end
plot(density(dat$log_LL_t, na.rm = TRUE), ylim = c(0, .9), xlab = "Plant size", main = "")
lines(IPM.true$meshpts, stable.z.dist.true/diff(IPM.true$meshpts)[1], col = "red")

# dev.copy2eps(file = "../../figures/c2/OenotheraSim.eps")

#### hand-calculated model for continuous-ized seedlings and discrete seedbank ####
# Density Independent 

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Section 1 - Define the demographic functions and parameters ####
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the true parameter vector, each parameter is given a name so the formulae, below, are easier to
## read. We'll use 'model.term' scheme to name elements of the vector

m.par.true <- c(## survival
  surv.int  =  coef(survMod_N)[1],
  surv.z    =   coef(survMod_N)[2],
  surv.N   = coef(survMod_N)[3],
  ## flowering
  flow.int  = coef(flwrMod_N)[1],
  flow.z    = coef(flwrMod_N)[2],
  flow.z.2  = coef(flwrMod_N)[3],
  flow.N   = coef(flwrMod_N)[4],
  ## growth
  grow.int  =   coef(sizeMod_N)[1],
  grow.z    =   coef(sizeMod_N)[2],
  grow.sd   =   summary(sizeMod_N)$sigma,
  grow.N    =   coef(sizeMod_all)[3],
  ## recruit size
  rcsz.mu  =   coef(recMod_all), #recruit size distribution 
  rcsz.sd  =   sd(residuals(recMod_all)),
  ## seed production by size
  seed.int  =   coef(seedMod_N)[1],
  seed.z    =   coef(seedMod_N)[2],
  ## Probability that a seed from the seedbank in year t will germinate to a seedling in year t+1 
  outSB.est <- outSB_all,
  ## Probability that a seed from the seedbank in year t will stay in the seedbank in year t+1 
  staySB.est <- staySB_all, 
  ## Probability that a seed produced by an adult plant in year t will enter the seedbank in year t+1 
  goSB.est <- goSB_all,
  ## Probability that a seed from a plant in year t will go directly to the continuous stage 
  goCont.est <- goCont_all
  )  
names(m.par.true) <- c("surv.int", "surv.z", "surv.N", "flow.int", "flow.z", "flow.z.2","flow.N", "grow.int", "grow.z", "grow.sd", "grow.N", "rcsz.int", "rcsz.sd", "seed.int", "seed.z", "outSB.est", "staySB.est", "goSB.est", "goCont.est")

## (Code below is from Maria Paniw)

# Empty list to save model coefficients 
paramCont=list(NULL)

# survival model is called 'survMod_all'
paramCont[[1]]=as.matrix(coef(survMod_all)) # save coefficients 

# growth model is called 'sizeMod_all'
paramCont[[2]]=cbind(as.matrix(coef(sizeMod_all)),sd(residuals(sizeMod_all))) # the third column is for the standard deviation of growth 

# seedling size distribution is a uniform distribution (of exp(size_2)) with a min of 0.1 and a max 0f 3
paramCont[[3]]= cbind(as.matrix(coef(recMod_all)), sd(residuals(sizeMod_all)))

# model for probability of flowering is flwrMod_all
paramCont[[4]]=as.matrix(coef(flwrMod_all))

# model for seed production per plant (if reproductive) is seedMod_all
paramCont[[5]]=as.matrix(coef(seedMod_all))

# name the paramCont list to keep track of coefficients
names(paramCont) <- c("survival", "growth", "recruitDist", "flowering", "seedProduction")

###########################################################################
### PART B - IPM and SIMULATIONS BASED ON KERNEL SELECTION
##########################################################################

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

## Define the survival kernel function
P_cont_fun <- function (zz, z) {
  
  return((1 - FL.fun(z)) * S.fun(z) * GR.fun(zz, z))
  
}

## Define the fecundity kernel function
F_cont_fun <- function (zz, z) {
  
  return(goCont * FL.fun(z) * SDP.fun(z) * SDS.fun(zz))
  
}

## Define the kernel funtion for the transition from continuous to seedbank
C_SB_cont_fun <- function (zz, z) {
  return(goSB * (1 - FL.fun(z)) * SDP.fun(z))
}

# Second, put together the kernels - four kernels for four years:
years <- 1 # for right now, treating everything as one year

# Define the lower and upper integration limit
L <-  1.2 * min(dat_all$log_LL_t, na.rm = TRUE) # minimum size
U <-  1.2 * max(dat_all$log_LL_t, na.rm = TRUE) # maximum size

n <-500 # bins

# These are the parameters for the discrete stages
# I usually only have seed banks (SB), but now I added a seedling stage

outSB <- 0.183 #SB to continuous stage
staySB <- 0.717 # staying in SB
goCont <- 0.119 # seeds become continuous right away (without going to the seed bank) 
goSB <- 0.466 # seeds go to the seedbank
surv.seeds <-  0.9 # survival of seeds

K <- array(0,c(n+1,n+1))

# I recommend you set i = 1, set n low to say 10 
  
  # Seeting up the kernels
  b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
  meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
  
  h=(U-L)/n # bin width 
  
  # Survival and growth 
  S <- diag(S.fun(meshp)) # Survival # put survival probabilities in the diagonal of the matrix
  G <- h*t(outer(meshp,meshp,GR.fun)) # Growth
  
  #Recruits distribution (seeds recruited from the seedbank into the continuous stage)
  c_o <- h*matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
  
  #Probability of flowering
  Pb = (FL.fun(meshp))
  
  #Number of seeds produced according to adult size
  b = (SDP.fun(meshp))
  
  FecALL= Pb * b
  
  # ammend the 'S' matrix by multiplying it by (1-Pb), since this is a monocarpic perennial
  S_new <- S * (1-Pb)
  
  # Control for eviction:
  # this is equivalent to redistributing evicted sizes evenly among existing size classes 
  G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  c_o <- c_o/matrix(as.vector(apply(c_o,2,sum)),nrow=n,ncol=n,byrow=TRUE)

  # make the continuous part of the P matrix
  Pkernel.cont <- as.matrix(G %*% S_new)
  
  # seedling (second column of your K)
  Pkernel.seedbank = c(staySB,outSB*c_o[,1]) # seeds survive and go to continuous
  
  # Make the full P kernel
  Pkernel <- cbind(Pkernel.seedbank,rbind(rep(0,length(meshp)),Pkernel.cont)) # discrete component
  
  ## make the F kernel
  Fkernel.cont <-  as.matrix(goCont * ( c_o %*% diag(FecALL))) # the size of seedlings that go into the seed bank from each 
  Fkernel.discr  <- matrix(c(0, goSB * FecALL), nrow = 1)
  Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))
  
  mat <-Pkernel+Fkernel
  

## 

library(popbio)

popbio::lambda(mat)

## the lambda is higher than when I use ipmr?? email Maria to see if it is correct? 