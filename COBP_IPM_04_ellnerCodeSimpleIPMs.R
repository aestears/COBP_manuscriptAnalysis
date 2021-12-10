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

