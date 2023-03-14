#//////////////////////////
# Integral Projection Models for Oenothera coloradensis
# Alice Stears
# 30 November 2021
#/////////////////////////

#### Load packages ####
library(ipmr)

#### load vital rate models from previous script ####
source("./analysis_scripts/01_VitalRateModels.R")

#### IPM A #### 
## deterministic, density-independent IPM using only continuous stages (no seedlings) and data from all transitions ##
## calculate starting population state vectors

## define a function to get predictions from a glm with a logit link
inv_logit <- function(lin.pred) {
  1/(1 + exp(-(lin.pred)))
}

# define starting domain bounds
L <- min(dat_all$log_LL_t, na.rm = TRUE) * 1.2 # lower bound (L)
U <- max(dat_all$log_LL_t, na.rm = TRUE) * 1.2 # upper bound (U)
n <- 500 

# calculate probability of establishment (skipping seedbank stage)
# calculate the number of seeds produced by each fruiting adult
ipm_A_seeds <- dat_all %>% 
  group_by(Plot_ID, Year) %>% 
  summarize("N_seeds_t" = sum(Num_seeds, na.rm = TRUE)) %>% 
  rename(Year_t = Year) %>% 
  mutate(Year_t = as.numeric(as.character(Year_t)))
# calculate the number of seedlings recruited in each year
ipm_A_recruits <- dat_all %>% 
  group_by(Plot_ID, Year) %>% 
  summarize("N_recruits_t" = length(seedling)) %>% 
  rename(Year_t = Year) %>% 
  mutate(Year_tplus1 = as.numeric(as.character(Year_t)) - 1)

#combine # seeds in previous year w/ number of seedlings
ipm_A_estabs <- ipm_A_seeds %>% 
  left_join(ipm_A_recruits, 
            by = c("Year_t" = "Year_tplus1",
                   "Plot_ID" = "Plot_ID")) %>% 
  dplyr::select(- Year_t.y) %>% 
  rename("N_recruits_tplus1" = "N_recruits_t") %>% 
  mutate("p_estab" = N_recruits_tplus1/N_seeds_t)
  
ipm_A_estabs[ipm_A_estabs$p_estab == Inf & 
                is.na(ipm_A_estabs$p_estab) == FALSE , 
              "p_estab"] <- NA

# the probability of establishment in the IPM is the average p(estab)
p.estab.ipm_A = mean(ipm_A_estabs$p_estab, na.rm = TRUE)

data_list <- list(
  g_int     = coef(sizeMod_all)[1],
  g_slope   = coef(sizeMod_all)[2],
  g_sd      = summary(sizeMod_all)$sigma,
  s_int     = coef(survMod_all)[1],
  s_slope   = coef(survMod_all)[2],
  p_b_int   = coef(flwrMod_all)[1], #probability of flowering
  p_b_slope = coef(flwrMod_all)[2],
  p_b_slope_2 = coef(flwrMod_all)[3],
  b_int   = coef(seedMod_all)[1], #seed production
  b_slope = coef(seedMod_all)[2],
  c_o_mu    = coef(recMod_all), #recruit size distribution
  c_o_sd    = summary(recMod_all)$sigma,
  p_estab = p.estab.ipm_A 
  )

# inital population state
init_size_state <- runif(500)

ipm_A <- init_ipm(sim_gen   = "simple", 
                       di_dd     = "di", 
                       det_stoch = "det") %>% 
  define_kernel(
    name          = "P",
    formula       =(1-p_b.) * s. * g.,
    
    s.            = 1/(1 + exp(-(s_int + s_slope * size_1))),
    g.            = dnorm(size_2, g_mu., g_sd), 
    g_mu.         = g_int + g_slope * size_1, 
    p_b.          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2)))),
    
    family        = "CC",
    data_list     = data_list,
    states        = list(c('size')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "g.")
    ) %>% 
  define_kernel(
    name          = "F", 
    formula       = p_b. * b. * p_estab. * c_o.,
    
    p_b.          = 1/(1 + exp(-(p_b_int + p_b_slope * size_1 + p_b_slope_2 * (size_1^2)))),
    b.            = exp(b_int + b_slope * size_1),
    c_o.          = dnorm(size_2, c_o_mu, c_o_sd),
    p_estab.       = p_estab,
    
    family        = "CC",
    data_list     = data_list,
    states        = list(c('size')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "c_o.")
    ) %>% 
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "F"), 
      int_rule = rep("midpoint", 2),
      state_start = rep("size", 2), 
      state_end = rep("size", 2)
      )
    ) %>% 
  define_domains(
    size = c(
      min(dat_all$log_LL_t, na.rm = TRUE) * .8, # lower bound (L)
      max(dat_all$log_LL_t, na.rm = TRUE) * 1.2, # upper bound (U)
      500 # number of mesh points
    )
    ) %>% 
  define_pop_state(
    n_size = runif(500)
      ) %>% 
  make_ipm(
    iterations = 1000
    )

lambda(ipm_A)

## estimate CIs using bootstrap resampling
## save the proto-ipm 
ipm_A_proto <- ipm_A$proto_ipm
# lists to hold outputs
ipm_A_CI_lambdas <- list()
ipm_A_CI_params <- list()
for(i in 1:1000) {
  ## sample continuous data
  sample_ind <- seq(1, nrow(dat_all), by = 1)
  
  boot_data_ind   <- sample(sample_ind, size = nrow(dat_all), replace = TRUE)
  
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
  seedMod_t_now <- glm(Num_seeds ~ log_LL_t , data = seedDat_now, family = poisson)
  ## Flowering probability ($p_b(z)$)
  flwrMod_t_now <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = boot_data, family = binomial)))
  ## Distribution of recruit size ($c_o(z')$)
  # subset the data
  recD_now <- boot_data[boot_data$seedling == 1 & is.na(boot_data$age) == FALSE,]
  # fit the model
  recMod_now <- lm(log_LL_t ~ 1, data = recD_now)
  
  # update the parameter list with parameters from this iteration
  data.list_now <- list(
    g_int     = coef(sizeMod_now)[1],
    g_slope   = coef(sizeMod_now)[2],
    g_sd      = summary(sizeMod_now)$sigma,
    s_int     = coef(survMod_now)[1],
    s_slope   = coef(survMod_now)[2],
    p_b_int   = coef(flwrMod_t_now)[1], #probability of flowering
    p_b_slope = coef(flwrMod_t_now)[2],
    p_b_slope_2 = coef(flwrMod_t_now)[3],
    b_int   = coef(seedMod_t_now)[1], #seed production
    b_slope = coef(seedMod_t_now)[2],
    c_o_mu    = coef(recMod_now), #recruit size distribution
    c_o_sd    = summary(recMod_now)$sigma, 
    p_estab = p.estab.ipm_A
  )
  
  # Insert the new vital rate models into the proto_ipm, then rebuild the IPM.
  
  parameters(ipm_A_proto) <- data.list_now
  
  boot_ipm <- ipm_A_proto %>%
    make_ipm(iterate = TRUE,
             normalize_pop_size = FALSE,
             iterations = 100)
  
  ipm_A_CI_lambdas[i] <- lambda(boot_ipm)
  
  ipm_A_CI_params[[i]] <- data.list_now
}

## calculate 95% CI for lambda
ipm_A_CI_lambdasVec <- log(sapply(ipm_A_CI_lambdas, unlist))
SE <- sd(ipm_A_CI_lambdasVec)/sqrt(1000)
mean <- mean(ipm_A_CI_lambdasVec)
ipm_A_CI_log <- c((mean - 1.96*SE),(mean + 1.96*SE))

## check for eviction from the model
#To check for eviction, we plot the survival model and the column sums of the survival/growth (P) matrix. Eviction occurs when the column sums are lower than the survival models suggests that they should be.
# define the x-axis values
meshpts <- seq(from = (min(dat_all$log_LL_t, na.rm = TRUE) * .8), to = (max(dat_all$log_LL_t, na.rm = TRUE) * 1.2) , length.out = 500)
# plot the model-predicted survival probs.
preds <- predict(object = survMod_all, newdata = data.frame("log_LL_t" = meshpts), type = 'response')
plot(x = meshpts, y = preds, ylab = "Survival Probability", type = "l")
# plot the survival values from the P matrix
points(meshpts,apply(ipm_A$sub_kernels$P,2,sum),col="red",lwd=3,cex=.1,pch=19)

#### IPM B ####
### Deterministic, density-independent IPM with all sites, all years, all continuous data + seedbank ###

data_list <- list(
  g_int     = coef(sizeMod_all)[1], # growth 
  g_slope   = coef(sizeMod_all)[2],
  g_sd      = sd(residuals(sizeMod_all)),
  s_int     = coef(survMod_all)[1], # survival
  s_slope   = coef(survMod_all)[2],
  p_b_int   = coef(flwrMod_all)[1], #probability of flowering
  p_b_slope = coef(flwrMod_all)[2],
  p_b_slope_2 = coef(flwrMod_all)[3],
  b_int   = coef(seedMod_all)[1], #seed production
  b_slope = coef(seedMod_all)[2],
  c_o_mu    = coef(recMod_all), #recruit size distribution
  c_o_sd    = sd(residuals(recMod_all)), 
  outSB  = outSB_all,
  staySB = staySB_all,
  goSB   = goSB_all, 
  goCont = goCont_all                  
)
# inital population state
init_size_state <- runif(500)

ipm_B <- init_ipm(sim_gen   = "general", 
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
    data_list     = data_list,
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
    data_list     = data_list,
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
    data_list     = data_list,
    states        = list(c('size')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions("norm", "c_o.")
  ) %>% define_kernel(
    name          = "seedbank_to_seedbank", 
    formula       = staySB.,
    
    staySB.       = staySB,
    
    family        = "DD",
    data_list     = data_list,
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
    data_list     = data_list,
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
    iterations = 1000
  )

ipmr::lambda(ipm_B)

## check for eviction from the model
#To check for eviction, we plot the survival model and the column sums of the survival/growth (P) matrix. Eviction occurs when the column sums are lower than the survival models suggests that they should be.
# define the x-axis values
meshpts <- seq(from = (min(dat_all$log_LL_t, na.rm = TRUE) * .8), to = (max(dat_all$log_LL_t, na.rm = TRUE) * 1.2) , length.out = 500)
# plot the model-predicted survival probs.
preds <- predict(object = survMod_all, newdata = data.frame("log_LL_t" = meshpts), type = 'response')
## plot the model predictions
plot(x = meshpts, y = preds, ylab = "Survival Probability", type = "l")
# plot the survival values from the P matrix (column sums)
points(meshpts,apply(ipm_B$sub_kernels$P,2,sum),col="red",lwd=3,cex=.1,pch=19)
# %% problem... it appears that there is eviction occuring? is it possible to have eviction in the center, rather than the edges? 


meshpts <- seq(from = (min(dat_all$log_LL_t, na.rm = TRUE) * .8), to = (max(dat_all$log_LL_t, na.rm = TRUE) * 1.2) , length.out = 500)
# plot the model-predicted survival probs.
preds <- predict(object = survMod_all, newdata = data.frame("log_LL_t" = meshpts), type = 'response')
plot(x = meshpts, y = preds, ylab = "Survival Probability", type = "l")
# plot the survival values from the P matrix
points(meshpts,apply(ipm_A$sub_kernels$P,2,sum),col="red",lwd=3,cex=.1,pch=19)

### Calculate IPM B by hand ###
# DAVE LOOK HERE# 
### hand-calculated model for continuous-ized seedlings and discrete seedbank (IPM B) ###
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
# define the vital rate functions
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

# Setting up the kernels
K <- array(0,c(n+1,n+1))
b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
h=(U-L)/n # bin width 

# Survival and growth 
S <- diag(S.fun(meshp, paramCont)) # Survival # put survival probabilities in the diagonal of the matrix
G <- h * t(outer(meshp,meshp,GR.fun, paramCont)) # Growth
# G <- t(outer(meshp,meshp,GR.fun)) # Growth

#Recruits distribution (seeds recruited from the seedbank into the continuous stage)
c_o <- h * matrix(rep(SDS.fun(meshp, paramCont),n),n,n,byrow=F)
# c_o <- matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)

#Probability of flowering
Pb = (FL.fun(meshp, paramCont))

#Number of seeds produced according to adult size
b_seed = (SDP.fun(meshp, paramCont))

FecALL= Pb * b_seed

# update the 'S' matrix by multiplying it by (1-Pb), since this is a monocarpic perennial
S_new <- S * (1-(Pb))

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

Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))

mat_all_DI <-Pkernel+Fkernel

eigenMat <- eigen(mat_all_DI)
# get the lambda
lam_all_DI <-eigenMat$values[1]

## calculate the uncertainty using bootstrapping
# make a vector to hold all of the lambdas from the resampling runs
allDI_lambdas <- numeric(1000L)
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
for(i in 1:1000) {
  ## sample continuous data
  sample_ind <- seq(1, nrow(dat_all), by = 1)
  
  boot_data_ind   <- sample(x = sample_ind, size = length(sample_ind), replace = TRUE)
  
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
  names(paramCont_now) <- c("survival", "growth", "recruitDist", "flowering", "seedProduction")
  
  # Insert the new vital rate models into the proto_ipm, then rebuild the IPM.
  K_now <- array(0,c(n+1,n+1))
  
  # I recommend you set i = 1, set n low to say 10 
  
  # Setting up the kernels
  b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
  meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
  
  h=(U-L)/n # bin width 
  
  # Survival and growth 
  S <- diag(S.fun(meshp, paramCont_now)) # Survival # put survival probabilities in the diagonal of the matrix
  G <- h * t(outer(meshp,meshp,GR.fun, paramCont_now)) # Growth
  # G <- t(outer(meshp,meshp,GR.fun)) # Growth
  
  #Recruits distribution (seeds recruited from the seedbank into the continuous stage)
  c_o <- h * matrix(rep(SDS.fun(meshp, paramCont_now),n),n,n,byrow=F)
  # c_o <- matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
  
  #Probability of flowering
  Pb = (FL.fun(meshp, paramCont_now))
  
  #Number of seeds produced according to adult size
  b_seed = (SDP.fun(meshp, paramCont_now))
  
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
for (i in 1:1000) {
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
abline(v = lam_all_DI, col = 'red', lwd = 2, lty = 2)
abline(v = mean(allDI_lambdas), col = "blue", lwd = 2, lty = 2)
## calculate 95% CI
SE <- sd(allDI_lambdas)/sqrt(1000)
mean <- mean(allDI_lambdas)
allDI_CI <- c((mean - 1.96*SE),(mean + 1.96*SE))

#### store the ipm results
save.image(file = "./analysis_scripts/ipmA_B_results.RData")
