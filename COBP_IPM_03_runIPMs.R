#//////////////////////////
# Integral Projection Models for Oenothera coloradensis
# Alice Stears
# 30 November 2021
#/////////////////////////

#### Load packages ####
library(ipmr)

#### load vital rate models from previous script ####
source("./analysis_scripts/COBP_IPM_02_VitalRateModels.R")

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
    n_b = 400, 
    
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

#### store the ipm results
save.image(file = "./analysis_scripts/ipmA_B_results.RData")
