#//////////////////////////
# Integral Projection Models for Oenothera coloradensis
# Alice Stears
# 30 November 2021
#/////////////////////////

#### Load packages ####
library(ipmr)

#### load vital rate models from previous script ####
source("./analysis_scripts/COBP_IPM_02_VitalRateModels.R")

#### deterministic, density-independent IPM using only continuous stages ####
# calculate probability of establishment (skipping seedling and seedbank stages)
simple_seeds <- dat %>% 
  group_by(Plot_ID, Year) %>% 
  summarize("N_seeds_t" = sum(Num_seeds, na.rm = TRUE)) %>% 
  rename(Year_t = Year) %>% 
  mutate(Year_t = as.numeric(as.character(Year_t)))
simple_recruits <- dat %>% 
  group_by(Plot_ID, Year) %>% 
  summarize("N_recruits_t" = length(recruit)) %>% 
  rename(Year_t = Year) %>% 
  mutate(Year_tplus1 = as.numeric(as.character(Year_t)) - 1)

simple_estabs <- simple_seeds %>% 
  left_join(simple_recruits, 
            by = c("Year_t" = "Year_tplus1",
                   "Plot_ID" = "Plot_ID")) %>% 
  dplyr::select(- Year_t.y) %>% 
  rename("N_recruits_tplus1" = "N_recruits_t") %>% 
  mutate("p_estab" = N_recruits_tplus1/N_seeds_t)
  
simple_estabs[simple_estabs$p_estab == Inf & 
                is.na(simple_estabs$p_estab) == FALSE , 
              "p_estab"] <- NA

p.estab.simple = mean(simple_estabs$p_estab, na.rm = TRUE)

data_list <- list(
  g_int     = coef(sizeMod)[1],
  g_slope   = coef(sizeMod)[2],
  g_sd      = summary(sizeMod)$sigma,
  s_int     = coef(survMod)[1],
  s_slope   = coef(survMod)[2],
  p_b_int   = coef(flwrMod_t)[1], #probability of flowering
  p_b_slope = coef(flwrMod_t)[2],
  p_b_slope_2 = coef(flwrMod_t)[3],
  b_int   = coef(seedMod_t)[1], #seed production
  b_slope = coef(seedMod_t)[2],
  c_o_mu    = coef(recMod), #recruit size distribution
  c_o_sd    = summary(recMod)$sigma,
  p_estab = p.estab.simple 
  )

# inital population state
init_size_state <- runif(500)

simple_ipm <- init_ipm(sim_gen   = "simple", 
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
      min(dat$log_LL_t, na.rm = TRUE) * .8, # lower bound (L)
      max(dat$log_LL_t, na.rm = TRUE) * 1.2, # upper bound (U)
      500 # number of mesh points
    )
    ) %>% 
  define_pop_state(
    n_size = runif(500)
      ) %>% 
  make_ipm(
    iterations = 100
    )

## check for eviction from the model
#To check for eviction, we plot the survival model and the column sums of the survival/growth (P) matrix. Eviction occurs when the column sums are lower than the survival models suggests that they should be.
# define the x-axis values
meshpts <- seq(from = (min(dat$log_LL_t, na.rm = TRUE) * .8), to = (max(dat$log_LL_t, na.rm = TRUE) * 1.2) , length.out = 500)
# plot the model-predicted survival probs.
preds <- predict(object = survMod, newdata = data.frame("log_LL_t" = meshpts), type = 'response')
plot(x = meshpts, y = preds, ylab = "Survival Probability", type = "l")
# plot the survival values from the P matrix
points(meshpts,apply(simple_ipm$sub_kernels$P,2,sum),col="red",lwd=3,cex=.1,pch=19)

#### deterministic, density-dependent IPM using only continuous stages ####
p.estab.simple = mean(simple_estabs$p_estab, na.rm = TRUE)

data_list <- list(
  g_int     = coef(sizeMod)[1],
  g_slope   = coef(sizeMod)[2],
  g_sd      = summary(sizeMod)$sigma,
  s_int     = coef(survMod)[1],
  s_slope   = coef(survMod)[2],
  p_b_int   = coef(flwrMod_t)[1], #probability of flowering
  p_b_slope = coef(flwrMod_t)[2],
  p_b_slope_2 = coef(flwrMod_t)[3],
  b_int   = coef(seedMod_t)[1], #seed production
  b_slope = coef(seedMod_t)[2],
  c_o_mu    = coef(recMod), #recruit size distribution
  c_o_sd    = summary(recMod)$sigma,
  p_estab = p.estab.simple 
)

# inital population state
init_size_state <- pnorm(rnorm(500, mean = mean(ht, na.rm = TRUE), sd = sd(ht, na.rm = TRUE)), mean = mean(ht, na.rm = TRUE), sd = sd(ht, na.rm = TRUE))

simple_ipm <- init_ipm(sim_gen   = "simple", 
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
      min(dat$log_LL_t, na.rm = TRUE) * .8, # lower bound (L)
      max(dat$log_LL_t, na.rm = TRUE) * 1.2, # upper bound (U)
      500 # number of mesh points
    )
  ) %>% 
  define_pop_state(
    n_size = runif(500)
  ) %>% 
  make_ipm(
    iterations = 100
  )
#### Deterministic, density-independent IPM for all data ####
# vital-rate model names: survMod, sizeMod, seedMod_t, flwrMod_t, recMod, p.estab.est, outSB.est, staySB.est, goSB.est, goSdlng.est 

### Implement the IPM 
# use ipmr to fit the IPM
## Set up the initial population conditions and parameters (example w/ only one discrete stage and dummy seedbank rates)
data_list <- list(
  g_int     = coef(sizeMod)[1],
  g_slope   = coef(sizeMod)[2],
  g_sd      = summary(sizeMod)$sigma,
  s_int     = coef(survMod)[1],
  s_slope   = coef(survMod)[2],
  p_b_int   = coef(flwrMod_t)[1], #probability of flowering
  p_b_slope = coef(flwrMod_t)[2],
  p_b_slope_2 = coef(flwrMod_t)[3],
  b_int   = coef(seedMod_t)[1], #seed production
  b_slope = coef(seedMod_t)[2],
  c_o_mu    = coef(recMod), #recruit size distribution
  c_o_sd    = summary(recMod)$sigma,
  goSdlng   = goSdlng.est, # Probability that non-seedbank seeds will germinate into seedlings in year t+1
  staySB = staySB.est, # Probability that a seed in the seedbank in year t will exit the seedbank in year t+1 
  goSB = goSB.est, # probability that a seed produced by an adult plant in year t will enter the seedbank
  outSB = outSB.est, # probability that a seedbank seed will germinate to a seedling in year t+1
  p_estab = p.estab.est # probability that a seedling will establish into a rosette in t+1
)

## Next, we set up two functions to pass into the model. These perform the inverse logit transformations for the probability of flowering model (r_r/𝑟𝑟(𝑧)).
# We'll set up some helper functions. The survival function in this model is a quadratic function, so we use an additional inverse logit function that can handle the quadratic term.

inv_logit <- function(int, slope, sv) {
  1/(1 + exp(-(int + slope * sv)))
}

inv_logit_2 <- function(int, slope, slope_2, sv) {
  1/(1 + exp(-(int + slope * sv + slope_2 * sv ^ 2)))
}

## Now, we’re ready to begin making the IPM kernels. We change the sim_gen argument of init_ipm() to "general".
det_ipm <- init_ipm(sim_gen = "general", # make a general IPM
                    di_dd = "di", # make it density independent
                    det_stoch = "det") %>% # make it deterministic
  define_kernel(
    name          = "P", # survival 
    # We add d_ht to formula to make sure integration is handled correctly.
    # This variable is generated internally by make_ipm(), so we don't need
    # to do anything else.
    formula       = (1-p_b.) * s. * g. * d_ht,
    family        = "CC",
    g.             = dnorm(ht_2, g_mu, g_sd),
    g_mu          = g_int + g_slope * ht_1,
    s.             = inv_logit(s_int, s_slope, ht_1),
    p_b.          = inv_logit_2(p_b_int, p_b_slope, p_b_slope_2, ht_1),
    data_list     = data_list,
    states        = list(c('ht')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions('norm', 'g.')
  ) %>%
  define_kernel(
    name          = "leave_seedlings", ## leave seedling stage and go to rosette stage
    formula       = p_estab. * c_o. * d_ht,
    family        = 'DC', # Note that now, family = "DC" because it denotes a discrete -> continuous transition
    p_estab.      = p_estab,
    c_o.          = dnorm(ht_2, c_o_mu, c_o_sd),
    data_list     = data_list,
    states        = list(c('ht', "s")),   # Note that here, we add "s" to our list in states, because this kernel uses seedlings 
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions('norm','c_o.')
  ) %>%
  define_kernel(
    name    = "repro_to_seedlings",
    formula       = (goSdlng.) * (p_b. * b. * d_ht),
    family        = "CD",
    goSdlng.      = goSdlng,
    p_b.          = inv_logit_2(p_b_int, p_b_slope, p_b_slope_2, ht_1),
    b.            = exp(b_int + b_slope * ht_1),
    data_list     = data_list,
    states        = list(c('ht', 's')),
    uses_par_sets = FALSE,
    evict_cor     = FALSE
  ) %>%
  define_kernel(
    name          = 'seedbank_to_seedlings',
    formula       = outSB.,
    family        = 'DD',
    outSB.        = outSB,
    data_list     = data_list,
    states        = list(c('b', 's')),
    uses_par_sets = FALSE,
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name    = "stay_seedbank",
    formula       = staySB.,
    family        = "DD",
    staySB.        = staySB,
    data_list     = data_list,
    states        = list(c('b')),
    uses_par_sets = FALSE,
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name          = 'repro_to_seedbank',
    formula       = (goSB.) * (p_b. * b. * d_ht),
    family        = 'CD',
    goSB.          = goSB, 
    p_b.          = inv_logit_2(p_b_int, p_b_slope, p_b_slope_2, ht_1),
    b.            = exp(b_int + b_slope * ht_1),
    data_list     = data_list,
    states        = list(c('b', 'ht')),
    uses_par_sets = FALSE,
    evict_cor = FALSE
  ) %>% ## define the starting and ending states for each kernel
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "leave_seedlings", "repro_to_seedlings", "seedbank_to_seedlings", "stay_seedbank", "repro_to_seedbank"),
      int_rule     = c(rep("midpoint", 6)),
      state_start    = c('ht', "s", "ht", "b", "b", "ht"),
      state_end      = c("ht", "ht", "s", "s", "b", "b")
    )
  ) %>%
  # actally run the IPM
  define_domains(
    # We can pass the variables we created above into define_domains
    ht = c(L, U, n)
  ) %>%
  define_pop_state(
    # We can also pass them into define_pop_state
    pop_vectors = list(
      n_ht = init_pop_vec,
      n_b  = init_seed_bank,
      n_s  = init_seedlings 
    )
  ) %>%
  make_ipm(iterations = 100,
           normalize_pop_size = FALSE,
           usr_funs = list(inv_logit   = inv_logit,
                           inv_logit_2 = inv_logit_2), return_main_env = TRUE )

## lambda is a generic function to compute per-capita growth rates. It has a number of different options depending on the type of model
lambda(det_ipm)

## If we are worried about whether or not the model converged to stable dynamics, we can use the exported utility is_conv_to_asymptotic. The default tolerance for convergence is 1e-10, but can be changed with the 'tol' argument.
is_conv_to_asymptotic(det_ipm, tol = 1e-10)
## additional calculations
lambda_ipmr <- lambda(det_ipm)
repro_value <- left_ev(det_ipm)
stable_dist <- right_ev(det_ipm)

## check for eviction from the model
#To check for eviction, we plot the survival model and the column sums of the survival/growth (P) matrix. Eviction occurs when the column sums are lower than the survival models suggests that they should be.
# define the x-axis values
meshpts <- seq(from = (min(dat$log_LL_t, na.rm = TRUE) * .8), to = (max(dat$log_LL_t, na.rm = TRUE) * 1.2) , length.out = 500)
# plot the model-predicted survival probs.
preds <- predict(object = survMod, newdata = data.frame("log_LL_t" = meshpts), type = 'response')
plot(x = meshpts, y = preds, ylab = "Survival Probability", type = "l")
# plot the survival values from the P matrix
points(meshpts,apply(det_ipm$sub_kernels$P,2,sum),col="red",lwd=3,cex=.1,pch=19)

### Visualize the IPM kernel
# first have to make a mega-kernel
mega_mat <- make_iter_kernel(ipm = det_ipm, 
                             mega_mat = c(stay_seedbank, 0, repro_to_seedbank, seedbank_to_seedlings, 0, repro_to_seedlings, 0, leave_seedlings, P))
# check to make sure I constructed the mega-kernel correctly (should be nearly equal values)
Re(eigen(mega_mat[[1]])$values[1]) - lambda(det_ipm)

## visualize the full kernel 
# define the meshpoints 
meshpts <- seq(from = L, to = U, length = 500)
# set up the paletteot the continuous part of the kernel (leave out first two rows and cols that correspond to the discrete stages)
## make the entire figure as a lattice plot
graphics::layout(mat = matrix(c(1,4, 7, 2, 5, 8, 3, 6, 9),nrow = 3, ncol = 3), widths = c(1.5,1.5,6),heights = c(6,1.5,1.5))
## B(t+1)
par(mar = c(3,3,3,1))
image(t(mega_mat$mega_matrix[1,3:502]^.1), xaxt = "n", yaxt = "n",
      main = "B(t+1)",
      col = pal[(round(min((mega_mat$mega_matrix[1,3:502])^.1),2)*100): 
                  (round(max((mega_mat$mega_matrix[1,3:502])^.1),2)*100)]) 
mtext(side = 1, text  = c("continous to \nseedbank"), line = -1, cex = .75)
## S(t+1)
par(mar = c(3,1,3,1))
image(t(mega_mat$mega_matrix[2,3:502]^.1), xaxt = "n", yaxt = "n",
      main = "S(t+1)",
      col = pal[(round(min((mega_mat$mega_matrix[2,3:502])^.1),2)*100): 
                  (round(max((mega_mat$mega_matrix[2,3:502])^.1),2)*100)]) 
mtext(side = 1, text  = c("continous to \nseedlings"), line = -1, cex = .75)
## K matrix
par(mar = c(3,3,3,3)
    ,mgp = c(1.75,.5,0)
)
image(x = meshpts, y = meshpts, t(mega_mat$mega_matrix[3:502,3:502])^.1,
      xlab = "n(z,t); log(cm)", ylab = "n(z',t+1); log(cm)", 
      main = "Continuous Stage (t+1)",
      col = pal[(round(min( t(mega_mat$mega_matrix[3:502,3:502])^.1),2)*100): (round(max( t(mega_mat$mega_matrix[3:502,3:502])^.1),2)*100)]
) ## get the correct values for the color ramp that correspond to the actual probabilities in the entire matrix
text(x = 4.8, y = 2.25, c("Continuous Stage (t)"), xpd = NA, srt = -90, cex = 1.25, font = 2)
abline(a = 0, b = 1, lty = 2)
contour(x = meshpts, y = meshpts, 
        t(mega_mat$mega_matrix[3:502,3:502]), 
        add = TRUE, drawlabels = TRUE, nlevels = 10, col = "grey30")
## seedlings to seedbank
par(mar = c(1,3,1,1))
image(as.matrix(0), xaxt = "n", yaxt = "n", col = "white") 
mtext(side = 1, text  = c("can't go from \nseedling \nto seedbank"), line = -1, cex = .75)
## seedlings to seedlings
par(mar = c(1,1,1,1))
image(as.matrix(0), xaxt = "n", yaxt = "n",  col = "white") 
mtext(side = 1, text  = c("can't stay in \nseedlings"), line = -1, cex = .75)
## S(t)
par(mar = c(1,3,1,3))
image(as.matrix(mega_mat$mega_matrix[3:502,2]^.1), yaxt = "n", xaxt = "n",
      col = pal[(round(min( t(mega_mat$mega_matrix[3:502,2])^.1),2)*100): 
                  (round(max( t(mega_mat$mega_matrix[3:502,2])^.1),2)*100)]) 
text(x = 1.05,y = .5, c("S(t)"), xpd = NA, srt = -90, cex = 1.25, font = 2)
mtext(side = 1, text  = c("seedling to continuous stage"), line = -1, cex = .75)
## plot the staySB probability
par(mar = c(3,3,1,1))
image(as.matrix(mega_mat$mega_matrix[1,1]^.1), xaxt = "n", yaxt = "n",
      col = pal[round(max(mega_mat$mega_matrix[1,1]^.1),2)*100])
mtext(side = 1, text  = c("stay in \nseedbank"), line = -1, cex = .75)
## plot the seedbank to seedlings probability
par(mar = c(3,1,1,1))
image(as.matrix(t(mega_mat$mega_matrix[2,1])^.1), xaxt = "n", yaxt = "n",
      col = pal[round(max(t(mega_mat$mega_matrix[2,1])^.1),2)*100])
mtext(side = 1, text = c("seedbank to \nseedlings"), line = -1, cex = .75)
## B(t)(is all zeros--can't transition to continuous stage from the seedbank)
par(mar = c(3,3,1,3))
image(t(mega_mat$mega_matrix[3:502,1]^.1), yaxt = "n", xaxt = "n",
      col = pal[(round(min( t(mega_mat$mega_matrix[3:502,1])^.1),2)*100): 
                  (round(max( t(mega_mat$mega_matrix[3:502,1])^.1),2)*100)]) 
text(x = 1.1,y = .5, c("B(t)"), xpd = NA, srt = -90, cex = 1.25, font = 2)
mtext(side = 1, text = c("can't go from seedbank to continuous stage"), line = -1, cex = .75)
dev.off()

#### determinisitic, density dependent IPM for all data -- FIRST YEAR ONLY ####
# use ipmr to fit the IPM
## Set up the initial population conditions and parameters (example w/ only one discrete stage and dummy seedbank rates)
data_list <- list(
  g_int     = coef(sizeMod_first)[1],
  g_slope   = coef(sizeMod_first)[2],
  g_sd      = summary(sizeMod_first)$sigma,
  s_int     = coef(survMod_first)[1],
  s_slope   = coef(survMod_first)[2],
  p_b_int   = coef(flwrMod_first)[1], #probability of flowering
  p_b_slope = coef(flwrMod_first)[2],
  p_b_slope_2 = coef(flwrMod_first)[3],
  b_int   = coef(seedMod_first)[1], #seed production
  b_slope = coef(seedMod_first)[2],
  c_o_mu    = coef(recMod_first), #recruit size distribution
  c_o_sd    = summary(recMod_first)$sigma,
  goSdlng   = goSdlng.est, # Probability that non-seedbank seeds will germinate into seedlings in year t+1
  staySB = staySB.est, # Probability that a seed in the seedbank in year t will exit the seedbank in year t+1 
  goSB = goSB.est, # probability that a seed produced by an adult plant in year t will enter the seedbank
  outSB = outSB.est, # probability that a seedbank seed will germinate to a seedling in year t+1
  p_estab = p.estab.est # probability that a seedling will establish into a rosette in t+1
)

## Now, we’re ready to begin making the IPM kernels. We change the sim_gen argument of init_ipm() to "general".
det_ipm_first <- init_ipm(sim_gen = "general", # make a general IPM
                          di_dd = "di", # make it density independent
                          det_stoch = "det") %>% # make it deterministic
  define_kernel(
    name          = "P", # survival 
    # We add d_ht to formula to make sure integration is handled correctly.
    # This variable is generated internally by make_ipm(), so we don't need
    # to do anything else.
    formula       = (1-p_b.) * s. * g. * d_ht,
    family        = "CC",
    g.             = dnorm(ht_2, g_mu, g_sd),
    g_mu          = g_int + g_slope * ht_1,
    s.             = inv_logit(s_int, s_slope, ht_1),
    p_b.          = inv_logit_2(p_b_int, p_b_slope, p_b_slope_2, ht_1),
    data_list     = data_list,
    states        = list(c('ht')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions('norm', 'g.')
  ) %>%
  define_kernel(
    name          = "leave_seedlings", ## leave seedling stage and go to rosette stage
    formula       = p_estab. * c_o. * d_ht,
    family        = 'DC', # Note that now, family = "DC" because it denotes a discrete -> continuous transition
    p_estab.      = p_estab,
    c_o.          = dnorm(ht_2, c_o_mu, c_o_sd),
    data_list     = data_list,
    states        = list(c('ht', "s")),   # Note that here, we add "s" to our list in states, because this kernel uses seedlings 
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions('norm','c_o.')
  ) %>%
  define_kernel(
    name    = "repro_to_seedlings",
    formula       = (goSdlng.) * (p_b. * b. * d_ht),
    family        = "CD",
    goSdlng.      = goSdlng,
    p_b.          = inv_logit_2(p_b_int, p_b_slope, p_b_slope_2, ht_1),
    b.            = exp(b_int + b_slope * ht_1),
    data_list     = data_list,
    states        = list(c('ht', 's')),
    uses_par_sets = FALSE,
    evict_cor     = FALSE
  ) %>%
  define_kernel(
    name          = 'seedbank_to_seedlings',
    formula       = outSB.,
    family        = 'DD',
    outSB.        = outSB,
    data_list     = data_list,
    states        = list(c('b', 's')),
    uses_par_sets = FALSE,
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name    = "stay_seedbank",
    formula       = staySB.,
    family        = "DD",
    staySB.        = staySB,
    data_list     = data_list,
    states        = list(c('b')),
    uses_par_sets = FALSE,
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name          = 'repro_to_seedbank',
    formula       = (goSB.) * (p_b. * b. * d_ht),
    family        = 'CD',
    goSB.          = goSB, 
    p_b.          = inv_logit_2(p_b_int, p_b_slope, p_b_slope_2, ht_1),
    b.            = exp(b_int + b_slope * ht_1),
    data_list     = data_list,
    states        = list(c('b', 'ht')),
    uses_par_sets = FALSE,
    evict_cor = FALSE
  )%>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "leave_seedlings", "repro_to_seedlings", "seedbank_to_seedlings", "stay_seedbank", "repro_to_seedbank"),
      int_rule     = c(rep("midpoint", 6)),
      state_start    = c('ht', "s", "ht", "b", "b", "ht"),
      state_end      = c("ht", "ht", "s", "s", "b", "b")
    )
  )  %>%
  define_domains(
    # We can pass the variables we created above into define_domains
    ht = c(L, U, n)
  ) %>%
  define_pop_state(
    # We can also pass them into define_pop_state
    pop_vectors = list(
      n_ht = init_pop_vec,
      n_b  = init_seed_bank,
      n_s  = init_seedlings 
    )
  ) %>%
  make_ipm(iterations = 100,
           normalize_pop_size = FALSE,
           usr_funs = list(inv_logit   = inv_logit,
                           inv_logit_2 = inv_logit_2), return_main_env = TRUE )

## lambda is a generic function to compute per-capita growth rates. It has a number of different options depending on the type of model
lambda(det_ipm_first)

## If we are worried about whether or not the model converged to stable dynamics, we can use the exported utility is_conv_to_asymptotic. The default tolerance for convergence is 1e-10, but can be changed with the 'tol' argument.
is_conv_to_asymptotic(det_ipm_first, tol = 1e-10)

#### determinisitic, density dependent IPM for all data -- SECOND YEAR ONLY ####
# use ipmr to fit the IPM
## Set up the initial population conditions and parameters (example w/ only one discrete stage and dummy seedbank rates)
data_list <- list(
  g_int     = coef(sizeMod_second)[1],
  g_slope   = coef(sizeMod_second)[2],
  g_sd      = summary(sizeMod_second)$sigma,
  s_int     = coef(survMod_second)[1],
  s_slope   = coef(survMod_second)[2],
  p_b_int   = coef(flwrMod_second)[1], #probability of flowering
  p_b_slope = coef(flwrMod_second)[2],
  p_b_slope_2 = coef(flwrMod_second)[3],
  b_int   = coef(seedMod_second)[1], #seed production
  b_slope = coef(seedMod_second)[2],
  c_o_mu    = coef(recMod_second), #recruit size distribution
  c_o_sd    = summary(recMod_second)$sigma,
  goSdlng   = goSdlng.est, # Probability that non-seedbank seeds will germinate into seedlings in year t+1
  staySB = staySB.est, # Probability that a seed in the seedbank in year t will exit the seedbank in year t+1 
  goSB = goSB.est, # probability that a seed produced by an adult plant in year t will enter the seedbank
  outSB = outSB.est, # probability that a seedbank seed will germinate to a seedling in year t+1
  p_estab = p.estab.est # probability that a seedling will establish into a rosette in t+1
)

## Now, we’re ready to begin making the IPM kernels. We change the sim_gen argument of init_ipm() to "general".
det_ipm_second <- init_ipm(sim_gen = "general", # make a general IPM
                           di_dd = "di", # make it density independent
                           det_stoch = "det") %>% # make it deterministic
  define_kernel(
    name          = "P", # survival 
    # We add d_ht to formula to make sure integration is handled correctly.
    # This variable is generated internally by make_ipm(), so we don't need
    # to do anything else.
    formula       = (1-p_b.) * s. * g. * d_ht,
    family        = "CC",
    g.             = dnorm(ht_2, g_mu, g_sd),
    g_mu          = g_int + s_slope, ht_1,
    p_b.          = inv_logit_2(p_b_int, p_b_slope, p_b_slope_2, ht_1),
    data_list     = data_list,
    states        = list(c('ht')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions('norm', 'g.')
  ) %>%
  define_kernel(
    name          = "leave_seedlings", ## leave seedling stage and go to rosette stage
    formula       = p_estab. * c_o. * d_ht,
    family        = 'DC', # Note that now, family = "DC" because it denotes a discrete -> continuous transition
    p_estab.      = p_estab,
    c_o.          = dnorm(ht_2, c_o_mu, c_o_sd),
    data_list     = data_list,
    states        = list(c('ht', "s")),   # Note that here, we add "s" to our list in states, because this kernel uses seedlings 
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions('norm','c_o.')
  ) %>%
  define_kernel(
    name    = "repro_to_seedlings",
    formula       = (goSdlng.) * (p_b. * b. * d_ht),
    family        = "CD",
    goSdlng.      = goSdlng,
    p_b.          = inv_logit_2(p_b_int, p_b_slope, p_b_slope_2, ht_1),
    b.            = exp(b_int + b_slope * ht_1),
    data_list     = data_list,
    states        = list(c('ht', 's')),
    uses_par_sets = FALSE,
    evict_cor     = FALSE
  ) %>%
  define_kernel(
    name          = 'seedbank_to_seedlings',
    formula       = outSB.,
    family        = 'DD',
    outSB.        = outSB,
    data_list     = data_list,
    states        = list(c('b', 's')),
    uses_par_sets = FALSE,
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name    = "stay_seedbank",
    formula       = staySB.,
    family        = "DD",
    staySB.        = staySB,
    data_list     = data_list,
    states        = list(c('b')),
    uses_par_sets = FALSE,
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name          = 'repro_to_seedbank',
    formula       = (goSB.) * (p_b. * b. * d_ht),
    family        = 'CD',
    goSB.          = goSB, 
    p_b.          = inv_logit_2(p_b_int, p_b_slope, p_b_slope_2, ht_1),
    b.            = exp(b_int + b_slope * ht_1),
    data_list     = data_list,
    states        = list(c('b', 'ht')),
    uses_par_sets = FALSE,
    evict_cor = FALSE
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "leave_seedlings", "repro_to_seedlings", "seedbank_to_seedlings", "stay_seedbank", "repro_to_seedbank"),
      int_rule     = c(rep("midpoint", 6)),
      state_start    = c('ht', "s", "ht", "b", "b", "ht"),
      state_end      = c("ht", "ht", "s", "s", "b", "b")
    )
  ) %>%
  define_domains(
    # We can pass the variables we created above into define_domains
    ht = c(L, U, n)
  ) %>%
  define_pop_state(
    # We can also pass them into define_pop_state
    pop_vectors = list(
      n_ht = init_pop_vec,
      n_b  = init_seed_bank,
      n_s  = init_seedlings 
    )
  ) %>%
  make_ipm(iterations = 100,
           normalize_pop_size = FALSE,
           usr_funs = list(inv_logit   = inv_logit,
                           inv_logit_2 = inv_logit_2), return_main_env = TRUE )

## lambda is a generic function to compute per-capita growth rates. It has a number of different options depending on the type of model
lambda(det_ipm_second)

## If we are worried about whether or not the model converged to stable dynamics, we can use the exported utility is_conv_to_asymptotic. The default tolerance for convergence is 1e-10, but can be changed with the 'tol' argument.
is_conv_to_asymptotic(det_ipm_second, tol = 1e-10)

#### deterministic, density dependent IPM for all data (no environmental covariates) ####
### Implement the IPM 
# use ipmr to fit the IPM
data_list <- list(
  g_int     = coef(sizeMod_dd)[1],
  g_slope   = coef(sizeMod_dd)[2],
  g_dd      = coef(sizeMod_dd)[3],
  g_sd      = summary(sizeMod_dd)$sigma,
  s_int     = coef(survMod_dd)[1],
  s_slope   = coef(survMod_dd)[2],
  s_dd      = coef(survMod_dd)[3],
  p_b_int   = coef(flwrMod_dd)[1], #probability of flowering
  p_b_slope = coef(flwrMod_dd)[2],
  p_b_slope_2 = coef(flwrMod_dd)[3],
  p_b_dd    = coef(flwrMod_dd)[4],
  b_int   = coef(seedMod_dd)[1], #seed production
  b_slope = coef(seedMod_dd)[2],
  b_dd     = coef(seedMod_dd)[3],
  c_o_int    = coef(recMod_dd)[1], #recruit size distribution
  c_o_dd      = coef(recMod_dd)[2],
  c_o_sd    = summary(recMod_dd)$sigma,
  goSdlng   = goSdlng.est, # Probability that non-seedbank seeds will germinate into seedlings in year t+1
  staySB = staySB.est, # Probability that a seed in the seedbank in year t will exit the seedbank in year t+1 
  goSB = goSB.est, # probability that a seed produced by an adult plant in year t will enter the seedbank
  outSB = outSB.est, # probability that a seedbank seed will germinate to a seedling in year t+1
  p_estab = p.estab.est # probability that a seedling will establish into a rosette in t+1
)

## Next, we set up two functions to pass into the model. These perform the inverse logit transformations for the probability of flowering model (r_r/𝑟𝑟(𝑧)).
# We'll set up some helper functions. The survival function in this model is a quadratic function, so we use an additional inverse logit function that can handle the quadratic term.

inv_logit <- function(lin.pred) {
  1/(1 + exp(-(lin.pred)))
}


## Now, we’re ready to begin making the IPM kernels. We change the sim_gen argument of init_ipm() to "general".
det_dd_ipm <- init_ipm(sim_gen = "general", # make a general IPM
                       di_dd = "dd", # make it density independent
                       det_stoch = "det") %>% # make it deterministic
  define_kernel(
    name          = "P", # survival 
    # We add d_ht to formula to make sure integration is handled correctly.
    # This variable is generated internally by make_ipm(), so we don't need
    # to do anything else.
    formula       = (1-p_b.) * s. * g. * d_ht,
    family        = "CC",
    g.             = dnorm(ht_2, g_mu., g_sd),
    g_mu.          = g_int + g_slope * ht_1 + g_dd * sum(n_ht_t),
    s.             = inv_logit(s_int + s_slope* ht_1 + s_dd * sum(n_ht_t)),
    p_b.          = inv_logit(p_b_int + p_b_slope* ht_1 + p_b_slope_2 * (ht_1^2) + p_b_dd * sum(n_ht_t)),
    data_list     = data_list,
    states        = list(c('ht')),
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions('norm', 'g.')
  ) %>%
  define_kernel(
    name          = "leave_seedlings", ## leave seedling stage and go to rosette stage
    formula       = p_estab. * c_o. * d_ht,
    family        = 'DC', # Note that now, family = "DC" because it denotes a discrete -> continuous transition
    p_estab.      = p_estab,
    c_o.          = dnorm(ht_2, c_o_mu., c_o_sd),
    c_o_mu.       = c_o_int + c_o_dd * sum(n_ht_t),
    data_list     = data_list,
    states        = list(c('ht', "s")),   # Note that here, we add "s" to our list in states, because this kernel uses seedlings 
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions('norm','c_o.')
  ) %>%
  define_kernel(
    name    = "repro_to_seedlings",
    formula       = (goSdlng.) * (p_b. * b. * d_ht),
    family        = "CD",
    goSdlng.      = goSdlng,
    p_b.          = inv_logit(p_b_int + p_b_slope * ht_1 + p_b_slope_2* (ht_1^2)  + p_b_dd * sum(n_ht_t)),
    b.            = exp(b_int + b_slope * ht_1 + b_dd * sum(n_ht_t)),
    data_list     = data_list,
    states        = list(c('ht', 's')),
    uses_par_sets = FALSE,
    evict_cor     = FALSE
  ) %>%
  define_kernel(
    name          = 'seedbank_to_seedlings',
    formula       = outSB.,
    family        = 'DD',
    outSB.        = outSB,
    data_list     = data_list,
    states        = list(c('b', 's')),
    uses_par_sets = FALSE,
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name    = "stay_seedbank",
    formula       = staySB.,
    family        = "DD",
    staySB.        = staySB,
    data_list     = data_list,
    states        = list(c('b')),
    uses_par_sets = FALSE,
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name          = 'repro_to_seedbank',
    formula       = (goSB.) * (p_b. * b. * d_ht),
    family        = 'CD',
    goSB.          = goSB, 
    p_b.          = inv_logit(p_b_int + p_b_slope * ht_1 + p_b_slope_2* (ht_1^2)  + p_b_dd * sum(n_ht_t)),
    b.            = exp(b_int + b_slope * ht_1 + b_dd * sum(n_ht_t)),
    data_list     = data_list,
    states        = list(c('b', 'ht')),
    uses_par_sets = FALSE,
    evict_cor = FALSE
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "leave_seedlings", "repro_to_seedlings", "seedbank_to_seedlings", "stay_seedbank", "repro_to_seedbank"),
      int_rule     = c(rep("midpoint", 6)),
      state_start    = c('ht', "s", "ht", "b", "b", "ht"),
      state_end      = c("ht", "ht", "s", "s", "b", "b")
    )
  ) %>%
  define_domains(
    # We can pass the variables we created above into define_domains
    ht = c(L, U, n)
  ) %>%
  define_pop_state(
    # We can also pass them into define_pop_state
    pop_vectors = list(
      n_ht = init_pop_vec,
      n_b  = init_seed_bank,
      n_s  = init_seedlings 
    )
  ) %>%
  make_ipm(iterations = 100,
           usr_funs = list(inv_logit   = inv_logit), return_main_env = TRUE )

## lambda is a generic function to compute per-capita growth rates. It has a number of different options depending on the type of model
lambda(det_dd_ipm)

## If we are worried about whether or not the model converged to stable dynamics, we can use the exported utility is_conv_to_asymptotic. The default tolerance for convergence is 1e-10, but can be changed with the 'tol' argument.
is_conv_to_asymptotic(det_dd_ipm, tol = 1e-10)
## additional calculations
lambda_ipmr <- lambda(det_ipm)
repro_value <- left_ev(det_ipm)
stable_dist <- right_ev(det_ipm)

### Visualize the IPM kernel
# first have to make a mega-kernel
mega_mat <- make_iter_kernel(ipm = det_ipm, 
                             mega_mat = c(stay_seedbank, 0, repro_to_seedbank, seedbank_to_seedlings, 0, repro_to_seedlings, 0, leave_seedlings, P),
)
# check to make sure I constructed the mega-kernel correctly (should be nearly equal values)
Re(eigen(mega_mat[[1]])$values[1]) - lambda(det_ipm)

## visualize the full kernel 
# define the meshpoints 
meshpts <- seq(from = L, to = U, length = 500)
# set up the palette
pal <- hcl.colors(n = 100, palette = "Heat 2", rev = TRUE)
# plot the continuous part of the kernel (leave out first two rows and cols that correspond to the discrete stages)
## make the entire figure as a lattice plot
graphics::layout(mat = matrix(c(1,4, 7, 2, 5, 8, 3, 6, 9),nrow = 3, ncol = 3), widths = c(1.5,1.5,6),heights = c(6,1.5,1.5))
## B(t+1)
par(mar = c(3,3,3,1))
image(t(mega_mat$mega_matrix[1,3:502]^.1), xaxt = "n", yaxt = "n",
      main = "B(t+1)",
      col = pal[(round(min((mega_mat$mega_matrix[1,3:502])^.1),2)*100): 
                  (round(max((mega_mat$mega_matrix[1,3:502])^.1),2)*100)]) 
mtext(side = 1, text  = c("continous to \nseedbank"), line = -1, cex = .75)
## S(t+1)
par(mar = c(3,1,3,1))
image(t(mega_mat$mega_matrix[2,3:502]^.1), xaxt = "n", yaxt = "n",
      main = "S(t+1)",
      col = pal[(round(min((mega_mat$mega_matrix[2,3:502])^.1),2)*100): 
                  (round(max((mega_mat$mega_matrix[2,3:502])^.1),2)*100)]) 
mtext(side = 1, text  = c("continous to \nseedlings"), line = -1, cex = .75)
## K matrix
par(mar = c(3,3,3,3)
    ,mgp = c(1.75,.5,0)
)
image(x = meshpts, y = meshpts, t(mega_mat$mega_matrix[3:502,3:502])^.1,
      xlab = "n(z,t); log(cm)", ylab = "n(z',t+1); log(cm)", 
      main = "Continuous Stage (t+1)",
      col = pal[(round(min( t(mega_mat$mega_matrix[3:502,3:502])^.1),2)*100): (round(max( t(mega_mat$mega_matrix[3:502,3:502])^.1),2)*100)]
) ## get the correct values for the color ramp that correspond to the actual probabilities in the entire matrix
text(x = 3.8, y = 2.25, c("Continuous Stage (t)"), xpd = NA, srt = -90, cex = 1.25, font = 2)
abline(a = 0, b = 1, lty = 2)
contour(x = meshpts, y = meshpts, 
        t(mega_mat$mega_matrix[3:502,3:502]), 
        add = TRUE, drawlabels = TRUE, nlevels = 10, col = "grey30")
## seedlings to seedbank
par(mar = c(1,3,1,1))
image(as.matrix(0), xaxt = "n", yaxt = "n", col = "white") 
mtext(side = 1, text  = c("can't go from \nseedling \nto seedbank"), line = -1, cex = .75)
## seedlings to seedlings
par(mar = c(1,1,1,1))
image(as.matrix(0), xaxt = "n", yaxt = "n",  col = "white") 
mtext(side = 1, text  = c("can't stay in \nseedlings"), line = -1, cex = .75)
## S(t)
par(mar = c(1,3,1,3))
image(as.matrix(mega_mat$mega_matrix[3:502,2]^.1), yaxt = "n", xaxt = "n",
      col = pal[(round(min( t(mega_mat$mega_matrix[3:502,2])^.1),2)*100): 
                  (round(max( t(mega_mat$mega_matrix[3:502,2])^.1),2)*100)]) 
text(x = 1.05,y = .5, c("S(t)"), xpd = NA, srt = -90, cex = 1.25, font = 2)
mtext(side = 1, text  = c("seedling to continuous stage"), line = -1, cex = .75)
## plot the staySB probability
par(mar = c(3,3,1,1))
image(as.matrix(mega_mat$mega_matrix[1,1]^.1), xaxt = "n", yaxt = "n",
      col = pal[round(max(mega_mat$mega_matrix[1,1]^.1),2)*100])
mtext(side = 1, text  = c("stay in \nseedbank"), line = -1, cex = .75)
## plot the seedbank to seedlings probability
par(mar = c(3,1,1,1))
image(as.matrix(t(mega_mat$mega_matrix[2,1])^.1), xaxt = "n", yaxt = "n",
      col = pal[round(max(t(mega_mat$mega_matrix[2,1])^.1),2)*100])
mtext(side = 1, text = c("seedbank to \nseedlings"), line = -1, cex = .75)
## B(t)(is all zeros--can't transition to continuous stage from the seedbank)
par(mar = c(3,3,1,3))
image(t(mega_mat$mega_matrix[3:502,1]^.1), yaxt = "n", xaxt = "n",
      col = pal[(round(min( t(mega_mat$mega_matrix[3:502,1])^.1),2)*100): 
                  (round(max( t(mega_mat$mega_matrix[3:502,1])^.1),2)*100)]) 
text(x = 1.1,y = .5, c("B(t)"), xpd = NA, srt = -90, cex = 1.25, font = 2)
mtext(side = 1, text = c("can't go from seedbank to continuous stage"), line = -1, cex = .75)


## visualize the P matrix
dev.off()
image(x = meshpts, 
      y = meshpts,
      t(make_iter_kernel(ipm = det_ipm, mega_mat = c(P))$mega_matrix),
      xlab = "n(z,t)", ylab = "n(z',t+1)", main = "P matrix")
contour(x = meshpts,
        y = meshpts,
        t(make_iter_kernel(ipm = det_ipm, mega_mat = c(P))$mega_matrix), 
        add = TRUE, drawlabels = TRUE, nlevels = 10, col = "grey30")

#### stochastic, density-independent IPM for all data with random effect for site and environmental covariates ####
### Vital Rate model names: survMod_env, sizeMod_env, seedMod_env, flwrMod_env, recMod_env, p.estab.est, outSB.est, staySB.est, goSB.est, goSdlng.est

### define the parameter list
fixed_list <- list(
  g_int     = fixef(sizeMod_env)[1], # growth model intercept
  g_slope   = fixef(sizeMod_env)[2], # growth model leaf size slope
  g_soilM   = fixef(sizeMod_env)[3], # growth model soil moisture coefficient
  g_temp    = fixef(sizeMod_env)[4], # growth model mean growing season temp coefficient
  g_precip  = fixef(sizeMod_env)[5], # growth model precip coefficient
  g_sd      = sd(resid(sizeMod_env)), # growth model sd
  s_int     = fixef(survMod_env)[1], # survival model intercept
  s_slope   = fixef(survMod_env)[2], # survival model leaf size slope 
  s_soilM   = fixef(survMod_env)[3], # survival model soil moisture coefficient
  s_soilT   = fixef(survMod_env)[4], # survival model soil temp coefficient
  s_temp    = fixef(survMod_env)[5], # survival model mean temp coefficient
  p_b_int   = fixef(flwrMod_env)[1], # probability of flowering
  p_b_slope = fixef(flwrMod_env)[2], # flowering model leaf size slope
  p_b_slope_2 = fixef(flwrMod_env)[3], # flowering model leaf size^2 slope
  p_b_soilM = fixef(flwrMod_env)[4], # flowering model soil moisture coefficient
  p_b_soilT = fixef(flwrMod_env)[5], # flowering model soil temp coefficient
  p_b_temp  = fixef(flwrMod_env)[6], # flowering model mean temp coefficient
  p_b_precip = fixef(flwrMod_env)[7], # flowering model precip coefficient
  b_int     = fixef(seedMod_env)[1], # seed production
  b_slope   = fixef(seedMod_env)[2], # seed model leaf size slope
  b_temp    = fixef(seedMod_env)[3], # seed model mean temp coefficient
  b_precip  = fixef(seedMod_env)[4], # seed model precip coefficient
  c_o_int   = fixef(recMod_env)[1], # recruit size model intercept
  c_o_soilM = fixef(recMod_env)[2], # recruit size model soil moisture coefficient
  c_o_sd    = sd(resid(recMod_env)), # recruit size model sd
  goSdlng   = goSdlng.est, # Probability that non-seedbank seeds will germinate into seedlings in year t+1
  staySB = staySB.est, # Probability that a seed in the seedbank in year t will exit the seedbank in year t+1
  goSB = goSB.est, # probability that a seed produced by an adult plant in year t will enter the seedbank
  outSB = outSB.est, # probability that a seedbank seed will germinate to a seedling in year t+1
  p_estab = p.estab.est # probability that a seedling will establish into a rosette in t+1
)

## now get a list of values for the random effects
# random values for the growth model
g_r_int <- unlist(ranef(sizeMod_env))
# random values for the survival model
s_r_int <- unlist(ranef(survMod_env))
# random values for the flowering model
p_b_r_int <- unlist(ranef(flwrMod_env))
# random values for the seed production model
b_r_int <- unlist(ranef(seedMod_env))
# random values for the recruit size model
c_o_r_int <- unlist(ranef(recMod_env))

# name the list of random effects
nms <- paste("r_", 1:18, sep = "")

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
all_params_list <- c(fixed_list, g_params, s_params, p_b_params, b_params, c_o_params)

### make environmental parameters
env_params <- list(
  soilM_mu = mean(dat$SoilMoisture_m3m3_s), # use a normal dist for soilM (?)
  soilM_sd = sd(dat$SoilMoisture_m3m3_s),
  soilT_mu = mean(dat$SoilTemp_grow_C_s, na.rm = TRUE), # use a norm dist for soilT (?)
  soilT_sd = sd(dat$SoilTemp_grow_C_s, na.rm = TRUE),
  temp_mu = mean(dat$tMean_grow_C_s),
  temp_sd = sd(dat$tMean_grow_C_s),
  precip_mean  = mean(dat$precipWaterYr_cm_s),
  precip_sd  = sd(dat$precipWaterYr_cm_s)
)

# define a wrapper function that samples from these distributions
sample_env <- function(env_params) {
  # We generate one value for each covariate per iteration, and return it 
  # as a named list
  soilM_now <- rnorm(1, mean = env_params$soilM_mu, sd = env_params$soilM_sd)
  soilT_now <- rnorm(1, mean = env_params$soilT_mu, sd = env_params$soilT_sd)
  temp_now  <- rnorm(1, mean = env_params$temp_mu, sd = env_params$temp_sd)
  precip_now <- rnorm(1, mean = env_params$precip_mean, sd  = env_params$precip_sd)
  
  out        <- list(soilM = soilM_now, soilT = soilT_now, temp = temp_now, precip = precip_now)
  return(out)
}

### define the IPM
stoch_env_DI_ipm <- init_ipm(sim_gen = "general", # make a general IPM
                             di_dd = "di", # make it density independent
                             det_stoch = "stoch",
                             kern_param = "param") %>% 
  define_kernel(
    name          = "P_plot", # survival 
    # We add d_ht to formula to make sure integration is handled correctly.
    # This variable is generated internally by make_ipm(), so we don't need
    # to do anything else.
    formula       = (1-p_b_plot) * s_plot * g_plot * d_ht,
    family        = "CC",
    g_plot        = dnorm(ht_2, mean = g_mu, sd = g_sd),
    g_mu          = g_int + g_slope * ht_1 + 
      g_soilM * soilM + g_precip * precip + g_temp * temp + #env covariates
      g_r_plot, # ranef of plot
    s_plot        = inv_logit_r(linPred = (s_int + s_slope*ht_1 + 
                                             s_soilM * soilM + s_soilT * soilT + s_temp * temp + # env covariates
                                             s_r_plot # ranef of plot
    )),
    p_b_plot      = inv_logit_r(linPred = (p_b_int + p_b_slope * ht_1 + p_b_slope_2 * I((ht_1)^2) + 
                                             p_b_soilM * soilM + p_b_soilT * soilT + p_b_temp * temp + p_b_precip * precip + # env covariates
                                             p_b_r_plot # ranef of plot 
    )), 
    data_list     = all_params_list,
    states        = list(c('ht')),
    uses_par_sets = TRUE,
    par_set_indices = list(plot = 1:18),
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions('norm', 'g_plot')
  ) %>%
  define_kernel(
    name          = "leave_seedlings_plot", ## leave seedling stage and go to rosette stage
    formula       = p_estab. * c_o_plot * d_ht,
    family        = 'DC', # Note that now, family = "DC" because it denotes a discrete -> continuous transition
    p_estab.      = p_estab,
    c_o_plot      = dnorm(ht_2, c_o_mu, c_o_sd),
    c_o_mu        = c_o_int + c_o_soilM * soilM  + # env covariate
      c_o_r_plot,
    data_list     = all_params_list,
    states        = list(c('ht', "s")),   # Note that here, we add "s" to our list in states, because this kernel uses seedlings 
    uses_par_sets = TRUE,
    par_set_indices = list(plot = 1:18),
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions('norm','c_o_plot')
  ) %>%
  define_kernel(
    name    = "repro_to_seedlings_plot",
    formula       = (goSdlng.) * (p_b_plot * b_plot * d_ht),
    family        = "CD",
    goSdlng.      = goSdlng,
    p_b_plot      = inv_logit_r(linPred = (p_b_int + p_b_slope * ht_1 + p_b_slope_2 *I(( ht_1)^2) + 
                                             p_b_soilM * soilM + p_b_soilT * soilT + p_b_temp * temp + p_b_precip * precip + # env covariates
                                             p_b_r_plot # ranef of plot 
    )),
    b_plot        = exp(b_int + b_slope * ht_1 + 
                          b_temp * temp + b_precip * precip + # env covariates
                          b_r_plot # ranef of plot
    ),
    data_list     = all_params_list,
    states        = list(c('ht', 's')),
    uses_par_sets = TRUE,
    par_set_indices = list(plot = 1:18),
    evict_cor     = FALSE
  ) %>%
  define_kernel(
    name          = 'seedbank_to_seedlings',
    formula       = outSB.,
    family        = 'DD',
    outSB.        = outSB,
    data_list     = all_params_list,
    states        = list(c('b', 's')),
    uses_par_sets = FALSE,
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name    = "stay_seedbank",
    formula       = staySB.,
    family        = "DD",
    staySB.        = staySB,
    data_list     = all_params_list,
    states        = list(c('b')),
    uses_par_sets = FALSE,
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name          = 'repro_to_seedbank_plot',
    formula       = (goSB.) * (p_b_plot * b_plot * d_ht),
    family        = 'CD',
    goSB.         = goSB, 
    p_b_plot      = inv_logit_r(linPred = (p_b_int + p_b_slope * ht_1 + p_b_slope_2 *I((ht_1)^2) + 
                                             p_b_soilM * soilM + p_b_soilT * soilT + p_b_temp * temp + p_b_precip * precip + # env covariates
                                             p_b_r_plot # ranef of plot 
    )),
    b_plot        = exp(b_int + b_slope * ht_1 + 
                          b_temp * temp + b_precip * precip + # env covariates
                          b_r_plot # ranef of plot
    ),
    data_list     = all_params_list,
    states        = list(c('b', 'ht')),
    uses_par_sets = TRUE,
    par_set_indices = list(plot = 1:18),
    evict_cor = FALSE
  ) 

## define the starting and ending states for each kernel
stoch_env_DI_ipm <- stoch_env_DI_ipm  %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P_plot", "leave_seedlings_plot", "repro_to_seedlings_plot", "seedbank_to_seedlings", "stay_seedbank", "repro_to_seedbank_plot"),
      int_rule     = c(rep("midpoint", 6)),
      state_start    = c('ht', "s", "ht", "b", "b", "ht"),
      state_end      = c("ht", "ht", "s", "s", "b", "b")
    )
  ) %>%
  define_domains(
    # We can pass the variables we created above into define_domains
    ht = c(L, U, n)
  ) %>%
  define_env_state(
    env_covs   = sample_env(env_params),
    data_list  = list(env_params = env_params,
                      sample_env = sample_env)
  ) %>% 
  define_pop_state(
    # We can also pass them into define_pop_state
    pop_vectors = list(
      n_ht = init_pop_vec,
      n_b  = init_seed_bank,
      n_s  = init_seedlings 
    )
  ) %>%
  make_ipm(iterate = TRUE, iterations = 100,
           normalize_pop_size = FALSE,
           usr_funs = my_funs, 
           kernel_seq = sample(1:18, 100, replace = TRUE)
  )

#### Stochastic, density-dependent IPM for all data with random effect for site and environmental covariates ####
### Vital Rate model names: survMod_dd_env, sizeMod_dd_env, seedMod_dd_env, flwrMod_dd_env, recMod_dd_env, p.estab.est, outSB.est, staySB.est, goSB.est, goSdlng.est
### define the parameter list
fixed_list <- list(
  g_int     = fixef(sizeMod_env_dd)[1], # growth model intercept
  g_slope   = fixef(sizeMod_env_dd)[2], # growth model leaf size slope
  g_soilM   = fixef(sizeMod_env_dd)[3], # growth model soil moisture coefficient
  g_temp    = fixef(sizeMod_env_dd)[4], # growth model mean growing season temp coefficient
  g_dd      = fixef(sizeMod_env_dd)[5], # growth model density dep. coefficient
  g_sd      = sd(resid(sizeMod_env_dd)), # growth model sd
  s_int     = fixef(survMod_env_dd)[1], # survival model intercept
  s_slope   = fixef(survMod_env_dd)[2], # survival model leaf size slope 
  s_temp    = fixef(survMod_env_dd)[3], # survival model mean temp coefficient
  s_dd      = fixef(survMod_env_dd)[4], # survival model density dep. coefficient
  p_b_int   = fixef(flwrMod_env_dd)[1], # probability of flowering
  p_b_slope = fixef(flwrMod_env_dd)[2], # flowering model leaf size slope
  p_b_slope_2 = fixef(flwrMod_env_dd)[3], # flowering model leaf size^2 slope
  #p_b_soilM = fixef(flwrMod_env_dd)[4], # flowering model soil moisture coefficient
  #p_b_soilT = fixef(flwrMod_env_dd)[5], # flowering model soil temp coefficient
  p_b_temp  = fixef(flwrMod_env_dd)[4], # flowering model mean temp coefficient
  #p_b_precip = fixef(flwrMod_env_dd)[7], # flowering model precip coefficient
  b_int     = fixef(seedMod_env_dd)[1], # seed production
  b_slope   = fixef(seedMod_env_dd)[2], # seed model leaf size slope
  b_temp    = fixef(seedMod_env_dd)[3], # seed model mean temp coefficient
  b_precip  = fixef(seedMod_env_dd)[4], # seed model precip coefficient
  c_o_int   = fixef(recMod_env_dd)[1], # recruit size model intercept
  c_o_soilM = fixef(recMod_env_dd)[2], # recruit size model soil moisture coefficient
  c_o_dd    = fixef(recMod_env_dd)[3], # recruit size model density dep. coefficient
  c_o_sd    = sd(resid(recMod_env_dd)), # recruit size model sd
  goSdlng   = goSdlng.est, # Probability that non-seedbank seeds will germinate into seedlings in year t+1
  staySB = staySB.est, # Probability that a seed in the seedbank in year t will exit the seedbank in year t+1
  goSB = goSB.est, # probability that a seed produced by an adult plant in year t will enter the seedbank
  outSB = outSB.est, # probability that a seedbank seed will germinate to a seedling in year t+1
  p_estab = p.estab.est # probability that a seedling will establish into a rosette in t+1
)

## now get a list of values for the random effects
# random values for the growth model
g_r_int <- unlist(ranef(sizeMod_env_dd))
# random values for the survival model
s_r_int <- unlist(ranef(survMod_env_dd))
# random values for the flowering model
p_b_r_int <- unlist(ranef(flwrMod_env_dd))
# random values for the seed production model
b_r_int <- unlist(ranef(seedMod_env_dd))
# random values for the recruit size model
c_o_r_int <- unlist(ranef(recMod_env_dd))

# name the list of random effects
nms <- paste("r_", 1:18, sep = "")

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
all_params_list <- c(fixed_list, g_params, s_params, p_b_params, b_params, c_o_params)

### make environmental parameters
env_params <- list(
  soilM_mu = mean(dat$SoilMoisture_m3m3_s), # use a normal dist for soilM (?)
  soilM_sd = sd(dat$SoilMoisture_m3m3_s),
  soilT_mu = mean(dat$SoilTemp_grow_C_s, na.rm = TRUE), # use a norm dist for soilT (?)
  soilT_sd = sd(dat$SoilTemp_grow_C_s, na.rm = TRUE),
  temp_mu = mean(dat$tMean_grow_C_s),
  temp_sd = sd(dat$tMean_grow_C_s),
  precip_mean  = mean(dat$precipWaterYr_cm_s),
  precip_sd  = sd(dat$precipWaterYr_cm_s)
)

# We define a wrapper function that samples from these distributions

sample_env <- function(env_params) {
  # We generate one value for each covariate per iteration, and return it 
  # as a named list. We can reference the names in this list in vital rate 
  # expressions.
  soilM_now <- rnorm(1, mean = env_params$soilM_mu, sd = env_params$soilM_sd)
  soilT_now <- rnorm(1, mean = env_params$soilT_mu, sd = env_params$soilT_sd)
  temp_now  <- rnorm(1, mean = env_params$temp_mu, sd = env_params$temp_sd)
  precip_now <- rnorm(1, mean = env_params$precip_mean, sd  = env_params$precip_sd)
  
  out        <- list(soilM = soilM_now, soilT = soilT_now, temp = temp_now, precip = precip_now)
  
  return(out)
  
}

### define the IPM
stoch_env_DD_ipm <- init_ipm(sim_gen = "general", # make a general IPM
                             di_dd = "dd", # make it density dependent
                             det_stoch = "stoch",
                             kern_param = "param") %>% 
  define_kernel(
    name          = "P_plot", # survival 
    # We add d_ht to formula to make sure integration is handled correctly.
    # This variable is generated internally by make_ipm(), so we don't need
    # to do anything else.
    formula       = (1-p_b_plot) * s_plot * g_plot * d_ht,
    family        = "CC",
    g_plot        = dnorm(ht_2, mean = g_mu, sd = g_sd),
    g_mu          = g_int + g_slope * ht_1 + 
      g_soilM * soilM  + g_temp * temp + #env covariates
      g_dd * sum(n_ht_t) + #density dependence
      g_r_plot, # ranef of plot
    s_plot        = inv_logit_r(linPred = (s_int + s_slope*ht_1 + 
                                             s_temp * temp + # env covariates
                                             s_dd * sum(n_ht_t) + # density dependence
                                             s_r_plot # ranef of plot
    )),
    p_b_plot      = inv_logit_r(linPred = (p_b_int + p_b_slope * ht_1 + p_b_slope_2 * I(( ht_1)^2) + 
                                             p_b_temp * temp + # env cov.
                                             p_b_r_plot # ranef of plot 
    )), 
    data_list     = all_params_list,
    states        = list(c('ht')),
    uses_par_sets = TRUE,
    par_set_indices = list(plot = 1:18),
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions('norm', 'g_plot')
  ) %>%
  define_kernel(
    name          = "leave_seedlings_plot", ## leave seedling stage and go to rosette stage
    formula       = p_estab. * c_o_plot * d_ht,
    family        = 'DC', # Note that now, family = "DC" because it denotes a discrete -> continuous transition
    p_estab.      = p_estab,
    c_o_plot      = dnorm(ht_2, c_o_mu, c_o_sd),
    c_o_mu        = c_o_int + c_o_soilM * soilM  + # env covariate
      c_o_dd * sum(n_ht_t) + # density dependence
      c_o_r_plot,
    data_list     = all_params_list,
    states        = list(c('ht', "s")),   # Note that here, we add "s" to our list in states, because this kernel uses seedlings 
    uses_par_sets = TRUE,
    par_set_indices = list(plot = 1:18),
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions('norm','c_o_plot')
  ) %>%
  define_kernel(
    name    = "repro_to_seedlings_plot",
    formula       = (goSdlng.) * (p_b_plot * b_plot * d_ht),
    family        = "CD",
    goSdlng.      = goSdlng,
    p_b_plot      = inv_logit_r(linPred = (p_b_int + p_b_slope * ht_1 + p_b_slope_2 * I(( ht_1)^2) + 
                                             p_b_temp * temp + # env cov.
                                             p_b_r_plot # ranef of plot 
    )),
    b_plot        = exp(b_int + b_slope * ht_1 + 
                          b_temp * temp + b_precip * precip + # env cov.
                          b_r_plot # ranef of plot
    ),
    data_list     = all_params_list,
    states        = list(c('ht', 's')),
    uses_par_sets = TRUE,
    par_set_indices = list(plot = 1:18),
    evict_cor     = FALSE
  ) %>%
  define_kernel(
    name          = 'seedbank_to_seedlings',
    formula       = outSB.,
    family        = 'DD',
    outSB.        = outSB,
    data_list     = all_params_list,
    states        = list(c('b', 's')),
    uses_par_sets = FALSE,
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name    = "stay_seedbank",
    formula       = staySB.,
    family        = "DD",
    staySB.        = staySB,
    data_list     = all_params_list,
    states        = list(c('b')),
    uses_par_sets = FALSE,
    evict_cor = FALSE
  ) %>%
  define_kernel(
    name          = 'repro_to_seedbank_plot',
    formula       = (goSB.) * (p_b_plot * b_plot * d_ht),
    family        = 'CD',
    goSB.         = goSB, 
    p_b_plot      = inv_logit_r(linPred = (p_b_int + p_b_slope * ht_1 + p_b_slope_2 * I(( ht_1)^2) + 
                                             p_b_temp * temp + # env cov.
                                             p_b_r_plot # ranef of plot 
    )),
    b_plot        = exp(b_int + b_slope * ht_1 + 
                          b_temp * temp + b_precip * precip + # env cov.
                          b_r_plot # ranef of plot
    ),
    data_list     = all_params_list,
    states        = list(c('b', 'ht')),
    uses_par_sets = TRUE,
    par_set_indices = list(plot = 1:18),
    evict_cor = FALSE
  ) %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P_plot", "leave_seedlings_plot", "repro_to_seedlings_plot", "seedbank_to_seedlings", "stay_seedbank", "repro_to_seedbank_plot"),
      int_rule     = c(rep("midpoint", 6)),
      state_start    = c('ht', "s", "ht", "b", "b", "ht"),
      state_end      = c("ht", "ht", "s", "s", "b", "b")
    )
  ) %>%
  define_domains(
    # pass the variables we created above into define_domains
    ht = c(L, U, n)
  ) %>%
  define_env_state(
    env_covs   = sample_env(env_params),
    data_list  = list(env_params = env_params,
                      sample_env = sample_env)
  ) %>% 
  define_pop_state(
    # We can also pass them into define_pop_state
    pop_vectors = list(
      n_ht = init_pop_vec,
      n_b  = init_seed_bank,
      n_s  = init_seedlings 
    )
  ) %>%
  make_ipm(iterate = TRUE, iterations = 100,
           normalize_pop_size = FALSE,
           usr_funs = my_funs, 
           kernel_seq = sample(1:18, 100, replace = TRUE)
  )

## lambda is a generic function to compute per-capita growth rates. It has a number of different options depending on the type of model
lambda(stoch_env_DD_ipm)
# L = 1.09; U = 3.66; n = 500; log(lambda) = -0.3565
# L = 0.9; U = 4.0; log(lambda) = -0.3389
# L = 0.8; U = 4.1; log(lambda) = -0.3348

#### deterministic, density independent IPM for each plot ####
## model parameters stored in the det_DI_mods list

### run IPMs inside a for-loop
sites <- unique(dat$Site)
for (i in 1:length(sites)) {
  # get the site name for this 'i'
  site_now <- sites[i]
  # get the model list for this site
  modList_now <- det_DI_mods[[which(names(det_DI_mods) == site_now)]]
  ## get the vital rate model parameters
  data_list <- list(
    g_int     = coef(modList_now$growth)[1],
    g_slope   = coef(modList_now$growth)[2],
    g_sd      = summary(modList_now$growth)$sigma,
    s_int     = coef(modList_now$surv)[1],
    s_slope   = coef(modList_now$surv)[2],
    p_b_int   = coef(modList_now$flowering)[1], #probability of flowering
    p_b_slope = coef(modList_now$flowering)[2],
    p_b_slope_2 = coef(modList_now$flowering)[3],
    b_int   = coef(modList_now$seedProduction)[1], #seed production
    b_slope = coef(modList_now$seedProduction)[2],
    c_o_mu    = coef(modList_now$recruitDist), #recruit size distribution
    c_o_sd    = summary(modList_now$recruitDist)$sigma,
    goSdlng   = goSdlng.est, # Probability that non-seedbank seeds will germinate into seedlings in year t+1
    staySB = staySB.est, # Probability that a seed in the seedbank in year t will exit the seedbank in year t+1 
    goSB = goSB.est, # probability that a seed produced by an adult plant in year t will enter the seedbank
    outSB = outSB.est, # probability that a seedbank seed will germinate to a seedling in year t+1
    p_estab = p.estab.est # probability that a seedling will establish into a rosette in t+1
  )
  ## make IPM kernels
  temp_ipm <- init_ipm(sim_gen = "general", # make a general IPM
                       di_dd = "di", # make it density independent
                       det_stoch = "det") %>% # make it deterministic
    define_kernel(
      name          = "P", # survival 
      # We add d_ht to formula to make sure integration is handled correctly.
      # This variable is generated internally by make_ipm(), so we don't need
      # to do anything else.
      formula       = (1-p_b.) * s. * g. * d_ht,
      family        = "CC",
      g.             = dnorm(ht_2, g_mu, g_sd),
      g_mu          = g_int + g_slope * ht_1,
      s.             = inv_logit(s_int, s_slope, ht_1),
      p_b.          = inv_logit_2(p_b_int, p_b_slope, p_b_slope_2, ht_1),
      data_list     = data_list,
      states        = list(c('ht')),
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm', 'g.')
    ) %>%
    define_kernel(
      name          = "leave_seedlings", ## leave seedling stage and go to rosette stage
      formula       = p_estab. * c_o. * d_ht,
      family        = 'DC', # Note that now, family = "DC" because it denotes a discrete -> continuous transition
      p_estab.      = p_estab,
      c_o.          = dnorm(ht_2, c_o_mu, c_o_sd),
      data_list     = data_list,
      states        = list(c('ht', "s")),   # Note that here, we add "s" to our list in states, because this kernel uses seedlings 
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm','c_o.')
    ) %>%
    define_kernel(
      name    = "repro_to_seedlings",
      formula       = (goSdlng.) * (p_b. * b. * d_ht),
      family        = "CD",
      goSdlng.      = goSdlng,
      p_b.          = inv_logit_2(p_b_int, p_b_slope, p_b_slope_2, ht_1),
      b.            = exp(b_int + b_slope * ht_1),
      data_list     = data_list,
      states        = list(c('ht', 's')),
      uses_par_sets = FALSE,
      evict_cor     = FALSE
    ) %>%
    define_kernel(
      name          = 'seedbank_to_seedlings',
      formula       = outSB.,
      family        = 'DD',
      outSB.        = outSB,
      data_list     = data_list,
      states        = list(c('b', 's')),
      uses_par_sets = FALSE,
      evict_cor = FALSE
    ) %>%
    define_kernel(
      name    = "stay_seedbank",
      formula       = staySB.,
      family        = "DD",
      staySB.        = staySB,
      data_list     = data_list,
      states        = list(c('b')),
      uses_par_sets = FALSE,
      evict_cor = FALSE
    ) %>%
    define_kernel(
      name          = 'repro_to_seedbank',
      formula       = (goSB.) * (p_b. * b. * d_ht),
      family        = 'CD',
      goSB.          = goSB, 
      p_b.          = inv_logit_2(p_b_int, p_b_slope, p_b_slope_2, ht_1),
      b.            = exp(b_int + b_slope * ht_1),
      data_list     = data_list,
      states        = list(c('b', 'ht')),
      uses_par_sets = FALSE,
      evict_cor = FALSE
    ) %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "leave_seedlings", "repro_to_seedlings", "seedbank_to_seedlings", "stay_seedbank", "repro_to_seedbank"),
        int_rule     = c(rep("midpoint", 6)),
        state_start    = c('ht', "s", "ht", "b", "b", "ht"),
        state_end      = c("ht", "ht", "s", "s", "b", "b")
      )
    )
  
  ## calculate starting population state vectors
  # starting size of the seedbank
  # seedbank estimate data stored for each site is in the seeds_est_site d.f.
  init_seed_bank <- deframe(round(seeds_est_site[seeds_est_site$Site == site_now,"seedbank_est"]))
  # Use seedling data to calculate the starting number of seedlings 
  # seedling data stored in the 'seedling_site' d.f.
  init_seedlings <- deframe(round(seedlings_site[seedlings_site$Site == site_now, "Seedlings_t"],0))
  # starting number of individuals in the continuous stage
  ht <- dat$log_LL_t
  init_pop_vec   <- pnorm(rnorm(500, mean = mean(ht, na.rm = TRUE), sd = sd(ht, na.rm = TRUE)), mean = mean(ht, na.rm = TRUE), sd = sd(ht, na.rm = TRUE)) #discretize_pop_vector(trait_values = ht, n_mesh = 500, pad_low = 0.8, pad_high = 1.2)$n_ht
  
  ## code to run the IPM
  temp_ipm <- temp_ipm %>%
    define_domains(
      # We can pass the variables we created above into define_domains
      ht = c(L, U, n)
    ) %>%
    define_pop_state(
      # We can also pass them into define_pop_state
      pop_vectors = list(
        n_ht = init_pop_vec,
        n_b  = init_seed_bank,
        n_s  = init_seedlings 
      )
    ) %>%
    make_ipm(iterations = 100,
             normalize_pop_size = FALSE,
             usr_funs = list(inv_logit   = inv_logit,
                             inv_logit_2 = inv_logit_2), return_main_env = TRUE )
  
  ## rename the ipm object to have the name of site
  assign(paste0(site_now,"__det_DI_ipm"), value = temp_ipm)
  
  ## Done!
}

# #### deterministic, density dependent IPM for each plot  ####
# #maybe cannot have deterministic, density dependent IPMs? 
# inv_logit_r <- function(linPred) {
#   1/(1 + exp(-(linPred)))
# }
# # run IPMs inside a for-loop
# for (i in 1:length(sites)) {
#   # get the site name for this 'i'
#   site_now <- sites[i]
#   # get the model list for this site
#   modList_now <- det_DD_mods[[which(names(det_DD_mods) == site_now)]]
#   ## get the vital rate model parameters
#   data_list <- list(
#     g_int     = coef(modList_now$growth)[1],
#     g_slope   = coef(modList_now$growth)[2],
#     g_dd      = coef(modList_now$growth)[3],
#     g_sd      = summary(modList_now$growth)$sigma,
#     s_int     = coef(modList_now$surv)[1],
#     s_slope   = coef(modList_now$surv)[2],
#     s_dd      = coef(modList_now$surv)[3],
#     p_b_int   = coef(modList_now$flowering)[1], #probability of flowering
#     p_b_slope = coef(modList_now$flowering)[2],
#     p_b_slope_2 = coef(modList_now$flowering)[3],
#     b_int   = coef(modList_now$seedProduction)[1], #seed production
#     b_slope = coef(modList_now$seedProduction)[2],
#     c_o_int   = coef(modList_now$recruitDist)[1], #recruit size distribution
#     c_o_dd    = coef(modList_now$recruitDist)[2],
#     c_o_sd    = summary(modList_now$recruitDist)$sigma,
#     goSdlng   = goSdlng.est, # Probability that non-seedbank seeds will germinate into seedlings in year t+1
#     staySB = staySB.est, # Probability that a seed in the seedbank in year t will exit the seedbank in year t+1 
#     goSB = goSB.est, # probability that a seed produced by an adult plant in year t will enter the seedbank
#     outSB = outSB.est, # probability that a seedbank seed will germinate to a seedling in year t+1
#     p_estab = p.estab.est # probability that a seedling will establish into a rosette in t+1
#   )
#   ## make IPM kernels
#   temp_ipm <- init_ipm(sim_gen = "general", # make a general IPM
#                        di_dd = "dd", # make it density independent
#                        det_stoch = "det") %>% # make it deterministic
#     define_kernel(
#       name          = "P", # survival 
#       # We add d_ht to formula to make sure integration is handled correctly.
#       # This variable is generated internally by make_ipm(), so we don't need
#       # to do anything else.
#       formula       = (1-p_b.) * s. * g. * d_ht,
#       family        = "CC",
#       g.             = dnorm(ht_2, g_mu, g_sd),
#       g_mu          = g_int + g_slope * ht_1 + g_dd * sum(n_ht_t),
#       s.             = inv_logit_r(s_int + s_slope * ht_1 + s_dd * sum(n_ht_t)),
#       p_b.          = inv_logit_2(p_b_int, p_b_slope, p_b_slope_2, ht_1),
#       data_list     = data_list,
#       states        = list(c('ht')),
#       uses_par_sets = FALSE,
#       evict_cor     = TRUE,
#       evict_fun     = truncated_distributions('norm', 'g.')
#     ) %>%
#     define_kernel(
#       name          = "leave_seedlings", ## leave seedling stage and go to rosette stage
#       formula       = p_estab. * c_o. * d_ht,
#       family        = 'DC', # Note that now, family = "DC" because it denotes a discrete -> continuous transition
#       p_estab.      = p_estab,
#       c_o.          = dnorm(ht_2, c_o_mu, c_o_sd),
#       c_o_mu        = c_o_int + c_o_dd * sum(n_ht_t),
#       data_list     = data_list,
#       states        = list(c('ht', "s")),   # Note that here, we add "s" to our list in states, because this kernel uses seedlings 
#       uses_par_sets = FALSE,
#       evict_cor     = TRUE,
#       evict_fun     = truncated_distributions('norm','c_o.')
#     ) %>%
#     define_kernel(
#       name    = "repro_to_seedlings",
#       formula       = (goSdlng.) * (p_b. * b. * d_ht),
#       family        = "CD",
#       goSdlng.      = goSdlng,
#       p_b.          = inv_logit_2(p_b_int, p_b_slope, p_b_slope_2, ht_1),
#       b.            = exp(b_int + b_slope * ht_1),
#       data_list     = data_list,
#       states        = list(c('ht', 's')),
#       uses_par_sets = FALSE,
#       evict_cor     = FALSE
#     ) %>%
#     define_kernel(
#       name          = 'seedbank_to_seedlings',
#       formula       = outSB.,
#       family        = 'DD',
#       outSB.        = outSB,
#       data_list     = data_list,
#       states        = list(c('b', 's')),
#       uses_par_sets = FALSE,
#       evict_cor = FALSE
#     ) %>%
#     define_kernel(
#       name    = "stay_seedbank",
#       formula       = staySB.,
#       family        = "DD",
#       staySB.        = staySB,
#       data_list     = data_list,
#       states        = list(c('b')),
#       uses_par_sets = FALSE,
#       evict_cor = FALSE
#     ) %>%
#     define_kernel(
#       name          = 'repro_to_seedbank',
#       formula       = (goSB.) * (p_b. * b. * d_ht),
#       family        = 'CD',
#       goSB.          = goSB, 
#       p_b.          = inv_logit_2(p_b_int, p_b_slope, p_b_slope_2, ht_1),
#       b.            = exp(b_int + b_slope * ht_1),
#       data_list     = data_list,
#       states        = list(c('b', 'ht')),
#       uses_par_sets = FALSE,
#       evict_cor = FALSE
#     ) %>%
#     define_impl(
#       make_impl_args_list(
#         kernel_names = c("P", "leave_seedlings", "repro_to_seedlings", "seedbank_to_seedlings", "stay_seedbank", "repro_to_seedbank"),
#         int_rule     = c(rep("midpoint", 6)),
#         state_start    = c('ht', "s", "ht", "b", "b", "ht"),
#         state_end      = c("ht", "ht", "s", "s", "b", "b")
#       )
#     )
#   
#   ## calculate starting population state vectors
#   # starting size of the seedbank
#   # seedbank estimate data stored for each site is in the seeds_est_site d.f.
#   init_seed_bank <- deframe(round(seeds_est_site[seeds_est_site$Site == site_now,"seedbank_est"]))
#   # Use seedling data to calculate the starting number of seedlings 
#   # seedling data stored in the 'seedling_site' d.f.
#   init_seedlings <- deframe(round(seedlings_site[seedlings_site$Site == site_now, "Seedlings_t"],0))
#   # starting number of individuals in the continuous stage
#   ht <- dat$log_LL_t
#   init_pop_vec   <- discretize_pop_vector(trait_values = ht, n_mesh = 500, pad_low = 0.8, pad_high = 1.2)$n_ht
#   
#   ## code to run the IPM
#   temp_ipm <- temp_ipm %>%
#     define_domains(
#       # We can pass the variables we created above into define_domains
#       ht = c(L, U, n)
#     ) %>%
#     define_pop_state(
#       # We can also pass them into define_pop_state
#       pop_vectors = list(
#         n_ht = init_pop_vec,
#         n_b  = init_seed_bank,
#         n_s  = init_seedlings 
#       )
#     ) %>%
#     make_ipm(iterations = 100,
#              usr_funs = list(inv_logit_r   = inv_logit_r,
#                              inv_logit_2 = inv_logit_2), return_main_env = TRUE )
#   
#   ## rename the ipm object to have the name of site
#   assign(paste0(site_now,"__det_DD_ipm"), value = temp_ipm)
#   
#   ## Done!
# }
#### stochastic, density independent IPM for each plot  ####
# define environmental covariates
env_params <- list(
  soilM_mu = mean(dat$SoilMoisture_m3m3_s), # use a normal dist for soilM (?)
  soilM_sd = sd(dat$SoilMoisture_m3m3_s),
  soilT_mu = mean(dat$SoilTemp_grow_C_s, na.rm = TRUE), # use a norm dist for soilT (?)
  soilT_sd = sd(dat$SoilTemp_grow_C_s, na.rm = TRUE),
  temp_mu = mean(dat$tMean_grow_C_s),
  temp_sd = sd(dat$tMean_grow_C_s),
  precip_mean  = mean(dat$precipWaterYr_cm_s),
  precip_sd  = sd(dat$precipWaterYr_cm_s)
)

# define a wrapper function that samples from these distributions
sample_env <- function(env_params) {
  # We generate one value for each covariate per iteration, and return it 
  # as a named list
  soilM_now <- rnorm(1, mean = env_params$soilM_mu, sd = env_params$soilM_sd)
  soilT_now <- rnorm(1, mean = env_params$soilT_mu, sd = env_params$soilT_sd)
  temp_now  <- rnorm(1, mean = env_params$temp_mu, sd = env_params$temp_sd)
  precip_now <- rnorm(1, mean = env_params$precip_mean, sd  = env_params$precip_sd)
  
  out        <- list(soilM = soilM_now, soilT = soilT_now, temp = temp_now, precip = precip_now)
  return(out)
}
# make for-loop
sites <- unique(dat$Site)
for (i in 1:length(sites)) {
  # get the site name for this 'i'
  site_now <- sites[i]
  # get the model list for this site
  modList_now <- stoch_DI_mods[[which(names(stoch_DI_mods) == site_now)]]
  ## get the vital rate model parameters
  data_list <- list(
    g_int     = coef(modList_now$growth)[1],
    g_slope   = coef(modList_now$growth)[2],
    g_soilM   = coef(modList_now$growth)[3], 
    g_temp    = coef(modList_now$growth)[4],
    g_sd      = summary(modList_now$growth)$sigma,
    s_int     = coef(modList_now$surv)[1],
    s_slope   = coef(modList_now$surv)[2],
    s_soilM   = coef(modList_now$surv)[3],
    s_temp    = coef(modList_now$surv)[4],
    s_soilT   = coef(modList_now$surv)[5],
    p_b_int   = coef(modList_now$flowering)[1], #probability of flowering
    p_b_slope = coef(modList_now$flowering)[2],
    p_b_slope_2 = coef(modList_now$flowering)[3],
    p_b_soilM = coef(modList_now$flowering)[4],
    p_b_temp  = coef(modList_now$flowering)[5],
    p_b_precip = coef(modList_now$flowering)[6],
    b_int   = coef(modList_now$seedProduction)[1], #seed production
    b_slope = coef(modList_now$seedProduction)[2],
    b_temp  = coef(modList_now$seedProduction)[3],
    c_o_int    = coef(modList_now$recruitDist)[1], #recruit size distribution
    c_o_soilM = coef(modList_now$recruitDist)[2],
    c_o_sd    = summary(modList_now$recruitDist)$sigma,
    goSdlng   = goSdlng.est, # Probability that non-seedbank seeds will germinate into seedlings in year t+1
    staySB = staySB.est, # Probability that a seed in the seedbank in year t will exit the seedbank in year t+1 
    goSB = goSB.est, # probability that a seed produced by an adult plant in year t will enter the seedbank
    outSB = outSB.est, # probability that a seedbank seed will germinate to a seedling in year t+1
    p_estab = p.estab.est # probability that a seedling will establish into a rosette in t+1
  )
  ## make IPM kernels
  temp_ipm <- init_ipm(sim_gen = "general", # make a general IPM
                       di_dd = "di", # make it density independent
                       det_stoch = "stoch", 
                       kern_param = "param") %>% # make it deterministic
    define_kernel(
      name          = "P", # survival 
      formula       = (1-p_b.) * s. * g. * d_ht,
      family        = "CC",
      g.            = dnorm(ht_2, mean = g_mu, sd = g_sd),
      g_mu          = g_int + g_slope * ht_1 + 
        g_soilM * soilM  + g_temp * temp,  #env covariates
      s.            = inv_logit_r(linPred = (s_int + s_slope*ht_1 + 
                                               s_soilM * soilM + s_soilT * soilT + s_temp * temp # env covariates
      )),
      p_b.      = inv_logit_r(linPred = (p_b_int + p_b_slope * ht_1 + p_b_slope_2 * (I(ht_1)^2) + 
                                           p_b_soilM * soilM  + p_b_temp * temp + p_b_precip * precip  # env cov.
      )), 
      data_list     = data_list,
      states        = list(c('ht')),
      uses_par_sets = FALSE,
      par_set_indices = list(plot = 1:18),
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm', 'g.')
    ) %>%
    define_kernel(
      name          = "leave_seedlings", ## leave seedling stage and go to rosette stage
      formula       = p_estab. * c_o. * d_ht,
      family        = 'DC', # Note that now, family = "DC" because it denotes a discrete -> continuous transition
      p_estab.      = p_estab,
      c_o.          = dnorm(ht_2, c_o_mu., c_o_sd),
      c_o_mu.       = c_o_int + c_o_soilM * soilM,
      data_list     = data_list,
      states        = list(c('ht', "s")),   # Note that here, we add "s" to our list in states, because this kernel uses seedlings 
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm','c_o.')
    ) %>%
    define_kernel(
      name    = "repro_to_seedlings",
      formula       = (goSdlng.) * (p_b. * b. * d_ht),
      family        = "CD",
      goSdlng.      = goSdlng,
      p_b.      = inv_logit_r(linPred = (p_b_int + p_b_slope * ht_1 + p_b_slope_2 * (I( ht_1)^2) + 
                                           p_b_soilM * soilM  + p_b_temp * temp + p_b_precip * precip  # env cov.
      )
      ),
      b.            = exp(b_int + b_slope * ht_1 + b_temp * temp ),
      data_list     = data_list,
      states        = list(c('ht', 's')),
      uses_par_sets = FALSE,
      evict_cor     = FALSE
    ) %>%
    define_kernel(
      name          = 'seedbank_to_seedlings',
      formula       = outSB.,
      family        = 'DD',
      outSB.        = outSB,
      data_list     = data_list,
      states        = list(c('b', 's')),
      uses_par_sets = FALSE,
      evict_cor = FALSE
    ) %>%
    define_kernel(
      name    = "stay_seedbank",
      formula       = staySB.,
      family        = "DD",
      staySB.        = staySB,
      data_list     = data_list,
      states        = list(c('b')),
      uses_par_sets = FALSE,
      evict_cor = FALSE
    ) %>%
    define_kernel(
      name          = 'repro_to_seedbank',
      formula       = (goSB.) * (p_b. * b. * d_ht),
      family        = 'CD',
      goSB.          = goSB, 
      p_b.      = inv_logit_r(linPred = (p_b_int + p_b_slope * ht_1 + p_b_slope_2 * (I( ht_1)^2) + 
                                           p_b_soilM * soilM +  p_b_temp * temp + p_b_precip * precip  # env cov.
      )
      ),
      b.            = exp(b_int + b_slope * ht_1 + b_temp * temp ),
      data_list     = data_list,
      states        = list(c('b', 'ht')),
      uses_par_sets = FALSE,
      evict_cor = FALSE
    ) %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "leave_seedlings", "repro_to_seedlings", "seedbank_to_seedlings", "stay_seedbank", "repro_to_seedbank"),
        int_rule     = c(rep("midpoint", 6)),
        state_start    = c('ht', "s", "ht", "b", "b", "ht"),
        state_end      = c("ht", "ht", "s", "s", "b", "b")
      )
    )
  
  ## calculate starting population state vectors
  # starting size of the seedbank
  # seedbank estimate data stored for each site is in the seeds_est_site d.f.
  init_seed_bank <- deframe(round(seeds_est_site[seeds_est_site$Site == site_now,"seedbank_est"]))
  # Use seedling data to calculate the starting number of seedlings 
  # seedling data stored in the 'seedling_site' d.f.
  init_seedlings <- deframe(round(seedlings_site[seedlings_site$Site == site_now, "Seedlings_t"],0))
  # starting number of individuals in the continuous stage
  ht <- dat$log_LL_t
  init_pop_vec   <- pnorm(rnorm(500, mean = mean(ht, na.rm = TRUE), sd = sd(ht, na.rm = TRUE)), mean = mean(ht, na.rm = TRUE), sd = sd(ht, na.rm = TRUE)) #discretize_pop_vector(trait_values = ht, n_mesh = 500, pad_low = 0.8, pad_high = 1.2)$n_ht
  
  ## code to run the IPM
  temp_ipm <- temp_ipm %>%
    define_domains(
      # We can pass the variables we created above into define_domains
      ht = c(L, U, n)
    ) %>%
    define_pop_state(
      # We can also pass them into define_pop_state
      pop_vectors = list(
        n_ht = init_pop_vec,
        n_b  = init_seed_bank,
        n_s  = init_seedlings 
      )
    ) %>%
    define_env_state(
      env_covs   = sample_env(env_params),
      data_list  = list(env_params = env_params,
                        sample_env = sample_env)) %>%
    make_ipm(iterations = 100,
             normalize_pop_size = FALSE,
             usr_funs = list(inv_logit_r   = inv_logit_r), return_main_env = TRUE )
  
  ## rename the ipm object to have the name of site
  assign(paste0(site_now,"__stoch_DI_ipm"), value = temp_ipm)
  
  ## Done!
}
#### stochastic, density dependent IPM for each plot  ####
# define environmental covariates
env_params <- list(
  soilM_mu = mean(dat$SoilMoisture_m3m3_s), # use a normal dist for soilM (?)
  soilM_sd = sd(dat$SoilMoisture_m3m3_s),
  soilT_mu = mean(dat$SoilTemp_grow_C_s, na.rm = TRUE), # use a norm dist for soilT (?)
  soilT_sd = sd(dat$SoilTemp_grow_C_s, na.rm = TRUE),
  temp_mu = mean(dat$tMean_grow_C_s),
  temp_sd = sd(dat$tMean_grow_C_s),
  precip_mean  = mean(dat$precipWaterYr_cm_s),
  precip_sd  = sd(dat$precipWaterYr_cm_s)
)

# define a wrapper function that samples from these distributions
sample_env <- function(env_params) {
  # We generate one value for each covariate per iteration, and return it 
  # as a named list
  soilM_now <- rnorm(1, mean = env_params$soilM_mu, sd = env_params$soilM_sd)
  soilT_now <- rnorm(1, mean = env_params$soilT_mu, sd = env_params$soilT_sd)
  temp_now  <- rnorm(1, mean = env_params$temp_mu, sd = env_params$temp_sd)
  precip_now <- rnorm(1, mean = env_params$precip_mean, sd  = env_params$precip_sd)
  
  out        <- list(soilM = soilM_now, soilT = soilT_now, temp = temp_now, precip = precip_now)
  return(out)
}
# make for-loop
sites <- unique(dat$Site)
for (i in 1:length(sites)) {
  # get the site name for this 'i'
  site_now <- sites[i]
  # get the model list for this site
  modList_now <- stoch_DD_mods[[which(names(stoch_DD_mods) == site_now)]]
  ## get the vital rate model parameters
  data_list <- list(
    g_int     = coef(modList_now$growth)[1],
    g_slope   = coef(modList_now$growth)[2],
    g_soilM   = coef(modList_now$growth)[3], 
    g_temp    = coef(modList_now$growth)[4],
    g_dd      = coef(modList_now$growth)[5],
    g_sd      = summary(modList_now$growth)$sigma,
    s_int     = coef(modList_now$surv)[1],
    s_slope   = coef(modList_now$surv)[2],
    s_soilM   = coef(modList_now$surv)[3],
    s_temp    = coef(modList_now$surv)[4],
    s_soilT   = coef(modList_now$surv)[5],
    s_dd      = coef(modList_now$surv)[6],
    p_b_int   = coef(modList_now$flowering)[1], #probability of flowering
    p_b_slope = coef(modList_now$flowering)[2],
    p_b_slope_2 = coef(modList_now$flowering)[3],
    p_b_soilM = coef(modList_now$flowering)[4],
    p_b_temp  = coef(modList_now$flowering)[5],
    p_b_precip = coef(modList_now$flowering)[6],
    b_int   = coef(modList_now$seedProduction)[1], #seed production
    b_slope = coef(modList_now$seedProduction)[2],
    b_temp  = coef(modList_now$seedProduction)[3],
    c_o_int    = coef(modList_now$recruitDist)[1], #recruit size distribution
    c_o_soilM = coef(modList_now$recruitDist)[2],
    c_o_dd    = coef(modList_now$recruitDist)[3],
    c_o_sd    = summary(modList_now$recruitDist)$sigma,
    goSdlng   = goSdlng.est, # Probability that non-seedbank seeds will germinate into seedlings in year t+1
    staySB = staySB.est, # Probability that a seed in the seedbank in year t will exit the seedbank in year t+1 
    goSB = goSB.est, # probability that a seed produced by an adult plant in year t will enter the seedbank
    outSB = outSB.est, # probability that a seedbank seed will germinate to a seedling in year t+1
    p_estab = p.estab.est # probability that a seedling will establish into a rosette in t+1
  )
  ## make IPM kernels
  temp_ipm <- init_ipm(sim_gen = "general", # make a general IPM
                       di_dd = "dd", # make it density independent
                       det_stoch = "stoch", 
                       kern_param = "param") %>% # make it deterministic
    define_kernel(
      name          = "P", # survival 
      # We add d_ht to formula to make sure integration is handled correctly.
      # This variable is generated internally by make_ipm(), so we don't need
      # to do anything else.
      formula       = (1-p_b.) * s. * g. * d_ht,
      family        = "CC",
      g.             = dnorm(ht_2, g_mu, g_sd), #0.90
      g_mu          = g_int + g_slope * ht_1 + g_soilM * soilM + g_temp * temp + g_dd * sum(n_ht_t),
      s.             = inv_logit_r(s_int + s_slope * ht_1 + s_soilM * soilM + s_soilT * soilT + s_dd * sum(n_ht_t)), # 0.62
      p_b.          = inv_logit_r(p_b_int + p_b_slope * ht_1 +  p_b_slope_2 * I(  ht_1)^2 + p_b_soilM * soilM + p_b_temp * temp + p_b_precip * precip), # 0.0029
      data_list     = data_list,
      states        = list(c('ht')),
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm', 'g.')
    ) %>%
    define_kernel(
      name          = "leave_seedlings", ## leave seedling stage and go to rosette stage
      formula       = p_estab. * c_o. * d_ht,
      family        = 'DC', # Note that now, family = "DC" because it denotes a discrete -> continuous transition
      p_estab.      = p_estab,
      c_o.          = dnorm(ht_2, c_o_mu., c_o_sd),
      c_o_mu.       = c_o_int + c_o_soilM * soilM + c_o_dd * sum(n_ht_t),
      data_list     = data_list,
      states        = list(c('ht', "s")),   # Note that here, we add "s" to our list in states, because this kernel uses seedlings 
      uses_par_sets = FALSE,
      evict_cor     = TRUE,
      evict_fun     = truncated_distributions('norm','c_o.')
    ) %>%
    define_kernel(
      name    = "repro_to_seedlings",
      formula       = (goSdlng.) * (p_b. * b. * d_ht),
      family        = "CD",
      goSdlng.      = goSdlng,
      p_b.          = inv_logit_r(p_b_int + p_b_slope * ht_1 +  p_b_slope_2 *I(ht_1)^2 + p_b_soilM * soilM + p_b_temp * temp + p_b_precip * precip),
      b.            = exp(b_int + b_slope * ht_1 + b_temp * temp ),
      data_list     = data_list,
      states        = list(c('ht', 's')),
      uses_par_sets = FALSE,
      evict_cor     = FALSE
    ) %>%
    define_kernel(
      name          = 'seedbank_to_seedlings',
      formula       = outSB.,
      family        = 'DD',
      outSB.        = outSB,
      data_list     = data_list,
      states        = list(c('b', 's')),
      uses_par_sets = FALSE,
      evict_cor = FALSE
    ) %>%
    define_kernel(
      name    = "stay_seedbank",
      formula       = staySB.,
      family        = "DD",
      staySB.        = staySB,
      data_list     = data_list,
      states        = list(c('b')),
      uses_par_sets = FALSE,
      evict_cor = FALSE
    ) %>%
    define_kernel(
      name          = 'repro_to_seedbank',
      formula       = (goSB.) * (p_b. * b. * d_ht),
      family        = 'CD',
      goSB.          = goSB, 
      p_b.          = inv_logit_r(p_b_int + p_b_slope * ht_1 +  p_b_slope_2 *I(ht_1)^2 + p_b_soilM * soilM + p_b_temp * temp + p_b_precip * precip),
      b.            = exp(b_int + b_slope * ht_1 + b_temp * temp ),
      data_list     = data_list,
      states        = list(c('b', 'ht')),
      uses_par_sets = FALSE,
      evict_cor = FALSE
    ) %>%
    define_impl(
      make_impl_args_list(
        kernel_names = c("P", "leave_seedlings", "repro_to_seedlings", "seedbank_to_seedlings", "stay_seedbank", "repro_to_seedbank"),
        int_rule     = c(rep("midpoint", 6)),
        state_start    = c('ht', "s", "ht", "b", "b", "ht"),
        state_end      = c("ht", "ht", "s", "s", "b", "b")
      )
    )
  
  ## calculate starting population state vectors
  # starting size of the seedbank
  # seedbank estimate data stored for each site is in the seeds_est_site d.f.
  init_seed_bank <- deframe(round(seeds_est_site[seeds_est_site$Site == site_now,"seedbank_est"]))
  # Use seedling data to calculate the starting number of seedlings 
  # seedling data stored in the 'seedling_site' d.f.
  init_seedlings <- deframe(round(seedlings_site[seedlings_site$Site == site_now, "Seedlings_t"],0))
  # starting number of individuals in the continuous stage
  ht <- dat$log_LL_t
  init_pop_vec   <- pnorm(rnorm(500, mean = mean(ht, na.rm = TRUE), sd = sd(ht, na.rm = TRUE)), mean = mean(ht, na.rm = TRUE), sd = sd(ht, na.rm = TRUE))
  #discretize_pop_vector(trait_values = ht, n_mesh = 500, pad_low = 1, pad_high = 1)$n_ht #, pad_low = 0.8, pad_high = 1.2)$n_ht
  
  ## code to run the IPM
  temp_ipm <- temp_ipm %>%
    define_domains(
      # We can pass the variables we created above into define_domains
      ht = c(L, U, n)
    ) %>%
    define_pop_state(
      # We can also pass them into define_pop_state
      pop_vectors = list(
        n_ht = init_pop_vec,
        n_b  = init_seed_bank,
        n_s  = init_seedlings 
      )
    ) %>%
    define_env_state(
      env_covs   = sample_env(env_params),
      data_list  = list(env_params = env_params,
                        sample_env = sample_env)) %>%
    make_ipm(iterate = TRUE, iterations = 100,
             normalize_pop_size = FALSE,
             usr_funs = list(inv_logit_r   = inv_logit_r), return_main_env = TRUE )
  
  ## rename the ipm object to have the name of site
  assign(paste0(site_now,"__stoch_DD_ipm"), value = temp_ipm)
  
  ## Done!
}

#### Simple, deterministic IPM with bootstrapped estimates of uncertainty ####
## initial IPM is called det_ipm

# make a vector to hold all of the lambdas from the resampling runs
all_lambdas <- numeric(200L)
# make a list to hold the model parameters
all_params <- list()
# get the proto-imp from the det_ipm model
use_proto <- det_ipm$proto_ipm
# Now, we refit the vital rate models, and use the parameters<- setter function to update the original proto_ipm object with the new vital rate models. This saves us from re-typing the whole model pipeline again. Normally, we would bootstrap more than 50 times, but for the sake of this example and saving time, we will only do 50.
for(i in 1:200) {
  ## sample continuous data
  sample_ind <- seq(1, nrow(dat), by = 1)
  
  boot_data_ind   <- sample(sample_ind, size = nrow(dat), replace = TRUE)
  
  boot_data <- dat[boot_data_ind,]
  
  ## sample discrete stage data
  sample_disc <- seq(1, nrow(discreteDat), by = 1)
  
  boot_disc_ind <- sample(sample_disc, size = nrow(discreteDat), replace = TRUE)
  
  boot_disc <- discreteDat[boot_disc_ind,]
  
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
  recD_now <- boot_data[boot_data$age == 0 & is.na(boot_data$age) == FALSE,]
  # fit the model
  recMod_now <- lm(log_LL_t ~ 1, data = recD_now)
  ## Probability of a seedling in year *t* establishing to a rosette in year *t+1* ($p_{estab}$)
  p.estab.est_now <- sum(boot_disc$Recruit_tplus1)/sum(boot_disc$Seedling_t)
  ## Probability that a seed from the seedbank in year t will germinate to a seedling in year t+1
  outSB.est_now <- sum(boot_disc$Seedling_tplus1)/sum(boot_disc$SeedBank_t)
  ## Probability that a seed from the seedbank in year t will stay in the seedbank in year t+1
  staySB.est_now <- sum(boot_disc$SeedBank_tplus1)/sum(boot_disc$SeedBank_t)
  ## Probability that a seed produced by an adult plant in year t will enter the seedbank in year t+1
  goSB.est_now <- sum(boot_disc[boot_disc$NewSeeds_t == 1, "SeedBank_tplus1"])/sum(boot_disc$NewSeeds_t)
  ## Probability that a seed from a plant in year t will go directly to the seedling stage
  goSdlng.est_now <- sum(boot_disc[boot_disc$NewSeeds_t == 1, "Seedling_tplus1"])/sum(boot_disc$NewSeeds_t)
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
    goSdlng   = goSdlng.est_now, 
    staySB = staySB.est_now, 
    goSB = goSB.est_now, 
    outSB = outSB.est_now, 
    p_estab = p.estab.est_now 
  )
  
  # Insert the new vital rate models into the proto_ipm, then rebuild the IPM.
  
  parameters(use_proto) <- data.list_now
  
  boot_ipm <- use_proto %>%
    make_ipm(iterate = TRUE,
             normalize_pop_size = FALSE,
             iterations = 100)
  
  all_lambdas[i] <- lambda(boot_ipm)
  
  all_params[[i]] <- data.list_now
}

# Plot the results
param_fig <- data.frame("param" = "g_int", 
                        "value" = sapply(X = all_params, FUN = function (x) x[[1]]))
param_fig <- rbind(param_fig,  data.frame("param" = "g_slope", 
                                          "value" = sapply(X = all_params,
                                                           FUN = function (x) x[[2]])))
param_fig <- rbind(param_fig, data.frame("param" = "g_sd", 
                                         "value" = sapply(X = all_params,
                                                          FUN = function (x) x[[3]])))
param_fig <- rbind(param_fig, data.frame("param" = "s_int", 
                                         "value" = sapply(X = all_params,
                                                          FUN = function (x) x[[4]])))
param_fig <- rbind(param_fig, data.frame("param" = "s_slope",
                                         "value" = sapply(X = all_params,
                                                          FUN = function (x) x[[5]])))
param_fig <- rbind(param_fig, data.frame("param" = "p_b_int", 
                                         "value" = sapply(X = all_params,
                                                          FUN = function (x) x[[6]])))
param_fig <- rbind(param_fig, data.frame("param" = "p_b_slope", 
                                         "value" = sapply(X = all_params,
                                                          FUN = function (x) x[[7]])))
param_fig <- rbind(param_fig, data.frame("param" = "p_b_slope_2", 
                                         "value" = sapply(X = all_params,
                                                          FUN = function (x) x[[8]])))
param_fig <- rbind(param_fig, data.frame("param" = "b_int", 
                                         "value" = sapply(X = all_params,
                                                          FUN = function (x) x[[9]])))
param_fig <- rbind(param_fig, data.frame("param" = "b_slope", 
                                         "value" = sapply(X = all_params,
                                                          FUN = function (x) x[[10]])))
param_fig <- rbind(param_fig, data.frame("param" = "c_o_mu", 
                                         "value" = sapply(X = all_params,
                                                          FUN = function (x) x[[11]])))
param_fig <- rbind(param_fig, data.frame("param" = "c_o_sd", 
                                         "value" = sapply(X = all_params,
                                                          FUN = function (x) x[[12]])))
param_fig <- rbind(param_fig, data.frame("param" = "goSdlng", 
                                         "value" = sapply(X = all_params,
                                                          FUN = function (x) x[[13]])))
param_fig <- rbind(param_fig, data.frame("param" = "staySB", 
                                         "value" = sapply(X = all_params,
                                                          FUN = function (x) x[[14]])))
param_fig <- rbind(param_fig, data.frame("param" = "goSB", 
                                         "value" = sapply(X = all_params,
                                                          FUN = function (x) x[[15]])))
param_fig <- rbind(param_fig, data.frame("param" = "outSB", 
                                         "value" = sapply(X = all_params,
                                                          FUN = function (x) x[[16]])))
param_fig <- rbind(param_fig, data.frame("param" = "p_estab", 
                                         "value" = sapply(X = all_params,
                                                          FUN = function (x) x[[17]])))
# calculate the mean values
test <- param_fig %>% 
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
plot(density(all_lambdas))
abline(v = lambda(det_ipm), col = 'red', lwd = 2, lty = 2)
abline(v = mean(all_lambdas), col = "blue", lwd = 2, lty = 2)

#### store the ipm results

save.image(file = "./analysis_scripts/ipm_results.RData")
