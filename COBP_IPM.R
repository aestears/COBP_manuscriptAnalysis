#//////////////////////////
# Integral Projection Models for Oenothera coloradensis
# Alice Stears
# 30 November 2021
#/////////////////////////

#### Load packages ####
library(tidyverse)
library(ipmr)

#### Load Data ####
## load spatialized data, which we'll need for neighborhood calculations
dat <- read.csv("../Raw Data/COBP_long_CURRENT.csv")

## transform the necessary variables
## make log-transformed size variables ####
dat$log_LL_t <- log(dat$LongestLeaf_cm)
dat$log_LL_tplus1 <- log(dat$longestLeaf_tplus1)
dat$log_LL_tminus1 <- log(dat$longestLeaf_tminus1)

## round capsule count numbers (that were modeled based on regression) to whole numbers 
dat$Num_capsules <- round(dat$Num_capsules, digits = 0)

## transform the capsule counts to seed counts (assume 4 seeds/capsule--from Burgess, 2005)
dat$Num_seeds <- dat$Num_capsules*4

## other potential covariates from the data
# Invertebrate leaf herbivory, stem herbivory, leaf spots (unknown source?)

#### Deterministic, non-density-dependent IPM for all data ####
### Make vital rate models 
## Survival ($s(z)$)
# subset the data to exclude flowering individuals
survDat <- dat[dat$flowering==0 | is.na(dat$flowering),]
# logistic glm with log-transformed size_t
survMod <- glm(survives_tplus1 ~ log_LL_t , data = survDat, family = binomial)
summary(survMod)
# plot model results 
plot(survives_tplus1 ~ log_LL_t, data = survDat)
newdata <- data.frame("log_LL_t" = seq(from = min(survDat$log_LL_t, na.rm = TRUE), 
                                       to = max(survDat$log_LL_t, na.rm = TRUE),
                                       length.out = 100))
lines(x = newdata$log_LL_t, y = predict(object = survMod, newdata =  newdata, type = "response"), col = "red")

## Growth ($G(z',z)$)
# lm w/ log-transformed size_t and size_t+1
sizeMod <- lm(log_LL_tplus1 ~ log_LL_t , data = dat)
summary(sizeMod)
# plot model results
plot(log_LL_tplus1 ~ log_LL_t, data = dat)
newdata <- data.frame("log_LL_t" = seq(from = min(dat$log_LL_t, na.rm = TRUE), 
                                       to = max(dat$log_LL_t, na.rm = TRUE),
                                       length.out = 100))
lines(x = newdata$log_LL_t, y = predict(object = sizeMod, newdata =  newdata), col = "red")


## Number of seeds produced, according to plant size ($b(z)$)
# using size in current year (no. of seeds/plant, for those that flowered ~ size_t)
# use only plants that flowered 
seedDat <- dat[dat$flowering==1,]
# fit poisson glm (for count data)
seedMod_t <- glm(Num_seeds ~ log_LL_t , data = seedDat, family = poisson)
summary(seedMod_t)
# plot model results
plot(Num_seeds ~ log_LL_t, data = seedDat)
newdata <- data.frame("log_LL_t" = seq(from = min(seedDat$log_LL_t, na.rm = TRUE), 
                                       to = max(seedDat$log_LL_t, na.rm = TRUE),
                                       length.out = 100))
lines(x = newdata$log_LL_t, y = predict(object = seedMod_t, newdata =  newdata, type = "response"), col = "red")

## Flowering probability ($p_b(z)$)
# using size in current year (w/ squared term)
# logistic glm with log-transformed size_t
flwrMod_t <- glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = dat, family = binomial)
summary(flwrMod_t)
# plot model results 
plot(flowering ~ log_LL_t, data = dat)
newdata <- data.frame("log_LL_t" = seq(from = min(dat$log_LL_t, na.rm = TRUE), 
                                       to = max(dat$log_LL_t, na.rm = TRUE),
                                       length.out = 100))
lines(x = newdata$log_LL_t, y = predict(object = flwrMod_t, newdata =  newdata, type = "response"), col = "red")

## Distribution of recruit size ($c_o(z')$)
# subset the data
recD <- dat[dat$age == 0 & is.na(dat$age) == FALSE,]
# fit the model
recMod <- lm(log_LL_t ~ 1, data = recD)
summary(recMod)
#plot the results
hist(recD$log_LL_t)
abline(v = recMod$coefficients, col = "blue", lwd = 2)

## Probability of a seedling in year *t* establishing to a rosette in year *t+1* ($p_{estab}$)
# make a column in 'dat' that labels 'recruits' to the rosette stage
dat$recruit <- 0
dat[dat$age==0 & is.na(dat$age)==FALSE,"recruit"] <- 1
# get the number of recruits/year 
estabTemp <- dat %>% 
  dplyr::select(Plot_ID, Year, recruit) %>% 
  group_by(Plot_ID, Year) %>% 
  summarize( recruits = sum(recruit)) %>% 
  rename( recruits_t = recruits)
# change the 2018 recruit values to NA, since we don't have a count of recruits to the rosette stage for that year
estabTemp[estabTemp$Year==2018, "recruits_t"] <- NA
# get the number of recruits in year t+1
estabTemp$Year <- estabTemp$Year - 1
names(estabTemp)[3] <- "recruits_tplus1"
# get the no.of seedlings/year
# load the seedling data
seedlings <- read.csv("../Raw Data/COBP_seedlings_8_23_21.csv") %>% 
  dplyr::select(Plot_ID, Seedlings_18, Seedlings_19, Seedlings_20) %>% 
  group_by(Plot_ID) %>% 
  summarize(Seedlings_18 = sum(Seedlings_18), Seedlings_19 = sum(Seedlings_19), Seedlings_20 = sum(Seedlings_20)) %>% 
  pivot_longer(cols = c(Seedlings_18, Seedlings_19, Seedlings_20), names_to = "Year", values_to = "Seedlings_t", names_pattern = "([[:digit:]]+)") %>% 
  mutate(Year = (as.numeric(Year) + 2000))
# combine seedling and rosette recruit data
estabs <- left_join(estabTemp, seedlings)
# calculate the probability of seedling in year t establishing to a rosette in year t+1
estabs$P_estab <- estabs$recruits_tplus1/estabs$Seedlings_t
# fix the 'inf' values
estabs[is.infinite(estabs$P_estab),"P_estab"] <- NA
# if the value is >1, round back to one (must have missed some seedlings)
estabs[estabs$P_estab > 1 & is.na(estabs$P_estab) == FALSE,"P_estab"] <- 1
estabs <- as.data.frame(estabs)
# calculate the probability of establishment value
p.estab.est <- sum(estabs$P_estab, na.rm = TRUE)/sum(is.na(estabs$P_estab)==FALSE)

## Calculate germination rate and rate of seed viability loss
# (data in SeedBagGreenhouseSeedlings.csv)--seeds from previous year
germ.rt.ours <- .03
#data from (Burgess, Hild & Shaw, 2005)--seedbank seed viability/germination rate doesn't seem to change much over time
germ.rt.Burgess <- mean(c(16.0, 13.0, 12, 8.3, 7.0, 5.3)/45)
germ.rt <- 0.08 #average together, but trending toward our estimate
# seed viability rate (from Burgess, Hild & Shaw, 2005)
viab.rt.Burgess <- mean(c(.79,.7,.63))
viab.rt <- viab.rt.Burgess

## Probability that a seed from the seedbank in year t will germinate to a seedling in year t+1 ($outSB$)
outSB.est <-germ.rt

## Probability that a seed from the seedbank in year t will stay in the seedbank in year t+1 ($staySB$)
# seed viability rate * (1-germination rate)
staySB.est <- (1-germ.rt)*viab.rt

## Probability that a seed produced by an adult plant in year t will enter the seedbank in year t+1 ($goSB$)
goSB.est <- (1 - germ.rt)*viab.rt

## Probability that a seed from a plant in year t will go directly to the seedling stage ($goSdlng$)
# use the germination rate, since it doesn't seem to change much with age (Burgess, Hild & Shaw, 2005)
goSdlng.est <- germ.rt

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
general_ipm <- init_ipm(sim_gen = "general", # make a general IPM
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
  ) 
## We’ve now defined all of the kernels, next are the implementation details. These also differ somewhat from simple IPMs. The key difference in the implementation arguments list lies in the state_start and state_end of each kernel, and is related to the family argument of each kernel. Kernels that begin with one state and end in a different state (e.g. moving from seed bank to a plant) will have different entries in the state_start and state_end slots. It is very important to get these correct, as ipmr uses this information to generate the model iteration procedure automatically (i.e. code corresponding to Equations 1-2).

## define the starting and ending states for each kernel
general_ipm <- general_ipm %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P", "leave_seedlings", "repro_to_seedlings", "seedbank_to_seedlings", "stay_seedbank", "repro_to_seedbank"),
      int_rule     = c(rep("midpoint", 6)),
      state_start    = c('ht', "s", "ht", "b", "b", "ht"),
      state_end      = c("ht", "ht", "s", "s", "b", "b")
    )
  )


##Actually run the IPM!
  
# lower limit of size
L <- 1.09
# upper limit of size
U <- 3.66
# number of break-points
n <- 500

set.seed(2312)

## Use data from seedbank to calculate the starting number of seeds (from COBPGreenhouseSeedbankStudyData.csv)
# the no. of seeds that germinated for sample of each plot ()
seedStudyGerms <- c(5,4,3,0.01,0.01,1,0.01,0.01,0.01,9,4,9,0,1,1,8,20,1)
# divide the seedcore germs by the germ rate to get the number of (seedbank + seed rain) seeds in each core
seedStudyGerms.1 <- seedStudyGerms/germ.rt
# convert the volume of the soil samples to the volume of the plots to get approx. number of (seedbank + seed rain) seeds in each plot
# cores are 3cm in diameter and 2.5cm deep %%% NEED TO ACTUALLY CHECK THIS; THIS IS JUST A PLACEHOLDER %%%, and there are 20 of them/plot
coreVol = ((pi*(3/2)^2)*2.5)*20 ## 753.98 cm^3 is the total volume of the soil cores
# calculate the volume of the first 3cm of soil in the plots
plotVol = 200*200*2.5
# estimate the number of seeds in the plot
seedsEst.all <- (plotVol*seedStudyGerms.1)/coreVol
# total number of seeds produced in each plot in each year
seedsTemp <- aggregate(x = dat[,c("Num_seeds")], by = list("Plot_ID" = dat$Plot_ID,
                                                           "Year" = dat$Year), FUN = sum, na.rm = TRUE)
seedRain.avg <- aggregate(x = seedsTemp$x, by = list("Plot_ID" = seedsTemp$Plot_ID), FUN = mean)
# re-order the plot seedRain data to be in the same order as the seed core data
seedRain.avg <- seedRain.avg[c(7:18,1:3,6,4:5),]
# subtract the seedRain from the seedcore seed est to get an estimate of the size of the seedbank state in year t (B(t))
# formula = (germs*(1-germ.rate) - mean no. of seeds produced in previous year)
seedBank_est <- mean(c(seedsEst.all - seedRain.avg$x))

# starting size of the seedbank
init_seed_bank <- seedBank_est
# Use seedling data to calculate the starting number of seedlings
init_seedlings <- round(sum(seedlings$Seedlings_t)/length(seedlings$Seedlings_t),0)
# starting number of individuals in the continuous stage
init_pop_vec   <- runif(500)

## code to run the IPM
general_ipm <- general_ipm %>%
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
           usr_funs = list(inv_logit   = inv_logit,
                           inv_logit_2 = inv_logit_2), return_main_env = TRUE )

## lambda is a generic function to compute per-capita growth rates. It has a number of different options depending on the type of model
lambda(general_ipm)

## If we are worried about whether or not the model converged to stable dynamics, we can use the exported utility is_conv_to_asymptotic. The default tolerance for convergence is 1e-10, but can be changed with the 'tol' argument.
is_conv_to_asymptotic(general_ipm, tol = 1e-10)
## additional calculations
lambda_ipmr <- lambda(general_ipm)
repro_value <- left_ev(general_ipm)
stable_dist <- right_ev(general_ipm)

### Visualize the IPM kernel
# first have to make a mega-kernel
mega_mat <- make_iter_kernel(ipm = general_ipm, 
                             mega_mat = c(stay_seedbank, 0, repro_to_seedbank, seedbank_to_seedlings, 0, repro_to_seedlings, 0, leave_seedlings, P),
)
# check to make sure I constructed the mega-kernel correctly (should be nearly equal values)
Re(eigen(mega_mat[[1]])$values[1]) - lambda(general_ipm)

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
      t(make_iter_kernel(ipm = general_ipm, mega_mat = c(P))$mega_matrix),
      xlab = "n(z,t)", ylab = "n(z',t+1)", main = "P matrix")
contour(x = meshpts,
        y = meshpts,
        t(make_iter_kernel(ipm = general_ipm, mega_mat = c(P))$mega_matrix), 
        add = TRUE, drawlabels = TRUE, nlevels = 10, col = "grey30")

#### Deterministic, density-dependent IPM for all data ####
## get the environmental data
source("./data_prep_scripts/prepping_Envi_Data.R")
# called "Climate", "soilTemp_plot", and "soilMoist"
## add the envi data to the 'dat' data.frame
# add soil moisture data (have that for every plot, from one time point)
dat <- left_join(dat, soilMoist[,c("Site", "SoilMoisture_m3m3")], by = c("Plot_ID" = "Site"))
# add soil temperature data (have mean values for every site, from one time point)
# reformat the data to have means for each site 
soilTemp <- soilTemp_plot %>% 
  group_by(Site, Location) %>% 
  summarize(SoilTemp_winter_C = mean(SoilTemp_winter_C, na.rm = TRUE),
            sd_soilTemp_winter = mean(sd_soilTemp_winter, na.rm = TRUE),
            SoilTemp_grow_C = mean(SoilTemp_grow_C, na.rm = TRUE), 
            sd_soilTemp_grow = mean(sd_soilTemp_grow, na.rm = TRUE))
dat <- left_join(dat, soilTemp)
# add the climate data (have mean values for each location)
dat <- left_join(dat, Climate)

# reformat variables
dat$Year <- as.factor(dat$Year)
dat$SoilMoisture_m3m3_s <- scale(dat$SoilMoisture_m3m3)
dat$SoilTemp_winter_C_s <- scale(dat$SoilTemp_winter_C)
dat$SoilTemp_grow_C_s <- scale(dat$SoilTemp_grow_C)
dat$tMean_grow_C_s <- scale(dat$tMean_grow_C)
dat$precipWaterYr_cm_s <- scale(dat$precipWaterYr_cm)
### fit the vital rate models including environmental variables

## Survival ($s(z)$)
# subset the data to exclude flowering individuals
survDat <- dat[dat$flowering==0 | is.na(dat$flowering),]
# logistic glm with log-transformed size_t
survMod_e <- glmer(survives_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s + SoilTemp_winter_C_s + SoilTemp_grow_C_s + tMean_grow_C_s + precipWaterYr_cm_s + (1|Plot_ID), data = survDat, family = binomial)
summary(survMod_e)
survMod_e_1 <- glmer(survives_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s + SoilTemp_grow_C_s + tMean_grow_C_s + (1|Plot_ID),  data = survDat, family = binomial)
summary(survMod_e_1)
## survMod_e_1 is the better fit, use this model
survMod_env <- survMod_e_1

## Growth ($G(z',z)$)
# lm w/ log-transformed size_t and size_t+1
sizeMod_e <- lmer(log_LL_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s + SoilTemp_winter_C_s + SoilTemp_grow_C_s + tMean_grow_C_s + precipWaterYr_cm_s + (1|Plot_ID), data = dat)
summary(sizeMod_e)
sizeMod_e_1 <- lmer(log_LL_tplus1 ~ log_LL_t + SoilMoisture_m3m3_s +  tMean_grow_C_s + precipWaterYr_cm_s + (1|Plot_ID), data = dat)
summary(sizeMod_e_1)
sizeMod_env <- sizeMod_e_1

## Number of seeds produced, according to plant size ($b(z)$)
# using size in current year (no. of seeds/plant, for those that flowered ~ size_t)
# use only plants that flowered 
seedDat <- dat[dat$flowering==1,]
# fit poisson glm (for count data)
seedMod_e <- glmer(Num_seeds ~ log_LL_t + SoilMoisture_m3m3_s + SoilTemp_winter_C_s + SoilTemp_grow_C_s + tMean_grow_C_s + precipWaterYr_cm_s + (1|Plot_ID), data = seedDat, family = poisson)
summary(seedMod_e)
seedMod_e_1 <- glmer(Num_seeds ~ log_LL_t + tMean_grow_C_s + precipWaterYr_cm_s + (1|Plot_ID) , data = seedDat, family = poisson)
summary(seedMod_e_1)
seedMod_env <- seedMod_e_1

## Flowering probability ($p_b(z)$)
# using size in current year (w/ squared term)
# logistic glm with log-transformed size_t
flwrMod_e <- glmer(flowering ~ log_LL_t + I(log_LL_t^2) + SoilMoisture_m3m3_s + SoilTemp_winter_C_s + SoilTemp_grow_C_s + tMean_grow_C_s + precipWaterYr_cm_s + (1|Plot_ID), data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(flwrMod_e)
flwrMod_e_1 <- glmer(flowering ~ log_LL_t + I(log_LL_t^2) + SoilMoisture_m3m3_s + SoilTemp_grow_C_s + tMean_grow_C_s + precipWaterYr_cm_s + (1|Plot_ID), data = dat, family = binomial, control = glmerControl(optimizer = "bobyqa"))
summary(flwrMod_e_1)
flwrMod_env <- flwrMod_e_1

## Distribution of recruit size ($c_o(z')$)
# subset the data
recD <- dat[dat$age == 0 & is.na(dat$age) == FALSE,]
# fit the model
recMod <- lm(log_LL_t ~ 1 , data = recD)
summary(recMod)
#plot the results

## Probability of a seedling in year *t* establishing to a rosette in year *t+1* ($p_{estab}$)
p.estab.est

## Probability that a seed from the seedbank in year t will germinate to a seedling in year t+1 ($outSB$)
outSB.est 

## Probability that a seed from the seedbank in year t will stay in the seedbank in year t+1 ($staySB$)
staySB.est 

## Probability that a seed produced by an adult plant in year t will enter the seedbank in year t+1 ($goSB$)
goSB.est

## Probability that a seed from a plant in year t will go directly to the seedling stage ($goSdlng$)
# use the germination rate, since it doesn't seem to change much with age (Burgess, Hild & Shaw, 2005)
goSdlng.est 

### define the parameter list
fixed_list <- list(
  g_int     = fixef(sizeMod_env)[1], # growth model intercept
  g_slope   = fixef(sizeMod_env)[2], # growth model slope
  g_sd      = sd(resid(sizeMod_env)), #growth model sd
  s_int     = fixef(survMod_env)[1], # survival model intercept
  s_slope   = fixef(survMod_env)[2], # survival model slope
  p_b_int   = fixef(flwrMod_env)[1], #probability of flowering
  p_b_slope = fixef(flwrMod_env)[2],
  p_b_slope_2 = fixef(flwrMod_env)[3],
  b_int   = fixef(seedMod_env)[1], #seed production
  b_slope = fixef(seedMod_env)[2],
  c_o_mu    = coef(recMod), #recruit size distribution
  c_o_sd    = summary(recMod)$sigma,
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

# name the list of random effects
nms <- paste("r_", 1:18, sep = "")

names(g_r_int)   <- paste('g_', nms, sep = "")
names(s_r_int)   <- paste('s_', nms, sep = "")
names(p_b_r_int) <- paste('p_b_', nms, sep = "")
names(b_r_int) <- paste('b_', nms, sep = "")

# Each set of parameters is converted to a named list. The names should match
# the variables referenced in each define_kernel() call.
g_params   <- as.list(g_r_int)
s_params   <- as.list(s_r_int)
p_b_params <- as.list(p_b_r_int)
b_r_params <- as.list(b_r_int)

# add them all together using c()
all_params_list <- c(fixed_list, g_params, s_params, p_b_params, b_r_params)

### define internal functions
inv_logit_r <- function(int, slope, sv, r_eff) {
  1/(1 + exp(-(int + slope * sv + r_eff)))
}

inv_logit_2_r <- function(int, slope, slope_2, sv, r_eff) {
  1/(1 + exp(-(int + slope * sv + slope_2 * sv ^ 2 + r_eff)))
}

my_funs <- list(inv_logit_r = inv_logit_r,
                inv_logit_2_r = inv_logit_2_r)

### define the IPM
better_ipm <- init_ipm(sim_gen = "general", # make a general IPM
                        di_dd = "di", # make it density independent
                        det_stoch = "stoch",
                        kern_param = "kern") %>% # make it deterministic
  define_kernel(
    name          = "P_plot", # survival 
    # We add d_ht to formula to make sure integration is handled correctly.
    # This variable is generated internally by make_ipm(), so we don't need
    # to do anything else.
    formula       = (1-p_b_plot) * s_plot * g_plot * d_ht,
    family        = "CC",
    g_plot        = dnorm(ht_2, mean = g_mu, sd = g_sd),
    g_mu          = g_int + g_slope * ht_1 + g_r_plot,
    s_plot        = inv_logit_r(int = s_int, slope = s_slope, sv = ht_1, r_eff = s_r_plot),
    p_b_plot      = inv_logit_2_r(int = p_b_int, slope = p_b_slope, 
                                  slope_2 = p_b_slope_2, sv = ht_1, r_eff = p_b_r_plot),
    data_list     = all_params_list,
    states        = list(c('ht')),
    uses_par_sets = TRUE,
    par_set_indices = list(plot = 1:18),
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions('norm', 'g_plot')
  ) %>%
  define_kernel(
    name          = "leave_seedlings", ## leave seedling stage and go to rosette stage
    formula       = p_estab. * c_o. * d_ht,
    family        = 'DC', # Note that now, family = "DC" because it denotes a discrete -> continuous transition
    p_estab.      = p_estab,
    c_o.          = dnorm(ht_2, c_o_mu, c_o_sd),
    data_list     = all_params_list,
    states        = list(c('ht', "s")),   # Note that here, we add "s" to our list in states, because this kernel uses seedlings 
    uses_par_sets = FALSE,
    evict_cor     = TRUE,
    evict_fun     = truncated_distributions('norm','c_o.')
  ) %>%
  define_kernel(
    name    = "repro_to_seedlings_plot",
    formula       = (goSdlng.) * (p_b_plot * b_plot * d_ht),
    family        = "CD",
    goSdlng.      = goSdlng,
    p_b_plot      = inv_logit_2_r(int = p_b_int, slope = p_b_slope, slope_2 = p_b_slope_2, 
                                  sv = ht_1, r_eff = p_b_r_plot),
    b_plot        = exp(b_int + b_slope * ht_1 + b_r_plot),
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
    p_b_plot      = inv_logit_2_r(int = p_b_int, slope = p_b_slope, slope_2 = p_b_slope_2, 
                                  sv = ht_1, r_eff = p_b_r_plot),
    b_plot        = exp(b_int + b_slope * ht_1 + b_r_plot),
    data_list     = all_params_list,
    states        = list(c('b', 'ht')),
    uses_par_sets = TRUE,
    par_set_indices = list(plot = 1:18),
    evict_cor = FALSE
  ) 

## define the starting and ending states for each kernel
better_ipm <- better_ipm %>%
  define_impl(
    make_impl_args_list(
      kernel_names = c("P_plot", "leave_seedlings", "repro_to_seedlings_plot", "seedbank_to_seedlings", "stay_seedbank", "repro_to_seedbank_plot"),
      int_rule     = c(rep("midpoint", 6)),
      state_start    = c('ht', "s", "ht", "b", "b", "ht"),
      state_end      = c("ht", "ht", "s", "s", "b", "b")
    )
  )


##Actually run the IPM!
# lower limit of size
L <- 1.09
# upper limit of size
U <- 3.66
# number of break-points
n <- 500
set.seed(2312)
## Use data from seedbank to calculate the starting number of seeds (From general_ipm)
#init_seed_bank 
# Use seedling data to calculate the starting number of seedlings
#init_seedlings 
# starting number of individuals in the continuous stage
#init_pop_vec   <- runif(500)

## code to run the IPM
better_ipm <- better_ipm %>%
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
  make_ipm(iterate = TRUE, iterations = 100,
           usr_funs = my_funs, 
           kernel_seq = sample(1:18, 100, replace = TRUE)
           )

## lambda is a generic function to compute per-capita growth rates. It has a number of different options depending on the type of model
lambda(general_ipm)

## If we are worried about whether or not the model converged to stable dynamics, we can use the exported utility is_conv_to_asymptotic. The default tolerance for convergence is 1e-10, but can be changed with the 'tol' argument.
is_conv_to_asymptotic(general_ipm, tol = 1e-10)
## additional calculations
lambda_ipmr <- lambda(general_ipm)
repro_value <- left_ev(general_ipm)
stable_dist <- right_ev(general_ipm)

### Visualize the IPM kernel
# first have to make a mega-kernel
mega_mat <- make_iter_kernel(ipm = general_ipm, 
                             mega_mat = c(stay_seedbank, 0, repro_to_seedbank, seedbank_to_seedlings, 0, repro_to_seedlings, 0, leave_seedlings, P),
)
# check to make sure I constructed the mega-kernel correctly (should be nearly equal values)
Re(eigen(mega_mat[[1]])$values[1]) - lambda(general_ipm)

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
image(x = meshpts, 
      y = meshpts,
      t(make_iter_kernel(ipm = general_ipm, mega_mat = c(P))$mega_matrix),
      xlab = "n(z,t)", ylab = "n(z',t+1)", main = "P matrix")
contour(x = meshpts,
        y = meshpts,
        t(make_iter_kernel(ipm = general_ipm, mega_mat = c(P))$mega_matrix), 
        add = TRUE, drawlabels = TRUE, nlevels = 10, col = "grey30")
