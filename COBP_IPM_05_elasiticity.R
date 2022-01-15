#//////////////////////////
# Integral Projection Models for Oenothera coloradensis: Elasticity Analysis
# Alice Stears
# 11 December 2021
#/////////////////////////

# load required packages
library(tidyverse)
library(ipmr)
#### load vital rate models from previous script ####
load(file = "./analysis_scripts/ipm_results.RData")
source(file = "./analysis_scripts/COBP_IPM_04_ellnerCodeSimpleIPMs.R")

# lower limit of size in the ipm
L = min(dat$log_LL_t, na.rm = TRUE) * .8
# upper limit of size in the ipm
U = max(dat$log_LL_t, na.rm = TRUE) * 1.2
# number of break-points
n <- 500
meshpts <- seq(from = L, to = U, length.out = 500)
h = meshpts[2] - meshpts[1]

#### for simple IPM from ellner code (called 'IPM.true') ####
(lam <- Re(eigen(IPM.true$K)$values[1])) # 1.417
w.eigen_temp <- Re(eigen(IPM.true$K)$vectors[,1]) 
stable.dist_ellner <- w.eigen_temp/sum(w.eigen_temp) 
v.eigen_temp <- Re(eigen(t(IPM.true$K))$vectors[,1]) 
repro.val_ellner <- v.eigen_temp/v.eigen_temp[1]

#### for simple IPM from impr code ####
simple_ipm_iterK <- make_iter_kernel(ipm = simple_ipm, name_ps = "P", f_forms = "F")$mega_matrix
# eigen things using ipmr functions
lambda(simple_ipm) # 1.417
w.eigen_simple <- ipmr::right_ev(ipm = simple_ipm)$size_w
v.eigen_temp2 <- ipmr::left_ev(ipm = simple_ipm)$size_v
v.eigen_simple <- v.eigen_temp2/v.eigen_temp2[1]
# eigen things using iteration matrix
Re(eigen(simple_ipm_iterK)$values[1]) # 1.417
w.eigen_simple.hand <- Re(eigen(simple_ipm_iterK)$vectors[,1])
stable.dist_simple.hand <- w.eigen_simple.hand/sum(w.eigen_simple.hand) 
v.eigen_simple.hand <- Re(eigen(t(simple_ipm_iterK))$vectors[,1])
repro.val_sipmle.hand=v.eigen_simple.hand/v.eigen_simple.hand[1] 

# compare stable age distributions
plot(x = meshpts, y = stable.dist_simple.hand, type = 'l', ylim = c(0,.015))
lines(x = meshpts, y =  w.eigen_simple, col = "red")
lines(x = meshpts, y = stable.dist_ellner, col = "orange")

# compare reproductive value
plot(x = meshpts, y = repro.val_sipmle.hand, type = "l", ylim = c(0,300))
lines(x = meshpts, y = v.eigen_simple, col = "red")
lines(x = meshpts, y = repro.val_ellner, col = "orange")

# The eigen-things can be combined to obtain the sensitivity and elasticity matrices. 
# calculate the dot*product of the left and right eigen vectors
v.dot.w <- sum(w.eigen_simple*v.eigen_simple)*h #(where h = step size; size btwn each mesh point)
# compute the sensitivity matrix
sens <- outer(v.eigen_simple,w.eigen_simple)/v.dot.w
# compute the elasticity matrix
elas <- matrix(as.vector(sens)*as.vector(simple_ipm_iterK)/lambda(simple_ipm),nrow=500)
image(sens)
image(elas)

# calculate the damping ratio (to determine if the asymptotic rates are good approximations of transient dynamics)
(lam_simple <- Re(eigen(simple_ipm_iterK)$values[1])) # asymptotic growth rate
(damp_simple=Re(eigen(simple_ipm_iterK)$values[1])/Re(eigen(simple_ipm_iterK)$values[2])) # damping ratio
# 5.65; indicating that the dominant eigenvalue (lambda) is dictating transient dynamics

#### general, deterministic, density-independent IPM ####
## compute the eigenvectors
# right eigenvector
w.eigen_det.di <- ipmr::right_ev(ipm = det_ipm, iterations = 100)
stable.size.dist_det.di <- w.eigen_det.di$ht_w/sum(w.eigen_det.di$ht_w) # convert to a probability density function
plot(x = meshpts, y =  stable.size.dist_det.di, type = 'l', xlab = "meshpts: (ln(longest leaf cm))", ylab = "probability")

# left eigenvector
v.eigen_det.di <- (ipmr::left_ev(ipm = det_ipm, iterations = 100))
repro.val_det.di <- v.eigen_det.di$ht_v/(v.eigen_det.di$ht_v)[1]
plot(x = meshpts, y = repro.val_det.di, type = 'l', xlab = "meshpts: (ln(longest leaf cm))", ylab = "probability")

# lambda
lambda_det.di <- lambda(det_ipm)

# The eigen-things can be combined to obtain the sensitivity and elasticity matrices.
# compute elasticity and sensitivity matrices
v.dot.w_det.di <- sum(stable.size.dist_det.di *repro.val_det.di)*h
mega.mat_det.di_proto <- format_mega_kernel(ipm = det_ipm, 
                                      mega_mat = c(stay_seedbank, 0, repro_to_seedbank, seedbank_to_seedlings, 0, repro_to_seedlings, 0, leave_seedlings, P))
mega.mat_det.di <- mega.mat_det.di_proto$mega_matrix[3:502,3:502]
# calculate the sensitivity function (whole kernel)
sens_det.di <- outer(repro.val_det.di,stable.size.dist_det.di, '*')/(v.dot.w_det.di)
# calculate the elasticity function (whole kernel)
elas_det.di <- matrix(as.vector(sens_det.di)*as.vector(mega.mat_det.di)/lambda_det.di,nrow=500)
# plot the sensitivity function of the entire continuous kernel
image(x = meshpts, y = meshpts, t(sens_det.di)^.1, 
      xlab = "ln(leaf) in year t", ylab = "ln(leaf) in year t+1")
contour(x = meshpts, y = meshpts, t(sens_det.di), add = TRUE)
# plot the elasticity function of the entire continuous kernel
image(x = meshpts, y = meshpts, t(elas_det.di)^.1, xlab = "ln(leaf) in year t", ylab = "ln(leaf) in year t+1")
contour(x = meshpts, y = meshpts, t(elas_det.di), add = TRUE)

## calculate elasticities of each of the vital rate functions and transition kernels
# calculate the kernel-level sensitivity function (from Ellner and Rees, pg. 96)
# calculate the elasticity of the entire kernel
# get the K matrix from the model
K <- mega.mat_det.di_proto$mega_matrix
# calculate the eigenvectors
w.z <- Re(eigen(K)$vectors[,1])
v.z1 <- Re(eigen(t(K))$vectors[,1])
# calculate the sensitivity matrix
K.sens_det.di <- outer(v.z1, w.z, "*")/sum(v.z1 * w.z * h)
# calculate the elasticity matrix
K.elas_det.di <- K.sens_det.di * (K / h) / (Re(eigen(K)$values[1]))

# calculate the elasticity of the P matrix
# get the K matrix from the model
P <- det_ipm$sub_kernels$P
# calculate the elasticity matrix
P.elas_det.di <- P * K.sens_det.di[3:502, 3:502] / (Re(eigen(K)$values[1]))
# get an actual number for the total elasticity of the P matrix
sum(P.elas_det.di) * h^2

#???? don't think the below stuff is correct
# s(z)
dK_by_ds <- outer(meshpts, meshpts, 
                  function(z1, z) {
                    (1 - 
                       (1/(1 + exp(-(data_list$p_b_int + data_list$p_b_slope * z + data_list$p_b_slope_2 * z^2))))) * # Pb(z)
                       (dnorm(x = z1, 
                              mean = (data_list$g_int + data_list$g_slope * z),
                              sd = data_list$g_sd)) # G(z',z)
                  }
                    )
s.sens.z <- apply(sens_det.di * dK_by_ds, 2, sum) * h
s.elas.z <- s.sens.z * (data_list$s_int + data_list$s_slope * meshpts)/lambda_det.di
plot(x = meshpts, y = s.sens.z, type = 'l', ylim = c(0,1), main = "s(z)", ylab = "sensitivity/elasticity", xlab = "ln(longest leaf_t)")
lines(x = meshpts, y  = s.elas.z, lty = 2)
legend('topleft', legend = c("sensitivity", "elasticity"), lty = c(1,2))

# Pb(z)
dK_by_dPb <- outer(meshpts, meshpts, 
                  function(z1, z) {
                    ((data_list$s_int + data_list$s_slope * z) * # s(z)
                      (dnorm(x = z1, # G(z',z)
                             mean = (data_list$g_int + data_list$g_slope * z),
                             sd = data_list$g_sd)))+
                      (data_list$goSdlng * 
                         (1/(1 + exp(-(data_list$b_int + data_list$b_slope * z ))))) + 
                      (data_list$goSB * 
                         (1/(1 + exp(-(data_list$b_int + data_list$b_slope * z )))))
                  }
)
Pb.sens.z <- apply(sens_det.di * dK_by_dPb, 2, sum) * h
Pb.elas.z <- Pb.sens.z * ((1/(1 + exp(-(data_list$p_b_int + data_list$p_b_slope * meshpts + data_list$p_b_slope_2 * meshpts^2))))/lambda_det.di)
plot(x = meshpts, y = Pb.sens.z, type = 'l', main = "Pb(z)", ylab = "sensitivity/elasticity", xlab = "ln(longest leaf_t)", ylim = c(0,2.8))
lines(x = meshpts, y  = Pb.elas.z, lty = 2)
legend('topleft', legend = c("sensitivity", "elasticity"), lty = c(1,2))

# pEstab
dK_by_dpEstab <- outer(meshpts, meshpts, 
                   function(z1, z) {
                     dnorm(z1, 
                           mean = data_list$c_o_mu, 
                           sd = data_list$c_o_sd
                           )
                   }
)
pEstab.sens.z <- apply(sens_det.di * dK_by_dpEstab, 2, sum) * h
pEstab.elas.z <- pEstab.sens.z * (data_list$p_estab/lambda_det.di)
plot(x = meshpts, y = pEstab.sens.z, type = 'l', main = "pEstab", ylab = "sensitivity/elasticity", xlab = "ln(longest leaf_t)")
lines(x = meshpts, y  = pEstab.elas.z, lty = 2)
legend('topleft', legend = c("sensitivity", "elasticity"), lty = c(1,2))

# b(z)
dK_by_db <- outer(meshpts, meshpts, 
                   function(z1, z) {
                     data_list$goSdlng * (1/(1 + exp(-(data_list$p_b_int + data_list$p_b_slope * z + data_list$p_b_slope_2 * z^2)))) # Pb(z) 
                     +
                       data_list$goSB * (1/(1 + exp(-(data_list$p_b_int + data_list$p_b_slope * z + data_list$p_b_slope_2 * z^2)))) # Pb(z)
                   }
)
b.sens.z <- apply(sens_det.di * dK_by_db, 2, sum) * h
b.elas.z <- b.sens.z * ((1/(1 + exp(-(data_list$b_int + data_list$b_slope * meshpts))))/lambda_det.di)
plot(x = meshpts, y = b.elas.z, type = 'l', main = "b(z)", ylab = "sensitivity/elasticity", xlab = "ln(longest leaf_t)", lty = 2)
lines(x = meshpts, y  = b.sens.z, lty = 1)
legend('topleft', legend = c("sensitivity", "elasticity"), lty = c(1,2))

# Co(z')
## Co(z') is only sensitive/elastic to changes in p_estab, and this sensitivity/elasticity is uniform across all sizes 

# outSB / staySB
## not sure how to calculate these... since they don't depend on size at all... maybe in Caswell book? 

# goSdlng / goSB
dK_by_dgoSdlng <- outer(meshpts, meshpts, 
                   function(z1, z) {
                     (1/(1 + exp(-(data_list$p_b_int + data_list$p_b_slope * z + data_list$p_b_slope_2 * z^2)))) *
                          (1/(1 + exp(-(data_list$b_int + data_list$b_slope * z )))) # b(z)
                   }
)
goSdlng.sens.z <- apply(sens_det.di * dK_by_dgoSdlng, 2, sum) * h
goSdlng.elas.z <- goSdlng.sens.z * data_list$goSdlng/lambda_det.di
plot(x = meshpts, y = goSdlng.sens.z, type = 'l', main = "goSdlng and goSB", ylab = "sensitivity/elasticity", xlab = "ln(longest leaf_t)")
lines(x = meshpts, y  = goSdlng.elas.z, lty = 2)
legend('topleft', legend = c("sensitivity", "elasticity"), lty = c(1,2))

# G(z',z)
dK_by_dG <- outer(meshpts, meshpts, 
                  function(z1, z) {
                    (1 - (1/(1 + exp(-(data_list$p_b_int + data_list$p_b_slope * z + data_list$p_b_slope_2 * z^2))))) * 
                       (data_list$s_int + data_list$s_slope * z)
                  }
)
G.sens.z <- sens_det.di * dK_by_dG
G.elas.z <- G.sens.z * outer(meshpts, meshpts, function(z1, z) {
  (dnorm(x = z1, 
         mean = (data_list$g_int + data_list$g_slope * z),
         sd = data_list$g_sd))
})/lambda_det.di
image(G.sens.z, main = "sensitivity of G(z',z)")
image(G.elas.z, main = "elasticity of G(z',z)")

## put into a mega-figure
pdf(file = "./images/det_di_IPM_elasticityPlots.pdf", onefile = TRUE, title = "Sensitivity and Elasticity of IPM parameters from a deterministic, density-independent IPM for O. coloradensis")
par(mfrow = c(5,2), mar = c(3,2,2,1))
# sensitivity function of the entire kernel
image(sens_det.di, main = "sensitivity of the continuous kernel")
contour(sens_det.di, add = TRUE)
# elasticity function of the entire kernel
image(elas_det.di, main = "elasticity of the continuous kernel")
contour(elas_det.di, add = TRUE)
#s(z)
plot(x = meshpts, y = s.sens.z, type = 'l', main = paste0("s(z); tot.elast. = ", round(sum(s.elas.z * h),2)), ylab = "sensitivity/elasticity", xlab = "ln(longest leaf_t)", ylim = c(0,2.8))
lines(x = meshpts, y  = s.elas.z, lty = 2)
legend('topleft', legend = c("sensitivity", "elasticity"), lty = c(1,2))
#b(z)
plot(x = meshpts, y = b.elas.z, type = 'l', main = paste0("b(z); tot.elast. = ", round(sum(b.elas.z * h),2)), ylab = "sensitivity/elasticity", xlab = "ln(longest leaf_t)", lty = 2, ylim = c(0,2.8),sub = paste0("tot. elasticity: ", round(sum(b.elas.z *h),2)))
lines(x = meshpts, y  = b.sens.z, lty = 1)
legend('topleft', legend = c("sensitivity", "elasticity"), lty = c(1,2))
#Pb(z)
plot(x = meshpts, y = Pb.sens.z, type = 'l',  main = paste0("Pb(z); tot.elast. = ", round(sum(Pb.elas.z * h),2)), ylab = "sensitivity/elasticity", xlab = "ln(longest leaf_t)", ylim = c(0,2.8))
lines(x = meshpts, y  = Pb.elas.z, lty = 2)
legend('topleft', legend = c("sensitivity", "elasticity"), lty = c(1,2))
#pEstab
plot(x = meshpts, y = pEstab.sens.z, type = 'l', main = paste0("pEstab; tot.elast. = ", round(sum(pEstab.elas.z * h),2)), ylab = "sensitivity/elasticity", xlab = "ln(longest leaf_t)", ylim = c(0,2.8))
lines(x = meshpts, y  = pEstab.elas.z, lty = 2)
legend('topleft', legend = c("sensitivity", "elasticity"), lty = c(1,2))
# goSdlng/goSB
plot(x = meshpts, y = goSdlng.sens.z, type = 'l',  main = paste0("goSdlng and go SB; tot.elast. = ", round(sum(goSdlng.elas.z * h),2)), ylab = "sensitivity/elasticity", xlab = "ln(longest leaf_t)", ylim = c(0,2.8))
lines(x = meshpts, y  = goSdlng.elas.z, lty = 2)
legend('topleft', legend = c("sensitivity", "elasticity"), lty = c(1,2))
#G(z',z)
# sensitivity
image(G.sens.z, main = "sensitivity of G(z',z)")
contour(G.sens.z, add = TRUE)
# elasticity
image(G.elas.z, main = paste0("elasticity of G(z',z); tot.elast. = ", round(sum(G.elas.z) * h^2,2)))
contour(G.elas.z, add = TRUE)
dev.off()
# calculate the damping ratio (to determine if the asymptotic rates are good approximations of transient dynamics)
(damp <- Re(eigen(mega.mat_det.di_proto$mega_matrix)$values[1])/Re(eigen(mega.mat_det.di_proto$mega_matrix)$values[2])) # damping ratio
# value is 1.3798, which is getting close to 1... meaning that the asymptotic rates/distributions aren't great estimates of the transient dynamics

## calculate elasticities for seedbank/seedling parameters by perturbing the ipm 
# original data list for the det_ipm model
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

## loop through w/ new parameters to update the ipm
# get the proto_ipm from the det_ipm model object
proto <- det_ipm$proto_ipm
# make a vector that contains the names of the parameters of interest
par_names <- c("goSdlng", "staySB", "goSB", "outSB", "p_estab")
# make a vector that contains the proportions by which the parameter will be changed
perc_changes <- c(.5,.6,.7,.75,.8,.85,.9,.95, 1.1, 1.15, 1.2, 1.25,1.3, 1.35, 1.4, 1.5)
rm(disc_perturbs)
for (i in 1:length(par_names)) {
  # get the name of the parameter
  par_now <- par_names[i]
  for (j in 1:length(perc_changes)) {
    # get the percentage by which it will be changed in this iteration
    perc_now <- perc_changes[j]
    # get a temporary data_list
    data_now <- data_list
    # replace the parameter in data_now with the changed value
    data_now[par_now] <- data_list[[par_now]] * perc_now
    
    # update the parameter list of the IPM
    parameters(proto) <- data_now
    
    # re-fit the ipm
    ipm_now <- proto %>% 
      make_ipm(iterate = TRUE, 
               iterations = 100)
    
    # get the lambda of this ipm
    lambda_now <- lambda(ipm_now)
    
    # save the data
    if (exists("disc_perturbs")) {
      disc_perturbs <- rbind(disc_perturbs, 
                             data.frame("param_name" = par_now, 
                                        "param_OGval" = data_list[[par_now]],
                                        "param_newVal" = data_list[[par_now]] * perc_now,
                                        "perturbation" = perc_now, 
                                        "lambda" = as.numeric(lambda_now)
                                        )
                             )
    } else {
      disc_perturbs <- data.frame(
        "param_name" = par_now, 
        "param_OGval" = data_list[[par_now]],
        "param_newVal" = data_list[[par_now]] * perc_now,
        "perturbation" = perc_now, 
        "lambda" = as.numeric(lambda_now)
        )
    }
  }
}
# calculate the sensitivity for each value
disc_perturbs$sens <- (disc_perturbs$lambda - as.numeric(lambda(det_ipm)))/(disc_perturbs$param_newVal - disc_perturbs$param_OGval)
# calculate elasticities
disc_perturbs$elas <- (disc_perturbs$param_OGval/as.numeric(lambda(det_ipm))) * disc_perturbs$sens
# average values for each parameter
mean_disc_perturbs <- disc_perturbs %>% 
  group_by(param_name) %>% 
  summarize(param_OGval = mean(param_OGval), sens_mean = mean(sens), elas_mean = mean(elas))
mean_disc_perturbs$sens_hand <- as.numeric(NA)
mean_disc_perturbs$elas_hand <- as.numeric(NA)


## try calculating sensitivities using eigen-things
# stay SB (is in position [1,1] of the matrix)
# calculate eigen-things
w.eigen_det.hand <- Re(eigen(mega.mat_det.di_proto$mega_matrix)$vectors[,1])
stable.dist_det.hand <- w.eigen_det.hand/sum(w.eigen_det.hand) 
v.eigen_det.hand <- Re(eigen(t(mega.mat_det.di_proto$mega_matrix))$vectors[,1])
repro.val_det.hand <- v.eigen_det.hand/v.eigen_det.hand[1] 

# calculate the sensitivity
staySB_sens <- (repro.val_det.hand[1]*stable.dist_det.hand[1])/(repro.val_det.hand %*% stable.dist_det.hand)
# calculate the elasticity
staySB_elas <- (data_list$staySB/as.numeric(lambda(det_ipm))) * staySB_sens
# store the values
mean_disc_perturbs[mean_disc_perturbs$param_name == "staySB", c("sens_hand", "elas_hand")] <- list(c(staySB_sens), c(staySB_elas))

# outSB (is in position [2,1] of the matrix)
# calculate the sensitivity
outSB_sens <- (repro.val_det.hand[2]*stable.dist_det.hand[1])/(repro.val_det.hand %*% stable.dist_det.hand)
# calculate the elasticity
outSB_elas <- (data_list$outSB/as.numeric(lambda(det_ipm))) * outSB_sens
# store the values
mean_disc_perturbs[mean_disc_perturbs$param_name == "outSB", c("sens_hand", "elas_hand")] <- list(c(outSB_sens), c(outSB_elas))

# goSB (is in position [2,1] of the matrix)
# calculate the sensitivity
goSB_sens <- (repro.val_det.hand[2]*stable.dist_det.hand[1])/(repro.val_det.hand %*% stable.dist_det.hand)
# calculate the elasticity
goSB_elas <- (data_list$goSB/as.numeric(lambda(det_ipm))) * goSB_sens
# store the values
mean_disc_perturbs[mean_disc_perturbs$param_name == "goSB", c("sens_hand", "elas_hand")] <- list(c(goSB_sens), c(goSB_elas))




#%%%AES%%% need to fix/work on this stuff below
## calculate the sensitivity/elasticity of the germination rate and viability rate parameters
# original data list for the det_ipm model
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
viab_OG <- 0.56
germ_OG <- 0.13

## loop through w/ new parameters to update the ipm
# get the proto_ipm from the det_ipm model object
proto <- det_ipm$proto_ipm

# make a d.f of the possible combinations of values
possible_vals <- data.frame("viab_rt" = c(rep_len(.05, length.out = 20),
                                          rep_len(.1, length.out = 20),
                                          rep_len(.15, length.out = 20),
                                          rep_len(.2, length.out = 20),
                                          rep_len(.25, length.out = 20),
                                          rep_len(.3, length.out = 20),
                                          rep_len(.35, length.out = 20),
                                          rep_len(.4, length.out = 20),
                                          rep_len(.45, length.out = 20),
                                          rep_len(.5, length.out = 20),
                                          rep_len(.55, length.out = 20),
                                          rep_len(.6, length.out = 20),
                                          rep_len(.65, length.out = 20),
                                          rep_len(.7, length.out = 20),
                                          rep_len(.75, length.out = 20),
                                          rep_len(.8, length.out = 20),
                                          rep_len(.85, length.out = 20),
                                          rep_len(.9, length.out = 20),
                                          rep_len(.95, length.out = 20),
                                          rep_len(1, length.out = 20)), 
                            "germ_rt" = rep(
                              seq(from = 0.05, to = 1, by = .05), 
                              times = 20)
                            )

rm(germ_viab_perturbs)
for (i in 1:nrow(possible_vals)) {
  germ_now <- possible_vals$germ_rt[i]
  viab_now <- possible_vals$viab_rt[i]
  if (viab_now > germ_now) {
    # calculate the parameter values using the current iteration of viab.rt and germ.rt
    outSB_now <- germ_now
    staySB_now <- (1 - germ_now) * .9
    goSB_now <- viab_now - germ_now
    goSdlng_now <- germ_now
    
    # get a temporary data_list
    data_now <- data_list
    # replace the parameter in data_now with the changed value
    data_now[c("outSB", "staySB", "goSdlng", "goSB")] <- c(outSB_now, staySB_now, goSdlng_now, goSB_now)
    
    # update the parameter list of the IPM
    parameters(proto) <- data_now
    
    # re-fit the ipm
    ipm_now <- proto %>% 
      make_ipm(iterate = TRUE, 
               iterations = 100)
    
    # get the lambda of this ipm
    lambda_now <- lambda(ipm_now)
    
    # save the data
    if (exists("germ_viab_perturbs")) {
      germ_viab_perturbs<- rbind(germ_viab_perturbs, 
                                 data.frame(
                                   "germ.rt" = germ_now,
                                   "viab.rt" = viab_now,
                                   "lambda" = as.numeric(lambda_now)
                                 )
      )
    } else {
      germ_viab_perturbs <- data.frame(
        "germ.rt" = germ_now,
        "viab.rt" = viab_now,
        "lambda" = as.numeric(lambda_now)
      )
    }
    }
  }


outSB_now <- germ_now
staySB_now <- (1 - germ_now) * .9
goSB_now <- viab_now - germ_now
goSdlng_now <- germ_now

# view lambda values
ggplot(data = germ_viab_perturbs) +
  geom_contour_filled(aes(x = germ.rt, y = viab.rt, z = lambda)) + 
  labs(title = "lambda") + 
  theme_classic()

# view seedbank parameter values
ggplot(data = germ_viab_perturbs) +
  geom_contour_filled(aes(x = germ.rt, y = viab.rt, z = outSB)) + 
  labs(title = "outSB") + 
  theme_classic()

ggplot(data = germ_viab_perturbs) +
  geom_contour_filled(aes(x = germ.rt, y = viab.rt, z = goSB)) + 
  labs(title = "goSB") + 
  theme_classic()

ggplot(data = germ_viab_perturbs) +
  geom_contour_filled(aes(x = germ.rt, y = viab.rt, z = staySB)) + 
  labs(title = "staySB") + 
  theme_classic()

ggplot(data = germ_viab_perturbs) +
  geom_contour_filled(aes(x = germ.rt, y = viab.rt, z = goSdlng)) + 
  labs(title = "goSdlng") + 
  theme_classic()

## LTRE analysis from Merow 2014 Appendix, pgs. 97-100
# compare IPM using the first transition to an IPM using the second transition




## make sure that the IPMs are not sensitive to change in size limits  
## 1. The size distribution and population growth rate (k) to which a population converges in the absence of perturbation (i.e. if the demographic transitions do not change), can be extracted directly from eigen analysis of the discretized kernel. The ‘stable size distribution’ is defined by the right eigenvector of the matrix, and the ‘asymptotic growth rate’ by the largest eigenvalue. The corresponding reproductive value, or contribution to long-term population size for each size, is defined by the left eigenvector. 
## 2. Asymptotic analyses may describe general characteristics of a population, but may poorly predict short-term dynamics. Transient dynamics are important when the population differs from the stable distribution (Williams et al. 2011). These changes can be simply quantified by projecting the population forward in time via matrix multiplication, starting from the initial size distribution and the discretized kernel (Appendix S1, A). One can determine whether transient dynamics are important using the damping ratio (q) (i.e. the ratio between the first and second eigenvalues) to describe the time-scale of transient dynamics (Caswell, 2001).
## 3. Using Markov chain theory,models of structured populations can be extended to incorporate stochastic changes in vital rates over time (Tuljapurkar 1990). These have proven powerful for exploring evolutionary dynamics (evolutionarily stable strategies; Appendix S1,G; Childs et al. 2004; Rees et al. 2006; Ellner &Rees 2007),management questions (risk of extinction; Fieberg & Ellner 2001; Morris & Doak 2002) and predicting environmental responses (Fieberg & Ellner 2001; Morris & Doak 2002). 
## 4. Sensitivity and elasticity (proportional sensitivity) analyses can determine how different parts of the kernel influence population statistics (de Kroon, van Groenendael & Ehrl?en 2000; Caswell 2001). These analyses can illustrate the relative importance of different transitions, showing in a continuous ‘landscape’ which vital rates and which size ranges contribute most to k or other population statistics (Appendix S1,A,B). Sensitivity and elasticity values can be used to estimate selection gradients in evolutionary studies (Caswell 2001) or to compare effects of different management options in conservation planning (Morris&Doak 2002). Methods exist to estimate sensitivity of an array of population characteristics in the context of transient dynamics (Caswell 2007; Haridas & Tuljapurkar 2007) and stochastic dynamics. 
## 5. Another strength of IPMs is the ability to explore vital rate parameter sensitivity (Appendix S1,B). Distinct from transition sensitivities (as above), for example, the sensitivity of k to growth regression parameters can be used to investigate the effects of changes in individual growth rate across all stages simultaneously (intercept), or in a manner that favours larger individuals over smaller individuals (slope). 
## 6. Many other population statistics are readily calculated from IPM matrices, such as passage times to life-history events (e.g. maturation; Fig. 3c), life expectancy (Fig. 3d), net reproductive rate (R0) or generation length (Appendix S1,A; Caswell 2001; Smallegange&Coulson 2013).