#//////////////////////////
# Integral Projection Models for Oenothera coloradensis: Elasticity Analysis
# Alice Stears
# 11 December 2021
#/////////////////////////

# load required packages
library(tidyverse)
library(ipmr)

#### load vital rate models from previous script ####
source("./COBP_manuscriptAnalysis/01_VitalRateModels.R")
load(file = "./intermediate_analysis_Data/ipmA_B_results.RData")

# lower limit of size in the ipm
L = min(dat_all$log_LL_t, na.rm = TRUE) * .8
# upper limit of size in the ipm
U = max(dat_all$log_LL_t, na.rm = TRUE) * 1.2
# number of break-points
n <- 500
meshpts <- seq(from = L, to = U, length.out = 500)
h = meshpts[2] - meshpts[1]

#### DI all dat_all IPM (IPM B)####
## IPM calculated by hand is called 'mat_all_DI' 
# calculate lambda
(lam_allDI <- Re(eigen(mat_all_DI)$values[1])) # 1.91

## calculate the stable size distribution
w.eigen_temp <-  Re(eigen(mat_all_DI)$vectors[,1]) 
stable.dist_allDI <- w.eigen_temp/sum(w.eigen_temp) 

## calculate the reproductive value distribution
v.eigen_temp <- Re(eigen(t(mat_all_DI))$vectors[,1]) 
repro.val_allDI <- v.eigen_temp/v.eigen_temp[1]

## calculate the sensitivity and elasticity matrices
# The eigen-things can be combined to obtain the sensitivity and elasticity matrices.
v.dot.w_allDI <- sum(stable.dist_allDI * repro.val_allDI)*h
# calculate the sensitivity function (whole kernel)
sens_allDI <- outer(repro.val_allDI,stable.dist_allDI, '*')/(v.dot.w_allDI)
# calculate the elasticity function (whole kernel)
elas_allDI <- matrix(as.vector(sens_allDI)*as.vector(mat_all_DI)/lam_allDI,nrow=501)
# plot the sensitivity function of the entire continuous kernel
image(x = meshpts, y = meshpts, t(sens_allDI[2:501,2:501])^.1, 
      xlab = "ln(leaf) in year t", ylab = "ln(leaf) in year t+1")
contour(x = meshpts, y = meshpts, t(sens_allDI[2:501,2:501]), add = TRUE)
# plot the elasticity function of the entire continuous kernel
image(x = meshpts, y = meshpts, t(elas_allDI[2:501,2:501])^.1, xlab = "ln(leaf) in year t", ylab = "ln(leaf) in year t+1")
contour(x = meshpts, y = meshpts, t(elas_allDI[2:501,2:501]), add = TRUE)

# make matrices using ipmr data to double-check
## calculate the stable size distribution
w.eigen_test <-  right_ev(ipm_B)
w.eigen_test <- c(w.eigen_test$b_w,w.eigen_test$size_w)
stable.dist_test <- w.eigen_test/sum(w.eigen_test) 

## calculate the reproductive value distribution
v.eigen_test <- left_ev(ipm_B)
v.eigen_test <- c(v.eigen_test$b_v, v.eigen_test$size_v)
repro.val_test <- v.eigen_test/v.eigen_test[1]

## make elas and sens matrices
v.dot.w_test <- sum(stable.dist_test * repro.val_test)*h
# calculate the sensitivity function (whole kernel)
sens_test <- outer(repro.val_test,stable.dist_test, '*')/(v.dot.w_test)
# calculate the elasticity function (whole kernel)
elas_test <- matrix(as.vector(sens_test)*as.vector(mat_all_DI)/lam_all_DI,nrow=501)
## use these instead...


## calculate sensitivity and elasticity of individual vital rates
## loop through w/ new parameters to update the ipm
# sens/elas of germination and viability rates
# make a vector that contains the names of the parameters of interest
par_names <- c("viab.rt", "germ.rt", "surv.rt")
# get info for all_DI IPM
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
#vital rate functions
S.fun <- function(z, paramCont) {
  mu.surv=paramCont[["survival"]]["(Intercept)",] + paramCont[["survival"]]["log_LL_t",]*z
  return(1/(1 + exp(-(mu.surv))))
}
GR.fun <- function(z,zz, paramCont){
  growth.mu = paramCont[["growth"]]["(Intercept)",1] + paramCont[["growth"]]["log_LL_t",1]*z
  return(dnorm(zz, mean = growth.mu, sd = paramCont[["growth"]][1,2]))
}
SDS.fun <- function(zz, paramCont){
  rec_mu <- paramCont[["recruitDist"]][1]
  rec_sd <- paramCont[["recruitDist"]][2]
  return(dnorm(zz, mean = rec_mu, sd = rec_sd))
}
FL.fun <- function(z, paramCont) {
  mu.fl = paramCont[["flowering"]][1,] + paramCont[["flowering"]][2,]*z +  paramCont[["flowering"]][3,]* (z^2)
  return(1/(1+ exp(-(mu.fl))))
}
SDP.fun <- function(z, paramCont) {
  mu.fps=exp(paramCont[["seedProduction"]][1,1] + paramCont[["seedProduction"]][2,1]*z)
  return(mu.fps)
}
# Define the lower and upper integration limit
L <-  1.2 * min(dat_all$log_LL_t, na.rm = TRUE) # minimum size
U <-  1.2 * max(dat_all$log_LL_t, na.rm = TRUE) # maximum size
n <-500 # bins

# previous SB vital rates (From script 01)
#germ.rt 
#viab.rt 
#surv.seeds

# make a vector that contains the proportions by which the parameter will be changed
perc_changes <- seq(0,1, by = .05)

for (i in 1:length(par_names)) {
  # get the name of the parameter
  par_now <- par_names[i]
  if (par_now == "viab.rt") {
    for (j in 1:length(perc_changes)) {
      # update the germ.rt
      viab_now <- perc_changes[j]
      # update the SB parameters
      outSB_now <- germ.rt * surv.seeds 
      staySB_now <- (1-germ.rt) * surv.seeds
      goSB_now <- viab_now * (1 - germ.rt)
      goCont_now <- viab_now * germ.rt
      # make the IPM
      K <- array(0,c(n+1,n+1))
      b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
      meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
      h=(U-L)/n # bin width 
      # Survival and growth 
      S <- diag(S.fun(meshp, paramCont)) # Survival # put survival probabilities in the diagonal of the matrix
      G <- h * t(outer(meshp,meshp,GR.fun, paramCont = paramCont)) # Growth
      #Recruits distribution 
      c_o <- h * matrix(rep(SDS.fun(meshp, paramCont),n),n,n,byrow=F)
      # c_o <- matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
      #Probability of flowering
      Pb = (FL.fun(meshp, paramCont))
      #Number of seeds produced according to adult size
      b_seed = (SDP.fun(meshp, paramCont))
      FecALL= Pb * b_seed
      # update the 'S' matrix by multiplying it by (1-Pb), since this is a monocarpic perennial
      S_new <- S * (1-Pb)
      # Control for eviction:
      G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
      c_o <- c_o/matrix(as.vector(apply(c_o,2,sum)),nrow=n,ncol=n,byrow=TRUE)
      # make the continuous part of the P matrix
      Pkernel.cont <- as.matrix(G %*% S_new)
      # seedbank (first column of your K)
      Pkernel.seedbank = c(staySB_now, outSB_now*c_o[,1]) # seeds survive and go to continuous
      # Make the full P kernel
      Pkernel <- cbind(Pkernel.seedbank,rbind(rep(0,length(meshp)),Pkernel.cont)) # discrete component
      ## make the F kernel
      Fkernel.cont <-  as.matrix(goCont_now * ((c_o) %*% diag(FecALL))) # the size of seedlings that go into the seed bank from each continuous size class
      Fkernel.discr  <- matrix(c(0, goSB_now * (FecALL)), nrow = 1)
      Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))
      # calculate the entire matrix
      mat_now <-Pkernel+Fkernel
      # calculate lambda
      lambda_now <- Re(eigen(mat_now)$values[1])
      # store the data
      if (i == 1 & j == 1) {
        disc_perturbs <- data.frame(
          "param_name" = par_now, 
          "param_OGval" = viab.rt,
          "param_newVal" = viab_now,
          "perturbation" = viab.rt - viab_now, 
          "lambda" = as.numeric(lambda_now))
      } else {
        disc_perturbs <- rbind(disc_perturbs, 
                               data.frame(
                                 "param_name" = par_now, 
                                 "param_OGval" = viab.rt,
                                 "param_newVal" = viab_now,
                                 "perturbation" = viab.rt - viab_now, 
                                 "lambda" = as.numeric(lambda_now)))
      }
    }
  } else if (par_now == "germ.rt") {
    for (j in 1:length(perc_changes)) {
      # update the germ.rt
      germ_now <- perc_changes[j]
      # update the SB parameters
      outSB_now <- germ_now * surv.seeds 
      staySB_now <- (1-germ_now) * surv.seeds
      goSB_now <- viab.rt * (1 - germ_now)
      goCont_now <- viab.rt * germ_now
      # make IPM
      K <- array(0,c(n+1,n+1))
      b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
      meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
      h=(U-L)/n # bin width 
      # Survival and growth 
      S <- diag(S.fun(meshp, paramCont = paramCont)) # Survival # put survival probabilities in the diagonal of the matrix
      G <- h * t(outer(meshp,meshp,GR.fun, paramCont = paramCont)) # Growth
      #Recruits distribution 
      c_o <- h * matrix(rep(SDS.fun(meshp, paramCont = paramCont),n),n,n,byrow=F)
      # c_o <- matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
      #Probability of flowering
      Pb = (FL.fun(meshp, paramCont = paramCont))
      #Number of seeds produced according to adult size
      b_seed = (SDP.fun(meshp, paramCont = paramCont))
      FecALL= Pb * b_seed
      # update the 'S' matrix by multiplying it by (1-Pb), since this is a monocarpic perennial
      S_new <- S * (1-Pb)
      # Control for eviction:
      G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
      c_o <- c_o/matrix(as.vector(apply(c_o,2,sum)),nrow=n,ncol=n,byrow=TRUE)
      # make the continuous part of the P matrix
      Pkernel.cont <- as.matrix(G %*% S_new)
      # seedbank (first column of your K)
      Pkernel.seedbank = c(staySB_now, outSB_now*c_o[,1]) # seeds survive and go to continuous
      # Make the full P kernel
      Pkernel <- cbind(Pkernel.seedbank,rbind(rep(0,length(meshp)),Pkernel.cont)) # discrete component
      ## make the F kernel
      Fkernel.cont <-  as.matrix(goCont_now * ((c_o) %*% diag(FecALL))) # the size of seedlings that go into the seed bank from each continuous size class
      Fkernel.discr  <- matrix(c(0, goSB_now * (FecALL)), nrow = 1)
      Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))
      # calculate the entire matrix
      mat_now <-Pkernel+Fkernel
      # calculate lambda
      lambda_now <- Re(eigen(mat_now)$values[1])
      # store the data
      if (i == 1 & j == 1) {
        disc_perturbs <- data.frame(
          "param_name" = par_now, 
          "param_OGval" = germ.rt,
          "param_newVal" = germ_now,
          "perturbation" = germ.rt - germ.nw, 
          "lambda" = as.numeric(lambda_now))
      } else {
        disc_perturbs <- rbind(disc_perturbs, 
                               data.frame("param_name" = par_now, 
                                          "param_OGval" = germ.rt,
                                          "param_newVal" = germ_now,
                                          "perturbation" = germ.rt - germ_now, 
                                          "lambda" = as.numeric(lambda_now)))
      }
    }
  } else if (par_now == "surv.rt"){
    for (j in 1:length(perc_changes)) {
      # update the germ.rt
      surv_now <- perc_changes[j]
      # update the SB parameters
      outSB_now <- germ.rt * surv_now 
      staySB_now <- (1-germ.rt) * surv_now
      goSB_now <- viab.rt * (1 - germ.rt)
      goCont_now <- viab.rt * germ.rt
      # make the IPM
      K <- array(0,c(n+1,n+1))
      b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
      meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
      h=(U-L)/n # bin width 
      # Survival and growth 
      S <- diag(S.fun(meshp, paramCont = paramCont)) # Survival # put survival probabilities in the diagonal of the matrix
      G <- h * t(outer(meshp,meshp,GR.fun, paramCont = paramCont)) # Growth
      #Recruits distribution 
      c_o <- h * matrix(rep(SDS.fun(meshp, paramCont = paramCont),n),n,n,byrow=F)
      # c_o <- matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
      #Probability of flowering
      Pb = (FL.fun(meshp, paramCont = paramCont))
      #Number of seeds produced according to adult size
      b_seed = (SDP.fun(meshp, paramCont = paramCont))
      FecALL= Pb * b_seed
      # update the 'S' matrix by multiplying it by (1-Pb), since this is a monocarpic perennial
      S_new <- S * (1-Pb)
      # Control for eviction:
      G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
      c_o <- c_o/matrix(as.vector(apply(c_o,2,sum)),nrow=n,ncol=n,byrow=TRUE)
      # make the continuous part of the P matrix
      Pkernel.cont <- as.matrix(G %*% S_new)
      # seedbank (first column of your K)
      Pkernel.seedbank = c(staySB_now, outSB_now*c_o[,1]) # seeds survive and go to continuous
      # Make the full P kernel
      Pkernel <- cbind(Pkernel.seedbank,rbind(rep(0,length(meshp)),Pkernel.cont)) # discrete component
      ## make the F kernel
      Fkernel.cont <-  as.matrix(goCont_now * ((c_o) %*% diag(FecALL))) # the size of seedlings that go into the seed bank from each continuous size class
      Fkernel.discr  <- matrix(c(0, goSB_now * (FecALL)), nrow = 1)
      Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))
      # calculate the entire matrix
      mat_now <-Pkernel+Fkernel
      # calculate lambda
      lambda_now <- Re(eigen(mat_now)$values[1])
      # store the data
      if (i == 1 & j == 1) {
        disc_perturbs <- data.frame(
          "param_name" = par_now, 
          "param_OGval" = surv.seeds,
          "param_newVal" = surv_now,
          "perturbation" = surv.seeds - surv_now, 
          "lambda" = as.numeric(lambda_now))
      } else {
        disc_perturbs <- rbind(disc_perturbs, 
                               data.frame("param_name" = par_now, 
                                          "param_OGval" = surv.seeds,
                                          "param_newVal" = surv_now,
                                          "perturbation" = surv.seeds - surv_now, 
                                          "lambda" = as.numeric(lambda_now)))
      }
    }
  }
}
# calculate the sensitivity for each value
disc_perturbs$sens <- (disc_perturbs$lambda - as.numeric((lam_allDI)))/(disc_perturbs$param_newVal - disc_perturbs$param_OGval)
# calculate elasticities
disc_perturbs$elas <- (disc_perturbs$param_OGval/as.numeric(lam_allDI)) * disc_perturbs$sens
# remove the 'inf' values, which result from having a difference of '0' between new and old values
disc_perturbs[is.infinite(disc_perturbs$sens),c("sens", "elas")] <- NA
# average values for each parameter
mean_disc_perturbs <- disc_perturbs %>% 
  group_by(param_name) %>% 
  summarize(param_OGval = mean(param_OGval, na.rm = TRUE), sens_mean = mean(sens, na.rm = TRUE), elas_mean = mean(elas, na.rm = TRUE))

#dir.create("./intermediate_analysis_Data/allSiteAllYears_noDDnoEnv")
saveRDS(mean_disc_perturbs, file = "./intermediate_analysis_Data/discreteParamElasticity.RDS")
## get elasticity/sensitivity for other vital rate parameters (model parameters in vital rate models)
# previous SB vital rates
# germ.rt 
# viab.rt 
# surv.seeds 
# 'actual' parameters are in the 'paramCont' list

# make a vector that contains the proportions by which the parameter will be changed
perc_changes <- seq(0.1,2, by = .05)

model_perturbs <-  data.frame( "param_name" = as.character(NA), 
                               "param_OGval" = NA,
                               "param_newVal" = NA, 
                               "lambda" = NA)
# loop through each of the vital rate processes
for (i in 1:length(names(paramCont))) {
  # get the name of the vital rate
  modelName <- names(paramCont)[i]
  if (ncol(paramCont[[modelName]]) == 1){
    for (j in 1: nrow(paramCont[[modelName]])) {
      for (k in 1:length(perc_changes)) {
        paramNow <- paramCont
        # replace the value in the parameter list
        paramNow[[modelName]][j,] <- paramCont[[modelName]][j]*perc_changes[k]
        # make the IPM
        K <- array(0,c(n+1,n+1))
        b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
        meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
        h=(U-L)/n # bin width 
        # Survival and growth 
        S <- diag(S.fun(meshp, paramNow)) # Survival # put survival probabilities in the diagonal of the matrix
        G <- h * t(outer(meshp,meshp,GR.fun, paramNow)) # Growth
        #Recruits distribution 
        c_o <- h * matrix(rep(SDS.fun(meshp, paramNow),n),n,n,byrow=F)
        # c_o <- matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
        #Probability of flowering
        Pb = (FL.fun(meshp, paramNow))
        #Number of seeds produced according to adult size
        b_seed = (SDP.fun(meshp, paramNow))
        FecALL= Pb * b_seed
        # update the 'S' matrix by multiplying it by (1-Pb), since this is a monocarpic perennial
        S_new <- S * (1-Pb)
        # Control for eviction:
        G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
        c_o <- c_o/matrix(as.vector(apply(c_o,2,sum)),nrow=n,ncol=n,byrow=TRUE)
        # make the continuous part of the P matrix
        Pkernel.cont <- as.matrix(G %*% S_new)
        # seedbank (first column of your K)
        Pkernel.seedbank = c(staySB_all, outSB_all*c_o[,1]) # seeds survive and go to continuous
        # Make the full P kernel
        Pkernel <- cbind(Pkernel.seedbank,rbind(rep(0,length(meshp)),Pkernel.cont)) # discrete component
        ## make the F kernel
        Fkernel.cont <-  as.matrix(goCont_all * ((c_o) %*% diag(FecALL))) # the size of seedlings that go into the seed bank from each continuous size class
        Fkernel.discr  <- matrix(c(0, goSB_all * (FecALL)), nrow = 1)
        Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))
        # calculate the entire matrix
        mat_now <-Pkernel+Fkernel
        # calculate lambda
        lambda_now <- Re(eigen(mat_now)$values[1])
        # store the data
        model_perturbs <- rbind(model_perturbs, 
                                data.frame(
                                  "param_name" = paste0(modelName,"_", names(paramNow[[modelName]][j,])), 
                                  "param_OGval" = paramCont[[modelName]][j],
                                  "param_newVal" = paramNow[[modelName]][j,], 
                                  "lambda" = as.numeric(lambda_now)))
      }
    }
  } else if (ncol(paramCont[[modelName]]) > 1) {
    for (j in 1: nrow(paramCont[[modelName]])) {
      for (k in 1:length(perc_changes)) {
        paramNow <- paramCont
        # replace the value in the parameter list
        paramNow[[modelName]][j,1] <- paramCont[[modelName]][j]*perc_changes[k]
        # make the IPM
        K <- array(0,c(n+1,n+1))
        b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
        meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
        h=(U-L)/n # bin width 
        # Survival and growth 
        S <- diag(S.fun(meshp, paramNow)) # Survival # put survival probabilities in the diagonal of the matrix
        G <- h * t(outer(meshp,meshp,GR.fun, paramNow)) # Growth
        #Recruits distribution 
        c_o <- h * matrix(rep(SDS.fun(meshp, paramNow),n),n,n,byrow=F)
        # c_o <- matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
        #Probability of flowering
        Pb = (FL.fun(meshp, paramNow))
        #Number of seeds produced according to adult size
        b_seed = (SDP.fun(meshp, paramNow))
        FecALL= Pb * b_seed
        # update the 'S' matrix by multiplying it by (1-Pb), since this is a monocarpic perennial
        S_new <- S * (1-Pb)
        # Control for eviction:
        G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
        c_o <- c_o/matrix(as.vector(apply(c_o,2,sum)),nrow=n,ncol=n,byrow=TRUE)
        # make the continuous part of the P matrix
        Pkernel.cont <- as.matrix(G %*% S_new)
        # seedbank (first column of your K)
        Pkernel.seedbank = c(staySB_all, outSB_all*c_o[,1]) # seeds survive and go to continuous
        # Make the full P kernel
        Pkernel <- cbind(Pkernel.seedbank,rbind(rep(0,length(meshp)),Pkernel.cont)) # discrete component
        ## make the F kernel
        Fkernel.cont <-  as.matrix(goCont_all * ((c_o) %*% diag(FecALL))) # the size of seedlings that go into the seed bank from each continuous size class
        Fkernel.discr  <- matrix(c(0, goSB_all * (FecALL)), nrow = 1)
        Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))
        # calculate the entire matrix
        mat_now <-Pkernel+Fkernel
        mat_now[which(is.nan(mat_now))] <- NA
        # calculate lambda
        lambda_now <- Re(eigen(mat_now)$values[1])
        # store the data
        model_perturbs <- rbind(model_perturbs, 
                                data.frame(
                                  "param_name" = paste0(modelName,"_", names(paramNow[[modelName]][j,1])), 
                                  "param_OGval" = paramCont[[modelName]][j],
                                  "param_newVal" = paramNow[[modelName]][j,], 
                                  "lambda" = as.numeric(lambda_now)))
      }
    }
    for (l in 2:ncol(paramCont[[modelName]])) {
      for (k in 1:length(perc_changes)) {
        paramNow <- paramCont
        # replace the value in the parameter list
        paramNow[[modelName]][1,l] <- paramCont[[modelName]][1,l]*perc_changes[k]
        # make the IPM
        K <- array(0,c(n+1,n+1))
        b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers 
        meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
        h=(U-L)/n # bin width 
        # Survival and growth 
        S <- diag(S.fun(meshp, paramNow)) # Survival # put survival probabilities in the diagonal of the matrix
        G <- h * t(outer(meshp,meshp,GR.fun, paramNow)) # Growth
        #Recruits distribution 
        c_o <- h * matrix(rep(SDS.fun(meshp, paramNow),n),n,n,byrow=F)
        # c_o <- matrix(rep(SDS.fun(meshp),n),n,n,byrow=F)
        #Probability of flowering
        Pb = (FL.fun(meshp, paramNow))
        #Number of seeds produced according to adult size
        b_seed = (SDP.fun(meshp, paramNow))
        FecALL= Pb * b_seed
        # update the 'S' matrix by multiplying it by (1-Pb), since this is a monocarpic perennial
        S_new <- S * (1-Pb)
        # Control for eviction:
        G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
        c_o <- c_o/matrix(as.vector(apply(c_o,2,sum)),nrow=n,ncol=n,byrow=TRUE)
        # make the continuous part of the P matrix
        Pkernel.cont <- as.matrix(G %*% S_new)
        # seedbank (first column of your K)
        Pkernel.seedbank = c(staySB_all, outSB_all*c_o[,1]) # seeds survive and go to continuous
        # Make the full P kernel
        Pkernel <- cbind(Pkernel.seedbank,rbind(rep(0,length(meshp)),Pkernel.cont)) # discrete component
        ## make the F kernel
        Fkernel.cont <-  as.matrix(goCont_all * ((c_o) %*% diag(FecALL))) # the size of seedlings that go into the seed bank from each continuous size class
        Fkernel.discr  <- matrix(c(0, goSB_all * (FecALL)), nrow = 1)
        Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))
        # calculate the entire matrix
        mat_now <-Pkernel+Fkernel
        mat_now[which(is.nan(mat_now))] <- NA
        # calculate lambda
        lambda_now <- Re(eigen(mat_now)$values[1])
        # store the data
        model_perturbs <- rbind(model_perturbs, 
                                data.frame(
                                  "param_name" = paste0(modelName,"_stndDev"), 
                                  "param_OGval" = paramCont[[modelName]][1,l],
                                  "param_newVal" = paramNow[[modelName]][1,l], 
                                  "lambda" = as.numeric(lambda_now)))
      }
    }
  }
}
# calculate the sensitivity for each value
model_perturbs$sens <- (model_perturbs$lambda - as.numeric((lam_allDI)))/(model_perturbs$param_newVal - model_perturbs$param_OGval)
# calculate elasticities
model_perturbs$elas <- (model_perturbs$param_OGval/as.numeric(lam_allDI)) * model_perturbs$sens
# remove the 'inf' values, which result from having a difference of '0' between new and old values
model_perturbs[is.infinite(model_perturbs$sens),c("sens", "elas")] <- NA
# average values for each parameter
mean_model_perturbs <- model_perturbs %>% 
  group_by(param_name) %>% 
  summarize(param_OGval = mean(param_OGval, na.rm = TRUE), sens_mean = mean(sens, na.rm = TRUE), elas_mean = mean(elas, na.rm = TRUE))

saveRDS(mean_model_perturbs, file = "./intermediate_analysis_Data/continuousParamElasticity.RDS")

