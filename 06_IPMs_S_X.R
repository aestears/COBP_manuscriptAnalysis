#//////////////////////////
# Integral Projection Models for Oenothera coloradensis: Code for IPMs S-X 
# (all transitions, each site, density dependent, w/ envi. covariates)
# Alice Stears
# 08 December 2021
#/////////////////////////

#### load vital rate models from previous script ####
source("./analysis_scripts/01_VitalRateModels.R")

#### fit IPMs ####
# use parameters in the 'dc_mods' d.f
# vital rate functions
# First define the elements that make up the IPM (vital rates):
# SURVIVAL:
S.fun <- function(z, N_all, tGrow, tWinter, param_list_N) {
  mu.surv=param_list_N$s_int + param_list_N$s_slope *z + param_list_N$s_dd * N_all + param_list_N$s_tGrow * tGrow
  return(1/(1 + exp(-(mu.surv))))
}
# GROWTH (we assume a constant variance)
GR.fun <- function(z,zz, N_all, tGrow, tWinter, param_list_N){
  growth.mu = param_list_N$g_int + param_list_N$g_slope *z + param_list_N$g_dd * N_all + param_list_N$g_tGrow * tGrow
  return(dnorm(zz, mean = growth.mu, sd = param_list_N$g_sd))
}
## SEEDLING SIZES (same approach as in growth function)
SDS.fun <- function(zz, N_all, tGrow, tWinter, param_list_N){
  rec_mu <- param_list_N$c_o_int + param_list_N$c_o_dd * N_all + param_list_N$c_o_tGrow * tGrow + param_list_N$c_o_tWinter * tWinter
  rec_sd <- param_list_N$c_o_sd
  return(dnorm(zz, mean = rec_mu, sd = rec_sd))
}
# PROBABILITY OF FLOWERING
FL.fun <- function(z, N_all, tGrow, tWinter, param_list_N) {
  mu.fl = param_list_N$p_b_int + param_list_N$p_b_slope*z +  param_list_N$p_b_slope_2 * (z^2) + param_list_N$p_b_dd * N_all + param_list_N$p_b_tGrow * tGrow + param_list_N$p_b_tWinter * tWinter
  return(1/(1+ exp(-(mu.fl))))
}
# SEED PRODUCTION
SDP.fun <- function(z, N_all, tGrow, tWinter, param_list_N) {
  mu.fps=exp(param_list_N$b_int + param_list_N$b_slope *z + param_list_N$b_dd * N_all + param_list_N$b_tGrow * tGrow + param_list_N$b_tWinter * tWinter)
  return(mu.fps)
}

outSB <- outSB_all #SB to continuous stage
staySB <- staySB_all # staying in SB
goCont <- goCont_all # seeds become continuous right away (without going to the seed bank)
goSB <- goSB_all # seeds go to the seedbank
surv.seeds <-  0.9 # survival of seeds

lambdas <- list()
IPMs_S_X <- vector(mode = "list", length = 6)
names(IPMs_S_X) <- unique(dc_mods$Site)
for (i in unique(dc_mods$Site)) {
  # get the data for this site
  dat_now <- dc_mods[dc_mods$Site == i,]
  # get the vital rates into the right format
  paramNow <- list(
    g_int = dat_now[dat_now$vitalRate == "Intercept","grow"],
    g_slope = dat_now[dat_now$vitalRate == "log_LL_t","grow"],
    g_dd  = dat_now[dat_now$vitalRate == "N_all_plot","grow"],
    g_tGrow = dat_now[dat_now$vitalRate == "tMean_grow_C","grow"],
    g_sd = 0.4837,
    s_int = dat_now[dat_now$vitalRate == "Intercept","surv"],
    s_slope = dat_now[dat_now$vitalRate == "log_LL_t","surv"],
    s_dd = dat_now[dat_now$vitalRate == "N_all_plot","surv"],
    s_tGrow = dat_now[dat_now$vitalRate == "tMean_grow_C","surv"],
    p_b_int = dat_now[dat_now$vitalRate == "Intercept","flwr"],
    p_b_slope = dat_now[dat_now$vitalRate == "log_LL_t","flwr"],
    p_b_slope_2 = dat_now[dat_now$vitalRate == "I(log_LL_t^2)","flwr"],
    p_b_dd = dat_now[dat_now$vitalRate == "N_all_plot","flwr"],
    p_b_tGrow = dat_now[dat_now$vitalRate == "tMean_grow_C","flwr"],
    p_b_tWinter = dat_now[dat_now$vitalRate == "tMean_winter_C","flwr"],
    b_int = dat_now[dat_now$vitalRate == "Intercept","seed"],
    b_slope = dat_now[dat_now$vitalRate == "log_LL_t","seed"],
    b_dd = dat_now[dat_now$vitalRate == "N_all_plot","seed"],
    b_tGrow = dat_now[dat_now$vitalRate == "tMean_grow_C","seed"],
    b_tWinter = dat_now[dat_now$vitalRate == "tMean_winter_C","seed"],
    c_o_int = dat_now[dat_now$vitalRate == "Intercept","rec"],
    c_o_dd = dat_now[dat_now$vitalRate == "N_all_plot","rec"],
    c_o_tGrow = dat_now[dat_now$vitalRate == "tMean_grow_C","rec"],
    c_o_tWinter = dat_now[dat_now$vitalRate == "tMean_winter_C","rec"],
    c_o_sd = 0.7684952
  )
  # set up IPM
  L <-  1.2 * min(dat_all$log_LL_t, na.rm = TRUE) # minimum size
  U <-  1.2 * max(dat_all$log_LL_t, na.rm = TRUE) # maximum size
  n <-500 # bins
  K <- array(0,c(n+1,n+1))
  b <- L+c(0:n)*(U-L)/n # interval that each cell of the matrix covers
  meshp <- 0.5*(b[1:n]+b[2:(n+1)]) # midpoint
  h=(U-L)/n # bin width
  S <- diag(S.fun(meshp, N_all = 500, tGrow = mean(unique(dat_all$tMean_grow_C), na.rm = TRUE), tWinter = mean(unique(dat_all$tMean_winter_C), na.rm = TRUE), param_list_N = paramNow)) # Survival # put survival probabilities in the diagonal of the matrix
  G <- h * t(outer(meshp,meshp,GR.fun, N_all = 500, tGrow = mean(unique(dat_all$tMean_grow_C), na.rm = TRUE), tWinter = mean(unique(dat_all$tMean_winter_C), na.rm = TRUE), param_list_N = paramNow)) # Growth
  c_o <- h * matrix(rep(SDS.fun(meshp, N_all = 500, tGrow = mean(unique(dat_all$tMean_grow_C), na.rm = TRUE), tWinter = mean(unique(dat_all$tMean_winter_C), na.rm = TRUE), param_list_N = paramNow),n),n,n,byrow=F)
  Pb <- (FL.fun(meshp, N_all = 500, tGrow = mean(unique(dat_all$tMean_grow_C), na.rm = TRUE), tWinter = mean(unique(dat_all$tMean_winter_C), na.rm = TRUE), param_list_N = paramNow))
  b_seed <- (SDP.fun(meshp, N_all = 500, tGrow = mean(unique(dat_all$tMean_grow_C), na.rm = TRUE), tWinter = mean(unique(dat_all$tMean_winter_C), na.rm = TRUE), param_list_N = paramNow))
  FecALL <- Pb * b_seed
  S_new <- S * (1-Pb)
  G <- G/matrix(as.vector(apply(G,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  c_o <- c_o/matrix(as.vector(apply(c_o,2,sum)),nrow=n,ncol=n,byrow=TRUE)
  Pkernel.cont <- as.matrix(G %*% S_new)
  Pkernel.seedbank <- c(staySB, outSB*c_o[,1])
  Pkernel <- cbind(Pkernel.seedbank,rbind(rep(0,length(meshp)),Pkernel.cont))
  Fkernel.cont <-  as.matrix(goCont * ((c_o) %*% diag(FecALL)))
  Fkernel.discr  <- matrix(c(0, goSB * (FecALL)), nrow = 1)
  Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))
  
  mat_all_DD <-Pkernel+Fkernel
  
  eigenMat <- base::eigen(mat_all_DD)
  # get the lambda
  lam_now <-  as.numeric(eigenMat$values[1])
  lambdas[[i]] <- lam_now
  
  IPMs_S_X[[i]] <- mat_all_DD
}

#### save the data to file ####
fileLoc <- "./intermediate_analysis_Data/"
## site-level DI IPM matrices
saveRDS(IPMs_S_X, file = paste0(fileLoc,"/IPMs_S_X.RDS"))
