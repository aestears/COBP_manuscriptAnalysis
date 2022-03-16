#//////////////////////////
# Integral Projection Models for Oenothera coloradensis: Analysis for spatial asynchrony
# Alice Stears
# 11 March 2022
#/////////////////////////

library(tidyverse)
library(rstatix)

# load data from script 1
dat_all <- read.csv(file = "/Users/astears/COBP_project/allDat_plus_contSeedlings.csv")

#### IPMs for each subpop over both transitions ####
#### DI IPM for each site--first half of data ####
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

## make an empty list to hold the IPM kernels 
site_IPMs_all_DI <- list()
for (i in c(1:length(unique(dat_all$Site)))) {
  ## get data for this 'current' site
  dat_now <- dat_all[dat_all$Site == unique(dat_all$Site)[i] # data for this site
                     ,]
  
  ## fit vital rate models
  ## Survival ($s(z)$)
  survDat_now <- dat_now[dat_now$flowering == 0 | is.na(dat_now$flowering),]
  survMod_now <- glm(survives_tplus1 ~ log_LL_t , data = survDat_now, family = binomial)
  ## Growth ($G(z',z)$)
  sizeMod_now <- lm(log_LL_tplus1 ~ log_LL_t , data = dat_now)
  ## Number of seeds produced, according to plant size ($b(z)$)
  seedDat_now <- dat_now[dat_now$flowering==1,]
  # fit poisson glm (for count data)
  seedMod_now <- MASS::glm.nb(Num_seeds ~ log_LL_t , data = seedDat_now)
  ## Flowering probability ($p_b(z)$)
  flwrMod_now <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = dat_now, family = binomial)))
  ## Distribution of recruit size ($c_o(z')$)
  # subset the data
  recD_now <- dat_all[dat_all$seedling == 1,]
  # fit the model
  recMod_now <- lm(log_LL_t ~ 1, data = recD_now)
  
  ## put in the parameter list (paramCont)
  paramCont <- list(
    g_int     = coef(sizeMod_now)[1], # growth 
    g_slope   = coef(sizeMod_now)[2],
    g_sd      = summary(sizeMod_now)$sigma,
    s_int     = coef(survMod_now)[1], # survival
    s_slope   = coef(survMod_now)[2],
    p_b_int   = coef(flwrMod_now)[1], #probability of flowering
    p_b_slope = coef(flwrMod_now)[2],
    p_b_slope_2 = coef(flwrMod_now)[3],
    b_int   = coef(seedMod_now)[1], #seed production
    b_slope = coef(seedMod_now)[2],
    c_o_mu    = coef(recMod_now), #recruit size distribution
    c_o_sd    = summary(recMod_now)$sigma,
    outSB  = outSB_all,
    staySB = staySB_all,
    goSB   = goSB_all, 
    goCont = goCont_all                  
  )
  # SURVIVAL:
  S.fun <- function(z) {
    mu.surv=paramCont$s_int + paramCont$s_slope *z
    return(1/(1 + exp(-(mu.surv))))
  }
  # GROWTH (we assume a constant variance)
  GR.fun <- function(z,zz){
    growth.mu = paramCont$g_int + paramCont$g_slope*z
    return(dnorm(zz, mean = growth.mu, sd = paramCont$g_sd))
  }
  ## SEEDLING SIZES (same approach as in growth function)
  SDS.fun <- function(zz){
    rec_mu <- paramCont$c_o_mu
    rec_sd <- paramCont$c_o_sd
    return(dnorm(zz, mean = rec_mu, sd = rec_sd))
  }
  # PROBABILITY OF FLOWERING 
  FL.fun <- function(z) {
    mu.fl = paramCont$p_b_int + paramCont$p_b_slope*z +  paramCont$p_b_slope_2 * (z^2)
    return(1/(1+ exp(-(mu.fl))))
  }
  # SEED PRODUCTION
  SDP.fun <- function(z) {
    mu.fps=exp(paramCont$b_int + paramCont$b_slope *z)
    return(mu.fps)
  }
  
  ## fit the IPM
  K <- array(0,c(n+1,n+1))
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
  # seedbank (first column of your K)
  Pkernel.seedbank = c(staySB, outSB*c_o[,1]) # seeds survive and go to continuous
  # Make the full P kernel
  Pkernel <- cbind(Pkernel.seedbank,rbind(rep(0,length(meshp)),Pkernel.cont)) # discrete component
  ## make the F kernel
  Fkernel.cont <-  as.matrix(goCont * ((c_o) %*% diag(FecALL))) # the size of seedlings that go into the seed bank from each continuous size class
  Fkernel.discr  <- matrix(c(0, goSB * (FecALL)), nrow = 1)
  # multiply the cont_to_disc distribution by the binwidth (h)
  Fkernel <- rbind(Fkernel.discr, cbind(rep(0, length.out = n),Fkernel.cont))
  
  mat <-Pkernel+Fkernel
  
  eigenMat <- eigen(mat)
  
  site_IPMs_all_DI[[i]] <- mat
}
names(site_IPMs_all_DI) <- paste0(unique(dat_all$Site))
lambdas_site_all_DI <- as.numeric(sapply(site_IPMs_all_DI, FUN = function(x) eigen(x)$values[1]))

#### load the bootstrap CI data ####
# for subpop-level IPMs w/ DI (all transitions, no env covariates)
subPop_bootCI_lambdas <- readRDS("~/Dropbox/Grad School/Research/Oenothera coloradensis project/COBP_analysis/intermediate_analysis_Data/site_level_IPMs_allYears/site_level_DI_bootCI_lambdas.RDS")
names(subPop_bootCI_lambdas) <- unique(dat_all$Site)
subPop_bootCI_params <- readRDS("~/Dropbox/Grad School/Research/Oenothera coloradensis project/COBP_analysis/intermediate_analysis_Data/site_level_IPMs_allYears/site_level_DI_bootCI_params.RDS")
names(subPop_bootCI_params) <- unique(dat_all$Site)

# get the lambdas into a matrix so I can make a correlation matrix
CI_lambdas <- (as.data.frame(sapply(subPop_bootCI_lambdas, function(x) unlist(as.numeric(x)))))
#reorder columns
CI_lambdas <- CI_lambdas[,c("Crow_Creek","Diamond_Creek",  "HQ5",  "HQ3",  "Meadow", "Unnamed_Creek")]
lambdaCor <- cor(log(CI_lambdas), method= "spearman")       

#### get spatial data for subPlot location ####
library(tidyverse)
library(sf)

# load dataset
setwd("/Users/Alice/Dropbox/Grad School/Research/Oenothera coloradensis project/Raw Data")
counts<- read.csv("./COBP_data_10_25_20.csv", stringsAsFactors = FALSE) #will have to update file name as it changes w/ most current version
sites <- read.csv("../Raw Data/COBP Plot Locations.csv", stringsAsFactors = FALSE)

#### get spatial plot locations ####
#order the sites in alphabetical order by site name 
sites$site_2 <- str_sub(sites$Site, start = str_locate(sites$Site, " ")[,1]+1, end = str_length(sites$Site)) #make a new variable that shows just the site name (without the location)
sites <- sites[order(sites$site_2),]
sites_NoSeed <- sites %>% 
  filter(site_2 != "Seed")
# make site data spatial 
sites_noSeed <- st_as_sf(sites_NoSeed, coords = c("Long_WGS84","Lat_WGS84"))
st_crs(sites_noSeed) <- 4326

# find the centroid of each subpopulation
for (i in unique(sites_noSeed$site_2)) {
  if (i == "Crow Creek"){
    site_subP <- data.frame("SubPop" = i,
                            geometry = st_centroid(
                              st_union(
                                sites_noSeed[sites_noSeed$site_2 == i,])))
  } else {
    site_subP <- rbind(site_subP,
                       data.frame("SubPop" = i,
                                  geometry = st_centroid(
                                    st_union(
                                      sites_noSeed[sites_noSeed$site_2 == i,])))
    )
  }
}
site_subP <- st_as_sf(site_subP)

## get a matrix of distances between subpopulations (in meters)
# rename subPops to match lambda matrix
site_subP$SubPop <- c("Crow_Creek", "Diamond_Creek", "HQ5", "HQ3", "Meadow", "Unnamed_Creek")
plotDist_mat <- st_distance(site_subP[match(site_subP$SubPop, names(CI_lambdas)),])
colnames(plotDist_mat) <- site_subP$SubPop
rownames(plotDist_mat) <- site_subP$SubPop

dim(plotDist_mat)
dim(lambdaCor)

#### perform a mantel test! ####
library(vegan)
vegan::mantel(plotDist_mat, lambdaCor, permutations = 1000, method = "spearman")

#### compare correlations just for each population ####
soapDistMat <- st_distance(site_subP[site_subP$SubPop %in% c("Meadow", "HQ3", "HQ5"),])
soapLamMat <- cor(CI_lambdas[,c("Meadow", "HQ3", "HQ5")], method = "spearman")
mantel(soapDistMat, soapLamMat)

baseDistMat <- st_distance(site_subP[site_subP$SubPop %in% c("Crow_Creek", "Diamond_Creek", "Unnamed_Creek"),])
baseLamMat <- cor(CI_lambdas[,c("Crow_Creek", "Diamond_Creek", "Unnamed_Creek")], method = "spearman")
mantel(baseDistMat, baseLamMat)


