#//////////////////////////
# Integral Projection Models for Oenothera coloradensis: Analysis for spatial asynchrony
# Alice Stears
# 11 March 2022
#/////////////////////////

library(tidyverse)
library(rstatix)

# load data from script 1
dat_all <- read.csv(file = "../Processed_Data/allDat_plus_contSeedlings.csv")

#### IPMs for each subpop over both transitions ####
## these are IPMs C-H
# load IPMs
IPMs_C_H <- readRDS(file = "./intermediate_analysis_Data/IPMs_C_H.RDS")
site_IPMs_all_DI <- IPMs_C_H
# calculate lambdas
lambdas_site_all_DI <- as.numeric(sapply(site_IPMs_all_DI, FUN = function(x) log(eigen(x)$values[1])))

#### load the bootstrap CI data ####
# for subpop-level IPMs w/ DI (all transitions, no env covariates)
subPop_bootCI_lambdas <- readRDS("./intermediate_analysis_Data/IPMs_C_H_bootCI_lambdas.RDS")
names(subPop_bootCI_lambdas) <- unique(dat_all$Site)
# convert to log-lambdas
subPop_bootCI_lambdas <- lapply(subPop_bootCI_lambdas, FUN = function(x) log(as.numeric(x)))

subPop_bootCI_params <- readRDS("./intermediate_analysis_Data/IPMs_C_H_bootCI_params.RDS")
names(subPop_bootCI_params) <- unique(dat_all$Site)

# get the lambdas into a matrix to then make a correlation matrix
CI_lambdas <- (as.data.frame(sapply(subPop_bootCI_lambdas, function(x) unlist(as.numeric(x)))))
#reorder columns
CI_lambdas <- CI_lambdas[,c("Crow_Creek","Diamond_Creek",  "HQ5",  "HQ3",  "Meadow", "Unnamed_Creek")]
lambdaCor <- cor(log(CI_lambdas), method= "spearman")       

#### get spatial data for subPlot location ####
library(tidyverse)
library(sf)

# load dataset
counts<- read.csv("../Raw Data/COBP_data_10_25_20.csv", stringsAsFactors = FALSE) #will have to update file name as it changes w/ most current version
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


