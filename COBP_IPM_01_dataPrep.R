#///////////////////////////////////////////////////////
# Integral Projection Models for Oenothera coloradensis
# Part 1: Data prep for IPMs
# Alice Stears
# 3 December 2021
#///////////////////////////////////////////////////////

#### Load packages ####
library(tidyverse)

#### Load Data ####
## load spatialized data, which we'll need for neighborhood calculations
dat <- read.csv("../Processed_Data/COBP_long_CURRENT.csv")

## get seedling data
seedlings <- read.csv("../Raw Data/COBP_seedlings_8_23_21.csv") %>% 
  dplyr::select(Plot_ID, Seedlings_18, Seedlings_19, Seedlings_20) %>% 
  group_by(Plot_ID) %>% 
  summarize(Seedlings_18 = sum(Seedlings_18), Seedlings_19 = sum(Seedlings_19), Seedlings_20 = sum(Seedlings_20)) %>% 
  pivot_longer(cols = c(Seedlings_18, Seedlings_19, Seedlings_20), names_to = "Year", values_to = "Seedlings_t", names_pattern = "([[:digit:]]+)") %>% 
  mutate(Year = (as.numeric(Year) + 2000))
# add the 'site' data 
seedlings$Site <- NA
seedlings[seedlings$Plot_ID %in% c("C4", "C5", "C8"), "Site"] <- "Crow_Creek"
seedlings[seedlings$Plot_ID %in% c("D10", "D11", "D7"), "Site"] <- "Diamond_Creek"
seedlings[seedlings$Plot_ID %in% c("U3", "U4", "U6"), "Site"] <- "Unnamed_Creek"
seedlings[seedlings$Plot_ID %in% c("S1", "S2", "S3"), "Site"] <- "HQ5"
seedlings[seedlings$Plot_ID %in% c("S4", "S5", "S6"), "Site"] <- "HQ3"
seedlings[seedlings$Plot_ID %in% c("S7", "S8", "S9"), "Site"] <- "Meadow"

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

#### transform variables ####
# scale environmental variables
dat$Year <- as.factor(dat$Year)
dat$SoilMoisture_m3m3_s <- scale(dat$SoilMoisture_m3m3)
dat$SoilTemp_winter_C_s <- scale(dat$SoilTemp_winter_C)
dat$SoilTemp_grow_C_s <- scale(dat$SoilTemp_grow_C)
dat$tMean_grow_C_s <- scale(dat$tMean_grow_C)
dat$precipWaterYr_cm_s <- scale(dat$precipWaterYr_cm)
# make log-transformed size variables ####
dat$log_LL_t <- log(dat$LongestLeaf_cm)
dat$log_LL_tplus1 <- log(dat$longestLeaf_tplus1)
dat$log_LL_tminus1 <- log(dat$longestLeaf_tminus1)
# round capsule count numbers (that were modeled based on regression) to whole numbers 
dat$Num_capsules <- round(dat$Num_capsules, digits = 0)
# remove spaces from the 'site' column
dat$Site <- str_replace(string = dat$Site, pattern = " ", replacement = "_")

#### calculate population size for each plot  ####
N_dat <- dat %>% 
  group_by(Plot_ID, Year) %>% 
  summarize("N_adults_t" = n())

# add seedling and adult N values together
N_dat <- seedlings %>% 
  rename(N_seedlings_t = Seedlings_t) %>% 
  mutate(Year = as.factor(Year)) %>% 
  left_join(N_dat) %>% 
  mutate(N_all_t = (N_seedlings_t + N_adults_t))
# add the 'N' data to the 'dat' data.frame
dat <- left_join(dat, N_dat)
# scale density dependence variables
dat$N_all_t_s <- scale(dat$N_all_t)

#### get establishment/recruit size data ####
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
estabTemp$Year <- as.numeric(as.character(estabTemp$Year)) - 1
names(estabTemp)[3] <- "recruits_tplus1"
# get the no.of seedlings/year
# load the seedling data ('seedlings' df)
# combine seedling and rosette recruit data
estabs <- left_join(estabTemp, seedlings)
# calculate the probability of seedling in year t establishing to a rosette in year t+1
estabs$P_estab <- estabs$recruits_tplus1/estabs$Seedlings_t
# fix the 'inf' values
estabs[is.infinite(estabs$P_estab),"P_estab"] <- NA
# if the value is >1, round back to one (must have missed some seedlings)
estabs[estabs$P_estab > 1 & is.na(estabs$P_estab) == FALSE,"P_estab"] <- 1
estabs <- as.data.frame(estabs)

#### Calculate germination rate and rate of seed viability loss ####
# (data in SeedBagGreenhouseSeedlings.csv)--seeds from previous year
germ.rt.ours <- .03
#data from (Burgess, Hild & Shaw, 2005)--seedbank seed viability/germination rate doesn't seem to change much over time
germ.rt.Burgess <- mean(c(16.0, 13.0, 12, 8.3, 7.0, 5.3)/45)
germ.rt <-.05  #average together, but trending toward our estimate
# seed viability rate (from Burgess, Hild & Shaw, 2005)
viab.rt.Burgess <- mean(c(.79,.7,.63))
# reduce by 20%
viab.rt <-viab.rt.Burgess*.8

#### calculate the starting number of seeds in the seedbank ####
# data from COBPGreenhouseSeedbankStudyData.csv
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
                                                           "Year" = dat$Year, 
                                                           "Site" = dat$Site), FUN = sum, na.rm = TRUE)
#get the mean seed rain in a year for each plot 
seedRain.avg <- aggregate(x = seedsTemp$x, by = list("Plot_ID" = seedsTemp$Plot_ID, 
                                                     "Site" = seedsTemp$Site), FUN = mean)
# re-order the plot seedRain data to be in the same order as the seed core data
seedRain.avg <- seedRain.avg[c(7:18,1:3,6,4:5),]
# put all of this data into one d.f
names(seedRain.avg)[3] <- "SeedRain_avg"
seedRain.avg$coreSeeds_est <- seedsEst.all
seeds_est <- seedRain.avg
# get an estimate of seedbank size for each plot (soil core seed est. - seedRain)
seeds_est$seedbank_est <- seeds_est$coreSeeds_est - seeds_est$SeedRain_avg
# subtract the seedRain from the soilcore seed est to get an estimate of the size of the seedbank state in year t (B(t))
# aggregate to the site level (site-wide totals, NOT means)
seeds_est_site <- seeds_est %>% 
  group_by(Site) %>% 
  summarize(SeedRain_avg = sum(SeedRain_avg), 
            coreSeeds_est = sum(coreSeeds_est), 
            seedbank_est = sum(seedbank_est))
# change the negative numbers to 100 %%% MAYBE SHOULD CHANGE THIS... JUST A TEMPORARY FIX%%%
seeds_est_site[seeds_est_site$seedbank_est < 0, "seedbank_est"] <- 100
# formula = (germs*(1-germ.rate) - mean no. of seeds produced in previous year)
seedBank_est <- mean(seeds_est$coreSeeds_est - seeds_est$SeedRain_avg)

## calculate total seedling no. by site (NOT mean seedling no.)--mean for each plot accross all years, then summed values for each site
seedlings_site <- seedlings %>% 
  group_by(Plot_ID, Site) %>% 
  summarise(Seedlings_t = mean(Seedlings_t)) %>% 
  group_by(Site) %>% 
  summarise(Seedlings_t = sum(Seedlings_t))


# change seedbank estimate that are negative to small positive numbers
seeds_est_faked <- seeds_est 
seeds_est_faked[seeds_est_faked$seedbank_est < 0, "seedbank_est"] <- 10

#### make d.f.s to illustrate the number of seeds/seedling/germs, etc. ####
plots <- unique(dat$Plot_ID)
years <- unique(dat$Year)
if (exists("seeds.out")) {
  remove("seeds.out")
}
for (i in 1:length(plots)) {
  plot_now <- plots[i]
  for (j in 1:length(years)) {
    year_now <- years[j]
    
    seeds_now <- data.frame("Year" = year_now, "Plot_ID" = plot_now, "NewSeeds_t" = 0, "SeedBank_t" = rep_len(1, length.out = round(seeds_est_faked[seeds_est_faked$Plot == plot_now , "seedbank_est"], 0)), "SeedBank_tplus1" = 0, "Seedling_tplus1" = 0, "Seedling_t" = 0, "Recruit_tplus1" = 0)
    # seeds that stay in the seedbank (staySB) (prob = (1-germ.rt)*viab.rt) 
    seeds_now[1:(round((1-germ.rt)*viab.rt * nrow(seeds_now),0)),"SeedBank_tplus1"] <- 1
    # seeds that leave the seedbank (outSB) (germ.rt)
    seeds_now[seeds_now$SeedBank_tplus1 == 0,][1:round(germ.rt * nrow(seeds_now),0), "Seedling_tplus1"] <- 1 
    # seeds that enter the seedbank from the continuous stage (goSB and goSdlng)
    n_newSeeds = sum(dat[dat$Plot_ID==plot_now & dat$Year == year_now,"Num_seeds"], 
                     na.rm = TRUE)
    n_newGoSB = round(((1-germ.rt)*viab.rt)*n_newSeeds, 0)
    n_newGoSdlng = round(n_newSeeds*germ.rt,0)
    if (n_newSeeds > 0) {
      seeds_now <- rbind(seeds_now, 
                         data.frame("Year" = year_now, "Plot_ID" = plot_now, 
                                    "NewSeeds_t" = rep_len(
                                      1, length.out = n_newSeeds), # number of seeds produced in this plot/year combo
                                    "SeedBank_t" = 0, # can't go into the seedbank in the current year 
                                    "SeedBank_tplus1" = c(rep_len(1, length.out = n_newGoSB),
                                                          rep_len(0, length.out = (n_newSeeds-n_newGoSB))), # could go into the seedbank (goSB) (1-germ.rt)*viab.rt
                                    "Seedling_tplus1" = c(rep_len(0, length.out = n_newGoSB), # 0 (went into seedbank)
                                                          rep_len(1, length.out = n_newGoSdlng), # 1 (became a seedling)
                                                          rep_len(0, length.out = (n_newSeeds - (n_newGoSB + n_newGoSdlng)))), #1 (died) could go directly to the seedling stage (goSdlng) (germ.rt)
                                    "Seedling_t" = 0, "Recruit_tplus1" = 0)) 
    }
    # probability of establishing into a rosette in year t+1 (do for each plot, not each plot/year combo)
    recruitRate_now <- mean(estabs[estabs$Plot_ID==plot_now,"P_estab"], na.rm = TRUE) # establishment rate
    if (is.nan(recruitRate_now)) {
      recruitRate_now <- NA
    }
    # get number of recruits
    n_recruits_tplus1 <- sum(dat[dat$Plot_ID == plot_now & dat$Year == as.numeric(as.character(year_now)) + 1,"recruit"], na.rm = TRUE)
    if (is.na(recruitRate_now) == TRUE & n_recruits_tplus1 > 0) {
     recruitRate_now <- .99
    }
    n_seedlings_t <- round(n_recruits_tplus1/recruitRate_now, 0)
    if (n_seedlings_t > 0 & is.na(n_seedlings_t) == FALSE) {
      seeds_now <- rbind(seeds_now, 
                         data.frame("Year" = year_now, "Plot_ID" = plot_now, 
                                    "NewSeeds_t" = 0,
                                    "SeedBank_t" = 0,
                                    "SeedBank_tplus1" = 0,
                                    "Seedling_tplus1" = 0,
                                    "Seedling_t" = rep_len(1, length.out = n_seedlings_t), 
                                    "Recruit_tplus1" = c(rep_len(1, length.out = n_recruits_tplus1),
                                                         rep_len(0, length.out = (n_seedlings_t - n_recruits_tplus1))
                                    ))) 
    } 
    
    ## add to output d.f
    if (exists("seeds.out")) {
      seeds.out <- rbind(seeds.out, seeds_now)
    } else {
      seeds.out <- seeds_now
    }
  }
}

# put in a 'discrete' stage d.f
discreteDat <- seeds.out

# write the discreteDat d.f to file
write.csv(x = discreteDat, 
          file = "../Processed_Data/discreteStageData.csv", row.names = FALSE)
