#///////////////////////////////////////////////////////
# Integral Projection Models for Oenothera coloradensis
# Part 1: Data prep for IPMs
# Alice Stears
# 3 December 2021
#///////////////////////////////////////////////////////

#### Load packages ####
library(tidyverse)

#### Load Data ####
# location of data files 
datFolder <- "/Users/alicestears/Dropbox/Grad School/Research/Oenothera coloradensis project"
## load continuous data
dat <- read.csv(paste0(datFolder, "/Processed_Data/COBP_long_CURRENT.csv"),
                sep = ",") 

## get seedling data
seedlings <- read.csv(paste0(datFolder, "/Raw Data/COBP_seedlings_8_23_21.csv"), sep = ",") %>% 
  dplyr::select(Plot_ID, Quadrant, Seedlings_18, Seedlings_19, Seedlings_20) %>% 
  group_by(Plot_ID, Quadrant) %>% 
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

#### make seedling data into 'continuous' data ####
# add the seedling data with each row representing a seedling
for (i in 1:nrow(seedlings)) {
  seedlingNow <- seedlings[i,]
  # make a d.f to rbind to the seedlings_long d.f. that contains the appropriate data and # of rows for the number of seedlings in the data
  if (seedlingNow$Seedlings_t > 0) {
    temp <- data.frame("Location" = NA, "Site" = rep(seedlingNow$Site, length.out = seedlingNow$Seedlings_t),
                       "Plot_ID" = seedlingNow$Plot_ID,
                       "Quadrant" = seedlingNow$Quadrant, "ID" = NA, "X_cm" = NA, "Y_cm" = NA, "Year" = seedlingNow$Year, "LongestLeaf_cm" = NA, "survives_t" = 1, "flowering" = 0, "Num_capsules" = NA, "Stem_Herb" = NA, "Invert_Herb" = NA, "LeafSpots" = NA, "survives_tplus1" = NA, "longestLeaf_tminus1" = NA, "longestLeaf_tplus1" = NA,  "age" = 0,  "seedling" = 1, "index" = NA, "Num_seeds" = NA, "log_LL_t" = NA, "log_LL_tplus1" = NA, "log_LL_tminus1" = NA, "N_seedlings_t" = seedlingNow$Seedlings_t, "N_adults_t" = NA, "N_all_t" = NA, "recruit" = NA)
    if (i == 1) {
      # add this data into the seedlings_long d.f
      seedlings_long <- temp
    } else {
      seedlings_long <- rbind(seedlings_long, temp)
    }
  }
}

# assign the seedlings a 'random' size from a uniform distribution (between .1 and 3 cm longest leaf size)
set.seed(12011993)
sizes <- runif(4425, min = 0.1, max = 3) # sample from a dist. skewed to the right, partially b/c I can't get the model to run without it, and partially b/c the seedlings actually were more likely to be larger (is multiplied by 3, since we want size to range from 0 to 3)
seedlings_long$LongestLeaf_cm = sizes

# # try setting size according to a normal distribution
# sizeTemp <- seq(from = 0.1, to = 3, by = .001)
# sizes <- rnorm(n = 4425, mean = 1.55, sd = 0.5)
# # truncate the sizes between .1 and 3 
# sizes[which(sizes<.1)] <- 0.1
# sizes[which(sizes>3)] <- 3
# 
# plot(density(rnorm(n = 4425, mean = 1.55, sd = 0.5)))
# seedlings_long$LongestLeaf_NORM = sizes

# try setting size according to a beta distributio
sizes <- 3 * rbeta(n = 4425, shape1 = 3, shape2 = 1)

seedlings_long$LongestLeaf_NORM = sizes

# because we know the number of individuals that were 'recruited' into the adult plant stage in each t+1, we can estimate the number of seedlings that survived from t in each quadrat. We will randomly assign the these 'new recruits' to a seedling 
# get the data for 'recruits' to the adult stage
recruitsTemp <- dat[dat$recruit == 1,]
# loop through by quadrat/year
quads <- unique(dat$Plot_ID)
quadrants <- unique(dat$Quadrant)
years <- unique(dat$Year)
for (i in 1:length(quads)) {
  # can't get recruit data for 2018,b/c we can't know how old the plants were in the first year of sampling
  for (k in 1:length(quadrants)) {
    for (j in 1:(length(years)-1)) {
      # get recruit data for this quad and the NEXT year
      recsNow <- recruitsTemp %>% 
        filter(Plot_ID == quads[i] & Quadrant == quadrants[k] & Year == years[j+1])
      # get seedling data for this quad and the CURRENT year
      seedsNow <- seedlings_long %>% 
        filter(Plot_ID == quads[i] & Quadrant == quadrants[k] & Year == years[j])
      # if # of seedlings is GREATER than # of recruits...
      if (nrow(seedsNow) >= nrow(recsNow)) {
        numSdlngsDead <- nrow(seedsNow) - nrow(recsNow)
        # add plant "ID" to seedlings that survived
        seedsNow$ID <- c(recsNow$ID,
                         rep(NA, length.out = numSdlngsDead))
        # add a '1' to "survives_tplus1" column for seedlings that survived
        seedsNow$survives_tplus1 <- c(rep(1, length.out = nrow(recsNow)),
                                      rep(0, length.out = numSdlngsDead))
        # add plant size in "longestLeaf_tplus1" to seedlings that survived
        seedsNow$longestLeaf_tplus1 <- c(recsNow$LongestLeaf_cm,
                                         rep(NA, length.out = numSdlngsDead))
        
      } else if (nrow(seedsNow) < nrow(recsNow)) { # if # of seedlings is LESS than # of recruits...
        # add a necessary number of rows to the seedling d.f to make it the same length as the recruit d.f
        numNewSdlngs <- nrow(recsNow) - nrow(seedsNow)
        
        # make a d.f of the new seedlings
        temp <- data.frame("Location" = NA, "Site" = rep(recsNow$Site, length.out = numNewSdlngs),
                           "Plot_ID" = quads[i],"Quadrant" = rep(recsNow$Quadrant, length.out = numNewSdlngs), "ID" = NA, "X_cm" = NA, 
                           "Y_cm" = NA, "Year" = years[j], 
                           "LongestLeaf_cm" = runif(numNewSdlngs, min = 0.1, max = 3),
                           "survives_t" = 1, "flowering" = 0, "Num_capsules" = NA, "Stem_Herb" = NA,
                           "Invert_Herb" = NA, "LeafSpots" = NA, "survives_tplus1" = NA,
                           "longestLeaf_tminus1" = NA, "longestLeaf_tplus1" = NA,  "age" = 0,  
                           "seedling" = 1, "index" = NA, "Num_seeds" = NA, "log_LL_t" = NA, 
                           "log_LL_tplus1" = NA, "log_LL_tminus1" = NA, 
                           "N_seedlings_t" = NA, "N_adults_t" = NA, 
                           "N_all_t" = NA, "recruit" = NA, 
                           "LongestLeaf_NORM" = rnorm(numNewSdlngs, mean = 1.55, sd = .5))
        
        # truncate the estimated seedling sizes from the Normal distribution to between .1 and 3 
        temp[temp$LongestLeaf_NORM<.1, "LongestLeaf_NORM"] <- 0.1
        temp[temp$LongestLeaf_NORM>3, "LongestLeaf_NORM"] <- 3
        # then add this to the seedsNow d.f
        if (nrow(seedsNow) == 0) {
          seedsNow <- temp
        } else {
          seedsNow <- rbind(seedsNow, temp)
        }
        # add plant "ID" to seedlings that survived
        seedsNow$ID <- recsNow$ID
        # add a '1' to "survives_tplus1" column for seedlings that survived
        seedsNow$survives_tplus1 <- 1
        # add plant size in "longestLeaf_tplus1" to seedlings that survived
        seedsNow$longestLeaf_tplus1 <- recsNow$LongestLeaf_cm
      }
      # save the results
      if (i == 1 & j == 1 & k == 1) {
        seedlings_cont <- seedsNow
      } else {
        seedlings_cont <- rbind(seedlings_cont, seedsNow)
      }
    }
  }
  }
  
## add the 2020 seedlings into the data (don't have survival, since we can't know whether they survive or not)
seedlings_2020 <- seedlings_long[seedlings_long$Year==2020,]
## add to the 'seedlings_cont' df
seedlings_cont <- rbind(seedlings_cont, seedlings_2020)

### TEMPORARILY change the "seedlings_cont$longestLeaf_cm" data to have the 
# seedling size information from a NORMAL DISTRIBUTION (rather than a UNIFORM DISTRIBUTION like it did before)
seedlings_cont <- seedlings_cont %>% 
  dplyr::select(-LongestLeaf_cm) %>% 
  rename(LongestLeaf_cm = LongestLeaf_NORM)

## "seedlings_cont" contains the continuous seedling data
# there are only values for 2018 and 2019, since we can't know whether seedlings from 2020 survived or not
# make sure both dfs have the same column order

datTemp <- dat[,names(seedlings_cont[,1:25])]
# add seedlings_cont to the dat dataframe
dat_all <- rbind(datTemp, seedlings_cont[,1:25])

## make sure that the "location" data in dat_all is correct
dat_all[dat_all$Site %in% c("Crow_Creek", "Diamond_Creek", "Unnamed_Creek"),"Location"] <- "FEWAFB"

dat_all[dat_all$Site %in% c("HQ5", "HQ3", "Meadow"),"Location"] <- "Soapstone"

## get the environmental data
source("./data_prep_scripts/prepping_Envi_Data.R")
# called "Climate", "soilTemp_plot", and "soilMoist"
## add the envi data to the 'dat' data.frame
# add soil moisture data (have that for every plot, from one time point)
dat_all <- left_join(dat_all, soilMoist[,c("Site", "SoilMoisture_m3m3")], by = c("Plot_ID" = "Site"))
# add soil temperature data (have mean values for every site, from one time point)
# reformat the data to have means for each site 
soilTemp <- soilTemp_plot %>% 
  group_by(Site, Location) %>% 
  summarize(SoilTemp_winter_C = mean(SoilTemp_winter_C, na.rm = TRUE),
            sd_soilTemp_winter = mean(sd_soilTemp_winter, na.rm = TRUE),
            SoilTemp_grow_C = mean(SoilTemp_grow_C, na.rm = TRUE), 
            sd_soilTemp_grow = mean(sd_soilTemp_grow, na.rm = TRUE))
dat_all <- left_join(dat_all, soilTemp)
# add the climate data (have mean values for each location)
dat_all$Year <- as.double(as.character(dat_all$Year))
dat_all <- left_join(dat_all, Climate)
# log-transform size variables
dat_all$log_LL_t <- log(dat_all$LongestLeaf_cm)
dat_all$log_LL_tplus1 <- log(dat_all$longestLeaf_tplus1)

# scale environmental variables
dat_all$Year <- as.factor(dat_all$Year)
dat_all$SoilMoisture_m3m3_s <- scale(dat_all$SoilMoisture_m3m3)
dat_all$SoilTemp_winter_C_s <- scale(dat_all$SoilTemp_winter_C)
dat_all$SoilTemp_grow_C_s <- scale(dat_all$SoilTemp_grow_C)
dat_all$tMean_grow_C_s <- scale(dat_all$tMean_grow_C)
dat_all$precipWaterYr_cm_s <- scale(dat_all$precipWaterYr_cm)

#### get the number of recruits/year ####
estabTemp <- dat %>% 
  dplyr::select(Plot_ID, Year, Quadrant, recruit) %>% 
  group_by(Plot_ID, Quadrant, Year) %>% 
  summarize( recruits = sum(recruit)) %>% 
  rename( recruits_t = recruits)
# change the 2018 recruit values to NA, since we don't have a count of recruits to the rosette stage for that year
estabTemp[estabTemp$Year==2018, "recruits_t"] <- NA
# get the number of recruits in year t+1
estabTemp$Year <- as.numeric(as.character(estabTemp$Year)) - 1
names(estabTemp)[4] <- "recruits_tplus1"
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
## data from Burgess, Hild & Shaw, 2005--NOT using data from MT place, only from FEWAFB
# Seeds per Capsule (total no., viable and inviable)
seed_per_cap <- mean(c(2.4, 1.0))
# Capsule Viability(%) (percentage of capsules that are viable--contain >=1 seed)
capsule_viab.rt <- mean(c(81, 61, 54))/100
# Seed Viability (%) (percentage of seeds in a viable capsule that are viable)
seed_viab.rt <- mean(c(1.9/2.4,  1/1))

## calculate the number of seeds based on this seed rate 
dat$Num_seeds <- round(dat$Num_capsules * seed_per_cap,0)

## calculate the rate at which a seed produced in a capsule in year t is viable (probability of a viable capsule * probability that a seed inside a viable capsule is viable) 
total_seed_viab.rt <- capsule_viab.rt * seed_viab.rt

# (data in SeedBagGreenhouseSeedlings.csv)--seeds from previous year
germ.rt.ours <- .03
#data from (Burgess, Hild & Shaw, 2005)--seedbank seed viability/germination rate doesn't seem to change much over time-- only from WAFB sites
germ.rt.Burgess <- mean(c(13.0, 12, 8.3, 7.0, 5.3)/(45 * seed_per_cap))

# germination rate from the Burgess paper incorporates both viability and germination. To isolate just the germination rate, divide the Burgess germination rate by the viability rate
germ.rt_temp <-germ.rt.Burgess/total_seed_viab.rt  
# because the Burgess germ.rate was estimated in a greenhouse, and we are confident that field rates are lower, multiply teh Burgess germ.rate by .80
germ.rt <- germ.rt_temp * .8

# the viability rate (proportion of seeds produced by an adult plant that are viable) is the 'total_seed_viab.rt' derived from the Burgess paper results
viab.rt <- total_seed_viab.rt

#### calculate the starting number of seeds in the seedbank ####
# data from COBPGreenhouseSeedbankStudyData.csv
# the no. of seeds that germinated for sample of each plot ()
seedStudyGerms <- c(5,4,3,0.01,0.01,1,0.01,0.01,0.01,9,4,9,0,1,1,8,20,1)
# divide the seedcore germs by the germ rate to get the number of (seedbank + seed rain) seeds in each core
seedStudyGerms.1 <- seedStudyGerms/germ.rt
# convert the volume of the soil samples to the volume of the plots to get approx. number of (seedbank + seed rain) seeds in each plot
# cores are 5.5cm in diameter and 3cm deep, and there are 20 of them/plot
coreVol = ((pi*(5.5/2)^2)*2.5)*20 ## 1187.915 cm^3 is the total volume of the soil cores for one plot
# calculate the volume of the first 3cm of soil in the plots
plotVol = 200*200*3
# estimate the number of seeds in the plot
seedsEst.all <- (plotVol*seedStudyGerms.1)/coreVol

# total number of seeds produced by adult plants in each plot in each year
seedsTemp <- aggregate(x = dat[,c("Num_seeds")], by = list("Plot_ID" = dat$Plot_ID,
                                                           "Year" = dat$Year, 
                                                           "Site" = dat$Site), FUN = sum, na.rm = TRUE)
names(seedsTemp)[4] <- "Num_seeds"

#get the mean seed rain in a year for each plot 
seedRain.avg <- aggregate(x = seedsTemp$Num_seeds, by = list("Plot_ID" = seedsTemp$Plot_ID, 
                                                     "Site" = seedsTemp$Site), FUN = mean)
# re-order the plot seedRain data to be in the same order as the seed core data
seedRain.avg <- seedRain.avg[c(7:18,1:3,6,4:5),]
# put all of this data into one d.f
names(seedRain.avg)[3] <- "SeedRain_avg"
# add seedbank estimates from the soil cores to this d.f
seedRain.avg$coreSeeds_est <- seedsEst.all
# put in a d.f called 'seeds_est'
seeds_est <- seedRain.avg
# subtract the seedRain from the soilcore seed est to get an estimate of the size of the seedbank state in year t (B(t))
seeds_est$seedbank_est <- seeds_est$coreSeeds_est - seeds_est$SeedRain_avg
# aggregate to the site level (site-wide totals, NOT means)
# get total number of seeds in the seedbank (estimated) for each site (only have one value, so need to extrapolate across years)
seeds_est_site <- seeds_est %>% 
  group_by(Site) %>% 
  summarize(SeedRain_avg = sum(SeedRain_avg), 
            coreSeeds_est = sum(coreSeeds_est), 
            seedbank_est = sum(seedbank_est))
# change the negative numbers to 100 %%% MAYBE SHOULD CHANGE THIS... JUST A TEMPORARY FIX%%%
seeds_est_site[seeds_est_site$seedbank_est < 0, "seedbank_est"] <- 100

# calculate a point-estimate of the mean seedbank size across all sites
# formula = (germs*(1-germ.rate) - mean no. of seeds produced in previous year)
seedBank_est <- mean(seeds_est$coreSeeds_est - seeds_est$SeedRain_avg)

## calculate average number of seedlings/year for each site (NOT mean seedling no.)--mean for each plot across all years, then summed values for each site
seedlings_site <- seedlings %>% 
  group_by(Plot_ID, Site) %>% 
  summarise(Seedlings_t = mean(Seedlings_t)) %>% 
  group_by(Site) %>% 
  summarise(Seedlings_t = sum(Seedlings_t))


# change seedbank estimate that are negative to small positive numbers
seeds_est_faked <- seeds_est 
seeds_est_faked[seeds_est_faked$seedbank_est < 0, "seedbank_est"] <- 10

#### make d.f to illustrate the number of seeds/seedling/germs, etc. ####
plots <- unique(dat$Plot_ID)
years <- unique(dat$Year)
if (exists("seeds.out")) {
  remove("seeds.out")
}
for (i in 1:length(plots)) {
  plot_now <- plots[i]
  for (j in 1:length(years)) {
    year_now <- years[j]
    
    ## number of seedbank seeds in year t
    seeds_now <- data.frame("Year" = year_now, "Plot_ID" = plot_now, "NewSeeds_t" = 0, "SeedBank_t" = rep_len(1, length.out = round(seeds_est_faked[seeds_est_faked$Plot == plot_now , "seedbank_est"], 0)), "SeedBank_tplus1" = 0, "Seedling_tplus1" = 0, "Seedling_t" = 0, "Recruit_tplus1" = 0)
    # seeds that stay in the seedbank (staySB) (1 - germ.rt) * 0.9
    seeds_now[1:(round((1-germ.rt)*0.9 * nrow(seeds_now),0)),"SeedBank_tplus1"] <- 1
    # seeds that leave the seedbank (outSB) (germ.rt)
    seeds_now[seeds_now$SeedBank_tplus1 == 0,][1:round((germ.rt*.9) * nrow(seeds_now),0), "Seedling_tplus1"] <- 1 
    
    ## seeds that enter the seedbank from the continuous stage (goSB and goSdlng)
    n_newSeeds = sum(dat[dat$Plot_ID==plot_now & dat$Year == year_now,"Num_seeds"], 
                     na.rm = TRUE)
    # get the number of seeds from the reproductive plants that go into the seedbank (viab.rt - germ.rt)
    n_newGoSB = round((viab.rt*(1 - germ.rt))*n_newSeeds, 0)
    # get the number of seeds from the reproductive plants that go into the seedling stage (germ.rt)
    n_newGoSdlng = round(n_newSeeds*viab.rt * germ.rt,0)
    
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
                                                          rep_len(0, length.out = (n_newSeeds - (n_newGoSB + n_newGoSdlng)))), #1 (died) 
                                    "Seedling_t" = 0, "Recruit_tplus1" = 0)) 
    }
    # probability of establishing into a rosette in year t+1 (do for each plot, not each plot/year combo)
    # get number of recruits in year t+1
    if (year_now %in% c(2018, 2019)) {
      n_recruits_tplus1 <- sum(dat[dat$Plot_ID == plot_now & dat$Year == as.numeric(as.character(year_now)) + 1,"recruit"], na.rm = TRUE)
    } else if (year_now == 2020) {
      n_recruits_tplus1 <- 0
    }
    # get the number of seedlings in year t
    n_seedlings_t <- sum(deframe(seedlings[seedlings$Plot_ID == plot_now & 
                                 seedlings$Year == year_now, "Seedlings_t"]))
    
    if (n_seedlings_t > n_recruits_tplus1) {
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
    } else if (n_seedlings_t <= n_recruits_tplus1) {
      # adjust the number of seedlings to be the same as the number of recruits
      n_seeds_adj <- (n_recruits_tplus1 + 1)
      seeds_now <- rbind(seeds_now, 
                         data.frame("Year" = year_now, "Plot_ID" = plot_now, 
                                    "NewSeeds_t" = 0,
                                    "SeedBank_t" = 0,
                                    "SeedBank_tplus1" = 0,
                                    "Seedling_tplus1" = 0,
                                    "Seedling_t" = rep_len(1, length.out = n_seeds_adj), 
                                    "Recruit_tplus1" = c(rep_len(1, length.out = n_recruits_tplus1),
                                                         rep_len(0, length.out = (n_seeds_adj - n_recruits_tplus1))
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
discDat <- seeds.out

####calculate the total population size for each quad/year combo ####
N_all <- dat_all %>% 
  group_by(Plot_ID, Year) %>% 
  summarize(N_all = n())
# add to the dat_all d.f
dat_all <- dat_all %>% left_join(N_all)
# make sure log(longest leaf) is calculated for all plants
dat_all$log_LL_t <- log(dat_all$LongestLeaf_cm)
# make sure pop size values are correct
N_all_dat <- dat_all %>% group_by(Site, Plot_ID, Year) %>% 
  summarize(N_all_t = n())
N_adultSeedling_dat <- dat_all %>% group_by(Site, Plot_ID, Year, seedling) %>% 
  summarize(N = n()) %>% 
  pivot_wider(id_cols = c(Plot_ID, Year), names_from = seedling, values_from = N) %>% 
  rename(N_adults_t = `0`, N_seedlings_t = `1`)

N_dat <- left_join(N_all_dat, N_adultSeedling_dat)
## get site-level N 
N_site <- dat_all %>% group_by(Site, Year) %>% 
  summarize(N_Site_t = n())
N_dat <- left_join(N_dat, N_site)

## add pop. size data to 'dat_all' df
dat_all <- dat_all  %>% 
  left_join(N_dat)

# # write the discreteDat d.f to file
#write.csv(x = discDat, file = "../Processed_Data/discreteStageData.csv", row.names = FALSE)
# # also write the continuous seedling d.f to file
#write.csv(x = dat_all, file = "../Processed_Data/allDat_plus_contSeedlings.csv", row.names = FALSE)
