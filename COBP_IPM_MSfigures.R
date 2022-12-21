#//////////////////////////
# Integral Projection Models for Oenothera coloradensis: Figures
# Alice Stears
# 11 March 2022
#/////////////////////////
library(tidyverse)
library(ggpubr)
library(leaflet)
library(sf)
library(mapedit)
library(ggmap)
library(ggspatial)
library(ggrepel)

#### figure of vital rates (all data, DI, all transitions) ####
datSoap <- dat_all[dat_all$Location=="Soapstone",]
survDat_soap <- datSoap[datSoap$flowering==0 | is.na(datSoap$flowering),]
survMod_soap <- glm(survives_tplus1 ~ log_LL_t , data = survDat_soap, family = binomial)
newdata <- data.frame("log_LL_t" = seq(from = min(survDat_all$log_LL_t, na.rm = TRUE), 
                                       to = max(survDat_all$log_LL_t, na.rm = TRUE),
                                       length.out = 100))
lines(x = newdata$log_LL_t, y = predict(object = survMod_all, newdata =  newdata, type = "response"), col = "red")

## Growth ($G(z',z)$)
# lm w/ log-transformed size_t and size_t+1
sizeMod_soap <- lm(log_LL_tplus1 ~ log_LL_t , data = datSoap)
## Number of seeds produced, according to plant size ($b(z)$)
seedDat_soap <- datSoap[datSoap$flowering == 1,]
# fit a negative binomial glm (poisson was overdispersed)
seedMod_soap <- MASS::glm.nb(Num_seeds ~ log_LL_t , data = seedDat_soap)
## Flowering probability ($p_b(z)$)
flwrMod_soap <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = datSoap, family = binomial)))
## Distribution of recruit size ($c_o(z')$)
recD_soap <- datSoap[datSoap$seedling == 1,]
recMod_soap <- lm(log_LL_t ~ 1, data = recD_soap)

datBase <- dat_all[dat_all$Location=="FEWAFB",]
survDat_Base <- datBase[datBase$flowering==0 | is.na(datBase$flowering),]
survMod_Base <- glm(survives_tplus1 ~ log_LL_t , data = survDat_Base, family = binomial)
newdata <- data.frame("log_LL_t" = seq(from = min(survDat_all$log_LL_t, na.rm = TRUE), 
                                       to = max(survDat_all$log_LL_t, na.rm = TRUE),
                                       length.out = 100))
lines(x = newdata$log_LL_t, y = predict(object = survMod_all, newdata =  newdata, type = "response"), col = "red")

## Growth ($G(z',z)$)
# lm w/ log-transformed size_t and size_t+1
sizeMod_Base <- lm(log_LL_tplus1 ~ log_LL_t , data = datBase)
## Number of seeds produced, according to plant size ($b(z)$)
seedDat_Base <- datBase[datBase$flowering == 1,]
# fit a negative binomial glm (poisson was overdispersed)
seedMod_Base <- MASS::glm.nb(Num_seeds ~ log_LL_t , data = seedDat_Base)
## Flowering probability ($p_b(z)$)
flwrMod_Base <- suppressWarnings((glm(flowering ~ log_LL_t + I(log_LL_t^2) , data = datBase, family = binomial)))
## Distribution of recruit size ($c_o(z')$)
recD_Base <- datBase[datBase$seedling == 1,]
recMod_Base <- lm(log_LL_t ~ 1, data = recD_Base)

# survival
survFig <- ggplot() + 
  geom_point(aes(x = dat_all$log_LL_t, y = dat_all$survives_tplus1), col = "grey40", alpha = .3) +
  geom_smooth(aes(x = log_LL_t, y = survives_tplus1), data = survDat_soap, method = "glm", method.args = list(family = "binomial"), col = '#00b159', fill = '#00b159',alpha = .4,  fullrange = TRUE) + 
  geom_smooth(aes(x = log_LL_t, y = survives_tplus1), data = survDat_Base, method = "glm", method.args = list(family = "binomial"), col = "#f37735", fill = "#f37735", alpha = .4, fullrange = TRUE)  +
  ylab(c("Probability(survival)")) +
  xlab(c("log(size_t)")) +
  ggtitle("Survival") +
  #geom_line(aes(x = newdata$log_LL_t, y = predict(object = survMod_soap, newdata =  newdata, type = "response")),  col = '#00b159', lwd = 1.5) +
  #geom_line(aes(x = newdata$log_LL_t, y = predict(object = survMod_Base, newdata =  newdata, type = "response")),  col = "#f37735", lwd = 1.5) +
  theme_classic()
# growth
sizeFig <- ggplot() +
  geom_point(aes(x = dat_all$log_LL_t, y = dat_all$log_LL_tplus1), col = "grey40", alpha = .3) +
  geom_smooth(aes(x = log_LL_t, y = log_LL_tplus1), data = datSoap, method = "lm", col = '#00b159', fill = '#00b159',alpha = .4,  fullrange = TRUE) + 
  geom_smooth(aes(x = log_LL_t, y =log_LL_tplus1), data = datBase, method = "lm", col = "#f37735", fill = "#f37735", alpha = .4, fullrange = TRUE) +
  ylab(c("log(size_t+1)")) +
  xlab(c("log(size_t)")) +
  ggtitle("Growth")+
  geom_abline(aes(intercept = 0, slope = 1), col = "grey20", lty = 2) +
  scale_colour_manual(values = c("#f69f71", '#4cc88a'),
                      guide_legend(title = "Population")) +
  #geom_line(aes(x = newdata$log_LL_t, y = predict(object = sizeMod_soap, newdata =  newdata, type = "response")),   col = '#00b159', lwd = 1.5) +
 # geom_line(aes(x = newdata$log_LL_t, y = predict(object = sizeMod_Base, newdata =  newdata, type = "response")),  col = "#f37735", lwd = 1.5) +
  theme_classic()

# # save legend
legend <- get_legend(ggplot() +
             geom_point(aes(x = dat_all$log_LL_t, y = dat_all$log_LL_tplus1, col = dat_all$Location)) +
             ylab(c("log(size_t+1)")) +
             xlab(c("log(size_t)")) +
             ggtitle("Growth")+
             scale_colour_manual(values = c("#f69f71", '#4cc88a'),
                                 guide_legend(title = "Population")) +
             guides(col = guide_legend(direction = "horizontal")) +
             geom_line(aes(x = newdata$log_LL_t, y = predict(object = sizeMod_soap, newdata =  newdata, type = "response")),
                       col = '#00b159', lwd = 1.5) +
             geom_line(aes(x = newdata$log_LL_t, y = predict(object = sizeMod_Base, newdata =  newdata, type = "response")),
                       col = "#f37735", lwd = 1.5) +
             theme_classic())

# flowering
flwrFig <- ggplot() + 
  geom_point(aes(x = dat_all$log_LL_t, y = dat_all$flowering), col = "grey40", alpha = .3) +
  geom_smooth(aes(x = log_LL_t, y = flowering), formula = (y ~ poly(x,2)), data = datSoap, method = "glm", 
              method.args = list(family = "binomial"), col = '#00b159', fill = '#00b159',alpha = .4,  fullrange = TRUE) + 
  geom_smooth(aes(x = log_LL_t, y =flowering), formula = (y ~ poly(x,2)), data = datBase, method = "glm", 
              method.args = list(family = "binomial"), col = "#f37735", fill = "#f37735", alpha = .4, fullrange = TRUE)  +
  ylab(c("Probability(flowering)")) +
  xlab(c("log(size_t)")) +
  ggtitle("Flowering ")+
  scale_colour_manual(values = c("#f69f71", '#4cc88a'), 
                      guide_legend(title = "Population")) +
  #geom_line(aes(x = newdata$log_LL_t, y = predict(object = flwrMod_soap, newdata =  newdata, type = "response")),  col = '#00b159', lwd = 1.5) +
  #geom_line(aes(x = newdata$log_LL_t, y = predict(object = flwrMod_Base, newdata =  newdata, type = "response")),  col = "#f37735", lwd = 1.5) +
  theme_classic()

# seeds
seedFig <- ggplot() +
  geom_point(aes(x = dat_all$log_LL_t, y = dat_all$Num_seeds), col = "grey40", alpha = .3) +
  geom_smooth(aes(x = log_LL_t, y = Num_seeds), 
              data =  seedDat_soap, method = MASS::glm.nb, 
              col = '#00b159', fill = '#00b159',alpha = .4,  fullrange = TRUE) + 
  geom_smooth(aes(x = log_LL_t, y = Num_seeds), 
              data = seedDat_Base, method = MASS::glm.nb, 
              col = "#f37735", fill = "#f37735", alpha = .4, fullrange = TRUE)  +
  ylab(c("Number of Seeds")) +
  xlab(c("log(size_t)")) +
  ggtitle("Seed Production")+
  #geom_line(aes(x = newdata$log_LL_t, y = predict(object = seedMod_soap, newdata =  newdata, type = "response")),  col = '#00b159', lwd = 1.5) +
  #geom_line(aes(x = newdata$log_LL_t, y = predict(object = seedMod_Base, newdata =  newdata, type = "response")),  col = "#f37735", lwd = 1.5) +
  theme_classic()

ggpubr::ggarrange(survFig, sizeFig, flwrFig, seedFig, align = "hv", 
                  legend = "bottom", legend.grob = legend,
                  labels = c("A", "B", "C", "D"))

#### figure of IPM matrix (all data, DI, all transitions) ####
# matrices are called 'mat_all_Soap' and 'mat_all_Base'
# "mat_all_DI" is from script 04
image(t(mat_all_DI)^.1)

#ipmr model 
contSeedlings_IPM$sub_kernels$continuous_to_seedbank

## visualize the full kernel 
# define the meshpoints 
meshpts <- seq(from = L, to = U, length = 500)
# set up the palette to the continuous part of the kernel (leave out first two rows and cols that correspond to the discrete stages)
## make the entire figure as a lattice plot
graphics::layout(mat = matrix(c(1,3,2,4), nrow = 2, ncol = 2), widths = c(1.5,6),heights = c(6,1.5))
## continuous to seedbank
pal <- hcl.colors(n = 100, palette = "Heat 2", rev = TRUE)
par(mar = c(3,3,3,1))
image((
  contSeedlings_IPM$sub_kernels$continuous_to_seedbank), xaxt = "n", yaxt = "n",
      #main = "B(t+1)",
      col = pal[(round(min((
        contSeedlings_IPM$sub_kernels$continuous_to_seedbank)^.2),2)*100): 
                  (round(max((
                    contSeedlings_IPM$sub_kernels$continuous_to_seedbank)^.2),2)*100)]) 
mtext(side = 3, text  = c("Seedbank"), line = 1, cex = 1, font = 2)
mtext(side = 2, c("Next Year (t+1)"), xpd = NA, cex = 1, font = 2, line = .75)

## K matrix
par(mar = c(3,3,3,3)
    ,mgp = c(1.75,.5,0)
)
image(x = meshpts, y = meshpts, t(contSeedlings_IPM$sub_kernels$P + contSeedlings_IPM$sub_kernels$F) ^.2,
      xlab = "n(z,t); log(cm)", ylab = "n(z',t+1); log(cm)", 
      main = "Continuous Stage",
      col = pal[(round(min( t(contSeedlings_IPM$sub_kernels$P + contSeedlings_IPM$sub_kernels$F)^.2),2)*100): 
                  (round(max(t(contSeedlings_IPM$sub_kernels$P + contSeedlings_IPM$sub_kernels$F)^.2),2)*100)]
) ## get the correct values for the color ramp that correspond to the actual probabilities in the entire matrix
text(x = 4.8, y = 1, c("Continuous Stage"), xpd = NA, srt = -90, cex = 1.25, font = 2)
abline(a = 0, b = 1, lty = 2, col = "grey30")
contour(x = meshpts, y = meshpts, 
        t(contSeedlings_IPM$sub_kernels$P + contSeedlings_IPM$sub_kernels$F), 
        add = TRUE, drawlabels = TRUE, nlevels = 10, col = "grey30")

## seedbank to seedbank
par(mar = c(2,3,1,1))
image(contSeedlings_IPM$sub_kernels$seedbank_to_seedbank^.2, xaxt = "n", yaxt = "n", 
      col = pal[(contSeedlings_IPM$sub_kernels$seedbank_to_seedbank^.2)*100]) 
mtext(side = 1, text  = c(0.71), line = -1, cex = .75)
## seedbank to continuous
par(mar = c(2,3,1,3))
image((
  contSeedlings_IPM$sub_kernels$seedbank_to_continuous), xaxt = "n", yaxt = "n",
  #main = "B(t+1)",
  col = pal[(round(min((
    contSeedlings_IPM$sub_kernels$seedbank_to_continuous)^.2),2)*100): 
      (round(max((
        contSeedlings_IPM$sub_kernels$seedbank_to_continuous)^.2),2)*100)]) 
#text(x = side  = 4, text = "Seedbank",las = 0, crt =180)
text(x = 1.05, y = .25, c("Seedbank"), xpd = NA, srt = -90, cex = 1.25, font = 2)
mtext(side = 1, c("Current Year (t)"), xpd = NA, cex = 1, font = 2, line = .75)
dev.off()

#### figures of elasticity and sensitivity matrices ####
# calculate elasticity and sensitivity matrices
## calculate the stable size distribution
w.eigen_test <-  right_ev(contSeedlings_IPM)
w.eigen_test <- c(w.eigen_test$b_w,w.eigen_test$size_w)
stable.dist_test <- w.eigen_test/sum(w.eigen_test) 

## calculate the reproductive value distribution
v.eigen_test <- left_ev(contSeedlings_IPM)
v.eigen_test <- c(v.eigen_test$b_v, v.eigen_test$size_v)
repro.val_test <- v.eigen_test/v.eigen_test[1]

## make elas and sens matrices
v.dot.w_test <- sum(stable.dist_test * repro.val_test)*h
# calculate the sensitivity function (whole kernel)
sens_test <- outer(repro.val_test,stable.dist_test, '*')/(v.dot.w_test)
# calculate the elasticity function (whole kernel)
elas_test <- matrix(as.vector(sens_test)*as.vector(mat_all_DI)/lambda(contSeedlings_IPM),nrow=501)

### elas figure
graphics::layout(mat = matrix(c(1,3,2,4), nrow = 2, ncol = 2), widths = c(1.5,6),heights = c(6,1.5))
## continuous to seedbank
pal <- hcl.colors(n = 100, palette = "Heat 2", rev = TRUE)
par(mar = c(3,3,3,1))
image((t(elas_test[1,2:501])^.1), xaxt = "n", yaxt = "n",
  #main = "B(t+1)",
  col = pal[(round(min(t(elas_test[1,2:501])^.1),2)*50): 
      (round(max(t(elas_test[1,2:501])^.2),2)*50)]) 
mtext(side = 3, text  = c("Seedbank"), line = 1, cex = 1, font = 2)
mtext(side = 2, c("Next Year (t+1)"), xpd = NA, cex = 1, font = 2, line = .75)

## K matrix
par(mar = c(3,3,3,3)
    ,mgp = c(1.75,.5,0)
)
image(x = meshpts, y = meshpts, t(elas_test[2:501,2:501])^.1,
      xlab = "n(z,t); log(cm)", ylab = "n(z',t+1); log(cm)", 
      main = "Continuous Stage",
      col = pal[(round(min(t(elas_test[2:501,2:501])^.1),2)*50): 
                  (round(max(t(elas_test[2:501,2:501])^.1), 2)*50)]
) ## get the correct values for the color ramp that correspond to the actual probabilities in the entire matrix
text(x = 4.8, y = 1, c("Continuous Stage"), xpd = NA, srt = -90, cex = 1.25, font = 2)
abline(a = 0, b = 1, lty = 2, col = "grey30")
contour(x = meshpts, y = meshpts, 
        t(elas_test[2:501,2:501]), 
        add = TRUE, drawlabels = TRUE, nlevels = 10, col = "grey30")

## seedbank to seedbank
par(mar = c(2,3,1,1))
image(as.matrix(elas_test[1,1]^.1), xaxt = "n", yaxt = "n", 
      col = pal[(elas_test[1,1]^.1)*50]) 
## seedbank to continuous
par(mar = c(2,3,1,3))
image(
  as.matrix(elas_test[2:501,1]^.1), xaxt = "n", yaxt = "n",
  #main = "B(t+1)",
  col = pal[(round(min(
    (elas_test[2:501,1]^.1)),2)*50): 
      (round(max((elas_test[2:501,1]^.1)),2)*50)]) 
#text(x = side  = 4, text = "Seedbank",las = 0, crt =180)
text(x = 1.05, y = .25, c("Seedbank"), xpd = NA, srt = -90, cex = 1.25, font = 2)
mtext(side = 1, c("Current Year (t)"), xpd = NA, cex = 1, font = 2, line = .75)
dev.off()

### sensitivity figure
graphics::layout(mat = matrix(c(1,3,2,4), nrow = 2, ncol = 2), widths = c(1.5,6),heights = c(6,1.5))
## continuous to seedbank
pal <- hcl.colors(n = 100, palette = "Heat 2", rev = TRUE)
par(mar = c(3,3,3,1))
image((t(sens_test[1,2:501])^.1), xaxt = "n", yaxt = "n",
      #main = "B(t+1)",
      col = pal[(round(min(t(sens_test[1,2:501])^.1),2)*25): 
                  (round(max(t(sens_test[1,2:501])^.2),2)*25)]) 
mtext(side = 3, text  = c("Seedbank"), line = 1, cex = 1, font = 2)
mtext(side = 2, c("Next Year (t+1)"), xpd = NA, cex = 1, font = 2, line = .75)

## K matrix
par(mar = c(3,3,3,3)
    ,mgp = c(1.75,.5,0)
)
image(x = meshpts, y = meshpts, t(sens_test[2:501,2:501])^.1,
      xlab = "n(z,t); log(cm)", ylab = "n(z',t+1); log(cm)", 
      main = "Continuous Stage",
      col = pal[(round(min(t(sens_test[2:501,2:501])^.1),2)*25): 
                  (round(max(t(sens_test[2:501,2:501])^.1), 2)*25)]
) ## get the correct values for the color ramp that correspond to the actual probabilities in the entire matrix
text(x = 4.8, y = 1, c("Continuous Stage"), xpd = NA, srt = -90, cex = 1.25, font = 2)
abline(a = 0, b = 1, lty = 2, col = "grey30")
contour(x = meshpts, y = meshpts, 
        t(sens_test[2:501,2:501]), 
        add = TRUE, drawlabels = TRUE, nlevels = 10, col = "grey30")

## seedbank to seedbank
par(mar = c(2,3,1,1))
image(as.matrix(sens_test[1,1]^.1), xaxt = "n", yaxt = "n", 
      col = pal[(sens_test[1,1]^.1)*25]) 
## seedbank to continuous
par(mar = c(2,3,1,3))
image(
  as.matrix(sens_test[2:501,1]^.1), xaxt = "n", yaxt = "n",
  #main = "B(t+1)",
  col = pal[(round(min(
    (sens_test[2:501,1]^.1)),2)*25): 
      (round(max((sens_test[2:501,1]^.1)),2)*25)]) 
#text(x = side  = 4, text = "Seedbank",las = 0, crt =180)
text(x = 1.05, y = .25, c("Seedbank"), xpd = NA, srt = -90, cex = 1.25, font = 2)
mtext(side = 1, c("Current Year (t)"), xpd = NA, cex = 1, font = 2, line = .75)

dev.off()


# density dependence figs -------------------------------------------------
siteDI_bootLambdas <- readRDS("~/Dropbox/Grad School/Research/Oenothera coloradensis project/COBP_analysis/intermediate_analysis_Data/site_level_IPMs_allYears/site_level_DI_bootCI_lambdas.RDS")
# get mean lambda for each subpop
names(siteDI_bootLambdas) <- unique(dat_all$Site)
# calculate CIs 
SE <- sapply(siteDI_bootLambdas, function(x) sd(x)/sqrt(1000)) 
mean <- sapply(siteDI_bootLambdas, mean)
allDI_CI <- 
  c((mean - 1.96*SE),(mean + 1.96*SE))

#log(sapply(siteDI_mats, function(x) eigen(x)$values[1]))

siteDD_bootLambdas <- readRDS("~/Dropbox/Grad School/Research/Oenothera coloradensis project/COBP_analysis/intermediate_analysis_Data/site_level_IPMs_allYears/site_level_DD_bootCI_lambdas.RDS")
names(siteDD_bootLambdas) <- unique(dat_all$Site)
# calculate CIs 
SE <- sapply(siteDD_bootLambdas, function(x) sd(x)/sqrt(1000)) 
mean <- sapply(siteDD_bootLambdas, mean)
allDI_CI <- 
  c((mean - 1.96*SE),(mean + 1.96*SE))
## something is very wierd about the CI's for DD Models... redo? double check? 

# log(sapply(siteDD_mats, function(x) eigen(x)$values[1]))


# figure of location of plots + species range ---------------------------------------------
# load data 
setwd("/Users/Alice/Dropbox/Grad School/Research/Oenothera coloradensis project/Raw Data")
sites <- read.csv("./COBP Plot Locations.csv", stringsAsFactors = FALSE)

#order the sites in alphabetical order by site name 
sites$site_2 <- str_sub(sites$Site, start = str_locate(sites$Site, " ")[,1]+1, end = str_length(sites$Site)) #make a new variable that shows just the site name (without the location)
sites <- sites[order(sites$site_2),]
sites_NoSeed <- sites %>% 
  filter(site_2 != "Seed")

sites_NoSeed <- st_as_sf(sites_NoSeed, coords = c("Long_WGS84","Lat_WGS84"))
# get centroid for each subpopulation
sites_NoSeed$population <- str_sub(sites_NoSeed$Site , start = 1, end = str_locate(sites_NoSeed$Site, pattern = " ")[,1]-1)
sites_NoSeed <- st_as_sf(sites_NoSeed)
st_crs(sites_NoSeed) <- 4326

for (i in unique(sites_NoSeed$site_2)) {
  if (i == "Crow Creek") {
    sites_centroids <- sites_NoSeed[sites_NoSeed$site_2 == i,][1,]
    sites_centroids$geometry <- 
      st_centroid(st_union(sites_NoSeed[sites_NoSeed$site_2 == i,]))
    } else {
      temp <- sites_NoSeed[sites_NoSeed$site_2 == i,][1,]
      temp$geometry <- 
        st_centroid(st_union(sites_NoSeed[sites_NoSeed$site_2 == i,]))
      sites_centroids <- rbind(sites_centroids, 
          temp)
  }
}

for (i in unique(sites_NoSeed$population)) {
  if (i == "FEWAFB") {
    sites_pops <- sites_NoSeed[sites_NoSeed$population == i,][1,]
    sites_pops$geometry <- 
      st_centroid(st_union(sites_NoSeed[sites_NoSeed$population == i,]))
  } else {
    temp <- sites_NoSeed[sites_NoSeed$population == i,][1,]
    temp$geometry <- 
      st_centroid(st_union(sites_NoSeed[sites_NoSeed$population == i,]))
    sites_pops <- rbind(sites_pops, 
                             temp)
  }
}
st_crs(sites_pops) <- 4326

## make a map of the distribution (present and historical) using mapedit package
# distMap <- mapview(sites_centroids)  %>% 
#   editMap()
# distMap_tempGeom <- distMap$finished
# distMap_tempTemp <- mapview(distMap_tempGeom) %>% 
#   editMap()
# distMap_final <- rbind(distMap_tempTemp$all, distMap_tempGeom)
# distMap_final$distType <- c("historical", "current", "historical", "current", "historical")
# # distMap_final <- st_make_valid(distMap_final)
# saveRDS(distMap_final, file = "../Raw Data/SpatialDistributionMap")

# define the palette for each subpopulation (plot)
pal <- colorFactor(
  palette = "Dark2",
  domain = sites$site_2)
# define the palette for each subpopulation 
pal_subPop <- colorFactor(
  palette = "Dark2",
  domain = sites_centroids$site_2)
pal_dist <- colorFactor(
  palette = "Spectral", 
  domain = distMap_final$distType
)

## map of species distributions and population locations
# leaflet()  %>% 
#   addProviderTiles(providers$Esri.WorldImagery) %>% 
#   addPolygons(data = distMap_final[distMap_final$distType == "historical",], color = "black", weight = 1, opacity = 1, fillOpacity = .5, fillColor = "#95b1c0") %>%
#   addPolygons(data = distMap_final[distMap_final$distType == "current",], color = "black", weight = 1, opacity = 1, fillOpacity = .7, fillColor = "#2B6482") %>%
#   #addCircles(data = sites_centroids, col = "black", radius = 600, opacity = .7, weight = 1, fillColor = ~pal_subPop(sites_centroids$site_2), fillOpacity = .7) %>%
#   addCircles(data = sites_pops, col = "black", radius = 2000, opacity = 1, weight = 1, fillColor = "black", fillOpacity = 1) %>% 
#   addLabelOnlyMarkers(lat = st_coordinates(sites_pops[1,])[2], 
#                       lng = st_coordinates(sites_pops[1,])[1],
#                       label = "FEWAFB", labelOptions = list(permanent = TRUE)) %>%
#   addLabelOnlyMarkers(lat = st_coordinates(sites_pops[2,])[2], 
#                       lng = st_coordinates(sites_pops[2,])[1],
#                       label = "Soapstone Prairie", labelOptions = list(permanent = TRUE, direction = "bottom")) %>% 
#   addLegend(position = "bottomright", colors = c("#2B6482", "#95b1c0"), labels = c("Current", "Historical"), title = "Species Distribution") %>% 
#   addScaleBar(position = "bottomright")

## ggmap version
dist_bbox <- st_bbox(st_buffer(st_union(distMap_final),10000))
names(dist_bbox) <- c("left", "bottom", "right", "top")
terrain_basemap <- get_map(
  location=dist_bbox, 
  zoom = 10,
  maptype = 'terrain-background', 
  source = 'stamen')

ggmap(terrain_basemap) + 
  geom_sf(data = distMap_final, 
          aes(fill = distType), 
          alpha = .8,
          inherit.aes = FALSE) +
  scale_fill_manual(values =  c("#2B6482", "#95b1c0"), 
                     name = "Species \n Distribution") + 
  geom_sf(data = sites_pops, aes(), inherit.aes = FALSE) + 
  geom_label_repel(inherit.aes = FALSE, 
                   data = sites_pops, 
                   aes(x = st_coordinates(sites_pops)[,1],
                                                              y = st_coordinates(sites_pops)[,2],
                                                              label = population)) + 
  xlab("") + 
  ylab("") + 
  annotation_scale(data = sites_NoSeed[sites_NoSeed$population == "FEWAFB",])


## map of plot locations (FEWAFB)
FEWAFB_bbox <- st_bbox(st_buffer(sites_NoSeed[sites_NoSeed$population == "FEWAFB",],1000))
names(FEWAFB_bbox) <- c("left", "bottom", "right", "top")
terrain_basemap <- get_map(
  location=FEWAFB_bbox, 
  zoom = 14,
  maptype = 'terrain-background', 
  source = 'stamen')

ggmap(terrain_basemap) + 
  geom_sf(data = sites_NoSeed[sites_NoSeed$population == "FEWAFB",], 
          aes(col = site_2), 
          inherit.aes = FALSE) +
  scale_color_manual(values = pal(unique(sites_NoSeed$site_2))[c(1,2,6)], 
                     name = "Subpopulation") + 
  xlim(c(-104.885, -104.869)) + 
  ylim(c(41.135, 41.156)) + 
  xlab("") + 
  ylab("") + 
  annotation_scale(data = sites_NoSeed[sites_NoSeed$population == "FEWAFB",]) + 
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank())

## map of plot locations (Soapstone)
Soapstone_bbox <- st_bbox(st_buffer(sites_NoSeed[sites_NoSeed$population == "Soapstone",],500))
names(Soapstone_bbox) <- c("left", "bottom", "right", "top")
terrain_basemap <- get_map(
  location=Soapstone_bbox, 
  zoom = 12,
  maptype = 'terrain-background', 
  source = 'stamen')

ggmap(terrain_basemap) + 
  geom_sf(data = sites_NoSeed[sites_NoSeed$population == "Soapstone",], 
          aes(col = site_2), 
          inherit.aes = FALSE) +
  scale_color_manual(values = pal(unique(sites_NoSeed$site_2))[c(3,4,5)], 
                     name = "Subpopulation") + 
  xlim(c(-105.023, -105.008)) + 
  ylim(c(40.985, 40.995)) + 
  xlab("") + 
  ylab("") + 
  annotation_scale(data = sites_NoSeed[sites_NoSeed$population == "FEWAFB",]) + 
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank())


# sample size calculations for Table S1 -------------------------------------
adult_seed_n <- dat_all %>% 
  group_by(Site, seedling, Year) %>% 
  summarize(n = sum(survives_t))

all_n <- dat_all %>% 
  group_by(Site,  Year) %>% 
  summarize(n = n())


seedlings_cont %>% 
  group_by(Site, Year) %>% 
  summarize(n = sum(seedling))
seedlings %>% 
  group_by(Site, Year) %>% 
  summarize(n = sum(Seedlings_t))
