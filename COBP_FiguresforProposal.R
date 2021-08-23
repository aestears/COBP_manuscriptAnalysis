#////////////////
# Oenothera coloradensis initial data analysis
# 24 May 2019
# Alice Stears
#////////////////

require(tidyverse)
require(sf)
require(leaflet)
require(wesanderson)
require(RColorBrewer)

# load dataset
setwd("/Users/Alice/Dropbox/Grad School/Research/Oenothera coloradensis project/Raw Data")
counts<- read.csv("./COBP_2018_data_7_14_20__all_REAL.csv", stringsAsFactors = FALSE) #will have to update file name as it changes w/ most current version
sites <- read.csv("./COBP Plot Locations.csv", stringsAsFactors = FALSE)

#make a column for year
counts$Date_2018 <- as.POSIXct(counts$Date_2018, tz = "UTC", format = "%m/%d/%Y")
#fix year
lubridate::year(counts$Date_2018) <- ifelse(is.na(lubridate::year(counts$Date_2018))==TRUE, NA, lubridate::year(counts$Date_2018) + 2000)
counts$Date_2019 <- as.POSIXct(counts$Date_2019, tz = "UTC", format = "%m/%d/%Y")
lubridate::year(counts$Date_2019) <- ifelse(is.na(lubridate::year(counts$Date_2019))==TRUE, NA, lubridate::year(counts$Date_2019) + 2000)

#make sure survival is coded as a factor
#counts$Alive_2018 <- as.factor(counts$Alive_2018)
#counts$Alive_2019 <- as.factor(counts$Alive_2019)

#### make a map of plot locations ####
#order the sites in alphabetical order by site name 
sites$site_2 <- str_sub(sites$Site, start = str_locate(sites$Site, " ")[,1]+1, end = str_length(sites$Site)) #make a new variable that shows just the site name (without the location)
sites <- sites[order(sites$site_2),]
sites_NoSeed <- sites %>% 
  filter(site_2 != "Seed")
# define the palette
pal <- colorFactor(
palette = "Dark2",
domain = sites$site_2)
sites_NoSeed <- st_as_sf(sites_NoSeed, coords = c("Long_WGS84","Lat_WGS84"))
st_crs(sites_NoSeed) <- 4326
leaflet(data = sites_NoSeed)  %>% 
  addProviderTiles(providers$OpenStreetMap.HOT) %>% 
  addCircles(col = "black", radius = 400, opacity = .7, weight = 1, fillColor = ~pal(sites_NoSeed$site_2), fillOpacity = .7) %>% 
  addLegend(position = "bottomright", colors = unique(pal(sites_NoSeed$site_2)), labels = unique(sites_NoSeed$site_2), title = "Sub-population")

#then ungeometrize the sites data.frame 
sites <- read.csv("./COBP Plot Locations.csv", stringsAsFactors = FALSE)
sites$site_2 <- str_sub(sites$Site, start = str_locate(sites$Site, " ")[,1]+1, end = str_length(sites$Site))

#### make a plot of plot density by plant type #### 
# for 2018 data
dens_2018 <- aggregate(counts$ID, by = list(counts$Plot_ID, counts$Bolting_2018), FUN = length)
names(dens_2018) <- c("Plot_ID", "Flowering", "Density")
dens_2018[dens_2018$Flowering==0,"Flowering"] <- "vegetative_18"
dens_2018[dens_2018$Flowering==1,"Flowering"] <- "reproductive_18"
names(dens_2018) <- c("Plot_ID", "DataType", "Density")
dens_2018$Year <- as.integer(2018) 
dens_2018 <- dens_2018[order(dens_2018$Plot_ID),]

# for 2019 data
dens_2019 <- aggregate(counts$ID, by = list(counts$Plot_ID, counts$Bolting_2019), FUN = length)
names(dens_2019) <- c("Plot_ID", "Flowering", "Density")
dens_2019[dens_2019$Flowering==0,"Flowering"] <- "vegetative_19"
dens_2019[dens_2019$Flowering==1,"Flowering"] <- "reproductive_19"
names(dens_2019) <- c("Plot_ID", "DataType", "Density")
dens_2019$Year <- as.integer(2019)
dens_2019 <- dens_2019[order(dens_2019$Plot_ID),]

#combine
dens <- rbind(dens_2018, dens_2019)

#read in seedling data
seeds <- read.csv("./COBP_seedlings_6_11_20.csv", stringsAsFactors = FALSE)
seedsTot <- aggregate(seeds[,c("Seedlings_18", "Seedlings_19")], by = list(seeds$Plot_ID), FUN = sum)
names(seedsTot) <- c("Plot_ID", "Density")
seedsTot$DataType <- "seedlings"
names(seedsTot) <- c("Plot_ID", "Seedlings_18", "Seedlings_19", "DataType")
#reshape seedling data frame
seedsTot <- seedsTot %>%  
  pivot_longer(c("Seedlings_18", "Seedlings_19"), names_to = "Type", values_to = "Density")
 
names(seedsTot) <- c("Plot_ID", "Type", "DataType", "Density")
seedsTot <- seedsTot[,c("Plot_ID","DataType", "Density")]
seedsTot$Year <- ifelse(str_detect(seedsTot$DataType, "_18"), 2018, 2019)

#join seedling and other plant data together
densTot <- rbind(dens, seedsTot)

#calculate total density (adult plants + seedlings)
densAll <- aggregate(dens[,c("Density")], by = list(dens$Plot_ID, dens$Year), FUN = sum)
names(densAll) <- c("Plot_ID", "Year", "Density")
densAll$DataType <- NA
densAll[densAll$Year==2018,"DataType"] <- "total_18"
densAll[densAll$Year==2019,"DataType"] <- "total_19"
densAll <- densAll[,c("Plot_ID", "DataType", "Density", "Year")]
densTot <- rbind(densTot, densAll)
densTot <- densTot[order(densTot$Plot_ID),]

#make a plot
densTot <- left_join(densTot, data.frame(Plot.Name = sites[,c("Plot.Name")]), by = c("Plot_ID" = "Plot.Name"))
densTot$SubPop <- NA
densTot[densTot$Plot_ID %in% c("C4","C5","C8"),"SubPop"] <- "Crow Creek"
densTot[densTot$Plot_ID %in% c("U3","U4","U6"),"SubPop"] <- "Unnamed"
densTot[densTot$Plot_ID %in% c("D10","D11","D7"),"SubPop"] <- "Diamond Creek"
densTot[densTot$Plot_ID %in% c("S1","S2","S3"),"SubPop"] <- "HQ5"
densTot[densTot$Plot_ID %in% c("S4","S5","S6"),"SubPop"] <- "HQ3"
densTot[densTot$Plot_ID %in% c("S7","S8","S9"),"SubPop"] <- "Meadow"

dat <- densTot[densTot$DataType %in% c("vegetative_18", "reproductive_18", "vegetative_19", "reproductive_19", "Seedlings_18", "Seedlings_19"),]
#add sub-population ID to data frame
dat <- left_join(dat, unique(counts[,c("Location","Site","Plot_ID")]), by = "Plot_ID")
#reorder dat dataframe
dat$Plot_ID <- factor(dat$Plot_ID, levels = c("C4", "C5", "C8", "D10", "D11", "D7", "U3", "U4", "U6", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9"))
dat <- dat[order(dat$Location),]

#get site-level sub-population size for 2018
SumTemp <- aggregate(dat[,c("Density")], by = list(dat$Site, dat$Year), FUN = sum) 
names(SumTemp) <- c("Site", "Year","Density")
SumTemp$Location <- NA
SumTemp[SumTemp$Site %in% c("Crow Creek", "Diamond Creek", "Unnamed Creek"), "Location"] <- "FEWAFB"
SumTemp[SumTemp$Site %in% c("HQ3", "HQ5", "Meadow"), "Location"] <- "Soapstone"

means<- aggregate(SumTemp, by = list(SumTemp$Location, SumTemp$Year), FUN = mean)
names(means) <- c("Location", "YearBad","SiteBad", "Year","MeanDensity", "LocationBad")

dat$SiteYear <- str_c(dat$Site, "_", dat$Year)

#make plot
palette <- wes_palette("Darjeeling1", type = "continuous")
palette

# no a very good plot ggplot(dat) +
  # geom_bar(aes(x = SiteYear, y = Density, fill = DataType),  width = .8, stat = "identity") + 
  # labs(x = "Sub Population", y = "Number of Individuals" , fill = "Life Stage") + 
  # scale_fill_manual(values = c("tomato", "tomato4", "turquoise", "turquoise4", "goldenrod", "goldenrod4")) + 
  # facet_wrap("Location", scale = "free") + 
  # theme_classic() +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  # geom_hline(data = means[means$Year==2018,], aes(yintercept=MeanDensity),
  #            linetype="dashed", color = "grey30") + 
  # geom_hline(data = means[means$Year==2019,], aes(yintercept=MeanDensity),
  #            linetype="dashed", color = "grey60") 

#### Make size-class distributions by plot ####
size <- select(counts, c("Location", "Site", "Plot_ID", "Quadrant", "ID", "LongestLeaf_cm_2018", "LongestLeaf_cm_2019", "Bolting_2018", "Bolting_2019", "No_Capsules_2018"))
mu_2018 <- plyr::ddply(size, "Location", summarise, grp.mean=mean(LongestLeaf_cm_2018, na.rm = TRUE), 
                  grp.sd = sd(LongestLeaf_cm_2018, na.rm = TRUE))
mu_2019 <- plyr::ddply(size, "Location", summarise, grp.mean=mean(LongestLeaf_cm_2019, na.rm = TRUE), 
                       grp.sd = sd(LongestLeaf_cm_2019, na.rm = TRUE))
# Less Good Plots (by year)
# #plot for 2018 values
# ggplot(counts, aes(x=LongestLeaf_cm_2018, color = Site)) + 
#   theme_classic()  +
#   labs(x = "Longest Leaf Length (cm)", y = "Density") +
#   geom_vline(data = mu_2018, aes(xintercept=grp.mean),
#              linetype="dashed", color = "grey40") + 
#   geom_density() +
#   facet_wrap(c("Location")) +
#   ylim(c(0,.20)) + 
#   xlim (c(1,40)) + 
#   ggtitle("2018 Size Class Distribution") + 
#   scale_colour_brewer(palette = "Dark2")
# 
# #plot for 2019 values
# ggplot(counts, aes(x=LongestLeaf_cm_2019, color = Site)) + 
#   theme_classic() +
#   theme(panel.grid.major = element_line(colour = "grey90")) +
#   labs(x = "Longest Leaf Length (cm)", y = "Density") +
#   geom_vline(data = mu_2019, aes(xintercept=grp.mean),
#              linetype="dashed", color = "grey40") + 
#   geom_density() +
#   facet_wrap(c("Location")) +
#   ylim(c(0,.20)) + 
#   xlim (c(1,40)) + 
#   ggtitle("2019 Size Class Distribution") + 
#   scale_colour_brewer(palette = "Dark2")

#plot by site
ggplot(counts, aes(color = Site)) + 
  theme_classic() +
  geom_density(aes(x=LongestLeaf_cm_2019), linetype = "dotted") +
  geom_density(aes(x=LongestLeaf_cm_2018)) +
  labs(x = "Longest Leaf Length (cm)", y = "Density") +
  facet_wrap(c("Site")) +
  ylim(c(0,.20)) + 
  xlim (c(1,40)) + 
  scale_colour_brewer(palette = "Dark2") + 
  guides(color = FALSE)

#make a better figure 
#define a function to add an alpha to a color
add.alpha <- function(col, alpha=1){
if(missing(col))
  stop("Please provide a vector of colours.")
apply(sapply(col, col2rgb)/255, 2, 
      function(x) 
        rgb(x[1], x[2], x[3], alpha=alpha))  
}
#define color palette
pal <- brewer.pal(6,"Dark2")

par(mfrow = c(2,3), mar = c(2.2,2.5,1.5,1), cex = 1)
plot(x = NA, y = NA,
     ylim = c(0,.2), 
     xlim = c(0,45),
     xlab = c(""),
     ylab = c(""),
     main = c("Crow Creek"))
abline(v = c(mean(counts[counts$Site=="Crow Creek","LongestLeaf_cm_2018"], na.rm = TRUE),mean(counts[counts$Site=="Crow Creek","LongestLeaf_cm_2019"], na.rm = TRUE)), lty = c(1,2), lwd = 2, col = "gray40")
#lines(density(counts[counts$Site=="Crow Creek","LongestLeaf_cm_2019"], na.rm = TRUE), lty = 2, lwd = 2, col = pal[1])
polygon(density(counts[counts$Site=="Crow Creek","LongestLeaf_cm_2018"], na.rm = TRUE), col = add.alpha(pal[1],.5), border = pal[1])
polygon(density(counts[counts$Site=="Crow Creek","LongestLeaf_cm_2019"], na.rm = TRUE), col = add.alpha(pal[1],.5), border = pal[1], lty = 2)
polygon( x = c(23,23,46.5,46.5), y = c(.07, .16, .16, .07), border = "transparent", col = "gray90")
legend(15,.19, legend = c("2018", "2019"), col = "gray40", lty = c(1,2), lwd = 2, text.width = 1, seg.len = 1, cex = .9, box.lty = 0, x.intersp = .2, bg = NA)


plot(x = NA, y = NA,
     ylim = c(0,.2),
     xlim = c(0,45), 
     xlab = c(""),
     ylab = c(""),
     main = c("Diamond Creek"))
abline(v = c(mean(counts[counts$Site=="Diamond Creek","LongestLeaf_cm_2018"], na.rm = TRUE),mean(counts[counts$Site=="Diamond Creek","LongestLeaf_cm_2019"], na.rm = TRUE)), lty = c(1,2), lwd = 2, col = "gray40")
#lines(density(counts[counts$Site=="Diamond Creek","LongestLeaf_cm_2019"], na.rm = TRUE), lty = 2, lwd = 2, col = "#D95F02")
polygon(density(counts[counts$Site=="Diamond Creek","LongestLeaf_cm_2018"], na.rm = TRUE), col = add.alpha(pal[2],.5), border = pal[2])
polygon(density(counts[counts$Site=="Diamond Creek","LongestLeaf_cm_2019"], na.rm = TRUE), col = add.alpha(pal[2],.5), border = pal[2], lty = 2)

plot(x = NA, y = NA,
     ylim = c(0,.2),
     xlab = c(""),
     xlim = c(0,45), 
     ylab = c(""),
     main = c("HQ3"))
abline(v = c(mean(counts[counts$Site=="HQ3","LongestLeaf_cm_2018"], na.rm = TRUE),mean(counts[counts$Site=="HQ3","LongestLeaf_cm_2019"], na.rm = TRUE)), lty = c(1,2), lwd = 2, col = "gray40")
# lines(density(counts[counts$Site=="HQ3","LongestLeaf_cm_2019"], na.rm = TRUE), lty = 2, lwd = 2, col = "#7570B3")
polygon(density(counts[counts$Site=="HQ3","LongestLeaf_cm_2018"], na.rm = TRUE), col = add.alpha(pal[3],.5), border = pal[3])
polygon(density(counts[counts$Site=="HQ3","LongestLeaf_cm_2019"], na.rm = TRUE), col = add.alpha(pal[3],.5), border = pal[3], lty = 2)


plot(x = NA, y = NA,
     ylim = c(0,.2),
     xlab = c(""),
     ylab = c(""),
     xlim = c(0,45), 
     main = c("HQ5"))
abline(v = c(mean(counts[counts$Site=="HQ5","LongestLeaf_cm_2018"], na.rm = TRUE),mean(counts[counts$Site=="HQ5","LongestLeaf_cm_2019"], na.rm = TRUE)), lty = c(1,2), lwd = 2, col = "gray40")
#lines(density(counts[counts$Site=="HQ5","LongestLeaf_cm_2019"], na.rm = TRUE), lty = 2, lwd = 2, col = "#E7298A")
polygon(density(counts[counts$Site=="HQ5","LongestLeaf_cm_2018"], na.rm = TRUE), col = add.alpha(pal[4],.5), border = pal[4])
polygon(density(counts[counts$Site=="HQ5","LongestLeaf_cm_2019"], na.rm = TRUE), col = add.alpha(pal[4],.5), border = pal[4], lty = 2)


plot(x = NA, y = NA,
     ylim = c(0,.2),
     xlab = c(""),
     ylab = c(""),
     xlim = c(0,45), 
     main = c("Meadow"))
abline(v = c(mean(counts[counts$Site=="Meadow","LongestLeaf_cm_2018"], na.rm = TRUE),mean(counts[counts$Site=="Meadow","LongestLeaf_cm_2019"], na.rm = TRUE)), lty = c(1,2), lwd = 2, col = "gray40")
#lines(density(counts[counts$Site=="Meadow","LongestLeaf_cm_2019"], na.rm = TRUE), lty = 2, lwd = 2, col = "#66A61E")
polygon(density(counts[counts$Site=="Meadow","LongestLeaf_cm_2018"], na.rm = TRUE), col = add.alpha(pal[5],.5), border = pal[5])
polygon(density(counts[counts$Site=="Meadow","LongestLeaf_cm_2019"], na.rm = TRUE), col = add.alpha(pal[5],.5), border = pal[5], lty = 2)
mtext("Plant Size (cm)", side = 1, line = 2.25, font = 2)


plot(x = NA, y = NA,
     ylim = c(0,.2),
     xlab = c(""),
     ylab = c(""),
     xlim = c(0,45), 
     main = c("Unnamed Creek"))
abline(v = c(mean(counts[counts$Site=="Unnamed Creek","LongestLeaf_cm_2018"], na.rm = TRUE),mean(counts[counts$Site=="Unnamed Creek","LongestLeaf_cm_2019"], na.rm = TRUE)), lty = c(1,2), lwd = 2, col = "gray40")
# lines(density(counts[counts$Site=="Unnamed Creek","LongestLeaf_cm_2019"], na.rm = TRUE), lty = 2, lwd = 2, col = "#E6AB02")
polygon(density(counts[counts$Site=="Unnamed Creek","LongestLeaf_cm_2018"], na.rm = TRUE), col = add.alpha(pal[6],.5), border = pal[6])
polygon(density(counts[counts$Site=="Unnamed Creek","LongestLeaf_cm_2019"], na.rm = TRUE), col = add.alpha(pal[6],.5), border = pal[6], lty = 2)


#### Make reproductive output-by-size plot ####
##2018
flowering <-counts[counts$Flowering_2018==1,]
flowering <- flowering[is.na(flowering$Location)==FALSE,]

ggplot(flowering, aes(x =LongestLeaf_cm_2018, y = No_Capsules_2018, color = Site)) +
  theme_classic() +
  theme(panel.grid.major = element_line(colour = "grey90")) +
  geom_point(aes(alpha = .8), show.legend = FALSE) +
  labs( x = "Longest Leaf Length (cm)", y = "Number of Seed Capsules") +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap("Location") + 
  scale_colour_brewer(palette = "Dark2")

ggplot(flowering, aes(x =LongestLeaf_cm_2018, y = No_Capsules_2018, color = Location)) +
  theme_classic()  +
  geom_point(aes(alpha = .8), show.legend = FALSE) +
  labs( x = "Longest Leaf Length (cm)", y = "Number of Seed Capsules") +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap("Location") + 
  scale_colour_brewer(palette = "Dark2")

dat <- flowering[flowering$Site=="Crow Creek",]
crowCreek<-lm(dat$No_Capsules_2018~dat$LongestLeaf_cm_2018, data=dat)
ccR <- 0.1848

dat <- flowering[flowering$Site=="Diamond Creek",]
diamondCreek <- lm(dat$No_Capsules_2018~dat$LongestLeaf_cm_2018, data=dat)
ddR <- 0.02329

dat <- flowering[flowering$Site=="Unnamed Creek",]
unnamedCreek <- lm(dat$No_Capsules_2018~dat$LongestLeaf_cm_2018, data=dat)
uuR <- 0.1478
  
dat <- flowering[flowering$Site=="HQ3",]
HQ3 <- lm(dat$No_Capsules_2018~dat$LongestLeaf_cm_2018, data=dat)
hq3R <- 0.2308

dat <- flowering[flowering$Site=="HQ5",]
HQ5 <- lm(dat$No_Capsules_2018~dat$LongestLeaf_cm_2018, data=dat)
hq5 <- 0.08747

dat <- flowering[flowering$Site=="Meadow",]
Meadow <- lm(dat$No_Capsules_2018~dat$LongestLeaf_cm_2018, data=dat)
meadowR <- 0.1801

#2018 and 2019
flowering_2_A <-counts %>% #get 2018 leaf size and reproductive output data
  filter(counts$Flowering_2018==1) %>% 
  select(Location,Plot_ID,ID, LongestLeaf_cm = LongestLeaf_cm_2018, No_Capsules = No_Capsules_2018) 
flowering_2_A$Year <- 2018

flowering_2_B <-counts %>% #get 2019 leaf size and reproductive output data
  filter(counts$Flowering_2019==1) %>% 
  select(Location,Plot_ID, ID, LongestLeaf_cm = LongestLeaf_cm_2019,No_Capsules = No_Capsules_2019)
flowering_2_B$Year <- 2019

#join into one data.frame
flowering_2 <- rbind(flowering_2_A, flowering_2_B)
flowering_2$Location_Year <- paste(flowering_2$Location,flowering_2$Year,sep = "_")

ggplot(flowering_2, aes(x =LongestLeaf_cm, y = No_Capsules, color = c(Location_Year))) +
  theme_classic()  +
  geom_point(aes(alpha = .8), show.legend = FALSE) +
  labs( x = "Longest Leaf Length (cm)", y = "Number of Seed Capsules") +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap("Location") + 
  scale_colour_brewer(palette = "Dark2")




#### reproductive output by size in previous year #### 
flowering_3_A <- counts %>% 
  filter(Flowering_2019==1) %>% 
  select(Location, Plot_ID, ID, LongestLeaf = LongestLeaf_cm_2018, No_Capsules_2019)
flowering_3_A$LeafYear <- 2018

flowering_3_B <- counts %>% 
  filter(Flowering_2019==1) %>% 
  select(Location, Plot_ID, ID, LongestLeaf = LongestLeaf_cm_2019, No_Capsules_2019)
flowering_3_B$LeafYear <- 2019

flowering_3 <- rbind(flowering_3_A, flowering_3_B)
  
ggplot(data = flowering_3, aes(x = LongestLeaf, y = No_Capsules_2019, col = Location)) +
  geom_point() + 
  facet_wrap(~LeafYear) + 
  geom_smooth(method = lm)

ggplot(data = counts, aes(x = LongestLeaf_cm_2018, y = No_Capsules_2019, col = Location)) + 
  geom_point() + 
  geom_smooth(method = lm)


dev.off()
par(mfrow = c(2,1))
plot(counts$No_Capsules_2019 ~ counts$LongestLeaf_cm_2018,
     xlab = "Longest Leaf, 2018",
     ylab = "No. Capsules, 2019",
     col = as.factor(counts$Location),
     pch = 16,
     xlim = c(0,40))
plot(counts$No_Capsules_2019 ~ counts$LongestLeaf_cm_2019,
     xlab = "Longest Leaf, 2019",
     ylab = "No. Capsules, 2019",
     col = as.factor(counts$Location),
     pch = 16,
     xlim = c(0,40))

#### identify zombie plants####
#identify those plants that were reproductive in 2018 and also alive in 2019!
zombie <- counts %>% 
  filter(Flowering_2018==1, Alive_2019==1)
#19 individuals that survived! (out of 317)--5.99% chance of surviving after flowering
sum(counts$Flowering_2018==1, na.rm = TRUE)


#### plot size in current year and probability of flowering ####
# subset for plants that are alive in both 2018 and 2019
temp2 <- counts %>% filter(counts$Alive_2018==1)
temp2$col <- ifelse(temp2$Location=="FEWAFB", "#f1511b", "#245f88")
temp2$log_LongestLeaf_2018 <- log(temp2$LongestLeaf_cm_2018)
plot(temp2$Bolting_2018 ~ temp2$log_LongestLeaf_2018, col = temp2$col, pch = 19,
     xlab = "2018 Leaf Length (cm)",
     ylab = "Bolting in 2018 (0-1)")

#make a model for FEWAFB
modAFB_2 <- glm(Bolting_2018 ~ log_LongestLeaf_2018, data = temp2[temp2$Location=="FEWAFB",], family = binomial)
summary(modAFB_2)
#predict new values
preddata <- with(temp2, data.frame(log_LongestLeaf_2018 = seq(min(temp2$log_LongestLeaf_2018), max(temp2$log_LongestLeaf_2018), length = 100)))
preds<- predict(modAFB_2, newdata = preddata, type = "link", se.fit = TRUE)
critval <- 1.96 ## approx 95% CI
upr<- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit

fit2 <- modAFB_2$family$linkinv(fit)
upr2 <- modAFB_2$family$linkinv(upr)
lwr2 <- modAFB_2$family$linkinv(lwr)

preddata$lwr <- lwr2
preddata$upr <- upr2 
preddata$preds <- fit2
#make a model for Soapstone
modSoap_2 <- glm(Bolting_2018 ~ log_LongestLeaf_2018, data = temp2[temp2$Location=="Soapstone",], family = binomial)
summary(modSoap_2)

#predict new values
preddata2 <- with(temp2, data.frame(log_LongestLeaf_2018 = seq(min(temp2$log_LongestLeaf_2018), max(temp2$log_LongestLeaf_2018), length = 100)))
preds<- predict(modSoap_2, newdata = preddata2, type = "link", se.fit = TRUE)
critval <- 1.96 ## approx 95% CI
upr<- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit

fit2 <- modSoap_2$family$linkinv(fit)
upr2 <- modSoap_2$family$linkinv(upr)
lwr2 <- modSoap_2$family$linkinv(lwr)

preddata2$lwr <- lwr2
preddata2$upr <- upr2 
preddata2$preds <- fit2
#define limits of CIs
x95_1 <- c(preddata$log_LongestLeaf_2018,rev(preddata$log_LongestLeaf_2018))
y95_1 <- c(preddata$lwr,rev(preddata$upr))
x95_2 <- c(preddata2$log_LongestLeaf_2018,rev(preddata2$log_LongestLeaf_2018))
y95_2 <- c(preddata2$lwr,rev(preddata2$upr))
#make a plot
 ggplot() + 
  theme_classic()+
  scale_colour_manual(values = c("#f1511b", "#245f88")) +
  geom_point(aes(x = temp2$log_LongestLeaf_2018, y = temp2$Bolting_2019, col = temp2$Location)) +
  xlab("log(Longest Leaf (cm) in time = t)")  + 
  ylab("Probability of Flowering in time = t") +
  xlim (1,4) +
  geom_polygon(aes(x = x95_1, y = y95_1), fill = "#f1511b", alpha = 0.3) +
  geom_polygon(aes(x = x95_2, y = y95_2), fill = "#245f88", alpha = 0.3) +
  geom_line(aes(y = preddata$preds, x = preddata$log_LongestLeaf_2018), col = "#f1511b") +
  geom_line(aes(y = preddata2$preds, x = preddata2$log_LongestLeaf_2018), col = "#245f88") +
  theme(legend.position = "none") 


# #make a better type of plot for changes in density
# dat$class <- NA
# dat[dat$DataType %in% c("vegetative_18","vegetative_19"),"class"] <- "vegetative"
# dat[dat$DataType %in% c("reproductive_18","reproductive_19"),"class"] <- "reproductive"
# dat[dat$DataType %in% c("Seedlings_18","Seedlings_19"),"class"] <- "seedling"
# 
# plot(dat[dat$class=="vegetative","Density"]~
#        dat[dat$class=="vegetative","Year"],
#      ylab = "Plot-level Density",
#      xlab = "Year",
#      main = "Changes in plot-level density of vegetative individuals over time")
# #crow creek
# lines(x = c(2018,2019), y = c(mean(dat[dat$Site=="Crow Creek" & dat$class=="vegetative" & dat$Year==2018,"Density"]), mean(dat[dat$Site=="Crow Creek" & dat$class=="vegetative" & dat$Year==2019,"Density"])), col = "green")
# #diamond creek
# lines(x = c(2018,2019), y = c(mean(dat[dat$Site=="Diamond Creek" & dat$class=="vegetative" & dat$Year==2018,"Density"]), mean(dat[dat$Site=="Diamond Creek" & dat$class=="vegetative" & dat$Year==2019,"Density"])), col = "blue")
# #unnamed creek
# lines(x = c(2018,2019), y = c(mean(dat[dat$Site=="Unnamed Creek" & dat$class=="vegetative" & dat$Year==2018,"Density"]), mean(dat[dat$Site=="Unnamed Creek" & dat$class=="vegetative" & dat$Year==2019,"Density"])), col = "red")
# #HQ5
# lines(x = c(2018,2019), y = c(mean(dat[dat$Site=="HQ5" & dat$class=="vegetative" & dat$Year==2018,"Density"]), mean(dat[dat$Site=="HQ5" & dat$class=="vegetative" & dat$Year==2019,"Density"])), col = "orange")
# #HQ3
# lines(x = c(2018,2019), y = c(mean(dat[dat$Site=="HQ3" & dat$class=="vegetative" & dat$Year==2018,"Density"]), mean(dat[dat$Site=="HQ3" & dat$class=="vegetative" & dat$Year==2019,"Density"])), col = "brown")
# #Meadow
# lines(x = c(2018,2019), y = c(mean(dat[dat$Site=="Meadow" & dat$class=="vegetative" & dat$Year==2018,"Density"]), mean(dat[dat$Site=="Meadow" & dat$class=="vegetative" & dat$Year==2019,"Density"])), col = "black")
 
#### Plot-level changes in density, accross type-classes ####
par(mfrow = c(2,3), mar = c(2.3,3,2.5,1), cex = 1)
dat$class <- NA
dat[dat$DataType %in% c("vegetative_18","vegetative_19"),"class"] <- "vegetative"
dat[dat$DataType %in% c("reproductive_18","reproductive_19"),"class"] <- "reproductive"
dat[dat$DataType %in% c("Seedlings_18","Seedlings_19"),"class"] <- "seedling"
#crow creek
plot(y = c(mean(dat[dat$Site=="Crow Creek" & dat$class=="vegetative" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="Crow Creek" & dat$class=="vegetative" & dat$Year==2019,"Density"])),
x = c(2018,2019),
xlab = "year", 
ylab = "",
main = "Crow Creek",
type = "l", 
col = "#1B9E77",
ylim = c(0,110),
xaxt = "n",
lty = 4,
cex.axis = 1.25,
cex.lab = 1.25)
text(2018.15, 64, "vegetative", col = "#1B9E77", cex = 1.25)
lines(y = c(mean(dat[dat$Site=="Crow Creek" & dat$class=="reproductive" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="Crow Creek" & dat$class=="reproductive" & dat$Year==2019,"Density"])),
x = c(2018,2019),
col = "#1B9E77",
lty = 2)
text(2018.2,15, "reproductive", col="#1B9E77", cex = 1.25)
lines(y = c(mean(dat[dat$Site=="Crow Creek" & dat$class=="seedling" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="Crow Creek" & dat$class=="seedling" & dat$Year==2019,"Density"])),
x = c(2018,2019),
col = "#1B9E77",
lty = 3,
lwd = 2)
text(2018.15,36,"seedling", col = "#1B9E77", cex = 1.25)
lines( y = c(mean(dat[dat$Site=="Crow Creek" & dat$class=="seedling" & dat$Year==2018,"Density"]) + 
mean(dat[dat$Site=="Crow Creek" & dat$class=="vegetative" & dat$Year==2018,"Density"]) + 
mean(dat[dat$Site=="Crow Creek" & dat$class=="reproductive" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="Crow Creek" & dat$class=="seedling" & dat$Year==2019,"Density"]) + 
mean(dat[dat$Site=="Crow Creek" & dat$class=="vegetative" & dat$Year==2019,"Density"]) + 
mean(dat[dat$Site=="Crow Creek" & dat$class=="reproductive" & dat$Year==2019,"Density"])),
x = c(2018,2019),
col = "#1B9E77",
lty = 1,
lwd = 2)
axis(1,at = c(2018,2019), labels = c("2018", "2019"), cex.axis = 1.25)
text(2018.08, 103, "total", col = "#1B9E77", cex = 1.25)

#unnamed creek
plot(y = c(mean(dat[dat$Site=="Unnamed Creek" & dat$class=="vegetative" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="Unnamed Creek" & dat$class=="vegetative" & dat$Year==2019,"Density"])),
x = c(2018,2019),
xlab = "year", 
ylab = "",
main = "Unnamed Creek",
type = "l", 
col = "#E6AB02",
ylim = c(0,450),
xaxt = "n",
lty = 4,
lwd = 2,
cex.axis = 1.25,
cex.lab = 1.25)
lines(y = c(mean(dat[dat$Site=="Unnamed Creek" & dat$class=="reproductive" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="Unnamed Creek" & dat$class=="reproductive" & dat$Year==2019,"Density"])),
x = c(2018,2019),
col = "#E6AB02",
lty = 2,
lwd = 2)
lines(y = c(mean(dat[dat$Site=="Unnamed Creek" & dat$class=="seedling" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="Unnamed Creek" & dat$class=="seedling" & dat$Year==2019,"Density"])),
x = c(2018,2019),
col = "#E6AB02",
lty =3,
lwd = 2)
lines( y = c(mean(dat[dat$Site=="Unnamed Creek" & dat$class=="seedling" & dat$Year==2018,"Density"]) + 
mean(dat[dat$Site=="Unnamed Creek" & dat$class=="vegetative" & dat$Year==2018,"Density"]) + 
mean(dat[dat$Site=="Unnamed Creek" & dat$class=="reproductive" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="Unnamed Creek" & dat$class=="seedling" & dat$Year==2019,"Density"]) + 
mean(dat[dat$Site=="Unnamed Creek" & dat$class=="vegetative" & dat$Year==2019,"Density"]) + 
mean(dat[dat$Site=="Unnamed Creek" & dat$class=="reproductive" & dat$Year==2019,"Density"])),
x = c(2018,2019),
lty = 1,
lwd = 2,
col = "#E6AB02")
axis(1,at = c(2018,2019), labels = c("2018", "2019"), cex.axis = 1.25)

#Diamond creek
plot(y = c(mean(dat[dat$Site=="Diamond Creek" & dat$class=="vegetative" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="Diamond Creek" & dat$class=="vegetative" & dat$Year==2019,"Density"])),
x = c(2018,2019),
xlab = "year", 
ylab = "",
main = "Diamond Creek",
type = "l", 
col = "#D95F02",
ylim = c(0,200),
xaxt = "n",
lty = 4,
lwd = 2,
cex.axis = 1.25,
cex.lab = 1.25)
lines(y = c(mean(dat[dat$Site=="Diamond Creek" & dat$class=="reproductive" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="Diamond Creek" & dat$class=="reproductive" & dat$Year==2019,"Density"])),
x = c(2018,2019),
col = "#D95F02",
lty = 2,
lwd = 2)
lines(y = c(mean(dat[dat$Site=="Diamond Creek" & dat$class=="seedling" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="Diamond Creek" & dat$class=="seedling" & dat$Year==2019,"Density"])),
x = c(2018,2019),
col = "#D95F02",
lty = 3,
lwd = 2)
lines( y = c(mean(dat[dat$Site=="Diamond Creek" & dat$class=="seedling" & dat$Year==2018,"Density"]) + 
mean(dat[dat$Site=="Diamond Creek" & dat$class=="vegetative" & dat$Year==2018,"Density"]) + 
mean(dat[dat$Site=="Diamond Creek" & dat$class=="reproductive" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="Diamond Creek" & dat$class=="seedling" & dat$Year==2019,"Density"]) + 
mean(dat[dat$Site=="Diamond Creek" & dat$class=="vegetative" & dat$Year==2019,"Density"]) + 
mean(dat[dat$Site=="Diamond Creek" & dat$class=="reproductive" & dat$Year==2019,"Density"])),
x = c(2018,2019),
col = "#D95F02",
lty = 1,
lwd = 2)
axis(1,at = c(2018,2019), labels = c("2018", "2019"), cex.axis = 1.25)

#HQ5
plot(y = c(mean(dat[dat$Site=="HQ5" & dat$class=="vegetative" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="HQ5" & dat$class=="vegetative" & dat$Year==2019,"Density"])),
x = c(2018,2019),
xlab = "year", 
ylab = "",
main = "HQ5",
type = "l", 
col = "#E7298A",
lty = 4,
ylim = c(0,510),
xaxt = "n",
lwd = 2,
cex.axis = 1.25,
cex.lab = 1.25)
lines(y = c(mean(dat[dat$Site=="HQ5" & dat$class=="reproductive" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="HQ5" & dat$class=="reproductive" & dat$Year==2019,"Density"])),
x = c(2018,2019),
col = "#E7298A",
lty = 2,
lwd = 2)
lines(y = c(mean(dat[dat$Site=="HQ5" & dat$class=="seedling" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="HQ5" & dat$class=="seedling" & dat$Year==2019,"Density"])),
x = c(2018,2019),
col = "#E7298A",
lty = 3,
lwd = 2)
lines( y = c(mean(dat[dat$Site=="HQ5" & dat$class=="seedling" & dat$Year==2018,"Density"]) + 
mean(dat[dat$Site=="HQ5" & dat$class=="vegetative" & dat$Year==2018,"Density"]) + 
mean(dat[dat$Site=="HQ5" & dat$class=="reproductive" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="HQ5" & dat$class=="seedling" & dat$Year==2019,"Density"]) + 
mean(dat[dat$Site=="HQ5" & dat$class=="vegetative" & dat$Year==2019,"Density"]) + 
mean(dat[dat$Site=="HQ5" & dat$class=="reproductive" & dat$Year==2019,"Density"])),
x = c(2018,2019),
col = "#E7298A",
lty = 1,
lwd = 2)
axis(1,at = c(2018,2019), labels = c("2018", "2019"), cex.axis = 1.25)

#HQ3
plot(y = c(mean(dat[dat$Site=="HQ3" & dat$class=="vegetative" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="HQ3" & dat$class=="vegetative" & dat$Year==2019,"Density"])),
x = c(2018,2019),
xlab = "year", 
ylab = "",
main = "HQ3",
type = "l", 
col = "#7570B3",
ylim = c(0,125),
xaxt = "n",
lty = 4,
lwd = 2,
cex.axis = 1.25,
cex.lab = 1.25)
lines(y = c(mean(dat[dat$Site=="HQ3" & dat$class=="reproductive" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="HQ3" & dat$class=="reproductive" & dat$Year==2019,"Density"])),
x = c(2018,2019),
col = "#7570B3",
lty = 2,
lwd = 2)
lines(y = c(mean(dat[dat$Site=="HQ3" & dat$class=="seedling" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="HQ3" & dat$class=="seedling" & dat$Year==2019,"Density"])),
x = c(2018,2019),
col = "#7570B3",
lty = 3,
lwd = 2)
lines( y = c(mean(dat[dat$Site=="HQ3" & dat$class=="seedling" & dat$Year==2018,"Density"]) + 
mean(dat[dat$Site=="HQ3" & dat$class=="vegetative" & dat$Year==2018,"Density"]) + 
mean(dat[dat$Site=="HQ3" & dat$class=="reproductive" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="HQ3" & dat$class=="seedling" & dat$Year==2019,"Density"]) + 
mean(dat[dat$Site=="HQ3" & dat$class=="vegetative" & dat$Year==2019,"Density"]) + 
mean(dat[dat$Site=="HQ3" & dat$class=="reproductive" & dat$Year==2019,"Density"])),
x = c(2018,2019),
lty = 1,
col = "#7570B3",
lwd = 2)
axis(1,at = c(2018,2019), labels = c("2018", "2019"), cex.axis = 1.25)

#Meadow
plot(y = c(mean(dat[dat$Site=="Meadow" & dat$class=="vegetative" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="Meadow" & dat$class=="vegetative" & dat$Year==2019,"Density"])),
x = c(2018,2019),
xlab = "year", 
ylab = "",
main = "Meadow",
type = "l", 
col = "#66A61E",
lty = 4,
ylim = c(0,30),
xaxt = "n",
lwd = 2,
cex.axis = 1.25,
cex.lab = 1.25)
lines(y = c(mean(dat[dat$Site=="Meadow" & dat$class=="reproductive" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="Meadow" & dat$class=="reproductive" & dat$Year==2019,"Density"])),
x = c(2018,2019),
col = "#66A61E",
lty = 2,
lwd = 2)
lines(y = c(mean(dat[dat$Site=="Meadow" & dat$class=="seedling" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="Meadow" & dat$class=="seedling" & dat$Year==2019,"Density"])),
x = c(2018,2019),
col = "#66A61E",
lty = 3,
lwd = 2)
lines( y = c(mean(dat[dat$Site=="Meadow" & dat$class=="seedling" & dat$Year==2018,"Density"]) + 
mean(dat[dat$Site=="Meadow" & dat$class=="vegetative" & dat$Year==2018,"Density"]) + 
mean(dat[dat$Site=="Meadow" & dat$class=="reproductive" & dat$Year==2018,"Density"]),
mean(dat[dat$Site=="Meadow" & dat$class=="seedling" & dat$Year==2019,"Density"]) + 
mean(dat[dat$Site=="Meadow" & dat$class=="vegetative" & dat$Year==2019,"Density"]) + 
mean(dat[dat$Site=="Meadow" & dat$class=="reproductive" & dat$Year==2019,"Density"])),
x = c(2018,2019),
col = "#66A61E",
lty = 1,
lwd = 2)
axis(1,at = c(2018,2019), labels = c("2018", "2019"), cex.axis = 1.25)
dev.off()
#### effect of herbivory on reproductive output
#subset dataset to exclude reproductive individuals

counts$herb_2018 <- as.integer(ifelse(counts$Stem_Herbivory_2018==1 | counts$Invert_Herbivory_2018==1, 1, 0))
counts$herb_2019 <- as.integer(ifelse(counts$Stem_Herbivory_2019==1 | counts$Invert_Herbivory_2019==1, 1, 0))
noRepro <- counts[counts$Flowering_2018!=1,]
Repro <- counts[counts$Flowering_2018==1,]
plot((Repro$No_Capsules_2018) ~ jitter(Repro$herb_2018),
     xlab = c("Herbivory in 2018"), 
     ylab = c("Reproductive Output in 2018"))
modHerb <- glm()

####compare size in 2018 to size in 2019 (show different colors/different regressions for browsed vs unbrowsed) ####
plot(counts$LongestLeaf_cm_2019 ~ counts$LongestLeaf_cm_2018, col = as.factor(counts$herb_2018), pch = 16)
plot(jitter(Alive_2019) ~ herb_2018, data = counts[counts$Alive_2018 & counts$Flowering_2018==0,])
modHerb<- glm(Alive_2019 ~ herb_2018, data = counts[counts$Alive_2018 & counts$Flowering_2018==0,], family = "binomial")
abline(modHerb)
#not really a good relationship (not enough data!)
### density dependence as another co-variate (easy metric would be plot level counts of plants)--effect on survival to the next year ####
#make a variable for plot-level counts of indivdiuals (including seedlings!) and each plot
#multiply number of seedlings by .5, since they probably have a lower effect than adult plants
densTot #data frame that has density within each plot
densAdj <- data.frame(Plot_ID = c(unique(densTot$Plot_ID), unique(densTot$Plot_ID))) #data frame for density adjusted for use in this model (*seedlings by .5)
densAdj$Plot_ID <- densAdj[order(densAdj$Plot_ID),]
densAdj$Year <- rep(c(2018,2019),18) 
densAdj$Tot <- NA
#make a function to create an adjusted total
adjSum <- function(Veg, Repro, Seedling) {
  adjTotal <- Veg + Repro + 0.5*Seedling
  return(adjTotal)
}
#make empty data frame for summed data
New <- data.frame()
#calculate the adjusted total for each plot
for (i in 1:length(unique(densAdj$Plot_ID))) {
  for (j in 1:length(unique(densAdj$Year))) {
    plot <- unique(densAdj$Plot_ID)[i]
    year <- unique(densAdj$Year)[j]
    Veg <- filter(densTot, Plot_ID==plot & Year == year & DataType == paste("vegetative_", str_sub(year,3,4) ,sep = ""))["Density"]
    Repro <- filter(densTot, Plot_ID==plot & Year == year & DataType == paste("reproductive_", str_sub(year,3,4) ,sep = ""))["Density"]
    Seedling <- filter(densTot, Plot_ID==plot & Year == year & DataType == paste("Seedlings_", str_sub(year,3,4) ,sep = ""))["Density"]
    adjTot <-  adjSum(Veg, Repro, Seedling)
    names(adjTot) <- "Tot"
    if(nrow(New)<1){
      New <- data.frame(Plot_ID = plot, Year = year, Tot = adjTot)
    }else{
      New <- rbind(New, data.frame(Plot_ID = plot, Year = year, Tot = adjTot)) 
    }
    }
}

#reformate New dataset
New <- spread(New, "Year", "Tot", convert = TRUE)
names(New) <- c("Plot_ID", "plotDens_2018", "plotDens_2019")
#join "New" data frame (which contains adjusted totals) to "Counts" data.frame
counts <- left_join(counts, New, by = "Plot_ID")

# make a figure comparing the number of individuals in a plot to the average survival of plants in the plot
#remove the plants that flowered in 2018
noFlr <- filter(counts, Alive_2018==1 & Bolting_2018==0)

plot(jitter(noFlr$Alive_2019)~ noFlr$plotDens_2018)
modDens <- glm(Alive_2019~plotDens_2018, data = noFlr, family = binomial)
plot(modDens)

preddataDens <- with(noFlr, data.frame(plotDens_2018 = seq(min(plotDens_2018), max(plotDens_2018), length = 100)))
preds <- predict(modDens, newdata = preddataDens, type = "link", se.fit = TRUE)
critval <- 1.96 ## approx 95% CI
upr <- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit

fit2 <- modDens$family$linkinv(fit)
upr2 <- modDens$family$linkinv(upr)
lwr2 <- modDens$family$linkinv(lwr)

preddataDens$lwr <- lwr2
preddataDens$upr <- upr2 
preddataDens$preds <- fit2

#plot the results of the regression, with 95% confidence bands
#define the boundaries of the polygons for the 95% confidence intervals 
#
x95_1 <- c(preddataDens$plotDens_2018,rev(preddataDens$plotDens_2018))
y95_1 <- c(preddataDens$lwr,rev(preddataDens$upr))

plot(jitter(noFlr$Alive_2019)~ noFlr$plotDens_2018, type = 'n')
lines(preddataDens$preds~preddataDens$plotDens_2018)
lines(preddataDens$upr ~ preddataDens$plotDens_2018, lty = 2)
lines(preddataDens$lwr ~ preddataDens$plotDens_2018, lty = 2)
points(jitter(noFlr$Alive_2019)~ noFlr$plotDens_2018)

#pretty inconclusive
#
#
##### probability of surviving from seedling to next year ####
sdNewVeg <-  filter(counts, Alive_2018==0 & Alive_2019==1 & Bolting_2019==0)
sdNewVegPlot <- aggregate(sdNewVeg[,"LongestLeaf_cm_2019"], by = list(sdNewVeg$Plot_ID), FUN = length)
names(sdNewVegPlot) <- c("Plot_ID", "Count_Veg_19")

sdNewRepro<-  filter(counts, Alive_2018==0 & Alive_2019==1 & Bolting_2019==1)
sdNewReproPlot <- aggregate(sdNewRepro[,"LongestLeaf_cm_2019"], by = list(sdNewRepro$Plot_ID), FUN = length)
names(sdNewReproPlot) <- c("Plot_ID", "Count_Repro_19")

sdNewAll <- left_join(sdNewVegPlot, sdNewReproPlot, by = "Plot_ID")

#get seedling data
sdNewSeeds <- aggregate(seeds$Seedlings_18, by = list(seeds$Plot_ID), FUN = sum)
names(sdNewSeeds) <- c("Plot_ID", "Seedlings_18")

sdNewAll <- left_join(sdNewAll, sdNewSeeds, by = "Plot_ID")

#calculate probability of surviving from seedling to vegetative next year
sdNewAll$Prob_Veg <- sdNewAll$Count_Veg_19/sdNewAll$Seedlings_18

#### plot reproductive output by plant size, incorporating 2019 data ####
repro_18 <- counts %>% 
  filter(Flowering_2018==1) %>% 
  dplyr::select(Location, Site, Plot_ID, Quadrant, ID, No_Capsules_2018, LongestLeaf_cm_2018) 
names(repro_18) <- c("Location", "Site", "Plot_ID", "Quadrant", "ID", "No_Capsules", "LongestLeaf") 
repro_18$Year <- 2018

repro_19 <- counts %>% 
  filter(Flowering_2019==1) %>% 
  dplyr::select(Location, Site, Plot_ID, Quadrant, ID, No_Capsules_2019, LongestLeaf_cm_2019)
names(repro_19) <- c("Location", "Site", "Plot_ID", "Quadrant", "ID", "No_Capsules", "LongestLeaf") 
repro_19$Year <- 2019
repro_19$No_Capsules <- as.numeric(repro_19$No_Capsules)

repro <- rbind(repro_18, repro_19)
repro$Year <- as.factor(repro$Year)
repro$No_Capsules <- as.numeric(repro$No_Capsules)

# size in current year with reproductive output in current year
plot(repro$No_Capsules ~ repro$LongestLeaf,
     ylim = c(0,600), col = as.factor(repro$Year), pch = 16,
     xlab = "Plant size in current year (cm)",
     ylab = "Number of Capsules Produced in current year")
abline(glm(No_Capsules ~ LongestLeaf, dat = repro_18[repro_18$Location=="Soapstone",], family =), col = "red", lwd = 1.5, lty = 2)
abline(glm(No_Capsules ~ LongestLeaf, dat = repro_18[repro_18$Location=="FEWAFB",]), col = "red", lwd = 1.5, lty = 3)
abline(glm(No_Capsules ~ LongestLeaf, dat = repro_19[repro_19$Location=="Soapstone",]), col = "black", lwd = 1.5, lty = 2)
abline(glm(No_Capsules ~ LongestLeaf, dat = repro_19[repro_19$Location=="FEWAFB",]), col = "black", lwd = 1.5, lty = 3)
abline(glm(No_Capsules ~ LongestLeaf, data = repro), lwd = 2)

ggplot(data = repro) + 
  geom_point(aes(x = LongestLeaf, y = No_Capsules, col = Year, pch = Location)) +
  theme_classic() + 
  geom_smooth(aes(x = LongestLeaf, y = No_Capsules), data = repro[repro$Year==2018 & repro$Location=="Soapstone",], method = lm, col = "red", lty = 2) + 
  geom_smooth(aes(x = LongestLeaf, y = No_Capsules), data = repro[repro$Year==2018 & repro$Location=="FEWAFB",], method = lm, col = "red", lty = 3) + 
  geom_smooth(aes(x = LongestLeaf, y = No_Capsules), data = repro[repro$Year==2019 & repro$Location=="Soapstone",], method = lm, col = "blue", lty = 2) + 
  geom_smooth(aes(x = LongestLeaf, y = No_Capsules), data = repro[repro$Year==2019 & repro$Location=="FEWAFB",], method = lm, col = "", lty = 3)

##try to fit poisson models to data
#using normal glms
modR_18_FEWAFB <- lm(No_Capsules ~ LongestLeaf, dat = repro_18[repro_18$Location=="FEWAFB",])
summary(modR_18_FEWAFB)
modR_18_Soapstone<- lm(No_Capsules ~ LongestLeaf, dat = repro_18[repro_18$Location=="Soapstone",])
summary(modR_18_Soapstone)

# modR_19_FEWAFB <- lm(No_Capsules ~ LongestLeaf, dat = repro_19[repro_19$Location=="FEWAFB",])
# summary(modR_19_FEWAFB) #no data for this model yet
modR_19_Soapstone<- lm(No_Capsules ~ LongestLeaf, dat = repro_19[repro_19$Location=="Soapstone",])
summary(modR_19_Soapstone)

modR_All <- lm(No_Capsules ~ LongestLeaf, dat = repro)
summary(modR_All)

#using poisson glms
#check that all of the capsule numbers are discrete counts
str(repro)
repro$No_Capsules <- round(as.numeric(repro$No_Capsules),0)
repro_18$No_Capsules <- round(as.numeric(repro_18$No_Capsules),0)
repro_19$No_Capsules <- round(as.numeric(repro_19$No_Capsules),0)

PmodR_18_FEWAFB <- glm(No_Capsules ~ LongestLeaf, dat = repro_18[repro_18$Location=="FEWAFB",], family = "poisson")
summary(PmodR_18_FEWAFB)
PmodR_18_Soapstone<- glm(No_Capsules ~ LongestLeaf, dat = repro_18[repro_18$Location=="Soapstone",], family = "poisson")
summary(PmodR_18_Soapstone)

# PmodR_19_FEWAFB <- glm(No_Capsules ~ LongestLeaf, dat = repro_19[repro_19$Location=="FEWAFB",], family = "poisson")
# summary(modR_19_FEWAFB) #no data for this model yet
PmodR_19_Soapstone<- glm(No_Capsules ~ LongestLeaf, dat = repro_19[repro_19$Location=="Soapstone",], family = "poisson")
summary(PmodR_19_Soapstone)

Pmod_All <- glm(No_Capsules ~ LongestLeaf, dat = repro, family = "poisson")
summary(Pmod_All)

##AIC indicates that poisson models are way worse than normal lms

##try and plot poisson regressions
Pmod_All$model$fitted <- predict(Pmod_All, type = "response")

ggplot()+
  geom_point(aes(x = LongestLeaf, y = No_Capsules, col = Location, pch = Year), data = repro) + 
  #geom_smooth(aes(repro$LongestLeaf, repro$No_Capsules), method = lm, col = "black") +
  geom_line(aes(y = Pmod_All$model$fitted, x = Pmod_All$model$LongestLeaf))
# geom_smooth(aes(LongestLeaf, No_Capsules), data = repro_18[repro_18$Location=="Soapstone",], method = lm, col = "blue", lty = 2) + 
# geom_smooth(aes(LongestLeaf, No_Capsules), data = repro_18[repro_18$Location=="FEWAFB",], method = lm, col = "red", lty = 2) + 
# geom_smooth(aes(LongestLeaf, No_Capsules), data = repro_19[repro_19$Location=="Soapstone",], method = lm, col = "blue", lty = 3) + 
# geom_smooth(aes(LongestLeaf, No_Capsules), data = repro_19[repro_19$Location=="FEWAFB",], method = lm, col = "red", lty = 3)


# size in previous year with reproductive output in next year
repro_2 <- counts %>% 
  filter(Flowering_2019==1) %>% 
  select(Location, Site, Plot_ID, ID, No_Capsules_2019, LongestLeaf_cm_2018, LongestLeaf_cm_2019)

plot(repro_2$No_Capsules_2019 ~ repro_2$LongestLeaf_cm_2018, pch = 16, col = "darkgreen")
points(repro_2$No_Capsules_2019 ~ repro_2$LongestLeaf_cm_2019, pch = 16, col = "blue")
abline(lm(repro_2$No_Capsules_2019 ~ repro_2$LongestLeaf_cm_2018), col = "darkgreen")
abline(lm(repro_2$No_Capsules_2019 ~ repro_2$LongestLeaf_cm_2019), col = "blue")

##not enough data yet to show any relationship

#### Make Vital Rate Function Plots ####
#### reproductive output by size accross two sites ####
flowering <-counts %>%
  filter(Flowering_2018==1, is.na(No_Capsules_2018)==FALSE, No_Capsules_2018!=0)
flowering$col <- ifelse(flowering$Location=="FEWAFB", "#f1511b", "#245f88")
flowering$log_LongestLeaf_2018 <- log(flowering$LongestLeaf_cm_2018)
flowering$log_Capsules_18 <- log(flowering$No_Capsules_2018)
plot(flowering$log_Capsules_18 ~ flowering$log_LongestLeaf_2018, col = flowering$col, pch = 16)

#make a model for FEWAFB
modAFB <- glm(No_Capsules_2018 ~ log_LongestLeaf_2018, data = flowering[flowering$Location=="FEWAFB",])
summary(modAFB)
#predict new values
preddata_D <- with(flowering, data.frame(log_LongestLeaf_2018 = seq(min(log_LongestLeaf_2018), max(log_LongestLeaf_2018), length = 100)))
preds <- predict(modAFB, newdata = preddata, se.fit = TRUE)
critval <- 1.96 ## approx 95% CI
upr <- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit

preddata_D$lwr <- lwr
preddata_D$upr <- upr 
preddata_D$preds <- fit

#make a model for Soapstone
modSoap <- glm(No_Capsules_2018 ~ log_LongestLeaf_2018, data = flowering[flowering$Location=="Soapstone",])
summary(modSoap)
#predict new values
preddata2_D <- with(flowering, data.frame(log_LongestLeaf_2018 = seq(min(log_LongestLeaf_2018), max(log_LongestLeaf_2018), length = 100)))
preds_2 <- predict(modSoap, newdata = preddata, se.fit = TRUE)
critval <- 1.96 ## approx 95% CI
upr_2 <- preds_2$fit + (critval * preds_2$se.fit)
lwr_2 <- preds_2$fit - (critval * preds_2$se.fit)
fit_2 <- preds_2$fit

preddata2_D$lwr <- lwr_2
preddata2_D$upr <- upr_2 
preddata2_D$preds <- fit_2

#plot the results of the regression, with 95% confidence bands
#define the boundaries of the polygons for the 95% confidence intervals 
#
x95_1_D <- c(preddata_D$log_LongestLeaf_2018,rev(preddata_D$log_LongestLeaf_2018))
y95_1_D <- c(preddata_D$lwr,rev(preddata_D$upr))
x95_2_D <- c(preddata2_D$log_LongestLeaf_2018,rev(preddata2_D$log_LongestLeaf_2018))
y95_2_D <- c(preddata2_D$lwr,rev(preddata2_D$upr))

plotD <- ggplot() + 
  theme_classic()+
  scale_colour_manual(values = c("#f1511b", "#245f88")) +
  geom_point(aes(x = flowering$log_LongestLeaf_2018, y = flowering$No_Capsules_2018, col = flowering$Location), alpha = 0.5) +
  xlab('log(Longest Leaf (cm) in  t)') + 
  ylab('No. of Capsules Produced in t') +
  geom_polygon(aes(x = x95_1_D, y = y95_1_D), fill = "#f1511b", alpha = 0.3) +
  geom_polygon(aes(x = x95_2_D, y = y95_2_D), fill = "#245f88", alpha = 0.3) +
  geom_line(aes(y = preddata_D$preds, x = preddata_D$log_LongestLeaf_2018), col = "#f1511b") +
  geom_line(aes(y = preddata2_D$preds, x = preddata2_D$log_LongestLeaf_2018), col = "#245f88") +
  theme(legend.position = "none") +
  theme(legend.position = "none")  

#### Growth from 2018 to 2019 ####
size <- counts %>% filter(counts$Alive_2018==1, counts$Alive_2019==1)
size$log_LongestLeaf_2018 <- log(size$LongestLeaf_cm_2018)
size$log_LongestLeaf_2019 <- log(size$LongestLeaf_cm_2019)

#make a model for FEWAFB
modAFB_C<- lm(log_LongestLeaf_2019 ~ log_LongestLeaf_2018, data = size[size$Location=="FEWAFB",])
summary(modAFB_C)
#predict new values
preddata_C <- with(size, data.frame(log_LongestLeaf_2018 = seq(min(size$log_LongestLeaf_2018), max(size$log_LongestLeaf_2018), length = 100)))
preds<- predict(modAFB_C, newdata = preddata_C, se.fit = TRUE)
critval <- 1.96 ## approx 95% CI
upr<- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit

preddata_C$lwr <- lwr
preddata_C$upr <- upr 
preddata_C$preds <- fit

#make a model for Soapstone
modSoap <- lm(log_LongestLeaf_2019 ~ log_LongestLeaf_2018, data = size[size$Location=="Soapstone",])
summary(modSoap)

#predict new values
preddata_C2 <- with(size, data.frame(log_LongestLeaf_2018 = seq(min(size$log_LongestLeaf_2018), max(size$log_LongestLeaf_2018), length = 100)))
preds<- predict(modSoap, newdata = preddata_C2, se.fit = TRUE)
critval <- 1.96 ## approx 95% CI
upr<- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit

preddata_C2$lwr <- lwr
preddata_C2$upr <- upr
preddata_C2$preds <- fit
#define limits of CIs
x95_C_1 <- c(preddata_C$log_LongestLeaf_2018,rev(preddata_C$log_LongestLeaf_2018))
y95_C_1 <- c(preddata_C$lwr,rev(preddata_C$upr))
x95_C_2 <- c(preddata_C2$log_LongestLeaf_2018,rev(preddata_C2$log_LongestLeaf_2018))
y95_C_2 <- c(preddata_C2$lwr,rev(preddata_C2$upr))
#make a plot
plotC <- ggplot() + 
  theme_classic()+
  scale_colour_manual(values = c("#f1511b", "#245f88")) +
  geom_point(aes(x = size$log_LongestLeaf_2018, y = size$log_LongestLeaf_2019, col = size$Location), alpha = .5) +
  xlab("log(Longest Leaf (cm) in t)")  + 
  ylab("log(Longest Leaf (cm) in t+1") +
  geom_polygon(aes(x = x95_C_1, y = y95_C_1), fill = "#f1511b", alpha = 0.3) +
  geom_polygon(aes(x = x95_C_2, y = y95_C_2), fill = "#245f88", alpha = 0.3) +
  geom_line(aes(y = preddata_C$preds, x = preddata_C$log_LongestLeaf_2018), col = "#f1511b") +
  geom_line(aes(y = preddata_C2$preds, x = preddata_C2$log_LongestLeaf_2018), col = "#245f88") +
  theme(legend.position = "none") 
#### plot size in previous year and reproductive probability ####
# subset for plants that are alive in both 2018 and 2019
temp2 <- counts %>% filter(counts$Alive_2019==1, counts$Alive_2018==1)
temp2$col <- ifelse(temp2$Location=="FEWAFB", "#f1511b", "#245f88")
temp2$log_LongestLeaf_2018 <- log(temp2$LongestLeaf_cm_2018)
plot(temp2$Bolting_2019 ~ temp2$log_LongestLeaf_2018, col = temp2$col, pch = 19,
     xlab = "2018 Leaf Length (cm)",
     ylab = "Bolting in 2019 (0-1)")

#make a model for FEWAFB
modAFB_B <- glm(Bolting_2019 ~ log_LongestLeaf_2018, data = temp2[temp2$Location=="FEWAFB",], family = binomial)
summary(modAFB_B)
#predict new values
preddata_B <- with(temp2, data.frame(log_LongestLeaf_2018 = seq(min(temp2$log_LongestLeaf_2018), max(temp2$log_LongestLeaf_2018), length = 100)))
preds<- predict(modAFB_B, newdata = preddata_B, type = "link", se.fit = TRUE)
critval <- 1.96 ## approx 95% CI
upr<- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit

fit2 <- modAFB_B$family$linkinv(fit)
upr2 <- modAFB_B$family$linkinv(upr)
lwr2 <- modAFB_B$family$linkinv(lwr)

preddata_B$lwr <- lwr2
preddata_B$upr <- upr2 
preddata_B$preds <- fit2
#make a model for Soapstone
modSoap_B <- glm(Bolting_2019 ~ log_LongestLeaf_2018, data = temp2[temp2$Location=="Soapstone",], family = binomial)
summary(modSoap_B)

#predict new values
preddata_B2 <- with(temp2, data.frame(log_LongestLeaf_2018 = seq(min(temp2$log_LongestLeaf_2018), max(temp2$log_LongestLeaf_2018), length = 100)))
preds<- predict(modSoap_B, newdata = preddata_B2, type = "link", se.fit = TRUE)
critval <- 1.96 ## approx 95% CI
upr<- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit

fit2 <- modSoap_B$family$linkinv(fit)
upr2 <- modSoap_B$family$linkinv(upr)
lwr2 <- modSoap_B$family$linkinv(lwr)

preddata_B2$lwr <- lwr2
preddata_B2$upr <- upr2 
preddata_B2$preds <- fit2
#define limits of CIs
x95_B_1 <- c(preddata_B$log_LongestLeaf_2018,rev(preddata_B$log_LongestLeaf_2018))
y95_B_1 <- c(preddata_B$lwr,rev(preddata_B$upr))
x95_B_2 <- c(preddata_B2$log_LongestLeaf_2018,rev(preddata_B2$log_LongestLeaf_2018))
y95_B_2 <- c(preddata_B2$lwr,rev(preddata_B2$upr))
#make a plot
plotB <- ggplot() + 
  theme_classic()+
  scale_colour_manual(values = c("#f1511b", "#245f88")) +
  geom_point(aes(x = temp2$log_LongestLeaf_2018, y = temp2$Bolting_2019, col = temp2$Location)) +
  xlab("log(Longest Leaf (cm) in t)")  + 
  ylab("Probability of Flowering in t+1") +
  xlim (1,4) +
  geom_polygon(aes(x = x95_B_1, y = y95_B_1), fill = "#f1511b", alpha = 0.3) +
  geom_polygon(aes(x = x95_B_2, y = y95_B_2), fill = "#245f88", alpha = 0.3) +
  geom_line(aes(y = preddata_B$preds, x = preddata_B$log_LongestLeaf_2018), col = "#f1511b") +
  geom_line(aes(y = preddata_B2$preds, x = preddata_B2$log_LongestLeaf_2018), col = "#245f88") + 
  theme(legend.position = "none") 

#### plot size in previous year by surival probability ####
# remove those plants that were reproductive in 2018
temp <- counts %>% filter(counts$Flowering_2018==0)
temp$col <- ifelse(temp$Location=="FEWAFB", "#f1511b", "#245f88")
temp$log_LongestLeaf_2018 <- log(temp$LongestLeaf_cm_2018)

par(mfrow=c(1,1), mar = c(5.1,4.1,4.1,2.1))
plot(jitter(temp$Alive_2019, .5) ~ temp$log_LongestLeaf_2018, col = temp$col, pch = 19,
     xlab = "2018 Leaf Length (cm)",
     ylab = "Alive in 2019 (0-1)")

#make a model for FEWAFB
modAFB_A <- glm(Alive_2019 ~ log_LongestLeaf_2018, data = temp[temp$Location=="FEWAFB",], family = binomial)
summary(modAFB_A)
#predict new values
preddata_A <- with(temp, data.frame(log_LongestLeaf_2018 = seq(min(log_LongestLeaf_2018), max(log_LongestLeaf_2018), length = 100)))
preds <- predict(modAFB_A, newdata = preddata_A, type = "link", se.fit = TRUE)
critval <- 1.96 ## approx 95% CI
upr <- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit

fit2 <- modAFB_A$family$linkinv(fit)
upr2 <- modAFB_A$family$linkinv(upr)
lwr2 <- modAFB_A$family$linkinv(lwr)

preddata_A$lwr <- lwr2 
preddata_A$upr <- upr2 
preddata_A$preds <- fit2

#plot the regression
plot(temp$Alive_2019 ~ log(temp$LongestLeaf_cm_2018), col = temp$col, pch = 19,
     xlab = "2018 Leaf Length (cm)",
     ylab = "Alive in 2019 (0-1)")
lines(probs~log(LongestLeaf_cm_2018), data = newdata, col = "tomato")

#make a model for Soapstone
modSoap_A <- glm(Alive_2019 ~ log_LongestLeaf_2018, data = temp[temp$Location=="Soapstone",], family = binomial)
summary(modSoap_A)
#predict new values
preddata_A2 <- with(temp, data.frame(log_LongestLeaf_2018 = seq(min(log_LongestLeaf_2018), max(log_LongestLeaf_2018), length = 100)))
preds_2 <- predict(modSoap_A, newdata = preddata_A, type = "link", se.fit = TRUE)
critval <- 1.96 ## approx 95% CI
upr_2 <- preds_2$fit + (critval * preds_2$se.fit)
lwr_2 <- preds_2$fit - (critval * preds_2$se.fit)
fit_2 <- preds_2$fit

fit2_2 <- modSoap_A$family$linkinv(fit_2)
upr2_2 <- modSoap_A$family$linkinv(upr_2)
lwr2_2 <- modSoap_A$family$linkinv(lwr_2)

preddata_A2$lwr <- lwr2_2
preddata_A2$upr <- upr2_2 
preddata_A2$preds <- fit2_2

#plot the results of the regression, with 95% confidence bands
#define the boundaries of the polygons for the 95% confidence intervals 
#
x95_A_1 <- c(preddata_A$log_LongestLeaf_2018,rev(preddata_A$log_LongestLeaf_2018))
y95_A_1 <- c(preddata_A$lwr,rev(preddata_A$upr))
x95_A_2 <- c(preddata_A2$log_LongestLeaf_2018,rev(preddata_A2$log_LongestLeaf_2018))
y95_A_2 <- c(preddata_A2$lwr,rev(preddata_A2$upr))

plotA <- ggplot() + 
  theme_classic()+
  scale_colour_manual(values = c("#f1511b", "#245f88")) +
  geom_point(data = temp, aes(x = log_LongestLeaf_2018, y = Alive_2019, col = Location)) +
  xlab("log(Longest Leaf in t (cm))")  + 
  ylab("Survived to t+1") +
  xlim (1,4) +
  geom_polygon(aes(x = x95_A_1, y = y95_A_1), fill = "#f1511b", alpha = 0.3) +
  geom_polygon(aes(x = x95_A_2, y = y95_A_2), fill = "#245f88", alpha = 0.3) +
  geom_line(aes(y = preddata_A$preds, x = preddata_A$log_LongestLeaf_2018), col = "#f1511b") +
  geom_line(aes(y = preddata_A2$preds, x = preddata_A2$log_LongestLeaf_2018), col = "#245f88") +
  theme(legend.position = "none") 



#### Make figure showing all preliminary vital rate plots in one panel ####
require(ggplot2)
require(ggpubr) #package used to put ggplots in one panel
theme_set(theme_pubr())

# plots are plotA, plotB, plotC, and plotD
VitalRates <- ggarrange(plotA, plotB, plotC, plotD,
                        labels = c("A", "B", "C", "D"),
                        common.legend = TRUE, legend = "right")

##---- VitalRatesFigure
VitalRates
