#-----------------------------------------------------------------
## Setup and data cleaning
#-----------------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lmer)
library(lmerTest)
source(here("analysis/", "Functions.R"))
library(car) 
source(here("analysis/", "density_analysis.R"))
#source(here("analysis", "Functions.R")) linked to exp analysis so we can extract model predicitons

#We are going to apply the models extracted from the size and density analysis to 9 LTER coastal sites. We are going to track urchin density through time and apply herbivory pressure estiamtes to each site based on our findings. We additionally are incoorperating observational data realting to detriatal supply to test its bearing on the amount of kelp biomass lost in a given area. 

betap <-as.vector(coef(lm1)[1])
betar<-as.vector(coef(lm1.r)[1])

purple.fun <- function(biomass){
  betap*biomass
} # this is the herbivory rate model prediction formed from the density analysis for purple urchins

red.fun <- function(biomass){
  betar*biomass} # this is the herbivory rate model prediction formed from the density analysis for red urchins


lt <- read.csv("data/survey_data/Annual_All_Species_Biomass_at_transect_20210108.csv", stringsAsFactors = F,na.strings ="-99999") %>% #LTER data density estimations collected from 50 transects across 11 sites between 2000-2018 in the Santa Barbara Channel.
  dplyr::select("YEAR", "MONTH", "SITE", "TRANSECT", "SP_CODE", "WM_GM2") %>%
  filter(SP_CODE %in% c("MAPY", "SFL", "SPL"))%>% #filtering the data to only include giant kelp, purple urchins, and red urchins
  filter(SITE != "SCTW", SITE != "SCDI") %>% 
  #Filtering out island sites. This study focuses on costal sites.
  rename_all(tolower) %>% 
  group_by(year, month, site, transect, sp_code) %>% 
  spread(sp_code, wm_gm2) %>%
  mutate(SFL.pred = red.fun(SFL),
         SPL.pred = purple.fun(SPL), 
         predicted.consumption = (SPL.pred + SFL.pred))


#-----------------------------------------------------------------
## Make the map

library(here)
library(rgdal)
library(rgeos)
library(sp)
library(maptools)
library(raster)
library(sf)
library(marmap)
library(grid)
library(gridBase)
library(rworldmap)
library(RColorBrewer)
library(maps)  
library(geosphere)
library(tidyverse)



sites <- read.csv(here("data/spatial", "lter_waypoints.csv")) 

mp <- lt %>%
  group_by(site, transect) %>%
  summarize(pc_mean = mean(predicted.consumption)) %>%
  left_join(sites) #Adding in the long lat data to the LTER density data

mp.mean <- lt %>%
  group_by(site) %>%
  summarize(pc_mean = mean(predicted.consumption), 
            se = sd(predicted.consumption)/n(), 
            sd = sd(predicted.consumption)) %>%
  left_join(sites %>%
              group_by(site) %>%
              summarize(lat = mean(lat), 
                        long = mean(long))) #Adding in the long lat data to the LTER density data

mp.mean$col <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6')

coordinates(mp) <- ~long + lat
proj4string(mp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
coordinates(mp.mean) <- ~long + lat
proj4string(mp.mean) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"


#shoreline and MPAS
us <- readOGR(here("data/spatial", "nos80k (1)"))
all_mpas <- readOGR(here("data/spatial", "state_mpas.shp"))
all_mpas <- all_mpas[all_mpas$NAME == "Naples SMCA" | all_mpas$NAME == "Campus Point SMCA (No-Take)", ]

proj <- proj4string(all_mpas)
us <- spTransform(us, proj)
clipper <- extent(-120.4955, -119.4544, 34.32426, 34.5)
shoreline <- crop(us, clipper)

# Kelp patches
p <- as(clipper, 'SpatialPolygons')
proj4string(p) <-  proj
p <- spTransform(p, "+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs")

patches <- readOGR(here("data/spatial", "patches.shp"))
proj4string(patches) <- "+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs"
polygone1 <- gBuffer(patches, byid=TRUE, width=0)
polygone2 <- gBuffer(p, byid=TRUE, width=0)
clip2 <- gIntersection(polygone1, polygone2, byid=TRUE)
patches <- spTransform(clip2, proj)


# Depth contours
contours <- readOGR(here("data/spatial", "contours.shp"))
contours <- spTransform(contours, proj)
contours <- crop(contours, clipper)

contours <- contours[contours$ELEV %in% c(seq(-10, -100, by = -10), -200, -300, -400, -500), ]


#c(bottom, left, top, right)

# Build the map
png(here::here("figures/", "fig3_pA.png"), width = 1200, height = 300)
d <- par(las = 1, mgp = c(3, 0.75, 0), mar = c(3, 6, 1, 1), tcl = -0.5)
plot(shoreline, col = "#e6c16c",  xlab = "", ylab = "", cex.axis = 1.5, axes = TRUE, lwd = 0.01, 
     xlim = c(-120.5, -119.4563), ylim = c(34.33406, 34.5)) # cal state
#plot(all_mpas, pch = 4, col = "#99ebff", add = T, border = "#99ebff") 
plot(contours, add = T, lwd = 0.25, ext = clipper, col = "gray80")
plot(patches, col = alpha("#00802b", 0.9), add = T, border = alpha("#00802b", 0))
plot(mp.mean, pch = 21, add = T, cex = (mp.mean$pc_mean + mp.mean$sd), bg = alpha("black", 0.25), col = NULL)
plot(mp.mean, pch = 21, add = T, lwd = 2, cex = mp.mean$pc_mean, bg = alpha(mp.mean$col, 0.6))
plot(mp.mean, pch = 16, add = T, cex = 0.5)
axis(tcl = -0.5, side = 1, labels = F)
axis(side = 2, tcl = -0.5, labels = F)
scalebar(xy = c(-120.45, 34.35), type = "bar", d = 8, divs = 4, lonlat = T, below = "kilometers")
par(d)
dev.off()




























































#--------------------------------------------------------------
## Old Code do not run
#--------------------------------------------------------------

# cal <- readOGR(here("data/spatial", "caloutline.shp"))
# all_mpas <- readOGR(here("data/spatial", "state_mpas.shp"))
# 
# d <- par(las = 1, mgp = c(3, 0.75, 0))
# plot(cal, col = "#FFEB9B", xlim = c(-120.5,-119.5), ylim = c(34.35,34.55), xlab = "", ylab = "", cex.axis = 1.5, axes = TRUE) # cal state
# plot(all_mpas, pch = 4, col = "#99ebff", add = T) 
# plot(mp, pch = 21, cex = mp$predited.consumption*0.25, add = T, lwd = 1.5,  bg = alpha("#e31a1c", .5))
# 
# par(d)
# 
# # Make a plot for each year...
# for(i in 2002:2018){
#   myfile <- file.path("figures/maps", paste("year", "_", i, ".png"))
#   png(myfile, width = 1000*3, height = 561*3, res = 300)
#   d <- par(las = 1, mgp = c(3, 0.75, 0), mar = c(4,5,3,1))
#   plot(cal, col = "#FFEB9B", xlim = c(-120.5,-119.5), ylim = c(34.35,34.55), xlab = "", ylab = "", cex.axis = 1.5, axes = TRUE) # cal state
#   plot(all_mpas, pch = 4, col = "#99ebff", add = T) 
#   plot(mp[mp$year == i,], pch = 21, cex = mp$predited.consumption[mp$year == i]*0.25, add = T, lwd = 1.5,  bg = alpha("#e31a1c", .25))
#   par(d)
#   dev.off()
# }
# 
# 
# # Make a time averaged plot...
# 
# av <- mp@data %>% 
#   mutate(group = cut(year, breaks = 4, labels = FALSE)) %>%
#   group_by(site, transect, group) %>%
#   summarize(predicted.consumption = mean(predited.consumption, na.rm = T)) %>%
#   left_join(sites)
# 
# coordinates(av) <- ~long + lat
# 
# proj4string(av) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
# 
# #c(bottom, left, top, right)
# id <- c("2000-2004", "2005-2009", "2010-2013", "2014-2018")
# 
# png(here("figures", "map.png"), width = 1000*1.2, 600*1.2)
# d <- par(mfrow = c(2, 2), las = 1, mgp = c(3, 0.75, 0), mar = c(3,5, 1.5, 0.5))
# for(i in 1:4){
#   plot(cal, col = "#FFEB9B", xlim = c(-120.5,-119.5), ylim = c(34.35,34.55), xlab = "", ylab = "", cex.axis = 1.5, axes = TRUE, main = paste(id[i])) # cal state
#   plot(all_mpas, pch = 4, col = "#99ebff", add = T) 
#   plot(av[av$group == i, ], pch = 21, cex = av$predicted.consumption[av$group == i], add = T, lwd = 1.5,  bg = alpha("#e31a1c", .5))
# }
# par(d)
# dev.off()

