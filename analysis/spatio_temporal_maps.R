## Setup and data cleaning
#-----------------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lmer)
library(lmerTest)
source(here("analysis", "Functions.R"))
library(car) 
source(here("analysis", "Functions.R"))
#source(here("analysis", "Functions.R")) linked to exp analysis so we can extract model predicitons

#We are going to make a map to visually respresent the predicted herbivory pressure of the 9 LTER coastal sites predicted by the models created int eh size and density analysis

lt <- read.csv("data/survey_data/Annual_All_Species_Biomass_at_transect.csv", stringsAsFactors = F,na.strings ="-99999") %>% #LTER data density estimations collected from 50 transects across 11 sites between 2000-2018 in the Santa Barbara Channel.
  dplyr::select("YEAR", "MONTH", "SITE", "TRANSECT", "SP_CODE", "PERCENT_COVER", "DENSITY", "WM_GM2", "DRY_GM2", "SCIENTIFIC_NAME", "COMMON_NAME", "GROUP", "MOBILITY", "GROWTH_MORPH", "COARSE_GROUPING" ) %>%
  mutate(id = paste(SITE, TRANSECT, sep = "")) %>%
  filter(COMMON_NAME == "Purple Urchin" | COMMON_NAME == "Red Urchin" | COMMON_NAME == "Giant Kelp")%>% #filtering the data to only include giant kelp, purple urchins, and red urchins
  filter(SITE != "SCTW", SITE != "SCDI") #Filtering out island sites. This study focuses on costal sites. 

names(lt) <- tolower(names(lt))

sites <- read.csv(here("data/spatial", "lter_waypoints.csv")) 

purple.fun <- function(biomass){
  0.023482*biomass
} # this is the herbivory rate model prediction formed from the density analysis for purple urchins

red.fun <- function(biomass){
  0.009042*biomass
} # this is the herbivory rate model prediction formed from the density analysis for red urchins

#-----------------------------------------------------------------
## Make the map

library(rgdal)
library(rgeos)
library(sp)
library(raster)

mp <- lt %>% as_tibble() %>%
  mutate(predited.consumption = ifelse(sp_code == "SPL", purple.fun(wm_gm2), 
                                       ifelse(sp_code == "SFL", red.fun(wm_gm2), NA))) %>%
  left_join(sites) #Adding in the long lat data to the LTER density data


coordinates(mp) <- ~long + lat

proj4string(mp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"



cal <- readOGR(here("data/spatial", "caloutline.shp"))
all_mpas <- readOGR(here("data/spatial", "state_mpas.shp"))

d <- par(las = 1, mgp = c(3, 0.75, 0))
plot(cal, col = "#FFEB9B", xlim = c(-120.5,-119.5), ylim = c(34.35,34.55), xlab = "", ylab = "", cex.axis = 1.5, axes = TRUE) # cal state
plot(all_mpas, pch = 4, col = "#99ebff", add = T) 
plot(mp, pch = 21, cex = mp$predited.consumption*0.25, add = T, lwd = 1.5,  bg = alpha("#e31a1c", .5))

par(d)

# Make a plot for each year...
for(i in 2002:2018){
  myfile <- file.path("figures/maps", paste("year", "_", i, ".png"))
  png(myfile, width = 1000*3, height = 561*3, res = 300)
  d <- par(las = 1, mgp = c(3, 0.75, 0), mar = c(4,5,3,1))
  plot(cal, col = "#FFEB9B", xlim = c(-120.5,-119.5), ylim = c(34.35,34.55), xlab = "", ylab = "", cex.axis = 1.5, axes = TRUE) # cal state
  plot(all_mpas, pch = 4, col = "#99ebff", add = T) 
  plot(mp[mp$year == i,], pch = 21, cex = mp$predited.consumption[mp$year == i]*0.25, add = T, lwd = 1.5,  bg = alpha("#e31a1c", .25))
  par(d)
  dev.off()
}


# Make a time averaged plot...

av <- mp@data %>% 
  mutate(group = cut(year, breaks = 4, labels = FALSE)) %>%
  group_by(site, transect, group) %>%
  summarize(predicted.consumption = mean(predited.consumption, na.rm = T)) %>%
  left_join(sites)

coordinates(av) <- ~long + lat

proj4string(av) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

#c(bottom, left, top, right)
id <- c("2000-2004", "2005-2009", "2010-2013", "2014-2018")

png(here("figures", "map.png"), width = 1000*1.2, 600*1.2)
d <- par(mfrow = c(2, 2), las = 1, mgp = c(3, 0.75, 0), mar = c(3,5, 1.5, 0.5))
for(i in 1:4){
  plot(cal, col = "#FFEB9B", xlim = c(-120.5,-119.5), ylim = c(34.35,34.55), xlab = "", ylab = "", cex.axis = 1.5, axes = TRUE, main = paste(id[i])) # cal state
  plot(all_mpas, pch = 4, col = "#99ebff", add = T) 
  plot(av[av$group == i, ], pch = 21, cex = av$predicted.consumption[av$group == i], add = T, lwd = 1.5,  bg = alpha("#e31a1c", .5))
}
par(d)
dev.off()
