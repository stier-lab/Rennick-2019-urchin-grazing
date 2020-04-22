#-----------------------------------------------------------------
## Setup and data cleaning
#-----------------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lmer)
library(lmerTest)
source(here("analysis", "Functions.R"))

lt <- read.csv("data/survey_data/Annual_All_Species_Biomass_at_transect.csv", stringsAsFactors = F,na.strings ="-99999") %>%
  dplyr::select("YEAR", "MONTH", "SITE", "TRANSECT", "SP_CODE", "PERCENT_COVER", "DENSITY", "WM_GM2", "DRY_GM2", "SCIENTIFIC_NAME", "COMMON_NAME", "GROUP", "MOBILITY", "GROWTH_MORPH", "COARSE_GROUPING" ) %>%
  mutate(id = paste(SITE, TRANSECT, sep = "")) %>%
  filter(COMMON_NAME == "Purple Urchin" | COMMON_NAME == "Red Urchin" | COMMON_NAME == "Giant Kelp")%>%
  filter(SITE != "SCTW", SITE != "SCDI")

names(lt) <- tolower(names(lt))

sites <- read.csv(here("data/spatial", "lter_waypoints.csv")) 

purple.fun <- function(biomass){
  0.0009784*biomass
}

red.fun <- function(biomass){
  0.0003767*biomass
}

mp <- lt %>% as_tibble() %>%
  mutate(predited.consumption = ifelse(sp_code == "SPL", purple.fun(wm_gm2), 
                                            ifelse(sp_code == "SFL", red.fun(wm_gm2), NA))) %>%
  left_join(sites)






#-----------------------------------------------------------------
## Make the map

library(rgdal)
library(rgeos)
library(sp)
library(raster)

coordinates(mp) <- ~long + lat

proj4string(mp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"



cal <- readOGR(here("data/spatial", "caloutline.shp"))
all_mpas <- readOGR(here("data/spatial", "state_mpas.shp"))

d <- par(las = 1, mgp = c(3, 0.75, 0))
plot(cal, col = "#FFEB9B", xlim = c(-120.5,-119.5), ylim = c(34.35,34.55), xlab = "", ylab = "", cex.axis = 1.5, axes = TRUE) # cal state
plot(all_mpas, pch = 4, col = "#99ebff", add = T) 
plot(mp, pch = 21, cex = mp$predited.consumption*5, add = T, lwd = 1.5,  bg = alpha("#e31a1c", .5))

par(d)

# Make a plot for each year...
for(i in 2002:2018){
  myfile <- file.path("figures/", paste("year", "_", i, ".png"))
  png(myfile, width = 1000*3, height = 561*3, res = 300)
  d <- par(las = 1, mgp = c(3, 0.75, 0), mar = c(4,5,3,1))
  plot(cal, col = "#FFEB9B", xlim = c(-120.5,-119.5), ylim = c(34.35,34.55), xlab = "", ylab = "", cex.axis = 1.5, axes = TRUE) # cal state
  plot(all_mpas, pch = 4, col = "#99ebff", add = T) 
  plot(mp[mp$year == i,], pch = 21, cex = mp$predited.consumption[mp$year == i]*5, add = T, lwd = 1.5,  bg = alpha("#e31a1c", .25))
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
plot(av[av$group == i, ], pch = 21, cex = av$predicted.consumption[av$group == i]*10, add = T, lwd = 1.5,  bg = alpha("#e31a1c", .5))
}
par(d)
dev.off()

# time series

p1 <- mp@data %>% 
  filter(sp_code == "MAPY") %>%
  group_by(year, site) %>%
  mutate(biomass = mean(wm_gm2, na.rm = T))%>%
  ggplot(aes(x = year, y = biomass))+
  geom_line(aes(color = site))+
  ggpubr::theme_pubr(legend = "right")

p3 <- mp@data %>% group_by(year, site, sp_code) %>%
  summarize(predicted.consumption = mean(predited.consumption)) %>%
  ungroup() %>%
  filter(sp_code != "MAPY") %>%
  group_by(year, site) %>%
  mutate(predicted.consumption = sum(predicted.consumption, na.rm = T)*24)%>% # convert to per day rather than per hour
  ggplot(aes(x = year, y = predicted.consumption))+
    geom_line(aes(color = site))+
    ggpubr::theme_pubr(legend = "right")

fig4 <- cowplot::plot_grid(p1, p3, align = "v", nrow = 2)

ggsave(here("figures", "timeseries.png"), fig4)


# Ok so idea: Part of what is complicated about urchin foraging is that they can be both detritovors and herbivores. if production of kelp detritus exceeds consumption rate then we wouldn't expect any change in standing stock of kelp biomass. BUT if consumption rate exceeds detritus production then we might see shifts in the standing stock of kelp. Based on Christie's paper (and some analysis of LTER data) we should be able to estimate summer production of kelp detritus (and spectifically the detritus that lands on the seafloor). We can then examine trends in the time series to look for periods when urchin consumption rate exceeds kelp detritus production to look for declines in kelp biomass.

# The NPP survey data includes an estimate of average % biomass lost on a NPP transect for each monthly survey as plants, fronds, exudates, cut tissue (i.e. prop cuts), and blade scenescence. I'm going to calculate a regional summer time average from the NPP data collected at MOHK, AQUE, and ABUR to estimate the fraction of kelp biomss lost as fronds and blades (mass per day). I will then use this estimate to estimate detrital production rates along each of the LTER core transects. 

# get NPP data


npp <- read.csv("data/survey_data/NPP_ALL_YEARS_seasonal_with_MC_stderr.csv", stringsAsFactors = F,na.strings ="-99999") %>%
  rename_all(tolower) %>%
  dplyr::select(year, season, site, plnt_dns, frnd_dns, dry_wt,f, p, d, c, b, l) %>%
  filter(season == "3-Summer") %>%
  mutate(det.r = f + b)

av <- npp %>% summarize(ave = mean(det.r, na.rm = T), 
                        sd = sd(det.r, na.rm = T)) # so approximately 2% of total biomass is lost as fronds and blades per day!
av <- av[1,1]
sd <- av[1,2]

lt <- lt %>% as_tibble() %>%
  dplyr::select(year, month, site, transect, sp_code, wm_gm2) %>%
  spread(sp_code, wm_gm2) %>%
  mutate(SFL.pred = red.fun(SFL),
         SPL.pred = purple.fun(SPL), 
         predicted.consumption = (SPL.pred + SFL.pred)*24, 
         detritus.production = av * MAPY)




ggplot(lt, aes(x = predicted.consumption, y = MAPY))+
  geom_point(aes(size = detritus.production), pch = 21)

lt <- lt %>% 
  drop_na(MAPY) %>%
  mutate(dummy = as.factor(ifelse(detritus.production >= predicted.consumption, "Detritus >= Consumption", "Detritus < Consumption")), 
         difference = detritus.production - predicted.consumption, 
         urc = SPL + SFL)


# comparison of kelp biomass across all sites/years with kelp biomass at sites with urchin biomass in the 90% percentile

1-(mean(lt$MAPY[lt$urc >= quantile(lt$urc, probs = .9)]) / mean(lt$MAPY[lt$urc < quantile(lt$urc, probs = .9)]))




lmer <- lmer(log(MAPY+1) ~ dummy + (1|site) + (1|year), data = lt)
summary(lmer)
library(car)
modelassump(lmer)


aov <- lm(log(MAPY +1) ~ dummy, data = lt)
summary(aov)

anova(lmer, aov)


sum <- lt %>% group_by(dummy) %>% summarize(mean = mean(MAPY), sd = sd(MAPY))

#percent difference for paper
(sum[1,2] - sum[2,2]) / sum[1,2]


ggplot(lt, aes(x = dummy, y = MAPY))+
  geom_jitter(pch = 21)+
  geom_boxplot(outlier.shape = NA, alpha = 0.75)+
  ggpubr::theme_pubclean()

sum <- lt %>% group_by(dummy) %>% summarize(mean = mean(predicted.consumption), sd = sd(predicted.consumption))
(sum[1,2] - sum[2,2]) / sum[1,2]

lmer2 <- lmer(log(predicted.consumption+1) ~ dummy + (1|site) + (1|year), data = lt)
summary(lmer2)
modelassump(lmer2)



glmer <- glmer((MAPY+1) ~ dummy + (1|site) + (1|year), data = lt, family = Gamma(link = "log"))
summary(glmer)

glmer <- glmer((MAPY) ~ dummy + (1|site) + (1|year), data = lt, family = gaussian(link = "inverse"))
summary(glmer)


glm <- glm(MAPY ~ dummy, data = lt, family = gaussian(link = log),  start = coef(lm(MAPY ~ dummy, data = lt)))




lmer2 <- lmer(MAPY ~ urc * dummy + (1|site) + (1|year), data = lt)
summary(lmer2)
newdat <- expand.grid(dummy = unique(lt$dummy), urc = seq(min(lt$urc), max(lt$urc), length.out = 1000))
newdat$y <- predict(lmer2, newdata = newdat, re.form = NA)

ggplot(lt, aes(x = urc, y = MAPY))+
  geom_jitter(aes(color = dummy), pch = 21)+
  geom_line(data = newdat, aes(x = urc, y = y, color = dummy)) +
  ggpubr::theme_pubclean()

ggplot()

# Summary stats

# Average and sd detrital production rate (proportion of kelp biomass lost per day)
av
sd

# percept of time point where consumption > detritus
summary(lt$dummy)[1] / (summary(lt$dummy)[1] + summary(lt$dummy)[2])



lt.sum <- lt %>%
  group_by(year, site) %>%
  summarize(kelp = mean(MAPY, na.rm = T), 
            spl = mean(SPL, na.rm = T), 
            sfl = mean(SFL, na.rm = T), 
            predicted.consumption = red.fun(spl) + purple.fun(sfl), 
            detritus.production = kelp * av, 
            difference = detritus.production - predicted.consumption, 
            prop = ((detritus.production+10^-1) - (predicted.consumption+10^-1)) / (detritus.production + 10^-1) )
  

p1 <- ggplot(lt.sum, aes(x = year, y = kelp))+
  geom_line(aes(color = site))+
  # geom_point(aes(x = date, y = MAPY), color = ifelse(lte$diff > 0, "gray", "red")) +
  labs(y = expression(paste("Giant kelp biomass (g m"^"-2"*")")), x = "")+
  ggpubr::theme_pubr()
p2 <- ggplot(lt.sum, aes(x = year, y = prop))+
  geom_line(aes(color = site))+
  labs(y = expression(paste("Kelp detritus biomass (g m"^"-2"*")")), x = "")+
  geom_hline(yintercept = 0, lty = "dashed") %>%
  ggpubr::theme_pubr(legend = "none")
p3 <- ggplot(lte, aes(x = date, y = predicted.consumption))+
  geom_line(aes(color = site))+
  labs(y = expression(paste("Predicted kelp consumption (g m"^"-2"*"d"^"-1"*")")), x = "")+
  ggpubr::theme_pubr(legend = "none")
# p4 <- ggplot(lte, aes(x = date, y = diff))+
#   geom_point(aes(color = site))+
#   ggpubr::theme_pubr(legend = "none")
cowplot::plot_grid(p1, p2, nrow = 2)

fig4 <- cowplot::plot_grid(p1, p2, p3, nrow = 3, align = "v")

ggsave(here("figures", "timeseries.png"), fig4, width = 7.5, height = 10)


# Ok so when I average the data to sites, very rarely does consumption > detritus. Therefore, I conclude that if there is an effect of urchin consumption it is going to be scale dependent, and likely occure at the transect level rather than the site level. There are 5 sites that had transects which experienced consumption > detritus. Lets plot the transect level kelp biomass with point color coded by the dummy variables for these sites. 

aoi <- unique(lt$site[lt$dummy == "Detritus < Consumption"])

labs <- c("Detritus < Consumption\nUrchin biomass > 667", 
          "Detritus < Consumption", 
          "Detritus >= Consumption")


temp <- lt

length(lt$site[lt$dummy == "Detritus < Consumption" & lt$urc >= 668])/
length(lt$site[lt$dummy == "Detritus < Consumption"])
length(lt$site[lt$dummy == "Detritus >= Consumption"])

lt$dummy2 <- ifelse(lt$dummy == "Detritus < Consumption" & lt$urc >= 668, "Urchin biomass >= 668", "Urchin biomass < 668")

s3 <- lt %>%
  filter(site %in% c("CARP", "NAPL", "IVEE")) %>%
  ggplot(aes(x = year, y = MAPY))+
  geom_line(aes(color = as.factor(transect)), show.legend = FALSE)+
  geom_point(aes(fill = dummy, shape = dummy2), alpha = 0.75, show.legend = F)+
  scale_shape_manual(values = c(21,24))+
  scale_fill_manual(values = c("red", "gray"))+
  labs(fill = "", y = expression(paste("Giant kelp biomass (g m"^"-2"*")")), x = "")+ 
  facet_wrap(~site, scales = "free_y")+
  ggpubr::theme_pubclean()

ggsave(here("figures", "kelpxyear_barrens.png"), s3, device = "png", width = 10, height = 3.33)

lt %>% 
  group_by(dummy2) %>%
  summarize(count = n())
##############################################################################################
## Scrap
##############################################################################################


# bring in the detritus data

dt <- read.csv(here("data/survey_data/", "LTE_Detritus_Biomass_All_Years.csv"), stringsAsFactors = F,na.strings ="-99999") %>%
  rename_all(tolower) %>%
  dplyr::select(year, month, date, site, transect, treatment, quad, side, sp_code, wet_wt, area) %>%
  filter(sp_code %in% c("MAPY", "MPJ", "MPS"), treatment == "CONTROL") %>%
  group_by(year, month, date, site, transect) %>%
  summarize(detritus = sum(wet_wt)/6) %>% #sampled in 6, 1 m2 quadrats
  ungroup() %>%
  mutate(date = lubridate::ymd(date))


ggplot(dt, aes(x = date, y = detritus))+
  geom_line(aes(color = site))


#How much detritus is produced per unit of kelp biomass per unit time? Where detritus is specifically fronds and blades, not whole plants. ???? Can I use Rodriquez et al 2013 (Ecology)?

#total biomass lost as fronds was four times higher than that lost as plants when averaged over all sites and months (5.6 +- 0.57 vs. 1.3 +- 0.23 g dry mass􏳘m􏳗2􏳘d􏳗1 for fronds and plants, respectively, mean 6 SE, Wilcoxon signed rank test, W 1⁄4 8692, P , 0.0001). 

# 5.6 +- 0.57 g dry mass per m2 per day

# dry mass to wet mass DM:WM = 0.095 +- 0.001 s.e from Yorke et al. 2019... from LTER data. Immediate source unknown.


5.6/0.095
# 58.94737 g wet mass per m2 per day...

# Ok so what if we looked for these patterns only at the LTE control transects. I don't have an estimate of detritus production (i.e. a rate) but we have data for urchin biomass, kelp biomass, and detritus biomass, sampled along the same 40 m transects at MOHK, CARP, IVEE, NAPLE, and AQUE... 

lte <- read.csv("data/survey_data/LTE_All_Species_Biomass_at_transect.csv", stringsAsFactors = F,na.strings ="-99999") %>%
  rename_all(tolower) %>%
  dplyr::select(year, month, date, site, transect, treatment, sp_code, wm_gm2) %>%
  filter(treatment == "CONTROL", sp_code %in% c("MAPY", "MPJ", "MPS", "SPL", "SFL")) %>%
  mutate(date = stringr::str_replace_all(date, "/", "-"),
         date = lubridate::mdy(date)) %>%
  spread(sp_code, wm_gm2) %>%
  mutate(SFL.pred = red.fun(SFL),
         SPL.pred = purple.fun(SPL), 
         predicted.consumption = (SPL.pred + SFL.pred)*24) %>%
  left_join(dt) %>%
  mutate(diff = detritus - predicted.consumption, 
         diff.prop = detritus/predicted.consumption)

p1 <- ggplot(lte, aes(x = date, y = MAPY))+
  geom_line(aes(color = site))+
  geom_point(aes(x = date, y = MAPY), color = ifelse(lte$diff > 0, "gray", "red")) +
  labs(y = expression(paste("Giant kelp biomass (g m"^"-2"*")")), x = "")+
  ggpubr::theme_pubr()
p2 <- ggplot(lte, aes(x = date, y = detritus))+
  geom_line(aes(color = site))+
  labs(y = expression(paste("Kelp detritus biomass (g m"^"-2"*")")), x = "")+
  ggpubr::theme_pubr(legend = "none")
p3 <- ggplot(lte, aes(x = date, y = predicted.consumption))+
  geom_line(aes(color = site))+
  labs(y = expression(paste("Predicted kelp consumption (g m"^"-2"*"d"^"-1"*")")), x = "")+
  ggpubr::theme_pubr(legend = "none")
# p4 <- ggplot(lte, aes(x = date, y = diff))+
#   geom_point(aes(color = site))+
#   ggpubr::theme_pubr(legend = "none")

fig4 <- cowplot::plot_grid(p1, p2, p3, nrow = 3, align = "v")

ggsave(here("figures", "timeseries.png"), fig4, width = 7.5, height = 10)

# ok so if the consumption exceeds detritus (assuming the biomass of detritus represents a constant supply rate) then the difference will be negative. I hypothesize that time points with negative difference in detritus - consumption will be correlated with lower kelp biomass. 

ggplot(lte, aes( x = date, y = MAPY))+
  geom_line()+
  geom_point(aes(x = date, y = MAPY), color = ifelse(lte$diff > 0, "white", "red")) +
  facet_wrap(~site)

lte %>% mutate(diff = as.factor(ifelse(diff > 0, 1, -1))) %>%
  drop_na(diff) %>%


ggplot(lte, aes(x = diff, y = MAPY))+
  geom_point(aes(color = site))


m <- lte %>% mutate(diff = as.factor(ifelse(diff > 0, 1, -1))) %>%
  drop_na(diff)

library(lme4)
library(lmerTest)
lm1 <- lmer(MAPY ~ diff + (1|site), lte)
summary(lm1)

ggplot(m)+
  geom_boxplot(aes(x = diff, y = MAPY))

ggplot(m)+
  geom_point(aes(x = detritus, y = MAPY))

# Ok so basically when consumption exceeds detritus, kelp biomass tends to be lower. This could be due to the fact that when kelp biomass is low detritus is low, but this relationship isn't that strong. So there is some evidence that urchins may be decreasing kelp biomass when detritus is limiting.





# by site

id <- unique(lte$site)
for(i in 1:5){
  
  myfile <- file.path("figures/", paste("site", "_", id[i], ".png"))
p1 <- ggplot(lte[lte$site == id[i], ], aes(x = date, y = MAPY))+
  geom_line()+
  geom_point(aes(x = date, y = MAPY), color = ifelse(lte$diff[lte$site == id[i]] > 0, "gray", "red")) +
  ggpubr::theme_pubr()
p2 <- ggplot(lte[lte$site == id[i], ], aes(x = date, y = detritus))+
  geom_line()+
  ggpubr::theme_pubr(legend = "none")
p3 <- ggplot(lte[lte$site == id[i], ], aes(x = date, y = predicted.consumption))+
  geom_line()+
  ggpubr::theme_pubr(legend = "none")

out <- cowplot::plot_grid(p1, p2, p3, nrow = 3, align = "v")
ggsave(myfile, out)


}


id <- unique(lte$site)
for(i in 1:5){
  
  myfile <- file.path("figures/", paste("site_v2", "_", id[i], ".png"))
  p1 <- ggplot(lte[lte$site == id[i], ], aes(x = date, y = MAPY))+
    geom_line()+
    geom_point(aes(x = date, y = MAPY), color = ifelse(lte$diff[lte$site == id[i]] > 0, "gray", "red")) +
    ggpubr::theme_pubr()
  p2 <- ggplot(lte[lte$site == id[i], ], aes(x = date, y = diff.prop))+
    geom_line()+
    ggpubr::theme_pubr(legend = "none")
  
  out <- cowplot::plot_grid(p1, p2, nrow = 2, align = "v")
  ggsave(myfile, out)
  
  
}


lmer <- lmer(MAPY ~ predicted.consumption * detritus + (1|site) + (1|year), data = lte)
summary(lmer)

lmer2 <- lmer(MAPY ~ lag(scale(MAPY)) + scale(predicted.consumption) * scale(detritus) + (1|site) + (1|year), data = lte)
summary(lmer2)

lmer4 <- lmer(MAPY ~ log(predicted.consumption+1) * log(detritus+1) + (1|site) + (1|year), data = lte)
summary(lmer4)

newdat <- expand.grid(detritus = seq(min(lte$detritus, na.rm = T), max(lte$detritus, na.rm = T), length.out = 100), predicted.consumption = seq(max(lte$predicted.consumption, na.rm = T), min(lte$predicted.consumption, na.rm = T), length.out = 100))

newdat$y <- predict(lmer, newdata = newdat, re.form = ~0)

ggplot(newdat, aes(x = predicted.consumption, y = detritus))+
  geom_raster(aes(fill = y))+
  geom_contour(aes(z = y))+
  geom_point(data = lte, aes(x = predicted.consumption, y = detritus, size = MAPY))


newdat <- expand.grid(detritus = c(1, 10, mean(lte$detritus, na.rm = T), mean(lte$detritus, na.rm = T) + sd(lte$detritus, na.rm = T)redicted.consumption = seq(max(lte$predicted.consumption, na.rm = T), min(lte$predicted.consumption, na.rm = T), length.out = 100))
newdat$y <- predict(lmer, newdata = newdat, re.form = ~0)

ggplot(lte, aes(x = predicted.consumption, y = MAPY))+
geom_point(aes(size = MAPYdetritusch = 21)+
geom_line(data = newdat, aes(x = predicted.consumption, y = y, color = as.factor(detritus)))
                                   
                                   
                                   
# univariate

lmer3 <- lmer(MAPY ~ lag(MAPY) + (1|site) + (1|year), data = lte)
summary(lmer3)

lmer5 <- lmer(MAPY ~ log(detritus+1) + (1|site) + (1|year), data = lte)
summary(lmer5)

lmer6 <- lmer(MAPY ~ predicted.consumption + (1|site) + (1|year), data = lte)
summary(lmer6)



# ok so I would predict that if we see an effect we would see an effect on the change in kelp biomass from the summer to the fall. In other words, if consumption exceeds detritus than urchins may attack live biomass. Therefore, the change in live kelp biomass should be negative when the difference between detritus and consumption is negative. 

lte <- read.csv("data/survey_data/LTE_All_Species_Biomass_at_transect.csv", stringsAsFactors = F,na.strings ="-99999") %>%
  as_tibble() %>%
  rename_all(tolower) %>%
  dplyr::select(year, month, date, site, transect, treatment, sp_code, wm_gm2) %>%
  filter(treatment == "CONTROL", sp_code %in% c("MAPY", "MPJ", "MPS", "SPL", "SFL")) %>%
  mutate(date = stringr::str_replace_all(date, "/", "-"),
         date = lubridate::mdy(date)) %>%
  spread(sp_code, wm_gm2) %>%
  left_join(dt) %>%
  mutate(quarter = ifelse(month >= 1 & month <= 3, "winter", 
                          ifelse(month >=4 & month <= 6, "spring", 
                                 ifelse(month >= 7 & month <= 9, "a.summer", 
                                        ifelse(month >= 10 & month <= 12, "b.fall", NA))))) %>%
  group_by(year, quarter, site) %>%
  summarize(kelp = mean(MAPY, na.rm = T), 
            sfl = mean(SFL, na.rm = T), 
            spl = mean(SPL, na.rm = T), 
            sfl.pred = red.fun(sfl), 
            spl.pred = purple.fun(spl), 
            predicted.consumption = (spl.pred + sfl.pred)*24, 
            detritus = mean(detritus, na.rm = T)) %>%
  filter(quarter %in% c("a.summer", "b.fall")) %>%
  arrange(site, year) %>%
  group_by(site, year) %>%
  mutate(delta.kelp = lead(kelp) - kelp, 
         difference = detritus - predicted.consumption, 
         diff.f = as.factor(ifelse(difference > 0, "detritus > consumption", "detritus < consumption"))) %>%
  filter(quarter == "a.summer")

ggplot(lte, aes(x = difference, y = delta.kelp))+
  geom_point(aes(color = site))

ggplot(lte, aes(x = diff.f, y = delta.kelp))+
  geom_boxplot()

lmer <- lmer(kelp ~ predicted.consumption * detritus + (1|site) + (1|year), data = lte)

summary(lmer)

lmer <- lmer(kelp ~ predicted.consumption + (1|site) + (1|year), data = lte)
summary(lmer)










