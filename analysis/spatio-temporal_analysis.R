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
library(car) 
source(here("analysis", "Functions.R"))
#source(here("analysis", "Functions.R")) linked to exp analysis so we can extract model predicitons

#We are going to apply the models extracted from the size and density analysis to 9 LTER coastal sites. We are going to track urchin density through time and apply herbivory pressure estiamtes to each site based on our findings. We additionally are incoorperating observational data realting to detriatal supply to test its bearing on the amount of kelp biomass lost in a given area. 

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


# -----------------------------------------------------
## Time series analysis
# -----------------------------------------------------

# time series

p1 <- mp@data %>% #LTER data 
  filter(sp_code == "MAPY") %>% #filtering for macrocystis pyrifera 
  group_by(year, site) %>%
  mutate(biomass = mean(wm_gm2, na.rm = T))%>% #adding a biomass column
  ggplot(aes(x = year, y = biomass))+
  geom_line(aes(color = site))+
  labs(y = expression(paste("Giant kelp biomass (g m"^"-2"*")")), x = "")+
  ggpubr::theme_pubr(legend = "right")

p3 <- mp@data %>% group_by(year, site, sp_code) %>%
  summarize(predicted.consumption = mean(predited.consumption)) %>%
  ungroup() %>%
  filter(sp_code != "MAPY") %>%
  group_by(year, site) %>%
  mutate(predicted.consumption = sum(predicted.consumption, na.rm = T))%>% # convert to per day rather than per hour? Why didnt you add = T)*24)%>% like in the spatial analysis? 
  ggplot(aes(x = year, y = predicted.consumption))+ #average consumption across sites
  geom_line(aes(color = site))+
  labs(y = expression(paste("Predicted kelp consumption (g m"^"-2"*"d"^"-1"*")")), x = "")+
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

av1 <- npp %>% summarize(ave = mean(det.r, na.rm = T), 
                        sd = sd(det.r, na.rm = T)) # so approximately 2% of total biomass is lost as fronds and blades per day!
av <- av1[1,1]
sd <- av1[1,2]

lt <- lt %>% as_tibble() %>%
  dplyr::select(year, month, site, transect, sp_code, wm_gm2) %>%
  spread(sp_code, wm_gm2) %>%
  mutate(SFL.pred = red.fun(SFL),
         SPL.pred = purple.fun(SPL), 
         predicted.consumption = (SPL.pred + SFL.pred), 
         detritus.production = av * MAPY) %>% 
  drop_na(MAPY) %>%
  mutate(dummy = as.factor(ifelse(detritus.production >= predicted.consumption, "Detritus >= Consumption", "Detritus < Consumption")), 
         difference = detritus.production - predicted.consumption, 
         urc = SPL + SFL) #if its negative that means that the urhins are eating more than is available in detritus by that much g/hr? 

# comparison of kelp biomass across all sites/years with kelp biomass at sites with urchin biomass in the 90% percentile

1-(mean(lt$MAPY[lt$urc >= quantile(lt$urc, probs = .9)]) / mean(lt$MAPY[lt$urc < quantile(lt$urc, probs = .9)]))

# Model kelp biomass as a function of a factor variable that desribes if detritus or consumption is greater... 

lmer <- lmer(log(MAPY+1) ~ dummy + (1|site) + (1|year), data = lt)
summary(lmer)
modelassump(lmer)


aov <- lm(log(MAPY +1) ~ dummy, data = lt)
summary(aov)

anova(lmer, aov) # so the mixed effects model does better...

# Plot it up!

      # summary stats for paper
      sum <- lt %>% group_by(dummy) %>% summarize(mean = mean(MAPY), sd = sd(MAPY))
      
      #percent difference in mean kelp biomass between consumption > detritus and visa versa
      (sum[1,2] - sum[2,2]) / sum[1,2]

ggplot(lt[lt$predicted.consumption > 0, ], aes(x = dummy, y = MAPY))+
  geom_jitter(aes(color = dummy), pch = 21, show.legend = F)+
  scale_color_manual(values = c("#8f4811", "#35753d"))+
  geom_boxplot(outlier.shape = NA, alpha = 0.75)+
  labs(x = "", y = expression(paste("Giant kelp biomass (g m"^"-2"*")")))+
  ggpubr::theme_pubclean()

# is this just due to low detritus/low kelp?
sum <- lt %>% group_by(dummy) %>% summarize(mean = mean(predicted.consumption), sd = sd(predicted.consumption))
(sum[1,2] - sum[2,2]) / sum[1,2] # consumption rates are 77% higher when consumption exceeds detritus. There NO! 

# model that for completeness
lmer2 <- lmer(log(predicted.consumption+1) ~ dummy + (1|site) + (1|year), data = lt)
summary(lmer2)
modelassump(lmer2)
    
    lmer2.1 <- lmer(log(predicted.consumption) ~ dummy + (1|site) + (1|year), data = lt[lt$predicted.consumption > 0, ])
    summary(lmer2.1)
    modelassump(lmer2.1)
    
    lmer2.2 <- glmer(predicted.consumption ~ dummy + (1|site) + (1|year), data = lt[lt$predicted.consumption > 0, ], family = Gamma(link = log))
    summary(lmer2.2)
    modelassump(lmer2.2)

# ok so we know that kelp biomass tends to be less when detritus < consumption and that this is likely not just due to low kelp biomass because consumption is significantly higher when detritus < consumption than when detritus >= consumption. What then is the relationship between kelp biomass and urchin biomass? 

lmer3 <- lmer(MAPY ~ predicted.consumption + (1|site) + (1|year), data = lt)
summary(lmer3)
modelassump(lmer3) # this model has considerable issues and describes the data poorly.. this needs to be improved.

# # 
# lmer3 <- lmer(MAPY ~ predicted.consumption + (1|site) + (1|year), data = lt[lt$predicted.consumption> 0, ])
# summary(lmer3)
# modelassump(lmer3)
# # 
# # 
# #   # refine model so that it fits the data more accuarately... TO BE CONTINUED
#       glmer3 <- glmer((MAPY+0.001) ~ log(predicted.consumption) + (1|site) + (1|year), data = lt[lt$predicted.consumption > 0, ], family = gaussian(link = "log"))
#       summary(glmer3)
#       modelassump(glmer3)

newdat <- data.frame(predicted.consumption = seq(min(lt$predicted.consumption), max(lt$predicted.consumption), length.out = 1000))
newdat$y <- predict(lmer3, newdata = newdat, re.form = NA)
newdat <- newdat[newdat$y > 0, ]
# newdat$y2 <- predict(glmer3, newdata = newdat, re.form = NA, type = "response")

#---------------------------------------------
## Figure 5
#---------------------------------------------

# 2 panel plot with (a) the relationship between kelp biomass and urchin biomass with points colored by if consumption > detritus, and fit with the model between kelp biomass and urchin biomass. Panel b is the boxplot showing that kelp biomass is less when detritus < consumption. 


p1 <- ggplot(lt, aes(x = predicted.consumption, y = MAPY))+
  geom_jitter(aes(color = dummy), pch = 21)+
  scale_color_manual(values = c("#8f4811", "#35753d"))+
  geom_line(data = newdat, aes(x = predicted.consumption , y = y)) +
  labs(x = expression(paste("Consumption rate (g m"^"-2"*"d"^"-1"*")")), y = expression(paste("Giant kelp biomass (g m"^"-2"*")")), color = "")+
  ggpubr::theme_pubclean()

p2 <- ggplot(lt, aes(x = dummy, y = MAPY))+
  geom_jitter(aes(color = dummy), pch = 21)+
  scale_color_manual(values = c("#8f4811", "#35753d"))+
  geom_boxplot(outlier.shape = NA, alpha = 0.75)+
  labs(x = "", y = "", color = "")+
  ggpubr::theme_pubclean()+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

fig5 <- cowplot::plot_grid(p1, p2, align = "h", rel_widths = c(1, 0.65) )
ggsave(here("figures", "kelpxurc.png"), fig5, device = "png", width = 10, height = 5)



#-----------------------------------------------------
## Supplemental figure 3
#-----------------------------------------------------

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

#-----------------------------------------------------
## Summary statistics
#-----------------------------------------------------

mean(lt$predicted.consumption, na.rm =T)
sd(lt$predicted.consumption, na.rm =T)
range(lt$predicted.consumption)

#temporal vs. spatial variation

lt %>% group_by(site) %>% summarize(cv.space = sd(predicted.consumption)/mean(predicted.consumption)) %>% summarize(mean(cv.space), sd(cv.space))
lt %>% group_by(year) %>% summarize(cv.time = sd(predicted.consumption)/mean(predicted.consumption)) %>% summarize(mean(cv.time), sd(cv.time))

