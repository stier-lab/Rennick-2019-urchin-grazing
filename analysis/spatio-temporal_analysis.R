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


lt <- read.csv("data/survey_data/Annual_All_Species_Biomass_at_transect.csv", stringsAsFactors = F,na.strings ="-99999") %>% #LTER data density estimations collected from 50 transects across 11 sites between 2000-2018 in the Santa Barbara Channel.
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



## Time series analysis
# -----------------------------------------------------

# time series

p1 <- lt %>% #LTER data 
  group_by(year, site) %>%
  summarize(biomass = mean(MAPY, na.rm = T))%>% #adding a biomass column
  ggplot(aes(x = year, y = biomass))+
  geom_line(aes(color = site))+
  labs(y = expression(paste("Giant kelp biomass (g m"^"-2"*")")), x = "")+
  ggpubr::theme_pubr(legend = "right")

p3 <- lt %>%
  group_by(year,site) %>%
  summarize(predicted.consumption = mean(predicted.consumption, na.rm = T)) %>%
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

lt <- lt %>%
  drop_na(MAPY) %>%
  mutate(detritus.production = av * MAPY, 
         dummy = as.factor(ifelse(detritus.production >= predicted.consumption, "Detritus >= Consumption", "Detritus < Consumption")), 
         urc = SPL + SFL)


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

# UPDATED 10/1/2020 BD: OK so our hypothesis is that urchin and living kelp biomass are decoupled under conditions when detrital supply is high. Therefore, we predict that the balance of consumption rate and detrial supply will explain kelp biomass dynamics better than urchin biomass alone. To test this we will build a model of kelp biomass ~ urchin biomass and compare that to a model of kelp biomass ~ urchin biomass +/* dummy factor

lt <- lt %>% mutate(urc.biomass = SFL + SPL)

lmer4 <- lmer(MAPY ~ urc.biomass * dummy + (1|site) + (1|year), data = lt)
summary(lmer4)
modelassump(lmer4)

lmer4.1 <- lmer(MAPY ~ urc.biomass + dummy + (1|site) + (1|year), data = lt)
summary(lmer4.1)
modelassump(lmer4.1)
    
    #The distribution of the response is relatively tricky. Zero-inflated to some extent... these are some models attempting to refine the fit.
    glm4.1.a <- glm(I(MAPY+1) ~ urc.biomass + dummy, data = lt, family = Gamma(link = "log"))
    summary(lmer4.1.a)
    modelassump(lmer4.1.a)
    
    lmer4.1.b <- lmer(log(MAPY+1) ~ urc.biomass + dummy + (1|site) + (1|year), data = lt)
    summary(lmer4.1.b)
    modelassump(lmer4.1.b)
    
    glmer4.1.c <- glm(MAPY ~ scale(urc.biomass) + dummy , family = gaussian(link = "log"), data = lt)
    summary(glmer4.1.c)
    modelassump(glmer4.1.c)
    

lmer4.2 <- lmer(MAPY ~ urc.biomass + (1|site) + (1|year), data = lt)
summary(lmer4.2)
modelassump(lmer4.2)

AIC(lmer4,lmer4.1, lmer4.2)
anova(lmer4,lmer4.1)
anova(lmer4.1, lmer4.2)

newdat <- data.frame(predicted.consumption = seq(min(lt$predicted.consumption), max(lt$predicted.consumption), length.out = 1000))
newdat$y <- predict(lmer3, newdata = newdat, re.form = NA)
newdat <- newdat[newdat$y > 0, ]
# newdat$y2 <- predict(glmer3, newdata = newdat, re.form = NA, type = "response")

#---------------------------------------------
## Figure 5
#---------------------------------------------

# 2 panel plot with (a) the relationship between kelp biomass and urchin biomass with points colored by if consumption > detritus, and fit with the model between kelp biomass and urchin biomass. Panel b is the boxplot showing that kelp biomass is less when detritus < consumption. 


p1 <- ggplot(lt, aes(x = SPL + SFL, y = MAPY))+
  geom_jitter(aes(color = dummy), pch = 21)+
  scale_color_manual(values = c("#8f4811", "#35753d"))+
  # geom_line(data = newdat, aes(x = predicted.consumption , y = y)) +
  labs(x = expression(paste("Combined urchin biomass (g m"^"-2"*")")), y = expression(paste("Giant kelp biomass (g m"^"-2"*")")), color = "")+
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

#------------------------------------------------------
## Ratio dependent plots 
#------------------------------------------------------

lt2 <- lt %>%
  mutate(ratio = predicted.consumption / detritus.production, 
         diff = predicted.consumption - detritus.production, 
         per.diff = (predicted.consumption - detritus.production) / (predicted.consumption + detritus.production)) %>%
  arrange(site, transect, year) %>%
  ungroup() %>%
  group_by(site, transect) %>%
  mutate(deltaK = MAPY - lag(MAPY, order_by = year))


lmer6 <- lmer(deltaK ~ per.diff + (1|site) + (1|year), lt2)
summary(lmer6)

newdat <- data.frame(per.diff = seq(min(lt2$per.diff, na.rm = T), max(lt2$per.diff, na.rm = T), length.out = 1000))
newdat$y <- predict(lmer6, newdata = newdat, re.form = NA)


mm<-model.matrix(~per.diff,data=newdat)
predFun<-function(.) mm%*%fixef(.)
bb<-bootMer(lmer6,FUN=predFun,nsim=200) #do this 200 timesCopy
#As we did this 200 times the 95% CI will be bordered by the 5th and 195th value:
  
bb_se<-apply(bb$t,2,function(x) x[order(x)][c(5,195)])
newdat$LC<-bb_se[1,]
newdat$UC<-bb_se[2,] 


lmer6.low <- lmer(deltaK ~ per.diff + (1|site) + (1|year), lt2[lt2$per.diff < 0, ])
summary(lmer6.low)

lmer6.high <- lmer(deltaK ~ per.diff + (1|site) + (1|year), lt2[lt2$per.diff >= 0, ])
summary(lmer6.high)

lt2 %>% filter(urc.biomass != 0) %>%
ggplot(aes(x = per.diff, y = deltaK))+
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu"))+
  geom_point(aes(size = urc.biomass, fill = deltaK), pch = 21)+
  geom_line(data = newdat, aes(x = per.diff, y = y))+
  geom_ribbon(data = newdat, aes(x = per.diff, ymin = LC, ymax = UC, y = y), alpha = .2) +
  # geom_line(data = newdat, aes(x = per.diff, y = LC), lty = 4)+
  # geom_line(data = newdat, aes(x = per.diff, y = UC), lty = 4)+
  geom_hline(yintercept = 0, lty = 3)+
  geom_vline(xintercept = 0, lty = 3)+
  labs(x = "Proportional difference between\nconsumption and detrial supply rate", 
       y = "delta Kelp biomass", 
       size = "Urchin biomass\ndensity")+
  ggpubr::theme_pubr(legend = "right")
                    


#------------------------------------------------------
## Lag plots
#------------------------------------------------------

ggplot(lt, aes(x = urc.biomass, y = lead(MAPY, n = 1)))+
  geom_jitter(aes(color = dummy), pch = 21)+
  scale_color_manual(values = c("#8f4811", "#35753d"))+
  # geom_line(data = newdat, aes(x = predicted.consumption , y = y)) +
  labs(x = expression(paste("Combined urchin biomass (g m"^"-2"*")")), y = expression(paste("Giant kelp biomass (g m"^"-2"*")")), color = "")+
  ggpubr::theme_pubclean()



