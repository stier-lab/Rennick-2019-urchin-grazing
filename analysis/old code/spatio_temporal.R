#-----------------------------------------------------------------
## Setup and data cleaning
#-----------------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
source(here("analysis/", "Functions.R"))
library(car) 
source(here("analysis/", "density_analysis.R"))

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
         urc.biomass = SPL + SFL) %>%
  filter(urc.biomass > 0)

ggplot(lt[lt$year > 2018, ], aes(x = site, y = MAPY))+
  geom_point(aes(color = site))

#-------------------------------------------------------
## Summary stats for paper
#-------------------------------------------------------

# comparison of kelp biomass across all sites/years with kelp biomass at sites with urchin biomass in the 90% percentile

1-(mean(lt$MAPY[lt$urc >= quantile(lt$urc, probs = .9)]) / mean(lt$MAPY[lt$urc < quantile(lt$urc, probs = .9)]))


# means and sd of kelp
sum <- lt %>% group_by(dummy) %>% summarize(mean = mean(MAPY), sd = sd(MAPY))

#percent difference in mean kelp biomass between consumption > detritus and visa versa
(sum[1,2] - sum[2,2]) / sum[1,2]

# is this just due to low detritus/low kelp?
sum <- lt %>% group_by(dummy) %>% summarize(mean = mean(predicted.consumption), sd = sd(predicted.consumption))
(sum[1,2] - sum[2,2]) / sum[1,2] # consumption rates are 77% higher when consumption exceeds detritus. There NO! 


# UPDATED 10/1/2020 BD: OK so our hypothesis is that urchin and living kelp biomass are decoupled under conditions when detrital supply is high. Therefore, we predict that the balance of consumption rate and detrial supply will explain kelp biomass dynamics better than urchin biomass alone. To test this we will build a model of kelp biomass ~ urchin biomass and compare that to a model of kelp biomass ~ urchin biomass +/* dummy factor

lmer4 <- lmer(sqrt(MAPY) ~ urc.biomass * dummy + (1|site) + (1|year), data = lt)
summary(lmer4)
modelassump(lmer4)

lmer4.1 <- lmer(sqrt(MAPY) ~ urc.biomass + dummy + (1|site) + (1|year), data = lt)
summary(lmer4.1)
modelassump(lmer4.1)

#The distribution of the response is relatively tricky. Zero-inflated to some extent... these are some models attempting to refine the fit.

lt$MAPY1 <- lt$MAPY+1
glm4.1.a <- glm(MAPY1 ~ urc.biomass + dummy, data = lt, family = Gamma(link = "log"))
summary(glm4.1.a)
modelassump(glm4.1.a)

lmer4.1.b <- glmer(MAPY1 ~ scale(urc.biomass) + dummy + (1|site) + (1|year), data = lt, family = gaussian(link = "log"))
summary(lmer4.1.b)
modelassump(lmer4.1.b)

glmer4.1.c <- glmer(MAPY1 ~ scale(urc.biomass) * dummy + (1|site) + (1|year) , family = Gamma(link = "log"), data = lt[lt$year < 2019, ])
summary(glmer4.1.c)
modelassump(glmer4.1.c)


glmerTMB4 <- glmmTMB::glmmTMB(MAPY1 ~ scale(urc.biomass) * dummy + (1|site) + (1|year) , family = Gamma(link = "log"), data = lt)
summary(glmerTMB4)

temp <- DHARMa::simulateResiduals(glmerTMB4)
plot(temp)
DHARMa::testZeroInflation(temp)
DHARMa::testDispersion(temp)
DHARMa::testTemporalAutocorrelation(simulationOutput = temp, time = lt$year)
modelassump(glmerTMB4)

  glmerTMB4.sqrt <- glmmTMB::glmmTMB(MAPY1 ~ scale(urc.biomass) * dummy + (1|site) , family = Gamma(link = "sqrt"), data = lt)
  summary(glmerTMB4.sqrt)

anova(glmerTMB4, glmerTMB4.sqrt)

  glmerTMB4.1.d <- glmmTMB::glmmTMB(MAPY1 ~ scale(urc.biomass) + dummy + (1|site) + (1|year) , family = Gamma(link = "log"), data = lt)
  summary(glmerTMB4.1.d)
  modelassump(glmerTMB4.1.d)
  
  glmerTMB4.1.e <- glmmTMB::glmmTMB(MAPY1 ~ urc.biomass + (1|site) + (1|year) , family = Gamma(link = "log"), data = lt)
  summary(glmerTMB4.1.e)
  modelassump(glmerTMB4.1.e)
  
AIC(glmerTMB4, glmerTMB4.1.d, glmerTMB4.1.e)
anova(glmerTMB4, glmerTMB4.1.d)

  
  
  
  
  
lmer4.2 <- lmer(sqrt(MAPY) ~ urc.biomass + (1|site) + (1|year), data = lt)
summary(lmer4.2)
modelassump(lmer4.2)

AIC(lmer4,lmer4.1, lmer4.2)
AIC(glm4.1.a, lmer4.1.b, glmer4.1.c)
AIC(glmer4.1.c, glmer4.1.d, glmer4.1.e)
anova(lmer4,lmer4.1)
anova(lmer4.1, lmer4.2)
anova(glm4.1.a, glmer4.1.c)




glmer4.1.c <- glmer(MAPY1 ~ scale(urc.biomass) * dummy + (1|site) + (1|year) , family = Gamma(link = "log"), data = lt)
summary(glmer4.1.c)
modelassump(glmer4.1.c)


newdat <- ggeffects::ggpredict(glmerTMB4, terms = c("urc.biomass", "dummy"))
plot(newdat)


#---------------------------------------------
## Figure 5
#---------------------------------------------

# 2 panel plot with (a) the relationship between kelp biomass and urchin biomass with points colored by if consumption > detritus, and fit with the model between kelp biomass and urchin biomass. Panel b is the boxplot showing that kelp biomass is less when detritus < consumption. 


p1 <- ggplot(lt, aes(x = urc.biomass, y = MAPY))+
  geom_jitter(colour="white",aes(fill=dummy, shape = dummy), alpha = 0.5, size = 2)+
  scale_shape_manual(values = c(24,21)) +
  scale_fill_manual(values = c("#272593", "#35753d"))+
  # geom_line(data = newdat, aes(x = predicted.consumption , y = y)) +
  labs(x = expression(paste("Combined urchin biomass (g m"^"-2"*")")), y = expression(paste("Giant kelp biomass (g m"^"-2"*")")), color = "")+
  # geom_line(data = newdat, aes(x = x, y = predicted, color = group))+
  # geom_ribbon(data = newdat, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, group = group), alpha = .1) +
  theme_classic()+
  theme(legend.position = c(0.75,0.9), 
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank())+
  coord_cartesian(ylim = c(0, 25000))

p2 <- ggplot(lt, aes(x = dummy, y = MAPY))+
  geom_jitter(colour="white",aes(fill=dummy, shape = dummy, alpha=0.5), size = 2)+
  scale_shape_manual(values = c(24,21)) +
  scale_fill_manual(values = c("#272593", "#35753d"))+
  geom_boxplot(outlier.shape = NA, alpha = 0.75,aes(fill=dummy))+
  labs(x = "", y = "", color = "")+
  scale_x_discrete(labels = c("Detritus <\n Consumption", "Detritus >=\n Consumption"))+
  scale_y_continuous(breaks = seq(0, 25000, by = 5000))+
  theme_classic()+
  theme(axis.text.y = element_blank())+
  coord_cartesian(ylim = c(0, 25000)) +
  theme(legend.position = "none")

fig5 <- cowplot::plot_grid(p1, p2, align = "h", rel_widths = c(1, 1), labels = "AUTO")
# fig5

ggsave(here("figures", "kelpxurc.png"), fig5, device = "png", width = 10, height = 5)
ggsave(here("figures", "kelpxurc.pdf"), fig5, device = "pdf", width = 8, height = 4, useDingbats = FALSE)




#-----------------------------------------------------
## Summary statistics
#-----------------------------------------------------

mean(lt$predicted.consumption, na.rm =T)
sd(lt$predicted.consumption, na.rm =T)
range(lt$predicted.consumption)

#temporal vs. spatial variation

lt %>% group_by(site) %>% 
  summarize(cv.space = sd(predicted.consumption)/mean(predicted.consumption)) %>%
  summarize(mean(cv.space), sd(cv.space))

lt %>% group_by(year) %>% 
  summarize(cv.time = sd(predicted.consumption)/mean(predicted.consumption)) %>%
  summarize(mean(cv.time), sd(cv.time))

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
  arrange(site,transect, year) %>%
  mutate(deltaK = MAPY - lag(MAPY, order_by = year)) %>% 
  filter(deltaK != 0 & MAPY != 0) %>%
  arrange(site,transect, year)

lt2$dummy2 <- ifelse(lt2$urc.biomass >= 668, "Urchin biomass >= 668", "Urchin biomass < 668")


lmer6 <- lmer(deltaK ~ per.diff + (1|site) + (1|year), lt2)
summary(lmer6)

newdat <- ggeffects::ggpredict(lmer6, terms = "per.diff", type = "fe")

f5 <- lt2 %>%
  ggplot(aes(x = per.diff, y = deltaK))+
  #scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu"))+
  geom_point(aes(size = urc.biomass, fill = MAPY, color = dummy2), pch = 21)+
  scale_color_manual(values = c("gray", "black"))+
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(n = 9, name = "Greens"), trans = "log10")+geom_line(data = newdat, aes(x = x, y = predicted))+
  geom_ribbon(data = newdat, aes(x = x, ymin = conf.low, ymax = conf.high, y = predicted), alpha = .2) +
  # geom_line(data = newdat, aes(x = per.diff, y = LC), lty = 4)+
  # geom_line(data = newdat, aes(x = per.diff, y = UC), lty = 4)+
  geom_hline(yintercept = 0, lty = 3)+
  geom_vline(xintercept = 0, lty = 3)+
  labs(x = "Proportional difference between\nconsumption and detrial supply rate", 
       y = "delta Kelp biomass", 
       size = "Urchin biomass\ndensity", 
       color = "", 
       fill = "Kelp biomass\nin year t")+
  annotate(x = c(-0.5, -0.5, 0.5, 0.5), y = c(15000, -15000, 15000, -15000), geom = "text", label = c("Kelp increases\nConsumption < detrital supply", "Kelp decreases\nConsumption < detrital supply", "Kelp increases\nConsumption > detrital supply", "Kelp decreases\nConsumption > detrial supply"))+
  ggpubr::theme_pubr(legend = "right")

ggsave(here("figures", "perdiffXdeltaK.png"), f5, device = "png", width = 10.63, height = 7.53)
ggsave(here("figures", "perdiffXdeltaK.pdf"), f5, device = "pdf", useDingbats = FALSE, width = 10.63, height = 7.53)


#------------------------------------------------------
## Time series plot
#------------------------------------------------------


# time series

t1 <- lt %>% #LTER data 
  group_by(year, site) %>%
  summarize(biomass = mean(MAPY, na.rm = T)) %>% #adding a biomass column
  ggplot(aes(x = year, y = biomass))+
  geom_line(aes(color = site), show.legend = F)+
  scale_color_manual(values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6'))+
  stat_summary(fun = mean, geom = "line", lwd = 2, alpha = 0.75)+
  labs(y = expression(paste("Giant kelp biomass (g m"^"-2"*")")), x = "")+
  theme_classic()

t2 <- lt %>%
  group_by(year,site) %>%
  summarize(predicted.consumption = mean(predicted.consumption, na.rm = T)) %>%
  ggplot(aes(x = year, y = predicted.consumption))+ #average consumption across sites
  geom_line(aes(color = site), show.legend = F)+
  scale_color_manual(values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6'))+
  stat_summary(fun = mean, geom = "line", lwd = 2, alpha = 0.75)+
  labs(y = expression(paste("Predicted kelp consumption (g m"^"-2"*"d"^"-1"*")")), x = "")+
  theme_classic()

fig4 <- cowplot::plot_grid(t1, t2, align = "v", nrow = 2)

ggsave(here("figures", "timeseries.pdf"), fig4, units = "in", width = 12, height = 6, useDingbats = FALSE)


#----------------------------------------------------
## Faceted time series plots
#----------------------------------------------------


lt %>% #LTER data 
  group_by(year, site) %>%
  summarize(kelp.biomass = mean(MAPY, na.rm = T), 
            urc.biomass = mean(urc.biomass, na.rm = T)) %>%
  pivot_longer(names_to = "species", values_to = "biomass", -c(year, site)) %>%
  ggplot(aes(x = year, y = biomass))+
  geom_line(aes(color = site), show.legend = F)+
  facet_grid(species ~ site , scales = "free")


means <- lt %>% #LTER data 
  group_by(year, site) %>%
  summarize(MAPY = mean(MAPY, na.rm = T), 
            urc.biomass = mean(urc.biomass, na.rm = T)) %>%
  pivot_longer(names_to = "species", values_to = "biomass", -c(year, site))

plot <- lt %>% #LTER data 
  group_by(year, site, transect) %>%
  select(year, site, transect, MAPY, urc.biomass) %>%
  pivot_longer(names_to = "species", values_to = "biomass", -c(year, site, transect)) %>%
  ggplot(aes(x = year, y = biomass))+
  geom_line(aes(group = transect), show.legend = F, alpha = 0.5)+
  geom_line(data = means, aes( x = year, y = biomass, color = site), lwd = 1.5, alpha = 0.75, show.legend = F)+
  facet_grid(species ~ site , scales = "free")+
  theme_bw()+
  theme(panel.grid = element_blank())


ggsave(filename = here::here("figures/", "facet_by_site_species.pdf"), plot = plot, device = "pdf", width = 14, height = 6)

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

plot2 <- lt %>% #LTER data 
  group_by(year, site) %>%
  summarize(kelp.biomass = mean(MAPY, na.rm = T),
            urc.biomass = mean(urc.biomass, na.rm = T)) %>%
  pivot_longer(names_to = "species", values_to = "biomass", -c(year, site)) %>%
  group_by(site, species) %>%
  mutate(biomass2 = as.numeric(scale2(biomass))) %>%
  ggplot(aes(x = year, y = biomass2))+
  geom_line(aes(color = species), lwd = 1)+
  geom_point(aes(color = species), lwd = 1.2)+
  scale_color_manual(values = c("#006d2c", "#810f7c"))+
  facet_wrap(~ site)+
  labs(x = "", y = "Biomass\n(z-scored by species across years within sites)")+
  theme_bw()+
  theme(panel.grid = element_blank())


ggsave(filename = here::here("figures/", "facet_by_site_zscored.pdf"), plot = plot2, device = "pdf")


forplot <- lt %>% #LTER data 
  group_by(year, site) %>%
  summarize(kelp.biomass = mean(MAPY, na.rm = T),
            urc.biomass = mean(urc.biomass, na.rm = T)) %>%
  pivot_longer(names_to = "species", values_to = "biomass", -c(year, site)) %>%
  group_by(site, species) %>%
  mutate(biomass2 = as.numeric(scale2(biomass)))

sites <- unique(forplot$site)

for(i in 1:length(sites)){
  forplot %>% filter(site == sites[i]) %>%
  ggplot(aes(x = year, y = biomass2))+
    geom_line(aes(color = species), lwd = 3, show.legend = F)+
    geom_point(aes(color = species), size = 4, show.legend = F, pch = 21, fill = "white")+
    scale_color_manual(values = c("#006d2c", "#810f7c"))+
    labs(x = "", y = "Biomass")+
    scale_x_continuous(breaks = c(2000, 2005, 2010, 2015, 2020))+
    scale_y_continuous(breaks = c(-2, -1, 0, 1, 2, 3))+
    coord_cartesian(ylim = c(-2, 3))+
    theme_bw()+
    theme(panel.grid = element_blank(), text = element_text(size = 30))+
  ggsave(filename = here::here("figures/", paste("facet_by_site_zscored", sites[i], ".pdf", sep = "")), device = "pdf", useDingbats = FALSE)
}












#-----------------------------------------------------------------------------------------
## OLD code
#-----------------------------------------------------------------------------------------








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


#------------------------------------------------------
## Lag plots
#------------------------------------------------------

ggplot(lt, aes(x = urc.biomass, y = lead(MAPY, n = 1)))+
  geom_jitter(aes(color = dummy), pch = 21)+
  scale_color_manual(values = c("#8f4811", "#35753d"))+
  # geom_line(data = newdat, aes(x = predicted.consumption , y = y)) +
  labs(x = expression(paste("Combined urchin biomass (g m"^"-2"*")")), y = expression(paste("Giant kelp biomass (g m"^"-2"*")")), color = "")+
  ggpubr::theme_pubclean()





#----------------------------------------------------
## Old Time series plots
#----------------------------------------------------

pd <- position_dodge(0.4)

p3b <- lt %>% #LTER data 
  mutate(id = paste(site,transect,"-")) %>%
  ggplot(aes(x = year, y = MAPY))+
  geom_line(aes(group = id), position = pd, color = alpha("black", 0.5))+
  geom_point(aes(group = id, color = dummy, fill = dummy, alpha = dummy, shape = dummy2), show.legend = F, position = pd, size = 3)+
  scale_shape_manual(values = c(21,24))+
  scale_color_manual(values = c("black", "black"))+
  scale_alpha_manual(values = c(0.75, 0))+
  scale_fill_manual(values = c("red", "black"))+
  stat_summary(fun = mean, geom = "line", lwd = 2, col = "forestgreen")+
  labs(fill = "", y = expression(paste("Giant kelp biomass (g m"^"-2"*")")), x = "")+
  ggpubr::theme_pubclean()


ggsave(here::here("figures/p3b.png"), p3b, device = "png", width = 12, height = 3)


p3c <- lt %>% #LTER data 
  mutate(id = paste(site,transect,"-")) %>%
  ggplot(aes(x = year, y = predicted.consumption))+
  geom_line(aes(group = id), position = pd, color = alpha("black", 0.5))+
  geom_point(aes(group = id, color = dummy, fill = dummy, alpha = dummy, shape = dummy2), show.legend = F, position = pd, size = 3)+
  scale_shape_manual(values = c(21,24))+
  scale_color_manual(values = c("black", "black"))+
  scale_alpha_manual(values = c(0.75, 0))+
  scale_fill_manual(values = c("red", "black"))+
  stat_summary(fun = mean, geom = "line", lwd = 2, col = "forestgreen")+
  labs(y = expression(paste("Predicted kelp consumption (g m"^"-2"*"d"^"-1"*")")), x = "")+
  ggpubr::theme_pubclean()

ggsave(here::here("figures/p3c.png"), p3c, device = "png", width = 12, height = 3)

cowplot::plot_grid(p3b,p3c, ncol = 1, align = "v")
ggsave(here::here("figures/p3bc.png"), device = "png", width = 12, height = 6)

