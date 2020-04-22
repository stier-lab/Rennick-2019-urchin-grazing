library(here)
library(tidyverse)

# So based on the Ocean tipping point guide and Cury et al. 2011, my intention was to identify a tipping point in the relationship between urchin biomass and kelp biomass. This is similar to the idea that Ling et al. 2015 (ProcB) used, except Ling used percent cover of macroalgea which I believe made from stronger relationships in the non-linearity because cover is bounded 0-100%. The idea was to use a model fit, and then run a change point analysis on the model prediction (sensu Cury et al. 2011). However, after plotting the data and fitting GAMs, LMs, and expentials, its apparent that the strongest relationship is actually linear... This suggests no non-linearity in the relationship which is a bit disconcerting. At this point, I'm going to move on without this analysis. NOTE: Ling et al. 2015 used wet biomass NOT dry biomass for urchin biomass density. We need to check and make sure we used wet biomass when we were estimating our densities! Ling et al. SOM used a conversion to go from stipe counts to %cover (%Cover=3.77*SC, which seems strange as it can be > 1...).  



lt <- read.csv("data/survey_data/Annual_All_Species_Biomass_at_transect.csv", stringsAsFactors = F,na.strings ="-99999") %>%
  dplyr::select("YEAR", "MONTH", "SITE", "TRANSECT", "SP_CODE", "PERCENT_COVER", "DENSITY", "WM_GM2", "DRY_GM2", "SCIENTIFIC_NAME", "COMMON_NAME", "GROUP", "MOBILITY", "GROWTH_MORPH", "COARSE_GROUPING" ) %>%
  mutate(id = paste(SITE, TRANSECT, sep = "")) %>%
  filter(COMMON_NAME == "Purple Urchin" | COMMON_NAME == "Red Urchin" | COMMON_NAME == "Giant Kelp")

names(lt) <- tolower(names(lt))

df <- lt %>%
  dplyr::select(year, month, site, transect, sp_code, dry_gm2) %>%
  spread(sp_code, dry_gm2) #THis is tracking 11 sites not 9, did you forget to take out the island sites? # filter(SITE != "SCTW", SITE != "SCDI") %>% 


p <- ggplot(df, aes(x = SPL, y = MAPY))+ #plotting the observed realtionship between purple urchins and kelp biomass
  geom_point(pch = 21, fill = alpha("black", 0.25))+
  geom_smooth(method = "gam", formula = y ~s(x)) 

ggExtra::ggMarginal(p, type = "histogram")

#Now we will try to determine the relationship between kelp biomass and urhcin biomass and develop a model that can predict kelp biomass as a function of urchin biomass. In order to determine this relationship we will build a mixed effect model that tracks this relationship through time across all sites in addition to looking at the global pattern across the entire dataset of kelp biomass as a result of urchin biomass. 

#fit linear mixed effect model tracking kelp biomass as a reult of total urchin biomass through time at each site in order to create model that represents the relationship, and can predict kelp biomass as a result of urchin biomass. 
lmer1 <- lmer(MAPY ~ I(SFL + SPL) + (1|site) + (1|year), data = df)
summary(lmer1) #what does REML criterion at convergence tell you? 
hist(residuals(lmer1))

#Is this significant? Where is the 'a' value? which value will tell you that? 

#testing a nonlinear exponential relationship between purple urchins and kelp biomass. 
#Why are we only using purples now? Becuase they have a greater impact overall? 
exp <- nls(MAPY ~ a * exp(-1 * b * SPL), data = df, start = list(a = 1000, b = 1/100))
summary(exp)

#testing a linear exponential relationship between purple urchins and kelp 
exp2 <- lm(MAPY ~ SPL + I(SPL^2), data = df)
summary(exp2)

#testing a nonlinear regression model on the realtionship between purple urchins and kelp 
hyp <- nls(MAPY ~ a / (b + SPL), data = df, start = list(a = 10, b = 1))
summary(hyp)

#based on these results the realtionship between kelp and urhcin biomass is most likely non linear? exp2 is not a statistically viable model for this relationship. How do you compare the models? Can you use AIC? 

newdat <- data.frame(SPL = seq(0, max(df$SPL), length.out = 1000))
newdat$y <- predict(exp, newdata = newdat)
newdat$y2 <- predict(exp2, newdata = newdat)
newdat$y3 <- predict(hyp, newdata = newdat)

p <- ggplot(df, aes(SPL, MAPY)) + 
  geom_point() + 
  geom_line(data = newdat, aes(x = SPL, y = y), color = "red")+
  #geom_line(data = newdat, aes(x = SPL, y = y2), color = "blue")+
  geom_line(data = newdat, aes(x = SPL, y = y3), color = "green")

ggExtra::ggMarginal(p, type = "histogram")

#All three models look realtively simillar. The standard deviations must be HUGE! 

library(mgcv)

df <- df %>% 
  mutate(mapy = scale(MAPY), 
         sfl = scale(SFL), 
         spl = scale(SPL))


#testing a Gaussian model to the realtionship between kelp and urhcin biomass
gam_y <- gam(mapy ~ s(spl), data = df, method = "REML")
summary(gam_y)  

newdat <- data.frame(spl = seq(min(df$spl), max(df$spl), length.out = 1000))
newdat$y <- predict(gam_y, newdata = newdat)

ggplot(df, aes(spl, mapy)) + 
  geom_point() + 
  geom_line(data = newdat, aes(x = spl, y = y))+
  geom_smooth(method = "lm", color = "red")


# use changepoint ctg.XXXX funciton for changepoint analysis



x <- 1:100
y <- 1/x
plot(y ~ x) #what is this? 


# Cross site and year mixed effect model of urchin x kelp biomass
#how is this different than what we already did up top? 

df <- df %>% as_tibble() %>% mutate(urc = SPL + SFL)

lmer <- lmer(MAPY ~ urc + (1|site) + (1|year), data = df)
summary(lmer)
source("analysis/Functions.R")
library(car)
modelassump(lmer)

#How is a generalized linear mixed effect model different from a normal lmer?
glmer <- glmer(MAPY ~ urc + (1|site) + (1|year), data = df, family = Gamma(link = "log"))
summary(glmer) # I get Error in eval(family$initialize, rho) : non-positive values not allowed for the 'gamma' family

df <- drop_na(df, MAPY)

lmer.ta <- lme(MAPY~urc,random=~1|site,data=df,correlation=corAR1())
summary(lmer.ta)

AIC(lmer.ta, lmer) #it says that models are not all fitted to the same number of observations

newdat <- data.frame(urc = seq(0, max(df$urc, na.rm = T), length.out = 1000))
newdat$y <- predict(lmer, newdata = newdat, re.form = ~0, se = T)

ggplot(df, aes(SPL, MAPY)) + 
  geom_point() + 
  geom_line(data = newdat, aes(x = urc, y = y), color = "red")

