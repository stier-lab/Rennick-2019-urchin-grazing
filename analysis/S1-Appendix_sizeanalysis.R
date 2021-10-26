####################################################################################################
## Size Analysis
####################################################################################################

## Mae Rennick and Bart DiFiore
## Urchin Size Data

#The following code analyzed foraging rates of urchins of different size classes based on trials conducted in summer 2018.
# --------------------------------------------------------------------------------------------------
## Setup and cleaning
# --------------------------------------------------------------------------------------------------

library(here)
library(tidyverse)
source(here("analysis", "Functions.R"))
library(car)
library(ggeffects)

df <- read.csv("data/size_data/combined_sizedata.csv") %>% 
  group_by(trial_id, sp, urchin_density, kelp_in, kelp_out, time_ran, size_class) %>%
  summarize(mean.testdiameter = mean(test_diameter), # average test diameter of urchins within a particular size class
            total.mass = sum(mass), # sum of all masses of urchins in a size class 
            mean.mass = mean(mass)) %>% # average mass of individual urchins in a size class
  mutate(herbivory_rate = ((kelp_in - kelp_out)/time_ran)*24, 
         percap.herb_rate = herbivory_rate / urchin_density, 
         sp = case_when(sp == "p" ~ "Purple Urchin", 
                        sp == "r" ~ "Red Urchin")) %>% 
  drop_na(total.mass)


# --------------------------------------------------------------------------------------------------
## Modeling and visualization
# --------------------------------------------------------------------------------------------------

df.p <- df %>% filter(sp == "Purple Urchin")
df.r <- df %>% filter(sp == "Red Urchin")

# Purple urchins only

# Linear model
lm1 <- lm(percap.herb_rate ~ mean.mass, df.p)
summary(lm1)
qqPlot(residuals(lm1))
hist(residuals(lm1))

# Log-log model (i.e. a power function) 
pwf <- lm(log(percap.herb_rate) ~ log(mean.mass), df.p)
summary(pwf)
qqPlot(residuals(pwf, type = "response"))

pwf.2 <- nls(percap.herb_rate ~ a*mean.mass^beta, df.p, start = list(a = exp(0.13527), beta = 0.4), )
summary(pwf.2)
qqPlot(residuals(pwf.2))
hist(residuals(pwf.2))

anova(lm1, pwf.2)
AIC(lm1, pwf.2)  


# Red urchins only

# Linear model
lm2 <- lm(percap.herb_rate ~ mean.mass, df.r)
summary(lm2)
qqPlot(residuals(lm2))
hist(residuals(lm2))

# Log-log model (i.e. a power function) 
pwf.r <- lm(log(percap.herb_rate+0.01) ~ log(mean.mass), df.r)
summary(pwf.r)
qqPlot(residuals(pwf.r, type = "response"))

pwf.r.2 <- nls(percap.herb_rate ~ a*mean.mass^beta, df.r, start = list(a = exp(0.13527), beta = 0.4), )
summary(pwf.r.2)
qqPlot(residuals(pwf.r.2))
hist(residuals(pwf.r.2))

# Exponential model
exp2 <- lm(percap.herb_rate ~ mean.mass + I(mean.mass^2), df.r)
summary(exp2)  
qqPlot(residuals(exp2))

AIC(lm2, pwf.r.2, exp2)

np <- data.frame(mean.mass = seq(min(df.p$mean.mass, na.rm = T), max(df.p$mean.mass, na.rm = T), length.out = 1000), sp = "Purple Urchin")
np$pred <- exp(predict(pwf, newdata = np))
np$pred.low <- np$pred - exp(predict(pwf, newdata = np, se.fit = T)$se.fit)
np$pred.high <- np$pred + exp(predict(pwf, newdata = np, se.fit = T)$se.fit)
np$pred.lm <- predict(lm1, newdata = np)
np$pred.lm.high <- np$pred.lm + predict(lm1, newdata = np, se.fit = T)$se.fit
np$pred.lm.low <- np$pred.lm - predict(lm1, newdata = np, se.fit = T)$se.fit
newdat <- np

# PLot it up


figs2 <- ggplot(df, aes(x = mean.mass, y = percap.herb_rate))+
  geom_jitter(aes(fill = sp), pch = 21, show.legend = F)+
  scale_fill_manual(values = c("#762a83", "#d73027"))+
  geom_line(data = newdat, aes(x = mean.mass, y = pred))+
  geom_ribbon(data = newdat, aes( x = mean.mass, ymax = pred.high, ymin = pred.low, y = pred), color = "gray", alpha = 0.1)+
  geom_line(data = newdat, aes(x = mean.mass, y = pred.lm), linetype= 3)+
  geom_ribbon(data = newdat, aes( x = mean.mass, ymax = pred.lm.high, ymin = pred.lm.low, y = pred.lm), color = "gray", alpha = 0.1)+
  facet_wrap(~sp, scales = "free_x")+
  cowplot::theme_cowplot()+
  #theme_bw()+
  labs(x = "Body mass (g)", y = expression(paste("Per capita herbivory rate (g m"^"-2"*"d"^"-1"*")")), color = "", linetype = "")+
  theme(strip.background = element_blank())

cowplot::save_plot("figures/herbivoryXsize.png", figs2, device = "png")


#--------------------------------------------------------------------------------------------
## Size-based prediction of density trials
#--------------------------------------------------------------------------------------------


biomass.trials <- read.csv("data/density_experiment/mesocosm_density_data.csv", stringsAsFactors = FALSE, fileEncoding="UTF-8-BOM") %>% 
  as_tibble() %>% 
  select(kelp_in, kelp_out, urchin_density, tank, date, trial_number, p_r, trial_id, total_time, urchin_size, urchin_mass, mortality) %>% 
  mutate(kelp_consumed=(kelp_in-kelp_out)) %>% 
  mutate (herbivory_rate = (kelp_consumed/total_time)*24, 
          abundance = urchin_density, 
          urchin_density = NULL) %>% 
  group_by(date, p_r, trial_number, trial_id, tank, total_time, kelp_in, kelp_out, mortality, kelp_consumed, abundance, herbivory_rate) %>%
  mutate(CB = ifelse(p_r == "p", exp(-1.18759)*urchin_mass^0.32019, 
                     urchin_mass)) %>% # estimate consumptive biomass to correct for size-specific differences
  summarize(biomass= sum(urchin_mass )/1.587, 
            predicted_CC_size = sum(CB)/1.587) %>%
  mutate(urchin_size = NULL, 
         urchin_mass  = NULL, 
         sp = ifelse(p_r == "p", "Purple urchin", "Red urchin")) %>% 
  pivot_longer(cols = c(herbivory_rate, predicted_CC_size))

figs3 <- ggplot(biomass.trials[biomass.trials$sp == "Purple urchin", ], aes(x = biomass, y = value, group = name))+
  geom_point(aes(color = name), show.legend = F)+
  scale_color_manual(values = c("black", "red"))+
  #geom_smooth(method = "lm", aes(linetype = name), show.legend = F)+
  labs(x = expression(paste("Biomass (g m"^"-2",")", sep = "")), y = expression(paste("Herbivory rate (g m"^"-2"*"d"^"-1"*")")))+
  cowplot::theme_cowplot()

cowplot::save_plot("figures/figure_S3.png", figs3)




# Scrap



ggplot(df, aes(x = mean.mass, y = percap.herb_rate/mean.mass))+
  geom_jitter(aes(fill = sp), pch = 21, show.legend = F)+
  scale_fill_manual(values = c("#762a83", "#d73027"))+
  facet_wrap(~sp, scales = "free_x")+
  cowplot::theme_cowplot()+
  #theme_bw()+
  #labs(x = "Body mass (g)", y = expression(paste("Per capita herbivory rate (g m"^"-2"*"d"^"-1"*")")), color = "", linetype = "")+
  theme(strip.background = element_blank())


summary(lm(log(I(percap.herb_rate/mean.mass)) ~ log(mean.mass), df.p))
summary(lm(log(I(percap.herb_rate/mean.mass + 0.001)) ~ log(mean.mass), df.r))









