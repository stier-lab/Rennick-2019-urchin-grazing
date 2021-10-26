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


# size and weight data for purple urhcins. 6 size bins in 10 mm increments from 10 mm-69 mm were tested over two weeks. 7 urchins per tank with 7 replicates per size bin. Total of 294 urchins in 42 trials.

pf <- read.csv("data/size_data/raw/purple_urchin_size_data_raw_c.csv") %>% 
  as_tibble() %>% #coercing the dataframe together
  dplyr::select(-c(date, days_starved, week_no, tank, time_out, side, use)) %>% #removing excess data
  mutate(abundance = urchin_dens, 
         urchin_dens = NULL, 
         sp = "p",
         round = NA, 
         time_ran = 48)%>%
  rename(size_class = urchin_size_cat)#cleaned purple urchin size trials dataset

pf$size_class <- as.character(pf$size_class)

df_r <- read.csv("data/size_data/raw/red_size_data_2.csv") %>% # size and weight data for red urhcins. Round 1. Red urchins were placed into 5 19mm size bins incrementally from 10-109 mm. Each size bin was replicated 3 times over a week, with a total of urchins 75 urchins.
  filter(rep != "control") %>% 
  mutate(round =1, 
         urchin_density=5, 
         rep = NULL, 
         weight= NA, 
         sp = "r", 
         use=1) %>% 
  filter(time_ran==96)

rf <-read.csv("data/size_data/raw/red_size_data_round2.csv") %>% #3 more replicates of each of the five size classes of red urchins with 5 urchins in each trial. #75 urchins
  rename(kelp_in = Kelp_in, kelp_out =Kelp_out, sp = urchin) %>%
  filter(time_ran ==96) %>% 
  filter(use==1) %>% 
  bind_rows(df_r) %>%
  dplyr::select(-c(tank, side, level)) %>%
  filter(size_class != "NA") %>% 
  rename(abundance = urchin_density)


col_order <- c("trial_id", "round", "sp", "abundance", "size_class", "size", "weight", "kelp_in", "kelp_out", "time_ran")
pf <- pf[, col_order]
rf <- rf[, col_order] #reordering the columns so that the red and purple dataframes match 

df <- bind_rows(pf, rf) %>%
  mutate(id = as.numeric(as.factor(paste(trial_id, round, sp, sep = "")))) %>%
  group_by(id, sp, abundance, kelp_in, kelp_out, round, time_ran, size_class) %>%
  summarize(mean.testdiameter = mean(size), # average test diameter of urchins within a particular size class
            total.mass = sum(weight), # sum of all masses of urchins in a size class 
            mean.mass = mean(weight)) %>% # average mass of individual urchins in a size class
  mutate(herbivory_rate = ((kelp_in - kelp_out)/time_ran)*24, 
         percap.herb_rate = herbivory_rate / abundance, 
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


biomass.trials <- read.csv("data/density_experiment/raw/urchin_density_data_raw_c.csv", stringsAsFactors = FALSE, fileEncoding="UTF-8-BOM") %>% #cdata from density dependent urchin herbivory lab trials for red and purple urhcins.  5 densities tested for red urchins with three replicates of each. 8 densities tested for purple urchins with four replicates of each. 
  as_tibble() %>% 
  select(kelp_in, kelp_out, urchin_density, tank, tank_size, date, trial_number, p_r, trial_id, total_time, urchin_size, urchin_mass, mortality) %>% 
  mutate(kelp_consumed=(kelp_in-kelp_out)) %>% 
  mutate (herbivory_rate = (kelp_consumed/total_time)*24, 
          abundance = urchin_density, 
          urchin_density = NULL) %>% 
  group_by(date, p_r, trial_number, trial_id, tank, tank_size, total_time, kelp_in, kelp_out, mortality, kelp_consumed, abundance, herbivory_rate) %>%
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









