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
         size_class2 = case_when(sp == "p" & size_class == "1" ~ "10-19", 
                                 sp == "p" & size_class == "2" ~ "20-29",
                                 sp == "p" & size_class == "3" ~ "30-39", 
                                 sp == "p" & size_class == "4" ~ "40-49", 
                                 sp == "p" & size_class == "5" ~ "50-59", 
                                 sp == "p" & size_class == "6" ~ "60-70", 
                                sp == "r" ~ size_class), 
         sp = case_when(sp == "p" ~ "Purple Urchin", 
                        sp == "r" ~ "Red Urchin"))


# --------------------------------------------------------------------------------------------------
## Modeling and visualization
# --------------------------------------------------------------------------------------------------

df.p <- df %>% filter(sp == "Purple Urchin")
df.r <- df %>% filter(sp == "Red Urchin", round == 2)

# Purple urchins only

  # Linear model
  lm1 <- lm(herbivory_rate ~ total.mass, df.p)
  summary(lm1)
  qqPlot(residuals(lm1))
  hist(residuals(lm1))

  # Log-log model (i.e. a power function) 
  pwf <- lm(log(herbivory_rate) ~ log(total.mass), df.p)
  summary(pwf)
  qqPlot(residuals(pwf, type = "response"))
  
  pwf.2 <- nls(herbivory_rate ~ a*total.mass^beta, df.p, start = list(a = exp(0.13527), beta = 0.4), )
  summary(pwf.2)
  qqPlot(residuals(pwf.2))
  hist(residuals(pwf.2))
  
  # Exponential model
  exp1 <- lm(herbivory_rate ~ total.mass + I(total.mass^2), df.p)
  summary(exp1)  
  qqPlot(residuals(exp1))

temp <- ggpredict(lm1)  
temp2 <- ggpredict(pwf)
temp3 <- ggpredict(pwf.2)
temp4 <- ggpredict(exp1)

plot(temp)
plot(temp2)
plot(temp3)
plot(temp4)  
  
anova(lm1, pwf.2)
AIC(lm1, pwf.2, exp1)  

ggplot(df.p, aes(x = size_class2, y = herbivory_rate))+
  geom_bar(stat = "identity")+
  facet_wrap(~sp)


ggplot(df, aes(x = mean.mass, y = herbivory_rate/abundance))+
  geom_point()+
  facet_wrap(~sp, scales = "free")


# Red urchins only

# Linear model
lm2 <- lm(herbivory_rate ~ total.mass, df.r)
summary(lm2)
qqPlot(residuals(lm2))
hist(residuals(lm2))

# Log-log model (i.e. a power function) 
pwf.r <- lm(log(herbivory_rate+0.01) ~ log(total.mass), df.r)
summary(pwf.r)
qqPlot(residuals(pwf.r, type = "response"))

pwf.r.2 <- nls(herbivory_rate ~ a*total.mass^beta, df.r, start = list(a = exp(0.13527), beta = 0.4), )
summary(pwf.r.2)
qqPlot(residuals(pwf.r.2))
hist(residuals(pwf.r.2))

# Exponential model
exp2 <- lm(herbivory_rate ~ total.mass + I(total.mass^2), df.r)
summary(exp2)  
qqPlot(residuals(exp2))

AIC(lm2, pwf.r.2, exp2)

np <- data.frame(total.mass = seq(min(df.p$total.mass, na.rm = T), max(df.p$total.mass, na.rm = T), length.out = 1000), sp = "Purple Urchin")
np$pred <- exp(predict(pwf, newdata = np))
np$pred.low <- np$pred - exp(predict(pwf, newdata = np, se.fit = T)$se.fit)
np$pred.high <- np$pred + exp(predict(pwf, newdata = np, se.fit = T)$se.fit)

# nr <- data.frame(total.mass = seq(min(df.r$total.mass, na.rm = T), max(df.r$total.mass, na.rm = T), length.out = 1000), sp = "r")
# nr$pred <- exp(predict(pwf.r, newdata = nr))
# nr$pred <- NA
# 
# newdat <- rbind(nr, np)

newdat <- np

# PLot it up


figs2 <- ggplot(df, aes(x = mean_mass, y = herbivory_rate))+
  geom_jitter(aes(fill = sp), pch = 21, show.legend = F)+
  scale_fill_manual(values = c("#762a83", "#d73027"))+
  geom_line(data = newdat, aes(x = total.mass, y = pred))+
  geom_ribbon(data = newdat, aes( x = total.mass, ymax = pred.high, ymin = pred.low, y = pred), color = "gray", alpha = 0.25)+
  facet_wrap(~sp, scales = "free_x")+
  cowplot::theme_cowplot()+
  #theme_bw()+
  labs(x = expression(paste("Biomass (g m"^"-2",")", sep = "")), y = expression(paste("Herbivory rate (g m"^"-2"*"d"^"-1"*")")), color = "", linetype = "")+
  theme(strip.background = element_blank())

cowplot::save_plot("figures/herbivoryXsize.png", figs2, device = "png")

#--------------------------------------------------------------------------------------------
## Test the assumption using consumptive biomass
#--------------------------------------------------------------------------------------------

biomass.trials <- read.csv("data/density_experiment/raw/urchin_density_data_raw_c.csv", stringsAsFactors = FALSE, fileEncoding="UTF-8-BOM") %>% #cdata from density dependent urchin herbivory lab trials for red and purple urhcins.  5 densities tested for red urchins with three replicates of each. 8 densities tested for purple urchins with four replicates of each. 
  as_tibble() %>% 
  select(kelp_in, kelp_out, urchin_density, tank, tank_size, date, trial_number, p_r, trial_id, total_time, urchin_size, urchin_mass, mortality) %>% 
  mutate(kelp_consumed=(kelp_in-kelp_out)) %>% 
  mutate (herbivory_rate = (kelp_consumed/total_time)*24, 
          abundance = urchin_density, 
          urchin_density = NULL) %>% 
  group_by(date, p_r, trial_number, trial_id, tank, tank_size, total_time, kelp_in, kelp_out, mortality, kelp_consumed, abundance, herbivory_rate) %>%
  mutate(CB = ifelse(p_r == "p", urchin_mass^0.32019, 
                     urchin_mass)) %>% # estimate consumptive biomass to correct for size-specific differences
  summarize(biomass= sum(urchin_mass )/1.587, 
            CB.tot = sum(CB)/1.587) %>%
  mutate(urchin_size = NULL, 
         urchin_mass  = NULL, 
         sp = ifelse(p_r == "p", "Purple urchin", "Red urchin"))

p1 <- biomass.trials %>% filter(sp == "Purple urchin") %>%
ggplot(aes(x = CB.tot, y = herbivory_rate))+
  geom_point()

p2 <- biomass.trials %>% filter(sp == "Purple urchin") %>%
  ggplot(aes(x = biomass, y = herbivory_rate))+
  geom_point()

cowplot::plot_grid(p1, p2)


lm.CB <- lm(herbivory_rate ~ 0 + CB.tot, biomass.trials[biomass.trials$sp == "Purple urchin", ])
summary(lm.CB)

lm.biomass <- lm(herbivory_rate ~ 0 + biomass, biomass.trials[biomass.trials$sp == "Purple urchin", ])
summary(lm.biomass)


# Run simulations to test the assumption

pop1 <- data.frame(urchin_mass = rnorm(1000, mean = mean(80), sd = 20))
pop1$CB = pop1$urchin_mass^0.32019 
urchin_biomass = sum(pop1$urchin_mass)
CB_biomass = sum(pop1$CB)

k <- 10

k*urchin_biomass*coef(lm.biomass)
k*CB_biomass*coef(lm.CB)


1^0.32019 
10^0.32019 
100^0.32019 



# temp



# Could bootstrap the main results to ensure that the main results are resistant 

# Could run the extropolated analysis using the size-based corrections, then use AIC to see which model best predicts kelp biomass ~ detrital supply + consumptive capacity. 






obs <- read.csv("data/survey_data/LTE_Urchin_All_Years_20210209.csv", stringsAsFactors = FALSE, na.strings = "-99999") %>% 
  janitor::clean_names() %>% 
  select(year, month, site, transect, treatment, common_name, size, count, area) %>% 
  filter(treatment == "CONTROL") %>% 
  group_by(year, month, site, transect, common_name, size, area) %>% 
  summarize(count = sum(count)) %>% #What does this change other than moving it? 
  group_by(year, month, site, transect, common_name, size, area) %>%
  complete(count = full_seq(1:count,1))%>% #transferring each count to one, and seperating each urchin into its own line in the data. 
  select(-count) %>% #removing the count column
  ungroup() %>%
  mutate(urchin_size = as.numeric(size)) #creating a size column and establishing it as a numeric value?? 


# Log-log model (i.e. a power function) 
pwf <- lm(log(herbivory_rate/abundance) ~ log(mean.mass), df.p)
summary(pwf)



a.p = 0.000592598
b.p = 2.872636198
se.p = 1.01
  
a.r = 0.000588783
b.r = 2.917071365
se.r = 1.02

temp <- obs %>% 
  mutate(urchin_td = urchin_size, 
         urchin_size = NULL, 
         urchin_mass = ifelse(common_name == "Purple Urchin", a.p*(urchin_td*10)^b.p+se.p, a.r*(urchin_td*10)^b.r+se.r), 
         urchin_cb = urchin_mass^0.32019, 
         CC_perurchin = exp(-1.18759)*urchin_mass^0.32019) %>%
  ungroup() %>%
  group_by(year, month, site, transect, common_name) %>%
  filter(common_name == "Purple Urchin") %>%
  summarize(CC_biomass = (sum(urchin_mass)/area)*coef(lm.biomass), 
            CC_size = (sum(urchin_cb)/area)*coef(lm.CB), 
            CC_sizedirect =(sum(CC_perurchin)/area) ) %>% 
  distinct()


lm.test <- lm(CC_size ~ CC_biomass, temp)
summary(lm.test)

lm.test2 <- lm(CC_sizedirect ~ CC_biomass, temp)
summary(lm.test2)

ggplot(temp, aes(x = CC_biomass, y = CC_sizedirect))+geom_point()+geom_abline(slope = 1, intercept = 0)

ggplot(temp, aes(x = CC_biomass, y = CC_size))+geom_point()+geom_abline(slope = 1, intercept = 0)


a.p*(obs$urchin_size[obs$common_name == "Purple Urchin"]*10)^b.p+se.p








temp <- obs %>% 
  mutate(urchin_td = urchin_size, 
         urchin_size = NULL, 
         urchin_mass = ifelse(common_name == "Purple Urchin", a.p*(urchin_td*10)^b.p+se.p, a.r*(urchin_td*10)^b.r+se.r)) %>% 
  group_by(common_name) %>% 
  summarize(mean = mean(urchin_mass), 
            sd = sd(urchin_mass))







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
         sp = ifelse(p_r == "p", "Purple urchin", "Red urchin"))

figs3 <- ggplot(biomass.trials[biomass.trials$sp == "Purple urchin", ], aes(x = biomass, y = herbivory_rate))+
  geom_point()+
  geom_point(aes(y = predicted_CC_size), color = "red")+
  labs(x = expression(paste("Biomass (g m"^"-2",")", sep = "")), y = expression(paste("Herbivory rate (g m"^"-2"*"d"^"-1"*")")))+
  cowplot::theme_cowplot()

cowplot::save_plot("figures/figure_S3.png", figs3)

lm.temp1 <- lm(predicted_CC_size ~ biomass, biomass.trials[biomass.trials$sp == "Purple urchin",])
summary(lm.temp1)

lm.temp2 <- lm(herbivory_rate ~ biomass, biomass.trials[biomass.trials$sp == "Purple urchin", ])
summary(lm.temp2)


ggplot(biomass.trials[biomass.trials$sp == "Purple urchin", ], aes(x = predicted_CC_size, y = herbivory_rate))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)


lm.temp3 <- lm(herbivory_rate ~ predicted_CC_size, biomass.trials[biomass.trials$sp == "Purple urchin",])
summary(lm.temp3)

ggplot(biomass.trials[biomass.trials$sp == "Purple urchin", ], aes(x = predicted_CC_size, y = herbivory_rate))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)+
  geom_smooth(method = "lm")


ggplot(biomass.trials[biomass.trials$sp == "Purple urchin", ], aes(x = predicted_CC_size/biomass, y = herbivory_rate/biomass))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)+
  geom_smooth(method = "lm")

#--------------------------------------------------------------------------------
## MTE scaling predictions - BART scrap
#--------------------------------------------------------------------------------

plot(I(herbivory_rate/total.mass) ~ total.mass, df[df$sp == "p", ])

mte <- lm(log(herbivory_rate/abundance) ~ log(mean.mass), df[df$sp == "p", ])
summary(mte)


mte.nls <- nls(I(herbivory_rate/abundance) ~ a*mean.mass^beta, data = df[df$sp == "p", ], start = list(a = 0.001, beta = 0.4) )
summary(mte.nls)

AIC(mte.nls, lm)

newdat <- data.frame(biomass = seq(min(p$biomass, na.rm = T), max(p$biomass, na.rm = T), length.out = 1000))
newdat$pred <- predict(mte, newdata = newdat, se.fit = T)[[1]]
newdat$pred.se <- predict(mte, newdata = newdat, se.fit = T)[[2]]
newdat$mte <- newdat$biomass^0.75

plot(log(herbivory_rate) ~ log(biomass), p)
lines(pred ~ log(biomass), newdat)

plot(herbivory_rate ~ biomass, p, ylab = expression(paste("Herbivory rate (g"['kelp']*"d"^"-1"*")")), xlab = "Total biomass")
lines(exp(pred) ~ biomass, newdat)
lines(exp(pred + pred.se) ~ biomass, newdat, lty =4, col = "gray")
lines(exp(pred - pred.se) ~ biomass, newdat, lty = 4, col = "gray")
abline(lm)

lm <- lm(herbivory_rate ~ biomass, p)

AIC(mte, lm) # looks like the power law function is much better!

plot(I(herbivory_rate/biomass) ~ biomass, p)

summary(lm(log(I(herbivory_rate/biomass)) ~ log(biomass), p))


df %>%
  filter(sp == "Purple urchin") %>%
  ggplot(aes(y = herbivory_rate, x = biomass))+
  geom_point()

df %>%
  filter(sp == "Red urchin") %>%
  ggplot(aes(y = herbivory_rate/mean.biomass, x = mean.biomass))+
  geom_point()

summary(lm(log(percap.herb_rate) ~ log(mean.biomass), df[df$sp == "p", ]))

