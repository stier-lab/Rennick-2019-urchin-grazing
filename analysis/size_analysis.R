#############################################################################################################
##Density Analysis
#############################################################################################################

## Mae Rennick and Bart DiFiore
## Urchin Size Data

#Here, we use the data from the size dependent herbivory trials in order to create a model that predicts the influence of red and purple urchin size on herbivory rate of giant kelp.519 red and purple urchins tested. 

# --------------------------------------------------------------------------------------------------
## Setup and cleaning
# --------------------------------------------------------------------------------------------------

library(here)
library(tidyverse)
source(here("analysis", "Functions.R"))
library(car)


# size and weight data for purple urhcins. 6 size bins in 10 mm increments from 10 mm-69 mm were tested over two weeks. 7 urchins per tank with 7 replicates per size bin. Total of 294 urchins in 42 trials.

pf <- read.csv("data/size_experiment/raw/purple_urchin_size_data_raw_c.csv") %>% 
  as_tibble() %>% #coercing the dataframe together
  dplyr::select(-c(date, days_starved, week_no, tank, time_out, side, use)) %>% #removing excess data
  mutate(abundance = urchin_dens, 
         urchin_dens = NULL, 
         sp = "p",
         round = NA, 
         time_ran = 48)%>%
  rename(size_class = urchin_size_cat)#cleaned purple urchin size trials dataset

pf$size_class <- as.character(pf$size_class)

df_r <- read.csv("data/size_experiment/raw/red_size_data_2.csv") %>% # size and weight data for red urhcins. Round 1. Red urchins were placed into 5 19mm size bins incrementally from 10-109 mm. Each size bin was replicated 3 times over a week, with a total of urchins 75 urchins.
  filter(rep != "control") %>% 
  mutate(round =1, 
         urchin_density=5, 
         rep = NULL, 
         weight= NA, 
         sp = "r", 
         use=1) %>% 
  filter(time_ran==96)

rf <-read.csv("data/size_experiment/raw/red_size_data_round2.csv") %>% #3 more replicates of each of the five size classes of red urchins with 5 urchins in each trial. #75 urchins
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
  group_by(id, sp, abundance, kelp_in, kelp_out, round, time_ran) %>%
  summarize(mean.size = mean(size), 
            biomass = sum(weight), 
            mean.biomass = mean(weight)) %>%
  mutate(herbivory_rate = ((kelp_in - kelp_out)/time_ran)*24, 
         percap.herb_rate = herbivory_rate / abundance)
  

# --------------------------------------------------------------------------------------------------
## Modeling and visulization
# --------------------------------------------------------------------------------------------------

ggplot(df, aes(x = biomass, y = herbivory_rate))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~sp, scales = "free_x") #my bad again for not measuring the mass of each of these urchins

ggplot(df, aes(x = mean.size, y = herbivory_rate))+
  geom_point(aes(color = as.factor(round)))+
  geom_smooth()+
  facet_wrap(~sp, scales = "free_x")

# would it be better to model it as biomass? or as test diameter? we only have biomass for one round of the red data, so it would limit the sample size.


# Test diameter

#testing different models to find the best fit model for predicting the effect of test diameter on herbivory rate

#purples
p <- df[df$sp == "p", ]
lm1 <- lm(herbivory_rate ~ mean.size, data = p)
summary(lm1)#linear regression model for the relationship between urchin test size and herbivory rate

exp1 <- lm(herbivory_rate ~ mean.size + I(mean.size^2), data = p)
summary(exp1) #exponential model for the relationship between purple urchin test size and herbivory rate

AIC(lm1, exp1) #linear model fits the data better

np <- data.frame(mean.size = seq(min(p$mean.size, na.rm = T), max(p$mean.size, na.rm = T), length.out = 1000), sp = "p")
np$pred <- predict(lm1, newdata = np) #developing linear model predictions for the purple data

#reds
r <- df[df$sp == "r", ]
lm2 <- lm(herbivory_rate ~ mean.size + round, data = r) #linear regression model for the relationship between red urchin test size and herbivory rate
summary(lm2) #Cannot confirm that the slop is different from zero according to the p-value. Low R2 value.

exp2 <- lm(herbivory_rate ~ (mean.size + I(mean.size^2)) * round, data = r) 
summary(exp2) #exponential model for the relationship between purple urchin test size and herbivory rate #negative intercept?

exp3 <- lm(herbivory_rate ~ mean.size + I(mean.size^2), data = r)
summary(exp3) #exponential model for the relationship between purple urchin test size and herbivory rate

AIC(lm2, exp2, exp3) #all relatively simillar fits
anova(exp3, exp2)

  # the exponential function predicts values less than zero, which is biologically impossible. So I"m going to test the fit of a ricker function. Ricker functions are commonly used as phenomenological models to describe ecological data that is non-negative, increase from zero to a peak and then decline to zero. 

# y ~ a*x*exp(-b*x) #Bart- can you explain to me what this equation means? 

ric <- nls(herbivory_rate ~ a * mean.size * exp(-1 * b * mean.size), data = r, start = list(a = 0.006, b = 1/50)) #where did you get the a be and x values from? 
summary(ric)

AIC(lm2, exp2, exp3, ric) # by AIC the exponential models seem to be better fits to the data, but they aren't biologically reasonable. So I'll plot the ricker function.

nr <- data.frame(mean.size = seq(min(r$mean.size, na.rm = T), max(r$mean.size, na.rm = T), length.out = 1000), sp = "r") #why are you calculating the mean size for all the reds? 
nr$pred <- predict(ric, newdata = nr)
nr$pred.exp <- predict(exp3, newdata = nr)

# So the ricker doesn't appear the fit the data very well. What about a power ricker? #how do you know htis other than the AIC values? Did you plot it and then delete it? 
pric <- nls2::nls2(herbivory_rate ~ a * mean.size ^ gamma * exp(-1 * b * mean.size), data = r[!is.na(r$mean.size), ], start = list(a = 0.01, b = 0.02, gamma = 0.001), algorithm = "plinear")
summary(pric) # can't get this to converge...
#lots of errors: is that due to nonconvergence? What does that mean? What do you need to converge and why? 


plot(herbivory_rate ~ mean.size, data = r)
lines(pred ~ mean.size, nr, col = "red") #exp3 model? 
lines(pred.exp ~ mean.size, nr, col = "darkgreen") #ricker model? 



nr <- data.frame(mean.size = seq(min(r$mean.size, na.rm = T), max(r$mean.size, na.rm = T), length.out = 1000), sp = "r") #how is this different from the nr above? Are you replacing it with new predictions? What happened? 
nr$pred <- predict(exp3, newdata = nr)
nr <- nr[nr$pred > 0, ]
newdat <- rbind(nr, np)

# PLot it up

df$sp <- ifelse(df$sp == "p", "Purple urchin", "Red urchin")
newdat$sp <- ifelse(newdat$sp == "p", "Purple urchin", "Red urchin")

fig3 <- ggplot(df, aes(x = mean.size, y = herbivory_rate))+
  geom_jitter(aes(fill = sp), pch = 21, show.legend = F)+
  scale_fill_manual(values = c("#762a83", "#d73027"))+
  geom_line(data = newdat, aes(x = mean.size, y = pred))+ #is it mapping 'pred' for both of them? , how did you get the red model on there?
  facet_wrap(~sp, scales = "free_x")+
  ggpubr::theme_pubclean()+
  labs(x = "Mean test diameter (mm)", y = expression(paste("Herbivory rate (g m"^"-2"*"d"^"-1"*")")), color = "", linetype = "")+
  theme(strip.background = element_blank())

ggsave("figures/herbivoryXsize_fig3.png", fig3, device = "png")



#--------------------------------------------------------------------------------
## MTE scaling predictions - BART scrap
#--------------------------------------------------------------------------------

plot(herbivory_rate ~ biomass, p)

mte <- lm(log(herbivory_rate) ~ log(biomass), p)
summary(mte)

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
  filter(sp == "p") %>%
  ggplot(aes(y = percap.herb_rate, x = mean.biomass))+
  geom_point()

df %>%
  filter(sp == "p") %>%
  ggplot(aes(y = herbivory_rate, x = biomass))+
  geom_point()

summary(lm(log(percap.herb_rate) ~ log(mean.biomass), df[df$sp == "p", ]))

