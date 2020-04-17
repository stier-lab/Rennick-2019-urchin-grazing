#############################################################################################################
##Density Analysis
#############################################################################################################

## Mae Rennick and Bart DiFiore
## Urchin Size Data

# --------------------------------------------------------------------------------------------------
## Setup and cleaning
# --------------------------------------------------------------------------------------------------

pf.s <- read.csv("data/size_experiment/raw/Purple_size_data.csv") # size and weight data for purple urhcins. 5 size bins in 10 mm increments from 10 mm-60 mm were tested over two weeks. 7 urchins per tank with 3-4 replicates per size bin each week. Total of 315 urchins in 45 trials.

pf <- read.csv("data/size_experiment/raw/Purple_consuption_data.csv") %>% 
  as_tibble() %>% #coercing the dataframe together
  filter(urchin_size_cat != 'control') %>% #removing the control
  left_join(pf.s) %>% #adding the urchin size and weight data to the dataframe
  dplyr::select(-c(date, days_starved,time_in, time_out, week_no, tank, side)) %>% #removing excess data
  mutate(abundance = urchin_dens, 
         urchin_dens = NULL, 
         sp = "p",
         round = NA, 
         time_ran = 48)%>%
  rename(size_class = urchin_size_cat) #cleaned purple urchin size trials dataset

df_r <- read.csv("data/size_experiment/raw/red_size_data_2.csv") %>% # size and weight data for red urhcins. Red urchins were placed into 5 19mm size bins incrementally from 10-109 mm. Each size bin was replicated three times, with a total of urchins in 45 .
  mutate(round =1, 
         urchin_density=5, 
         rep = NULL, 
         weight= NA, 
         sp = "r", 
         use=1)

rf <-read.csv("data/size_experiment/raw/red_size_data_round2.csv") %>%
  rename(kelp_in = Kelp_in, kelp_out =Kelp_out, sp = urchin) %>%
  bind_rows(df_r) %>%
  dplyr::select(-c(tank, side, level)) %>%
  rename(abundance = urchin_density)

col_order <- c("trial_id", "round", "sp", "abundance", "size_class", "size", "weight", "kelp_in", "kelp_out", "time_ran", "use")
pf <- pf[, col_order]
rf <- rf[, col_order] #reordering the columns so that the red and purple dataframes match - is this too much info? It seems like overkill. 

df <- bind_rows(pf, rf) %>% #combining red and purple datasets
  mutate(id = as.numeric(as.factor(paste(trial_id, round, sp, sep = "")))) %>%
  group_by(id, sp, abundance, kelp_in, kelp_out, use, round, time_ran) %>%
  summarize(mean.size = mean(size), 
            biomass = sum(weight)) %>%
  mutate(herbivory_rate = ((kelp_in - kelp_out)/time_ran)*24) %>%
  filter(time_ran == 48, use == 1)  #what do you do about all of the NA 'rounds'?

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

exp2 <- lm(herbivory_rate ~ (mean.size + I(mean.size^2)) * round, data = r) #what does 'round' do to it? 
summary(exp2) #exponential model for the relationship between purple urchin test size and herbivory rate #negative intercept?

exp3 <- lm(herbivory_rate ~ mean.size + I(mean.size^2), data = r) #no round
summary(exp3) #exponential model for the relationship between purple urchin test size and herbivory rate

AIC(lm2, exp2, exp3) #all relatively simillar fits
anova(exp3, exp2)

  # the exponential function predicts values less than zero, which is biologically impossible. So I"m going to test the fit of a ricker function. Ricker functions are commonly used as phenomenological models to describe ecological data that is non-negative, increase from zero to a peak and then decline to zero. 

# y ~ a*x*exp(-b*x) #can you explain to me what this equation means? 

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

