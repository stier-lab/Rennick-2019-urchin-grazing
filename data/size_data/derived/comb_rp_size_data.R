## Mae Rennick and Bart DiFiore
## Urchin Size Data

# --------------------------------------------------------------------------------------------------
## Setup and cleaning
# --------------------------------------------------------------------------------------------------

pf.s <- read.csv("data/size_experiment/raw/Purple_size_data.csv")

pf <- read.csv("data/size_experiment/raw/Purple_consuption_data.csv") %>% 
  as_tibble() %>%
  filter(urchin_size_cat != 'control') %>% 
  left_join(pf.s) %>%
  dplyr::select(-c(date, days_starved,time_in, time_out, week_no, tank, side)) %>%
  mutate(abundance = urchin_dens, 
         urchin_dens = NULL, 
         sp = "p",
         round = NA, 
         time_ran = 48)%>%
  rename(size_class = urchin_size_cat)

df_r <- read.csv("data/size_experiment/raw/red_size_data_2.csv") %>% 
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
rf <- rf[, col_order]

df <- bind_rows(pf, rf) %>% 
  mutate(id = as.numeric(as.factor(paste(trial_id, round, sp, sep = "")))) %>%
  group_by(id, sp, abundance, kelp_in, kelp_out, use, round, time_ran) %>%
  summarize(mean.size = mean(size), 
            biomass = sum(weight)) %>%
  mutate(herbivory_rate = (kelp_in - kelp_out)/time_ran) %>%
  filter(time_ran == 48)

# --------------------------------------------------------------------------------------------------
## Modeling and visulization
# --------------------------------------------------------------------------------------------------

ggplot(df, aes(x = biomass, y = herbivory_rate))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~sp, scales = "free_x")

ggplot(df, aes(x = mean.size, y = herbivory_rate))+
  geom_point(aes(color = as.factor(round)))+
  geom_smooth()+
  facet_wrap(~sp, scales = "free_x")

# would it be better to model it as biomass? or as test diameter? we only have biomass for one round of the red data, so it would limit the sample size.


# Test diameter

#purples
p <- df[df$sp == "p", ]
lm1 <- lm(herbivory_rate ~ mean.size, data = p)
summary(lm1)

exp1 <- lm(herbivory_rate ~ mean.size + I(mean.size^2), data = p)
summary(exp1)

AIC(lm1, exp1)

np <- data.frame(mean.size = seq(min(p$mean.size, na.rm = T), max(p$mean.size, na.rm = T), length.out = 1000), sp = "p")
np$pred <- predict(lm1, newdata = np)

#reds
r <- df[df$sp == "r", ]
lm2 <- lm(herbivory_rate ~ mean.size + round, data = r)
summary(lm2)

exp2 <- lm(herbivory_rate ~ mean.size + I(mean.size^2) * round, data = r)
summary(exp2)

exp3 <- lm(herbivory_rate ~ mean.size + I(mean.size^2), data = r)
summary(exp3)

AIC(lm2, exp2, exp3)
anova(exp3, exp2)

  # the exponential function predicts values less than zero, which is biologically impossible. So I"m going to test the fit of a ricker function. Ricker functions are commonly used as phenomenological models to describe ecological data that is non-negative, increase from zero to a peak and then decline to zero. 

# y ~ a*x*exp(-b*x)

ric <- nls(herbivory_rate ~ a * mean.size * exp(-1 * b * mean.size), data = r, start = list(a = 0.006, b = 1/50))
summary(ric)

AIC(lm2, exp2, exp3, ric) # by AIC the exponential models seem to be better fits to the data, but they aren't biologically reasonable. So I'll plot the ricker function.

nr <- data.frame(mean.size = seq(min(r$mean.size, na.rm = T), max(r$mean.size, na.rm = T), length.out = 1000), sp = "r")
nr$pred <- predict(ric, newdata = nr)
nr$pred.exp <- predict(exp3, newdata = nr)

# So the ricker doesn't appear the fit the data very well. What about a power ricker? 
pric <- nls2::nls2(herbivory_rate ~ a * mean.size ^ gamma * exp(-1 * b * mean.size), data = r[!is.na(r$mean.size), ], start = list(a = 0.01, b = 0.02, gamma = 0.001), algorithm = "plinear")
summary(pric) # can't get this to converge...


plot(herbivory_rate ~ mean.size, data = r)
lines(pred ~ mean.size, nr, col = "red")
lines(pred.exp ~ mean.size, nr, col = "darkgreen")



nr <- data.frame(mean.size = seq(min(r$mean.size, na.rm = T), max(r$mean.size, na.rm = T), length.out = 1000), sp = "r")
nr$pred <- predict(exp3, newdata = nr)
nr <- nr[nr$pred > 0, ]
newdat <- rbind(nr, np)
# PLot it up

df$sp <- ifelse(df$sp == "p", "Purple urchin", "Red urchin")
newdat$sp <- ifelse(newdat$sp == "p", "Purple urchin", "Red urchin")

fig3 <- ggplot(df, aes(x = mean.size, y = herbivory_rate))+
  geom_jitter(aes(fill = sp), pch = 21, show.legend = F)+
  scale_fill_manual(values = c("#762a83", "#d73027"))+
  geom_line(data = newdat, aes(x = mean.size, y = pred))+
  facet_wrap(~sp, scales = "free_x")+
  ggpubr::theme_pubclean()+
  labs(x = "Mean test diameter (mm)", y = expression(paste("Herbivory rate (g h"^"-1"*")")), color = "", linetype = "")+
  theme(strip.background = element_blank())

ggsave("figures/herbivoryXsize_fig3.png", fig3, device = "png")






















