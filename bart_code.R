#######################################################
## Plot it up
#######################################################
library(tidyverse)

df <- read.csv("data/urchin_density_data_all.csv") %>% 
  select(trial = trial_id, den = urchin_density, urc_size = urchin_size, urc_mass = urchin_mass, kelp_in = kelp_in_total, kelp_out = Kelp_out_total) %>%
  mutate(kelp_con = kelp_in - kelp_out, 
         trial_id = trial) %>%
  separate(trial, 
           into = c("species", "junk"), 
           sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  select(-junk) %>%
  group_by(trial_id) %>%
  mutate(biomass = sum(urc_mass)) %>%
  drop_na(den)

df %>% select(-urc_size, -urc_mass) %>%
  filter(species == "P") %>%
  distinct() %>%
  ggplot(aes(x = den, y = kelp_con))+
    geom_jitter(pch = 1)+
    geom_smooth()+
    geom_smooth(method = "lm", color = "red")+
    labs(x = "Number of individuals", y = "Total kelp consumption")

df %>% select(-urc_size, -urc_mass) %>%
  filter(species == "P") %>%
  distinct() %>%
  mutate(percap.con = kelp_con/den) %>%
  ggplot(aes(x = den, y = percap.con))+
    geom_jitter(pch = 1, aes(color = as.factor(den)))+
    labs(x = "Number of individuals", y = "Per capita kelp consumption")+
    guides(color=FALSE)

df %>% select(-urc_size, -urc_mass) %>%
  filter(species == "P") %>%
  distinct() %>%
  ggplot(aes(x = biomass, y = kelp_con))+
  geom_jitter(pch = 1)+
  geom_smooth()

df %>% select(-urc_size, -urc_mass) %>%
  filter(species == "P") %>%
  distinct() %>%
  mutate(kelp_con = ifelse(kelp_con <0, 0, kelp_con)) %>%
  ggplot(aes(x = biomass, y = kelp_con/biomass))+
  geom_jitter(pch = 1)+
  geom_smooth()


mf <- df %>% ungroup() %>%
  select(-urc_size, -urc_mass, -trial_id) %>%
  filter(species == "P") %>%
  distinct() %>%
  mutate(kelp_con = ifelse(kelp_con <0, 0, kelp_con)) %>%
  mutate(week = ifelse(den > 24 & den != 44, 2, 1)) %>%
  as.data.frame()

mf$den <- as.numeric(mf$den)

lm1 <- lm(kelp_con ~ 0+ den, data = as.data.frame(mf))
summary(lm1)
hist(residuals(lm1))
plot(residuals(lm1) ~ fitted(lm1))

lm2 <- lm(kelp_con ~ 0 + I(den^2) + den, data = mf)
summary(lm2)
AIC(lm1,lm2)

plot(kelp_con ~ den, data = mf)
abline(lm1, col = "red")
lines(seq(min(mf$den), max(mf$den), length.out = 1000), predict(lm2, newdata = list(den = seq(min(mf$den), max(mf$den), length.out = 1000))), lwd = 1.5, col = "blue")


typeIII <- function(den, a, h, q){
  (a*den^(1+q)) / (1 + a*h*den^(1+q))
}

starts <- expand.grid(a = seq(0.01, 1, length.out = 10), h = seq(0,10, length.out = 10), q = seq(1,1.9, length.out = 10))

fit <- nls2::nls2(kelp_con ~ typeIII(den = den, a = a, h = h, q = q), data = mf, start = starts, algorithm = "brute-force")

fit2 <- nls(kelp_con ~ typeIII(den = den, a = a, h = h, q = q), data = mf, start = coef(fit))

plot(kelp_con ~ den, data = mf, xlab = "Number of urchins", ylab = "Total kelp consumption (g)", col = mf$week)
abline(lm1, col = "red")
lines(seq(min(mf$den), max(mf$den), length.out = 1000), predict(lm2, newdata = list(den = seq(min(mf$den), max(mf$den), length.out = 1000))), lwd = 1.5, col = "blue")
curve(typeIII(den = x, a = coef(fit2)[1], h = coef(fit2)[2], q = coef(fit2)[3]),seq(min(mf$den), max(mf$den), length.out = 1000), add = T, col = "green" )
abline(v=22,lty =4)

AIC(lm1, lm2, fit2)

lm3 <- lm(I(kelp_con/biomass) ~ biomass, data = mf)
summary(lm3)

lm.break <- segmented(lm1, seg.Z = ~ den)

