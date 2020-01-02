library(tidyverse)
# ------------------------------------------------------------------------------------------------
## Set up and visualization
# ------------------------------------------------------------------------------------------------

df <-read.csv("data/density_experiment/derived/urchin_density_data_cleaned.csv") %>% 
  mutate(kelp_consumed = kelp_in_total-Kelp_out_total) %>% 
  mutate (herbivory_rate = kelp_consumed/total_time, 
          abundance = urchin_density, 
          urchin_density = NULL) %>% 
  group_by(trial_id) %>%
  mutate(biomass= sum(urchin_mass), 
         urchin_size = NULL, 
         urchin_mass = NULL) %>%
  distinct()

ggplot(df, aes(x = biomass, y = herbivory_rate))+
  geom_point(pch = 21)+
  geom_smooth()+
  facet_wrap(~type)+
  ggpubr::theme_pubclean()

# ------------------------------------------------------------------------------------------------
## Purple Analaysis
# ------------------------------------------------------------------------------------------------

pf <- df[df$type == "p", ]

lm1 <- lm(herbivory_rate ~ biomass, pf)
summary(lm1)
exp1 <- lm(herbivory_rate ~ biomass + I(biomass^2), pf)
summary(exp1)

pred <- data.frame(biomass = seq(min(pf$biomass), max(pf$biomass), length.out = 1000))
pred$lm <- predict(lm1, newdata = pred)
pred$exp1 <- predict(exp1, newdata = pred)

sig <- nls2::nls2(herbivory_rate ~ (biomass^(1+q) * a) / (1 + a*h*biomass^(1+q)), data = pf, start = list(a = 0.1, h = 1.5, q = 0), algorithm = "port") # coulnd't get this parameterization to work
summary(sig) 

sig <- nls(herbivory_rate ~ (a * biomass^2) / (b^2 + biomass^2), data = pf, start = list(a = 10, b = 1000)) # this is an alternative parameterization based on Bolker et al. 2008, bestiary of functions p. 22. "a" is similar to handling time, and b is similar to attack rate in this parameterization.

summary(sig)
pred$sig <- predict(sig, newdata = pred)

AIC(lm1, exp1, sig)
# So it seems that there is no evidence for any differences between curves.

ggplot(pf, aes(x = biomass, y = herbivory_rate))+
  geom_jitter(pch = 21, width =30)+
  geom_line(data = pred, aes(x = biomass, y = lm), color = "red")+
  geom_line(data = pred, aes(x = biomass, y = exp1), color = "blue")+
  geom_line(data = pred, aes(x = biomass, y = sig), color = "green")+
  geom_vline(xintercept = coef(sig)[2], lty = "dashed")+
  ggpubr::theme_pubclean()


# does the same pattern apply with abundance (not biomass)

lm2 <- lm(herbivory_rate ~ abundance, pf)
summary(lm2)
exp2 <- lm(herbivory_rate ~ abundance + I(abundance^2), pf)
summary(exp2)

pred2 <- data.frame(abundance = seq(min(pf$abundance), max(pf$abundance), length.out = 1000))
pred2$lm2 <- predict(lm2, newdata = pred2)
pred2$exp2 <- predict(exp2, newdata = pred2)

sig2 <- nls(herbivory_rate ~ (a * abundance^2) / (b^2 + abundance^2), data = pf, start = list(a = 10, b = 22)) # this is an alternative parameterization based on Bolker et al. 2008, bestiary of functions p. 22. "a" is similar to handling time, and b is similar to attack rate in this parameterization.

summary(sig2)
pred2$sig2 <- predict(sig2, newdata = pred2)

AIC(lm2, exp2, sig2)

ggplot(pf, aes(x = abundance, y = herbivory_rate))+
  geom_jitter()+
  geom_line(data = pred2, aes(x = abundance, y = lm2), color = "red")+
  geom_line(data = pred2, aes(x = abundance, y = exp2), color = "blue")+
  geom_line(data = pred2, aes(x = abundance, y = sig2), color = "green")+
  geom_vline(xintercept = coef(sig2)[2], lty = "dashed")+
  ggpubr::theme_pubclean()


#-----------------------------------------------------------------------------------------
## Breakpoint analysis

# package strucchange can only can only be used for time series and ordered data

# package segments seems like it might be a good option


lm1 <- lm(herbivory_rate ~ biomass, pf)
summary(lm1)

# my.seg <- segmented(lm1, 
#                     seg.Z = ~ biomass, 
#                     psi = NA)

my.seg <- segmented(lm1, 
                    seg.Z = ~ biomass, 
                    psi = list(biomass = c(1200)))

summary(my.seg)
plot(herbivory_rate ~ biomass, pf)
lines(fitted(my.seg) ~ biomass, pf)

# get the fitted data
my.fitted <- fitted(my.seg)
my.model <- data.frame(biomass = pf$biomass, herb = my.fitted)

# plot the fitted model
ggplot(my.model, aes(x = biomass, y = herb)) + geom_line()+
  geom_point(data = pf, aes(x = biomass, y = herbivory_rate))

# Doesnt seem like there is much evidence for any breakpoints given the variance in the data


# ------------------------------------------------------------------------------------------------
## Red analysis
# ------------------------------------------------------------------------------------------------ 

rf <- df[df$type == "r", ]


lm1 <- lm(herbivory_rate ~ biomass, rf)
summary(lm1)
exp1 <- lm(herbivory_rate ~ biomass + I(biomass^2), rf)
summary(exp1)

pred3 <- data.frame(biomass = seq(min(rf$biomass), max(rf$biomass), length.out = 1000))
pred3$lm <- predict(lm1, newdata = pred3)
pred3$exp1 <- predict(exp1, newdata = pred3)

sig <- nls2::nls2(herbivory_rate ~ (biomass^(1+q) * a) / (1 + a*h*biomass^(1+q)), data = rf, start = list(a = 0.1, h = 1.5, q = 0), algorithm = "port") # coulnd't get this parameterization to work
summary(sig) 

sig <- nls(herbivory_rate ~ (a * biomass^2) / (b^2 + biomass^2), data = rf, start = list(a = 10, b = 1000)) # this is an alternative parameterization based on Bolker et al. 2008, bestiary of functions p. 22. "a" is similar to handling time, and b is similar to attack rate in this parameterization.

summary(sig)
pred3$sig <- predict(sig, newdata = pred3)

AIC(lm1, exp1, sig)
# So it seems that there is no evidence for any differences between curves.

ggplot(rf, aes(x = biomass, y = herbivory_rate))+
  geom_jitter(pch = 21, width =30)+
  geom_line(data = pred3, aes(x = biomass, y = lm), color = "red")+
  geom_line(data = pred3, aes(x = biomass, y = exp1), color = "blue")+
  geom_line(data = pred3, aes(x = biomass, y = sig), color = "green")+
  geom_vline(xintercept = coef(sig)[2], lty = "dashed")+
  ggpubr::theme_pubclean()

# ------------------------------------------------------------------------------------------------
## Figure 2
# ------------------------------------------------------------------------------------------------ 

temp1 <- pred %>% gather(model, prediction, -biomass) %>% mutate(sp = "Purple urchin")
temp2 <- pred3 %>% gather(model, prediction, -biomass) %>% mutate(sp = "Red urchin")

gg <- bind_rows(temp1, temp2) %>% filter(model != "exp1")
gg$Model <- ifelse(gg$model == "lm", "Linear", "Sigmoid")
df$sp <- ifelse(df$type == "p", "Purple urchin", "Red urchin")


fig2 <- ggplot(df, aes(x = biomass, y = herbivory_rate))+
  geom_jitter(color = "white", pch = 21, show.legend = F)+
  geom_rect(aes(xmin= 668 - 115, xmax=668 + 115, ymin=0, ymax=Inf), color = "gray90", fill = "gray90")+
  geom_vline(xintercept = 668, linetype = 4)+
  geom_jitter(aes(fill = sp), pch = 21, show.legend = F)+
  scale_fill_manual(values = c("#762a83", "#d73027"))+
  geom_line(data = gg, aes(x = biomass, y = prediction, color = Model, linetype = Model))+
  #scale_color_manual(c("#0B61A4", "#FF9200"))+
  facet_wrap(~sp)+
  ggpubr::theme_pubclean()+
  labs(x = "Urchin biomass (g)", y = expression(paste("Herbivory rate (g h"^"-1"*")")), color = "", linetype = "")+
  theme(strip.background = element_blank())
  
ggsave("figures/herbivoryXdensity_fig2.png", fig2, device = "png")











