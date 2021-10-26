#############################################################################################################
##Density Analysis
#############################################################################################################

## Mae Rennick and Bart DiFiore
## Urchin Density Data

library(here)
library(tidyverse)
source(here("analysis", "Functions.R"))
library(car)

#Here, we use the data from the density dependent herbivory trials in order to create a model that predicts the influence of red and purple urchin denisty on herbivory rate on giant kelp. 

# ------------------------------------------------------------------------------------------------
## Set up and visualization and clean data
# ------------------------------------------------------------------------------------------------


df <- read.csv("data/density_experiment/mesocosm_density_data.csv") %>%
  as_tibble() %>% 
  select(kelp_in, kelp_out, urchin_density, tank, date, trial_number, p_r, trial_id, total_time, urchin_size, urchin_mass, mortality) %>% 
  mutate(kelp_consumed=(kelp_in-kelp_out)) %>% 
  mutate (herbivory_rate = (kelp_consumed/total_time)*24, 
          abundance = urchin_density, 
          urchin_density = NULL) %>% 
  group_by(date, p_r, trial_number, trial_id, tank, total_time, kelp_in, kelp_out, mortality, kelp_consumed, abundance, herbivory_rate) %>%
  summarize(biomass= sum(urchin_mass )/1.587) %>%
  mutate(urchin_size = NULL, 
         urchin_mass  = NULL, 
         sp = ifelse(p_r == "p", "Purple urchin", "Red urchin"))

###summary statistics 
df %>% group_by(sp) %>% 
  summarize(mean = mean(herbivory_rate), 
            sd = sd(herbivory_rate), 
            se = sd(herbivory_rate)/n(),
            min = min(herbivory_rate), 
            max = max(herbivory_rate))

# ------------------------------------------------------------------------------------------------
## Purple Analaysis
# ------------------------------------------------------------------------------------------------
#here we are testing different models to determine the best fit model for the red density data

pf <- df[df$p_r == "p", ] #limiting our dataset to only purple urchins

lm1 <- lm(herbivory_rate ~ 0 + biomass, pf) #linear model representing herbivory rate as a function of biomass.
summary(lm1)
modelassump(lm1)

exp1 <- lm(herbivory_rate ~  0 + biomass + I(biomass^2), pf) #testing the exponential relationship between the biomass of purple urchins and their herbivory rate. 
summary(exp1)

pow1 <- lm(log(herbivory_rate+1) ~ 0 + log(biomass), pf) # this is fine but need to fit untransformed so that I can use AIC

summary(pow1)

pred <- data.frame(biomass = seq(min(pf$biomass), max(pf$biomass), length.out = 1000))
pred$lm <- predict(lm1, newdata = pred)#model prediction for a linear function 
pred$exp1 <- predict(exp1, newdata = pred)#model prediction for an exponential function 
pred$pow1 <- predict(pow1, newdat = pred)#model prediction for a power law function 

pow2 <- nls(herbivory_rate ~ a*(biomass^b), data = pf, start = list(a = exp(-1.377), b = 0.29))
summary(pow2) #This is a nonlinear regression model for the realtionship between biomass and herbivory rate.

sig <- nls(herbivory_rate ~ (a * biomass^2) / (b^2 + biomass^2), data = pf, start = list(a = 10, b = 1000)) # this is an alternative parameterization based on Bolker et al. 2008, bestiary of functions p. 22. "a" is similar to handling time, and b is similar to attack rate in this parameterization.
summary(sig)
pred$sig <- predict(sig, newdata = pred) 


AIC(lm1, exp1, pow2, sig) #comparing models
# So it seems that there is no evidence for any differences between curves.

model_compare <- ggplot(pf, aes(x = biomass, y = herbivory_rate))+
  geom_jitter(pch = 21, width =30)+
  geom_line(data = pred, aes(x = biomass, y = lm), color = "red")+
  geom_line(data = pred, aes(x = biomass, y = exp1), color = "blue")+
  geom_line(data = pred, aes(x = biomass, y = sig), color = "green")+
  #geom_line(data = pred2, aes(x = biomass, y = fit), color = "black")+
  geom_vline(xintercept = coef(sig)[2], lty = "dashed")+ 
  ggpubr::theme_pubclean()


ggplot(pf, aes(x = biomass, y = herbivory_rate/biomass))+
  geom_point()+
  geom_smooth(method = "lm") #linear model of biomass v. herbivory

lin_mod <- lm(I(herbivory_rate/biomass) ~ biomass, pf) #Linear model tracking per capita herbivory rate?
summary(lin_mod) #High p value and low R2 value. Cannot prove that the slope is different from zero which means there is no evidence for a nonlinearity. 

exp_mod <- lm(I(herbivory_rate/biomass) ~ biomass + I(biomass^2), pf)
summary(exp_mod)

# does the same pattern apply with abundance (not biomass)

lm2 <- lm(herbivory_rate ~ abundance, pf) #linear regression of abundance v herbivory rate
summary(lm2) #Why does the intercept have a p-value higher than 0.05? 
exp2 <- lm(herbivory_rate ~ abundance + I(abundance^2), pf) #exponential regression of herbivory rate as a function of abundance.
summary(exp2) 


pred2 <- data.frame(abundance = seq(min(pf$abundance), max(pf$abundance), length.out = 1000))
pred2$lm2 <- predict(lm2, newdata = pred2) #linear model predictions for herbivory rate as a function of biomass
pred2$exp2 <- predict(exp2, newdata = pred2) #exponential model predictions for herbivory rate as a function of biomass

sig2 <- nls(herbivory_rate ~ (a * abundance^2) / (b^2 + abundance^2), data = pf, start = list(a = 10, b = 22)) # this is an alternative parameterization based on Bolker et al. 2008, bestiary of functions p. 22. "a" is similar to handling time, and b is similar to attack rate in this parameterization. 

summary(sig2)
pred2$sig2 <- predict(sig2, newdata = pred2)

AIC(lm2, exp2, sig2) #all of the model predictions are simillar

ggplot(pf, aes(x = abundance, y = herbivory_rate))+
  geom_jitter()+
  geom_line(data = pred2, aes(x = abundance, y = lm2), color = "red")+
  geom_line(data = pred2, aes(x = abundance, y = exp2), color = "blue")+
  geom_line(data = pred2, aes(x = abundance, y = sig2), color = "green")+
  geom_vline(xintercept = coef(sig2)[2], lty = "dashed")+ 
  ggpubr::theme_pubclean() 



# ------------------------------------------------------------------------------------------------
## Red analysis
# ------------------------------------------------------------------------------------------------ 
#here we are testing different models to determine the best fit model for the red ensity data

rf <- df[df$p_r == "r", ] #limiting our dataset to only purple urchins


lm1.r <- lm(herbivory_rate ~ 0 + biomass, rf)#linear model representing herbivory rate as a function of biomass.
summary(lm1.r)
modelassump(lm1.r)

exp1.r <- lm(herbivory_rate ~ biomass + I(biomass^2), rf)
summary(exp1) #exponential regression of herbivory rate as a function of biomass

pred3 <- data.frame(biomass = seq(min(rf$biomass), max(rf$biomass), length.out = 1000))
pred3$lm1.r <- predict(lm1.r, newdata = pred3) #linear model predictions.
pred3$exp1.r <- predict(exp1.r, newdata = pred3) #exponential model predictions

sig.r <- nls(herbivory_rate ~ (a * biomass^2) / (b^2 + biomass^2), data = rf, start = list(a = 10, b = 1000)) # this is an alternative parameterization based on Bolker et al. 2008, bestiary of functions p. 22. "a" is similar to handling time, and b is similar to attack rate in this parameterization.
summary(sig.r)

pow2.r <- nls(herbivory_rate ~ (a * biomass^b), data = rf, start = list(a = 10, b = 1)) 
summary(pow2.r)

pred3$sig.r <- predict(sig.r, newdata = pred3)

AIC(lm1.r, pow2.r, sig.r)
# So it seems that there is no evidence for any differences between linear and sigmoidal curves. 

model_compare3 <- ggplot(rf, aes(x = biomass, y = herbivory_rate))+
  geom_jitter(pch = 21, width =30)+
  geom_line(data = pred3, aes(x = biomass, y = lm), color = "red")+
  geom_line(data = pred3, aes(x = biomass, y = exp1), color = "blue")+
  geom_line(data = pred3, aes(x = biomass, y = sig), color = "green")+
  geom_vline(xintercept = coef(sig)[2], lty = "dashed")+
  ggpubr::theme_pubclean()

# ------------------------------------------------------------------------------------------------
## Figure 2: The relationship between biomass and herbivory rate
# ------------------------------------------------------------------------------------------------ 
pred <- data.frame(biomass = seq(min(pf$biomass), max(pf$biomass), length.out = 1000))
pred$herbivory_rate <- predict(lm1, newdata = pred, interval = "confidence")[,1]
pred$low <- predict(lm1, newdata = pred, interval = "confidence")[,2]
pred$high <- predict(lm1, newdata = pred, interval = "confidence")[,3]
pred$sp <- "Purple urchin"

pred.r <- data.frame(biomass = seq(min(rf$biomass), max(rf$biomass), length.out = 1000))
pred.r$herbivory_rate <- predict(lm1.r, newdata = pred, interval = "confidence")[,1]
pred.r$low <- predict(lm1.r, newdata = pred, interval = "confidence")[,2]
pred.r$high <- predict(lm1.r, newdata = pred, interval = "confidence")[,3]
pred.r$sp <- "Red urchin"

pred <- bind_rows(pred, pred.r)


fig2 <- ggplot(df, aes(x = biomass, y = herbivory_rate))+
  geom_rect(aes(xmin= 668, xmax=1246, ymin=-Inf, ymax=Inf), fill = "gray90", alpha = 0.1)+
  #geom_vline(xintercept = 668, linetype = 4)+
  geom_point(aes(fill = sp), pch = 21, size=3)+
  scale_fill_manual(values = c("#550f7a", "#E3493B"))+
  geom_line(data = pred, aes( x= biomass, y = herbivory_rate), size = 1, show.legend = F)+
  geom_ribbon(data = pred, aes( ymin = low, ymax = high,fill=sp), alpha = 0.3, show.legend = F)+
  facet_wrap(~sp)+
  theme_classic()+
  theme(strip.text = element_text(size = 10))+
  labs(x = expression(paste("Urchin biomass (g m"^"-2"*")")), y = expression(paste("Herbivory rate (g m"^"-2"*"d"^"-1"*")")), color = "", linetype = "")+
  theme(strip.background = element_blank())+
  theme(legend.position = c(.85,.90))+
  theme(legend.title=element_blank())+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )+
  theme(axis.title.x= element_text(color= "black", size=14),
        axis.title.y= element_text(color= "black", size=14))+
  theme(legend.background = element_blank(), legend.box.background = element_blank()  )+
  theme(axis.text = element_text(size = 12))+
  theme(legend.text=element_text(size=12))
  
fig2

ggsave("figures/herbivoryXdensity_fig2.png", fig2, device = "png",width=7,height=3.5)
ggsave("figures/herbivoryXdensity_fig2.pdf", fig2, device = "pdf", useDingbats = FALSE)

tiff(filename="figures/Fig2.tif",height=5600/2,width=5200,units="px",res=800,compression="lzw")
fig2
dev.off()


########################
## Summary stats 
########################

predict(lm1, newdata = list(biomass = 668), se.fit = T) #Using the model to predict herbivory rate at the transition density cited in Ling et al. 2016

#------------------------------------------------
# Supplemental figure 1
#-----------------------------------------------

pf$con_per_g_biomass <- pf$herbivory_rate / pf$biomass #adding con_per_g_biomass to the purple urchin dataset

s1 <- lm(con_per_g_biomass ~ biomass, pf) #linear model of per capita consumption for purple urchins. what is con? concentration? concentration of what? This is percapita consumption? 
summary(s1)
modelassump(s1) #Is this heterodacictic too? 


rf$con_per_g_biomass <- rf$herbivory_rate / rf$biomass #adding con_per_g_biomass to the red urchin dataset


s2 <- lm(con_per_g_biomass ~ biomass, rf) #linear model of per capita consumption of red urchins.
summary(s2)
modelassump(s2) #This looks less conclusive. P value over .05 meaning we cannot confirm that the slope is not zero therefore there is a positve correlation between biomass and consumption? 

plot(con_per_g_biomass ~ biomass, pf)
plot(con_per_g_biomass ~ biomass, rf) #This one looks really spread out

gg <- data.frame(biomass = seq(0, max(df$biomass, na.rm = T), length.out = 1000)) #biomass values from the dataset
gg$`Purple urchin` <- predict(s1, newdata = gg)
gg$`Red urchin` <- predict(s2, newdata = gg) #merging the linear model predictions of per capita consumption for purple and red urchins in the gg dataset.

gg <- gg %>% gather(sp, prediction, -biomass) #merging the purple and red predictions

df$con_per_g_biomass <- df$herbivory_rate / df$biomass

S1 <- ggplot(df, aes(x = biomass, y = con_per_g_biomass))+
  geom_jitter(aes(fill = sp), pch = 21, show.legend = F)+
  scale_fill_manual(values = c("#762a83", "#d73027"))+
  facet_wrap(~sp, scales = "free")+
  scale_x_continuous(limits = c(0, 2000))+
  scale_y_continuous(limits = c(-0.02, 0.06))+
  geom_hline(yintercept = 0, lty = "dashed", color = "gray")+
  labs(x = expression(paste("Urchin biomass density (g m"^"-2"*")")), 
       y = expression(paste("Herbivory rate (g"["kelp"]*"g"["urc"]^"-1"*"m"^"-2"*"d"^"-1"*")")), 
       color = "", linetype = "")+
  theme_classic()+
  theme(strip.text = element_text(size = 10))+
  theme(strip.background = element_blank())+
  theme(legend.position = c(.90,.90))+
  theme(legend.title=element_blank())+
  theme(axis.title.x= element_text(color= "black", size=20),
        axis.title.y= element_text(color= "black", size=20))+
  theme(legend.text=element_text(size=10))+
  theme(legend.background = element_rect( 
    size=0.5, linetype ="solid"))+
  theme(axis.text = element_text(size = 15))+
  theme(legend.text=element_text(size=15))

ggsave(here("figures", "percapconsumptionxbiomass.png"), S1, device = "png", width = 8.5, height = 5)


#----------------------------------------
# deltaAIC table
#----------------------------------------

pl <- c(
  `(Intercept)` = "Intercept"
)

sjPlot::tab_model(lm1, pow2, sig, 
                  show.aic = T, 
                  show.icc = F, 
                  show.loglik = T, 
                  show.ngroups = F, 
                  pred.labels = pl, 
                  dv.labels = c("Linear", "Power-law", "Sigmoid"), file = here::here("figures/", "AICtablePurple.html"))

sjPlot::tab_model(lm1.r, pow2.r, sig.r, 
                  show.aic = T, 
                  show.icc = F, 
                  show.loglik = T, 
                  show.ngroups = F, 
                  pred.labels = pl, 
                  dv.labels = c("Linear", "Power-law", "Sigmoid"), file = here::here("figures/", "AICtableRed.html"))




#-----------------------------------------------------------------------------------------
## Test the effects of the tank dividers
#-----------------------------------------------------------------------------------------

bd <- df %>% separate(tank, into = c("tank", "side"), sep = "[-]")

aov <- aov(herbivory_rate ~ side, bd)
summary(aov)

TukeyHSD(aov)

out <- list()
tanks <- unique(bd$tank)
tanks <- tanks[tanks != 9]
for(i in 1:length(unique(bd$tank))){
  out[[i]] <- summary(aov(herbivory_rate ~ side, bd[bd$tank == tanks[i], ]))
}

table(bd$tank, bd$side)

ggplot(bd, aes(x = tank, group = side, y = herbivory_rate))+
  geom_bar(aes(fill = side), stat = "identity", position = "dodge")
