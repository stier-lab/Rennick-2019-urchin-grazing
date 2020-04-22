Mae Rennick
Urchin Density/Herbivory Rate Data 
07/10/2019

urchin_density_data_2 <-read.csv("urchin_density_data.csv", header= TRUE)

red_urchins_2 <- urchin_density_data_2 %>% 
  filter(trial_number == "1-R" | trial_number == "R-1")

ggplot(red_urchins_2) +
  geom_point(aes(urchin_density, kelp_consumed))+
  geom_smooth(aes(urchin_density, kelp_consumed), method="lm")

ggplot(red_urchins_2)+
  geom_point(aes(tank_biomass, kelp_consumed))+
  geom_smooth(aes(tank_biomass, kelp_consumed), method="lm")
