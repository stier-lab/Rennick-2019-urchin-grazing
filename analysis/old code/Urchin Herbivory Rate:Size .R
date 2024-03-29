library(here)

consumption_data <- read.csv(here("data/urchin_kelp_consumption_data 2019_5_20.csv"))
size_data <- read.csv(here("data/urchin_size_data 2019_5_20.csv"))
library(esquisse)

grouped_trial_id <- size_data %>% 
  group_by(trial_id) %>% 
  summarize(avg_mass= mean(weight), avg_size= (mean(size))) %>%
  signif(digits=3)

combined_data <- merge(grouped_trial_id, consumption_data, by="trial_id")

trial_time <- combined_data %>% 
  mutate(trial_time=time_in- time_out +48)

combined_data_2 <- combined_data %>% 
  filter(use=="1")

herbivory_rate <- combined_data_2 %>% 
  mutate(herbivory_rates = kelp_consumed/336)

esquisse::esquisser(data = herbivory_rate)
library(ggplot2)

ggplot(data = herbivory_rate) +
  aes(x = avg_size, y = avg_mass) +
  geom_point(color = '#0c4c8a') +
  geom_smooth(span = 0.75) +
  labs(title = 'Relationship Between Urchin Size and Mass ',
    x = 'Urchin Size (mm)',
    y = 'Urchin Mass (g)') +
  theme_minimal()
library(ggplot2)

ggplot(data = herbivory_rate) +
  aes(x = avg_size, y = herbivory_rates) +
  geom_point(color = '#0c4c8a') +
  geom_smooth(span = 0.5) +
  labs(title = 'Size Dependent Herbivory Rates of Purple Urchins on Giant Kelp',
    x = 'Urchin Size (mm)',
    y = 'Herbivory Rate (kelp consumed/hour)') +
  theme_stata()

esquisse::esquisser(data = herbivory_rate)
library(ggplot2)

ggplot(data = herbivory_rate) +
  aes(x = avg_mass, y = herbivory_rates) +
  geom_point(color = '#0c4c8a') +
  geom_smooth(method = "lm")+
  labs(title = 'Mass Dependent Herbivory Rates of Pruple Urchins on Giant Kelp ',
    x = 'Urchin Mass (g)',
    y = 'Herbivory Rate (kelp consumed (g)/ hour)') +
  theme_minimal()

read.csv()



