#Mae
#Red urchin size 

red_size <- read.csv("data/red_size_data.csv")

red_size$size <- as.numeric(as.character(red_size$size))

grouped_rs <- red_size %>% 
  group_by(avg_size)

ggplot(data = grouped_rs) +
  aes(x = avg_size, y =kelp_consumed ) +
  geom_point(color = '#0c4c8a') +
  geom_smooth() +
  labs(title = 'Red Urchin Kelp Consumption as a Function of Size',
       x = 'Urchin Size (mm)',
       y = 'Kelp Consuption (g)') +
  theme_minimal()

#######################
## Per unit size plot
#######################

red_size$per_size <- (red_size$kelp_consumed/5)/red_size$avg_size

red_size %>% 
  group_by(avg_size) %>%
  mutate(per_size_tidy= kelp_consumed/5/avg_size) %>% 
  ggplot() +
  aes(x = avg_size, y =per_size_tidy ) +
  geom_point(color = '#0c4c8a') +
  geom_smooth() +
  labs(title = 'Per Unit Size Consumption vs Size',
       x = 'Urchin Size (mm)',
       y = 'Kelp Consuption (g)') +
  theme_minimal()


#####################
#herbivory rate plot

red_size %>% 
  group_by(avg_size) %>%
  mutate(per_size_tidy= kelp_consumed/5/avg_size) %>% 
  ggplot() +
  aes(x = avg_size, y =herb_rate ) +
  geom_point(color = '#0c4c8a') +
  geom_smooth() +
  labs(title = 'Per Unit Size Consumption vs Size',
       x = 'Urchin Size (mm)',
       y = 'herbivory rate') +
  theme_minimal()


 library(ggplot2)