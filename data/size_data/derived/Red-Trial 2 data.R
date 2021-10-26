##Mae Rennick
## Red Urchin Size trials round 2

df_r2 <-read.csv("data/size_experiment/raw/red_size_data_round2.csv")

df_r2_2 <-df_r2 %>%
  filter(use!= 0) %>% 
  filter(time_ran != 48) %>% 
  group_by(trial_id) %>% 
  mutate(avg_size= mean(size)) %>% 
  mutate(avg_weight =mean(weight)) %>% 
  mutate (kelp_consumed =Kelp_in -Kelp_out)

ggplot(df_r2_2, aes(avg_size, kelp_consumed)) +
  geom_point() +
  geom_smooth()

ggplot(df_r2_2, aes(avg_weight, kelp_consumed)) +
  geom_point() +
  geom_smooth()
