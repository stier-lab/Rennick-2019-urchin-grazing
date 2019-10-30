## Mae Rennick
## Urchin Size Data

df_p1 <- read.csv("data/size_experiment/raw/Purple_size_data.csv")

df_p1 <- df_p1 %>% 
  mutate(urchin = "p")


str(df_p)
head(df_p)

df_p2 <- read.csv("data/size_experiment/raw/Purple_consuption_data.csv")

df_p3 <-df_p2 %>%
  filter(urchin_size_cat != 'control') %>% 
  mutate(urchin = "p")%>% 
  left_join(df_p1) %>% 
  select(-date) %>% 
  select (-days_starved) %>% 
  select (-time_in) %>% 
  select (-time_out) %>% 
  select (-kelp_per_urchin) %>%
  select (-kelp_consumed) %>% 
  mutate(level= "B") %>% 
  mutate (time_ran = 48) %>% 
  mutate( round= week_no) %>% 
  select(-week_no) %>% 
  mutate(size_class= urchin_size_cat) %>% 
  select (-urchin_size_cat) %>% 
  mutate (urchin_density =urchin_dens) %>% 
  select (-urchin_dens)

df_r <- read.csv("data/size_experiment/raw/red_size_data_2.csv")

df_r2 <-df_r %>% 
  mutate(round =1) %>% 
  mutate(urchin_density=5) %>% 
  select(-rep) %>% 
  mutate (weight= NA) %>% 
  mutate (urchin= "r") %>% 
  mutate (use=1)

df_rp_size <- rbind(df_r2,df_p3)

df_rs_2 <-read.csv("data/size_experiment/raw/red_size_data_round2.csv")

df_rs_3 <- df_rs_2 %>%
  rename(kelp_in = Kelp_in) %>% 
  rename (kelp_out =Kelp_out)

df_rp_size_2 <- rbind(df_rs_3, df_rp_size)

#################################

df_rp_size_plot <- df_rp_size_2 %>% 
  filter(use!= 0) %>% 
  filter(time_ran != 48) %>%
  filter(urchin != "p") %>% 
  group_by(trial_id) %>% 
  mutate(avg_size= mean(size)) %>% 
  mutate (kelp_consumed =kelp_in -kelp_out)

ggplot(df_rp_size_plot , aes(avg_size, kelp_consumed)) +
  geom_point() +
  geom_smooth()
