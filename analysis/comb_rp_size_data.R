## Mae Rennick
## Urchin Size Data

df_p1 <- read.csv("data/size_experiment/raw/Purple_size_data.csv") #raw size and mass data from the size trials

df_p1 <- df_p1 %>% 
  mutate(urchin = "p") 


str(df_p)
head(df_p)

df_p2 <- read.csv("data/size_experiment/raw/Purple_consuption_data.csv") #raw experimental data from the herbivory trials (purple)

df_p3 <-df_p2 %>% #combined data set for purple urchins
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

df_rp <- rbind(df_r2,df_p3) #combined dataset for red and purple urchins size trials
