# Purple Urchin Size/ Herbivory Rate
#Mae
#08/06/2019

df<-read.csv("data/urchin_kelp_consumption_data 2019_5_20.csv") %>% 
  select(-notes) %>%
  mutate(urchin_size_cat = as.numeric(as.character(urchin_size_cat)))

temp <- df %>%
  group_by(urchin_size_cat, week_no, tank, side) %>%
  mutate(percap = kelp_consumed/urchin_dens/urchin_size_cat)

ggplot(data = temp, aes(x = urchin_size_cat, y = percap))+ 
  geom_point()+
  geom_smooth()

ggplot(data= df, aes(x=urchin_size_cat, y=kelp_consumed)) +
  geom_point()+
  geom_smooth()
