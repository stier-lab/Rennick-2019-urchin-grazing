## MAE Density and Size Herbivory Figures 
## 10-08-19

###############DensityTrials 

df_d<-read.csv("data/density_experiment/derived/urchin_density_data_cleaned.csv") %>% 
  mutate(kelp_consumed = kelp_in_total-Kelp_out_total) %>% 
  mutate (herbivory_rate = kelp_consumed/total_time) %>% 
  group_by(trial_id) %>%
  mutate(total_mass= sum(urchin_mass)) %>% 
  ungroup()
  

############## Red Denisty 

df_rdd <- df_d %>%
  filter(type == 'r')
df_rd <- within(df_rdd, herbivory_rate[herbivory_rate<0] <- 0)


plot(df_rd$herbivory_rate ~ df_rd$urchin_density)
mod <- lm(herbivory_rate ~ urchin_density, data = df_rd)
summary(mod)
plot(herbivory_rate ~ urchin_density, df_rd)
abline(mod)


############## Red Density Plot
    ##NUMBER

ggplot(df_rd, aes(urchin_density,herbivory_rate))+
  theme_classic()+
  geom_point(color="black")+
  geom_jitter(color="black")+
  geom_smooth(method="lm", color="orange")+
  ggtitle("Density Dependent Herbivory Rates of Red Urchins")+
  xlab("Urchin Density")+
  ylab("Kelp Consumed (g) /time")

    ## BIOMASS

ggplot(df_rd, aes(total_mass,herbivory_rate))+
  theme_classic()+
  geom_point(color="black")+
  geom_jitter(color="black")+
  geom_smooth(method="lm", color="orange")+
  ggtitle("Density Dependent Herbivory Rates of Red Urchins")+
  xlab("Urchin Biomass (g)")+
  ylab("Kelp Consumed (g) /time")

############## Purple Density

df_pd <- df_d %>% 
  filter(type == 'p')

plot(df_pd$herbivory_rate ~ df_pd$urchin_density)
mod <- lm(herbivory_rate ~ urchin_density, data = df_pd)
summary(mod)
plot(herbivory_rate ~ urchin_density, df_pd)
abline(mod)


############# Purple Density Plot

    ####NUMBER

ggplot(df_pd, aes(urchin_density,herbivory_rate))+
  theme_classic()+
  geom_point(color="black")+
  geom_jitter(color="black")+
  geom_smooth(method="lm", color="orange")+
  ggtitle("Density Dependent Herbivory Rates of Purple Urchins")+
  xlab("Urchin Density")+
  ylab("Kelp Consumed (g) /time")


    ####BIOMASS

ggplot(df_pd, aes(total_mass,herbivory_rate))+
  theme_classic()+
  geom_point(color="black")+
  geom_jitter(color="black")+
  geom_smooth(method="lm", color="orange")+
  ggtitle("Density Dependent Herbivory Rates of Purple Urchins")+
  xlab("Urchin Biomass (g)")+
  ylab("Kelp Consumed (g) /time")

################# R and P Density Plot

ggplot(df_d, aes(urchin_density,herbivory_rate, color=type))+
  theme_classic()+
  geom_point()+
  geom_jitter()+
  scale_color_manual(breaks= c("r", "p"), values= c("orange", "black"), labels=c("Red", "Purple"))+
  geom_smooth(method="lm", color="black")+
  ggtitle("Density Dependent Herbivory Rates of Urchins")+
  xlab("Urchin Density")+
  ylab("Kelp Consumed (g) /time")+
  theme(legend.position="bottom")+
  theme(legend.title = element_blank())

plot(df_d$herbivory_rate ~ df_d$urchin_density)
mod_d <- lm(herbivory_rate ~ urchin_density, data = df_d)
summary(mod_d)


############### Size Trials 

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

df_rs_2 <-read.csv("data/size_experiment/raw/red_size_data_round2-3.csv")

df_rs_3 <- df_rs_2 %>%
  rename(kelp_in = Kelp_in) %>% 
  rename (kelp_out =Kelp_out)

df_rp_size_2 <- rbind(df_rs_3, df_rp_size)

df_s <- df_rp_size_2 %>% 
  group_by(trial_id) %>%
  mutate(avg_size = mean(size)) %>%
  mutate(avg_mass=mean(weight)) %>% 
  ungroup() %>% 
  filter(avg_size != "NA") %>% 
  mutate(kelp_consumed = kelp_in-kelp_out) %>% 
  mutate (herbivory_rate = kelp_consumed/time_ran)

######## Red Size

df_rs_4 <- df_s %>% 
  filter(urchin == 'r') %>% 
  filter(time_ran=='96')

plot(df_rs_4$herbivory_rate ~ df_rs_4$avg_size)
mod_rs <- lm(herbivory_rate ~ avg_size, data = df_rs_4)
summary(mod_rs)
plot(herbivory_rate ~ avg_size, df_rs_4, main="Size Dependent Herbivory Rates of Red Urchins", xlab="Average Size (mm)", ylab="Herbivory Rate (Kelp Consumed (g)/ Time)")
abline(mod_rs)


############## Red Size Plot

###### SIZE 

ggplot(df_rs_4, aes(avg_size,herbivory_rate))+
  theme_classic()+
  geom_point(color="black")+
  geom_jitter(color="black")+
  ggtitle("Size Dependent Herbivory Rates of Red Urchins")+
  xlab("Urchin Size (mm)")+
  ylab("Kelp Consumed (g) /time")


### BIOMASS

bms_rs <- df_s %>% 
  filter(urchin == 'r') %>% 
  filter(avg_mass != "NA")

ggplot(bms_rs, aes(avg_mass, herbivory_rate))+
  theme_classic()+
  geom_point(color="black")+
  geom_jitter(color="black")+
  geom_smooth(method="lm", color="dark green")+
  ggtitle("Size Dependent Herbivory Rates of Red Urchins")+
  xlab("Average Mass (g)")+
  ylab("Kelp Consumed (g) /time")


######## Purple Size 

df_ps<- df_s %>% 
  filter(urchin == 'p')

plot(df_ps$herbivory_rate ~ df_ps$avg_size, main="Size Dependent Herbivory Rates of Purple Urchins", xlab="Average Size (mm)", ylab="Herbivory Rate (Kelp Consumed (g)/ Time)")

mod_ps <- lm(herbivory_rate ~ avg_size, data = df_ps)
summary(mod_ps)
plot(herbivory_rate ~ avg_size, df_ps, main="Size Dependent Herbivory Rates of Purple Urchins", xlab="Average Size (mm)", ylab="Herbivory Rate (Kelp Consumed (g)/ Time)")
abline(mod_ps)

################ Purple Size Plot

###SIZE (mm)

ggplot(df_ps, aes(avg_size,herbivory_rate))+
  theme_classic()+
  geom_point(color="black")+
  geom_jitter(color="black")+
  geom_smooth(method="lm", color="dark green")+
  ggtitle("Size Dependent Herbivory Rates of Purple Urchins")+
  xlab("Urchin Size (mm)")+
  ylab("Kelp Consumed (g) /time")

### BIOMASS

bms_ps <- df_s %>% 
  filter(urchin == 'p') %>% 
  filter(avg_mass != "NA")

ggplot(bms_ps, aes(avg_mass, herbivory_rate))+
  theme_classic()+
  geom_point(color="black")+
  geom_jitter(color="black")+
  geom_smooth(method="lm", color="dark green")+
  ggtitle("Size Dependent Herbivory Rates of Purple Urchins")+
  xlab("Average Mass (g)")+
  ylab("Kelp Consumed (g) /time")

####### Red and Purple Size 

df_rps <-df_s %>% 
  filter(time_ran=="48")
mod_s <- lm(herbivory_rate ~ avg_size, data = df_rps) 
summary(mod_s)
plot(herbivory_rate ~ avg_size, df_rps, main="Size Dependent Herbivory Rates of Urchins", xlab="Average Size (mm)", ylab="Herbivory Rate (Kelp Consumed (g)/ Time)")
abline(mod_s)

ggplot(df_rps, aes(avg_size,herbivory_rate, color=urchin))+
  theme_classic()+
  geom_point()+
  geom_jitter()+
  scale_color_manual(breaks= c("r", "p"), values= c("orange", "black"), labels=c("Red", "Purple"))+
  geom_smooth(method="lm", color="black",data=df_rps, aes(avg_size,herbivory_rate))+
  ggtitle("Size Dependent Herbivory Rates of Urchins")+
  xlab("Average Size (mm)")+
  ylab("Kelp Consumed (g) /time")+
  theme(legend.position="bottom")+
  theme(legend.title = element_blank())

