####################################################################################
## Size and Density Frequency Distributions
####################################################################################

library(tidyverse)
library(ggpubr)

# Here we will be creating histograms that explore the natural distribution of size and abundance of purple and red urhcins across (5????) LTER coastal sites in the Santa Barbara Channel. These realtionships will be used to create size bins for experimentaion and exemplify the natural distribution seen in our density focused herbivory trials. Becuase the abundance and distribution of herbivorous urchins will no doubt have an affect on kelp biomass, it is important to understand and map these distributions in these coastal sites where we will later be testing our experimental models. We will also be looking at our local density in realtionship to the transition density cited in Ling et al. 2015 that predicts the urhcin density required to incite a forward transition from a kelp dominated state to an urhcin dominated state. 



# Get size data

urc <- read.csv("data/survey_data/LTE_Urchin_All_Years_20190611.csv", header = T) %>% # LTER dataset: urchin size and frequency data collected from 5 sites between 2008-2019 in the Santa Barbara Channel. #They only take this data once at each site? 
  filter(TREATMENT == "CONTROL") %>% select(YEAR, MONTH, DATE, SITE, TRANSECT, SIZE, COUNT, COMMON_NAME) %>% 
  rename_all(tolower) %>% 
  group_by(year, month, date, site, transect, common_name, size) %>% 
  summarize(count = sum(count)) %>% #What does this change other than moving it? 
  group_by(year, month, date, site, transect, common_name, size) %>%
  complete(count = full_seq(1:count,1))%>% #transferring each count to one, and seperating each urchin into its own line in the data. 
  select(-count) %>% #removing the count column
  ungroup() %>%
  mutate(size = as.numeric(size)) #creating a size column and establishing it as a numeric value?? 


dh <- ggplot(urc, aes(x = size))+
  geom_density(aes(fill = common_name), alpha = 0.75, adjust = 2)+
  scale_fill_manual(values = c("#762a83", "#d73027"))+
  scale_x_continuous(breaks = seq(4,14, by = 2), labels = seq(4,14,by=2))+
  labs(x = "Urchin size (cm)", y = "Density", fill = "")+
  theme_pubclean()+
  theme(legend.position = "right") #histogram depicting the relationship between density and urchin sizes across sites, seperately in red and urhcin populations. 

#more purple than red urchins but red urchins get much bigger than purple urchins.

ggsave("figures/urchinsize_densityhisto.png", dh, device = "png")



# Get density data

lt <- read.csv("data/survey_data/Annual_All_Species_Biomass_at_transect.csv", stringsAsFactors = F,na.strings ="-99999") %>% #LTER dataset: urchin desnity data collected from 50 transects across 11 sites between 2000-2018 in the Santa Barbara Channel. (what are strings?)
  dplyr::select("YEAR", "MONTH", "SITE", "TRANSECT", "SP_CODE", "WM_GM2") %>%
  filter(SP_CODE == "SPL" | SP_CODE == "SFL") %>% #filtering for only pruple and red urchins
  filter(SITE != "SCTW", SITE != "SCDI") %>% #Removing two island sites. We are only working with coastal sites. 
  group_by(YEAR, MONTH, SITE, TRANSECT) %>%
  spread(SP_CODE, WM_GM2) %>% #dividing purple and red urhcin densities into their own respective columns.
  rename_all(tolower) %>%
  mutate(urc.biomass = spl + sfl) #red+purple urchin densities combined into a total urchin biomass column.

dh.d <- ggplot(lt, aes(x = urc.biomass))+
  geom_density(fill = NA, color = NA )+
  geom_rect(aes(xmin= 668 - 115, xmax=668 + 115, ymin=0, ymax=Inf), color = "gray90", fill = "gray90")+ # 688gm^-2 is the transition denisty cited in Ling et. al 2016 as the biomass of urchins required to incite a forward transition from a kelp dominated to an urhcin domianted state, with an error range of plus or minus 115gm^-2.
  geom_vline(xintercept = 668, linetype = 4)+ 
  geom_density(fill = "#a8325e", alpha = 0.8)+
  labs(x = expression(paste("Combined urchin biomass (g m"^"-2"*")")), y = "Density" )+
  theme_pubclean() #histogram looking at the relationship of urchin density as a product of combined urchin biomass. 

ggsave("figures/urchindensity_densityhisto.png", dh.d, device = "png")

p1 <- cowplot::plot_grid(dh.d, dh+labs(y = ""), align = "h", rel_widths = c(0.75, 1 )) #combining histograms to display the relationship between red and purple urchin size (individually) and density, in addition to combined urhcin biomass and density. 

ggsave("figures/urchinhistos.png", p1, device = "png", width = 10, height = 4)

#----------------------------------
## Summary stats for paper
#----------------------------------

# mean and range of combined urchin biomass 

mean(lt$urc.biomass, na.rm = T)
sd(lt$urc.biomass)
sd(lt$urc.biomass, na.rm = T) / sqrt(length(lt$urc.biomass)) 
range(lt$urc.biomass)

# purple sea urchins composed XX  x % of urchin biomass

lt %>% ungroup() %>% summarize(sum(spl)/sum(urc.biomass)) #Purple urchins make up 63.1% of the total biomass(How do you know the error range??)

# density comparison
read.csv("data/survey_data/Annual_All_Species_Biomass_at_transect.csv", stringsAsFactors = F,na.strings ="-99999") %>% #LTER dataset: urchin desnity data collected from 50 transects across 11 sites between 2000-2018 in the Santa Barbara Channel. 
  dplyr::select("YEAR", "MONTH", "SITE", "TRANSECT", "SP_CODE", "DENSITY") %>%
  filter(SP_CODE == "SPL" | SP_CODE == "SFL") %>%
  filter(SITE != "SCTW", SITE != "SCDI") %>%
  group_by(YEAR, MONTH, SITE, TRANSECT) %>%
  rename_all(tolower) %>%
  ungroup() %>%
  group_by(sp_code) %>%
  summarize(mean(density), 
            sd(density))


# means of urchin size by species

urc %>% group_by(common_name) %>%
  summarize(mean(size), 
            sd(size))
