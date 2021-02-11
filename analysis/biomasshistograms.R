library(tidyverse)
library(scales)
library(wesanderson)

pal <- wes_palette("Darjeeling1", 9, type = "continuous")


df <- read.csv("data/survey_data/Annual_All_Species_Biomass_at_transect_20210108.csv", stringsAsFactors = F,na.strings ="-99999") %>% #LTER dataset: urchin desnity data collected from 50 transects across 11 sites between 2000-2018 in the Santa Barbara Channel. (what are strings?)
  filter(SP_CODE == "SPL" | SP_CODE == "SFL" ) %>%
  dplyr::select("YEAR", "MONTH", "SITE", "TRANSECT", "SP_CODE", "WM_GM2", "DENSITY") %>%
  filter(SITE != "SCTW", SITE != "SCDI") %>% #Removing two island sites. We are only working with coastal sites. 
  group_by(YEAR, MONTH, SITE, TRANSECT) %>%
  rename_all(tolower) %>%
  mutate(biomass = wm_gm2, 
         sp_code = case_when(sp_code == "SPL" ~ "Purple urchin",
                   sp_code == "SFL" ~ "Red urchin"))


df$site<-as.factor(df$site)

# for size data:

sz <- read.csv("data/survey_data/Annual_Quad_Swath_All_Years_20200108.csv", stringsAsFactors = F,na.strings ="-99999") %>% #LTER dataset: urchin desnity data collected from 50 transects across 11 sites between 2000-2018 in the Santa Barbara Channel. (what are strings?)
  filter(SP_CODE == "SPL" | SP_CODE == "SFL" ) %>%
  filter(SITE != "SCTW", SITE != "SCDI") %>%
  dplyr::select("YEAR", "MONTH", "SITE", "TRANSECT", "QUAD", "SIDE", "SIZE", "SP_CODE") %>%
  rename_all(tolower)

# citation: Reed, D. 2020. SBC LTER: Reef: Kelp Forest Community Dynamics: Invertebrate and algal density ver 26. Environmental Data Initiative. https://doi.org/10.6073/pasta/cd4009368880fcc32d67ab0542e6b048. Accessed 2020-12-15.

#------------------------------------------------------------------
## Stats for text
#------------------------------------------------------------------

# red urchin vs purple urchin biomass

(aov1 <- aov(biomass ~ sp_code, df))
summary(aov1)

TukeyHSD(aov1)

# red urchin vs purple urchin density

(aov2 <- aov(density ~ sp_code, df))
summary(aov2)

TukeyHSD(aov2)

# red urchin vs. purple urchin size

(aov3 <- aov(size ~ sp_code, sz))
summary(aov3)
TukeyHSD(aov3)


#---------------------------------------------------------------------
## Plotting
#---------------------------------------------------------------------


p1 <- ggplot(df, aes( x= sp_code, y= biomass))+
  geom_point(alpha = 0, show.legend = F)+
  geom_rect(aes(ymin= 668 - 115, ymax=668 + 115, xmin=0, xmax=Inf), fill = "gray90", alpha = 1/60)+ # 688gm^-2 is the transition denisty cited in Ling et. al 2016 as the biomass of urchins required to incite a forward transition from a kelp dominated to an urhcin domianted state, with an error range of plus or minus 115gm^-2.
  geom_hline(yintercept = 668, linetype = 4)+ 
  geom_boxplot(aes(fill = sp_code), alpha = 0.5, show.legend = F, outlier.shape = NA)+
  geom_jitter(aes(color = sp_code,alpha=0.5), width = 0.3, show.legend = F)+
  scale_color_manual(values = c("#550f7a", "#E3493B"))+
  scale_fill_manual(values = c("#550f7a", "#E3493B"))+
  scale_y_log10(breaks = c(1, 10, 100, 1000), labels = c(1, 10, 100, 1000))+
  labs(y = expression(paste("Urchin biomass (g m"^"-2"*")")), x = "")+
  coord_cartesian(ylim = c(0.8, 3000))+
  theme_classic()+
  theme(axis.text=element_text(size=12))+
  theme(legend.title=element_blank())+
  theme(legend.position = c(0.2,0.85))+
  theme(axis.title.x= element_text(color= "black", size=16),
        axis.title.y= element_text(color= "black", size=16))+
  theme(legend.text=element_text(size=10))+
  expand_limits(y = 0.006)+
  theme(axis.text = element_text(size = 12))


p2 <- ggplot(df, aes(x = biomass))+
  geom_density(aes(fill = sp_code), alpha = 0.75, adjust = 2)+
  scale_fill_manual(values = c("#550f7a", "#E3493B"))+
  scale_x_log10()+
  labs(x = expression(paste("Urchin Biomass (g m"^"-2"*")")), y = "Density" )+
  theme_classic()+
  theme(axis.text=element_text(size=12))+
  theme(legend.title=element_blank())+
  theme(legend.position = c(0.2,0.85))+
  theme(axis.title.x= element_text(color= "black", size=20),
        axis.title.y= element_text(color= "black", size=20))+
  theme(legend.text=element_text(size=10))+
  expand_limits(y = 0.006)+
  theme(axis.text = element_text(size = 15))

show_col(colorblind_pal()(8))

p2_alt <- df %>% #LTER data 
  group_by(year, site, transect) %>%
  summarize(biomass = sum(biomass, na.rm = T)) %>% #adding a biomass column
  group_by(year, site) %>%
  summarize(biomass = mean(biomass)) %>%
  ggplot(aes(x = year, y = biomass))+
  geom_line(aes(color = site,lty=site), size=.75,show.legend = F)+
  scale_colour_manual(values=pal)+
  # scale_color_manual(values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6'))+
  stat_summary(fun = mean, geom = "line", lwd = 2, alpha = 0.75)+
  labs(y = "", x = "")+
  theme_classic()+
  theme(axis.text=element_text(size=12))+
  theme(legend.title=element_blank())+
  theme(legend.position = c(0.2,0.85))+
  theme(axis.title.x= element_text(color= "black", size=16),
        axis.title.y= element_text(color= "black", size=16))+
  theme(legend.text=element_text(size=10))+
  expand_limits(y = 0.006)


p2_alt

fig1 <- cowplot::plot_grid(p1,p2_alt, align = "h")
ggsave("figures/biomasshisto.png", fig1, device = "png", width = 8.5, height = 5)
ggsave("figures/biomasshisto.pdf", fig1, device = "pdf", width = 7, height = 3, useDingbats = FALSE)



