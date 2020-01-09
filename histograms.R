####################################################################################
## Size Frequency Distributions
####################################################################################

library(tidyverse)
library(ggpubr)


# Get size data

urc <- read.csv("data/survey_data/LTE_Urchin_All_Years_20190611.csv", header = T) %>% 
  filter(TREATMENT == "CONTROL") %>% select(YEAR, MONTH, DATE, SITE, TRANSECT, SIZE, COUNT, COMMON_NAME) %>% rename_all(tolower) %>% group_by(year, month, date, site, transect, common_name, size) %>%
  summarize(count = sum(count)) %>%
  group_by(year, month, date, site, transect, common_name, size) %>%
  complete(count = full_seq(1:count,1))%>%
  select(-count) %>%
  ungroup() %>%
  mutate(size = as.numeric(size))

dh <- ggplot(urc, aes(x = size))+
  geom_density(aes(fill = common_name), alpha = 0.75, adjust = 2)+
  scale_fill_manual(values = c("#762a83", "#d73027"))+
  scale_x_continuous(breaks = seq(4,14, by = 2), labels = seq(4,14,by=2))+
  labs(x = "Urchin size (cm)", y = "Density", fill = "")+
  theme_pubclean()+
  theme(legend.position = "right")

ggsave("figures/urchinsize_densityhisto.png", dh, device = "png")



# Get density data

lt <- read.csv("data/survey_data/Annual_All_Species_Biomass_at_transect.csv", stringsAsFactors = F,na.strings ="-99999") %>%
  dplyr::select("YEAR", "MONTH", "SITE", "TRANSECT", "SP_CODE", "WM_GM2") %>%
  filter(SP_CODE == "SPL" | SP_CODE == "SFL") %>%
  filter(SITE != "SCTW", SITE != "SCDI") %>%
  group_by(YEAR, MONTH, SITE, TRANSECT) %>%
  spread(SP_CODE, WM_GM2) %>%
  rename_all(tolower) %>%
  mutate(urc.biomass = spl + sfl)

dh.d <- ggplot(lt, aes(x = urc.biomass))+
  geom_density(fill = NA, color = NA )+
  geom_rect(aes(xmin= 668 - 115, xmax=668 + 115, ymin=0, ymax=Inf), color = "gray90", fill = "gray90")+
  geom_vline(xintercept = 668, linetype = 4)+
  geom_density(fill = "#a8325e", alpha = 0.8)+
  labs(x = expression(paste("Combined urchin biomass (g m"^"-2"*")")), y = "Density" )+
  theme_pubclean()

ggsave("figures/urchindensity_densityhisto.png", dh.d, device = "png")

p1 <- cowplot::plot_grid(dh.d, dh+labs(y = ""), align = "h", rel_widths = c(0.75, 1 ))

ggsave("figures/urchinhistos.png", p1, device = "png", width = 10, height = 4)

#----------------------------------
## Summary stats for paper
#----------------------------------

# mean and range of combined urchin biomass 

mean(lt$urc.biomass, na.rm = T)
sd(lt$urc.biomass)
sd(lt$urc.biomass, na.rm = T) / sqrt(length(lt$urc.biomass)) # se
range(lt$urc.biomass)

# purple sea urchins composed XX  x % of urchin biomass

lt %>% ungroup() %>% summarize(sum(spl)/sum(urc.biomass))

# density comparison
read.csv("data/survey_data/Annual_All_Species_Biomass_at_transect.csv", stringsAsFactors = F,na.strings ="-99999") %>%
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




























