####################################################################################
## Size Frequency Distributions
####################################################################################

library(tidyverse)
library(ggpubr)


# Get size data

urc <- read.csv("data/survey_data/LTE_Urchin_All_Years_20190611.csv", header = T) %>% 
  filter(TREATMENT == "CONTROL") %>% select(YEAR, MONTH, DATE, SITE, TRANSECT, SIZE, COUNT, COMMON_NAME) 

names(urc) <- tolower(names(urc))

urc <- urc %>% group_by(year, month, date, site, transect, common_name, size) %>%
  summarize(count = sum(count)) %>%
  spread(size, count) %>%
  gather(size, count, -c(year, month, date, site, transect, common_name)) %>%
  group_by(year, month, date, site, transect, common_name, size) %>%
  summarize(count = sum(count)) %>%
  group_by(year, month, date, site, transect, common_name, size) %>%
  drop_na(count)%>%
  complete(count = full_seq(1:count, 1))%>%
  select(-count) %>%
  ungroup() %>%
  mutate(size = as.numeric(size)) %>%
  mutate(protection = ifelse(site == "IVEE" | site == "NAPL", "MPA", "FISHED"))

bysite <- ggplot(urc, aes(x = size))+
  geom_density(aes(fill = site), alpha = 0.5, adjust = 2)+
  facet_wrap(~common_name)

ggsave("figures/urchinsize_bysite.png", bysite, device = "png")


dh <- ggplot(urc, aes(x = size))+
  geom_density(aes(fill = common_name), alpha = 0.75, adjust = 1.75)+
  scale_fill_manual(values = c("#762a83", "#d73027"))+
  scale_x_continuous(breaks = seq(4,14, by = 2), labels = seq(4,14,by=2))+
  labs(x = "Urchin size (cm)", y = "Density", fill = "", title = "Urchin size-frequency" )+
  theme_pubclean()

ggsave("figures/urchinsize_densityhisto.png", dh, device = "png")


ggplot(urc, aes(x = size))+
  geom_histogram(aes(fill = common_name), position = "dodge", binwidth = 0.5)

ggplot(urc, aes(x = size))+
  geom_histogram(binwidth = 0.5)+
  facet_wrap(~common_name)

ggplot(urc, aes(x = size))+
  geom_histogram(binwidth = 0.5)+
  facet_wrap(~common_name)



# Get density data

lt <- read.csv("data/survey_data/Annual_All_Species_Biomass_at_transect.csv", stringsAsFactors = F,na.strings ="-99999") %>%
  dplyr::select("YEAR", "MONTH", "SITE", "TRANSECT", "SP_CODE", "PERCENT_COVER", "DENSITY", "DRY_GM2", "SCIENTIFIC_NAME", "COMMON_NAME", "GROUP", "MOBILITY", "GROWTH_MORPH", "COARSE_GROUPING" ) %>%
  mutate(id = paste(SITE, TRANSECT, sep = "")) %>%
  filter(COMMON_NAME == "Purple Urchin" | COMMON_NAME == "Red Urchin")

names(lt) <- tolower(names(lt))

dh.d <- ggplot(lt, aes(x = (density+1)))+
  geom_density(aes(fill = common_name), alpha = 0.75)+
  scale_x_continuous(trans=log10_trans())+
  scale_fill_manual(values = c("#762a83", "#d73027"))+
  labs(x = "Urchin density (ind./m2)", y = "Density", fill = "", title = "Urchin density" )+
  theme_pubclean()

ggsave("figures/urchindensity_densityhisto.png", dh.d, device = "png")

p1 <- cowplot::plot_grid(dh, dh.d+labs(y = ""))

ggsave("figures/urchinhistos.png", p1, device = "png", width = 10, height = 6)
