library(here)
library(tidyverse)

npp <- read.csv("data/survey_data/NPP_ALL_YEARS_seasonal_with_MC_stderr.csv", stringsAsFactors = F,na.strings ="-99999") %>%
  rename_all(tolower) %>%
  dplyr::select(year, season, site, plnt_dns, frnd_dns, dry_wt,f, p, d, c, b, l) %>% 
  filter(season == "3-Summer") %>%
  mutate(det.r = f + b)

av1 <- npp %>% summarize(ave = mean(det.r, na.rm = T), 
                         sd = sd(det.r, na.rm = T)) # so approximately 2% of total biomass is lost as fronds and blades per day!
av <- av1[1,1]
sd <- av1[1,2]

dt <- read.csv("data/survey_data/LTE_Detritus_Biomass_All_Years.csv", stringsAsFactors = F,na.strings ="-99999") %>%
  rename_all(tolower) %>% 
  select(year, month, date, site, transect, treatment, sp_code, wet_wt, area) %>%
  filter(sp_code == "MAPY", treatment == "CONTROL", month %in% c(7, 8, 9)) %>% 
  group_by(year, month, site, transect, area) %>% 
  summarize(detritus = sum(wet_wt, na.rm = T), 
            area = sum(area), 
            detritus.gm2 = detritus/area)

lte <- read.csv("data/survey_data/LTE_All_Species_Biomass_at_transect_20200605.csv", stringsAsFactors = F,na.strings ="-99999") %>%
  rename_all(tolower) %>% 
  select(year, month, site, transect, treatment, sp_code, wm_gm2) %>%
  filter(sp_code == "MAPY", treatment == "CONTROL", month %in% c(7, 8, 9)) %>% 
  mutate(predicted.detritus = av * wm_gm2) %>%
  select(year, month, site, transect, predicted.detritus, wm_gm2) %>%
  group_by(year, month, site, transect)

dt <- left_join(dt, lte)

ggplot(dt, aes(x = detritus.gm2, y = predicted.detritus))+
  geom_smooth(method = "lm", formula = y ~ 0+x)+
  geom_abline(aes(slope = 1, intercept = 0), lty = "dashed", col = "gray")+
  geom_point()+
  labs(x = expression(paste("Observed kelp detritus (g ",m^-2,")")), 
       y = expression(paste("Predicted detrital supply rate (g ",m^-2,d^-1,")")))+
  theme_classic()+
  theme(text = element_text(size = 22))

ggsave(filename = here::here("figures", "appendixS2-figS1.png"), device = "png")

lm1 <- lm(predicted.detritus ~ 0 + detritus.gm2, dt)
summary(lm1)

lm2 <- lm(detritus ~ wm_gm2, dt)
summary(lm2)


