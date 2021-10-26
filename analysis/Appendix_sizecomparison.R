#---------------------------------------------------------------------
## Size-frequency comparison
#---------------------------------------------------------------------

# The following code compares the distribution of urchin sizes used in the density trials with the observed size-frequency distribution of urchins in the SBC LTER dataset. 

# Get experimental data at the individual level


exp <- read.csv("data/density_experiment/raw/urchin_density_data_raw_c.csv", stringsAsFactors = FALSE, fileEncoding="UTF-8-BOM") %>% #cdata from density dependent urchin herbivory lab trials for red and purple urhcins.  5 densities tested for red urchins with three replicates of each. 8 densities tested for purple urchins with four replicates of each. 
  as_tibble() %>% 
  select(urchin_density, p_r, urchin_size, urchin_mass) %>% 
  mutate (abundance = urchin_density, 
          urchin_density = NULL, 
          common_name = case_when(p_r == "p" ~ "Purple Urchin", 
                          p_r == "r" ~ "Red Urchin"))


# Get observational data from SBC LTER

obs <- read.csv("data/survey_data/LTE_Urchin_All_Years_20210209.csv", stringsAsFactors = FALSE, na.strings = "-99999") %>% 
  janitor::clean_names() %>% 
  select(year, month, site, transect, treatment, common_name, size, count, area) %>% 
  filter(treatment == "CONTROL") %>% 
  group_by(year, month, site, transect, common_name, size) %>% 
  #summarize(count = sum(count)) %>% #What does this change other than moving it? 
  group_by(year, month, site, transect, common_name, size) %>%
  complete(count = full_seq(1:count,1))%>% #transferring each count to one, and seperating each urchin into its own line in the data. 
  select(-count) %>% #removing the count column
  ungroup() %>%
  mutate(urchin_size = as.numeric(size)) #creating a size column and establishing it as a numeric value?? 

# 
# df <- exp %>% select(common_name, urchin_size) %>% 
#   mutate(technique = "experimental") %>%
#   bind_rows(select(obs, common_name, urchin_size) %>% mutate(urchin_size = 10*urchin_size, technique = "observation")) %>% 
#   group_by(technique, common_name)

den.fun <- function(x, adjust = 2){
  x.vec = density(x, adjust = adjust)[1]$x
  y.vec = density(x, adjust = adjust )[2]$y
  data.frame(x = x.vec, y = y.vec)
}

sum <- obs %>% 
  group_by(common_name) %>%
  summarize(quant.low = quantile(urchin_size, probs = c(0.025, 0.97)))


obs.den <- obs %>% 
  select(common_name, urchin_size) %>% 
  mutate(urchin_size = 10*urchin_size, 
         technique = "observation") %>% 
  group_by(common_name, technique) %>%
  summarize(den.fun(urchin_size)) %>%
  group_by(common_name) %>%
  mutate(quant = as.factor(ifelse((common_name == "Purple Urchin" & x >=25 & x <= 70) |
                                           (common_name == "Red Urchin" & x >= 30 & x <=100), "yes", "no")))

exp.den <- exp %>%
  select(common_name, urchin_size) %>% 
  mutate( technique = "experiment") %>% 
  group_by(common_name, technique) %>%
  summarize(den.fun(urchin_size, adjust = 1))


appendixX_s1 <- ggplot(obs.den, aes(x = x, y = y))+
  geom_line(aes(group = technique), color = "gray", alpha = 0.25)+
  geom_ribbon(aes(ymax = y, fill = quant), ymin = 0)+ 
  scale_fill_brewer(guide="none")+
  facet_wrap(~ common_name, scales = "free")+
  geom_line(data = exp.den, aes(x = x, y = y, group = technique), alpha = 0.5)+
  labs(x = "Urchin test diameter (mm)", y = "Density")+
  cowplot::theme_cowplot()

cowplot::save_plot("figures/size_comparison.png", appendixX_s1)  








