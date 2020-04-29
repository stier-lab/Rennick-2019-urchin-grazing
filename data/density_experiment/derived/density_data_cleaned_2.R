#### Density Data Cleaning 
###Mae Rennick 

library(here)
library(tidyverse)
source(here("analysis", "Functions.R"))
library(car)

df <- read.csv("data/density_experiment/raw/urchin_density_data_raw_2.csv") %>% 
  filter(trial_number != 't') %>% 
  select(-kelp_in1, -kelp_in2, -kelp_in3, -kelp_out1, -kelp_out2, -kelp_out_R, -mortality)
