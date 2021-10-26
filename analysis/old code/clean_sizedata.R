####################################################################################################
## Size Analysis
####################################################################################################

## Mae Rennick and Bart DiFiore
## Urchin Size Data

#The following code analyzed foraging rates of urchins of different size classes based on trials conducted in summer 2018.
# --------------------------------------------------------------------------------------------------
## Setup and cleaning
# --------------------------------------------------------------------------------------------------

library(here)
library(tidyverse)
source(here("analysis", "Functions.R"))
library(car)
library(ggeffects)


# size and weight data for purple urhcins. 6 size bins in 10 mm increments from 10 mm-69 mm were tested over two weeks. 7 urchins per tank with 7 replicates per size bin. Total of 294 urchins in 42 trials.

pf <- read.csv("data/size_data/raw/purple_urchin_size_data_raw_c.csv") %>% 
  as_tibble() %>% #coercing the dataframe together
  dplyr::select(-c(date, days_starved, week_no, tank, time_out, side, use)) %>% #removing excess data
  mutate(abundance = urchin_dens, 
         urchin_dens = NULL, 
         sp = "p",
         round = NA, 
         time_ran = 48)%>%
  rename(size_class = urchin_size_cat)#cleaned purple urchin size trials dataset

pf$size_class <- as.character(pf$size_class)

df_r <- read.csv("data/size_data/raw/red_size_data_2.csv") %>% # size and weight data for red urhcins. Round 1. Red urchins were placed into 5 19mm size bins incrementally from 10-109 mm. Each size bin was replicated 3 times over a week, with a total of urchins 75 urchins.
  filter(rep != "control") %>% 
  mutate(round =1, 
         urchin_density=5, 
         rep = NULL, 
         weight= NA, 
         sp = "r", 
         use=1) %>% 
  filter(time_ran==96)

rf <-read.csv("data/size_data/raw/red_size_data_round2.csv") %>% #3 more replicates of each of the five size classes of red urchins with 5 urchins in each trial. #75 urchins
  rename(kelp_in = Kelp_in, kelp_out =Kelp_out, sp = urchin) %>%
  filter(time_ran ==96) %>% 
  filter(use==1) %>% 
  bind_rows(df_r) %>%
  dplyr::select(-c(tank, side, level)) %>%
  filter(size_class != "NA") %>% 
  rename(abundance = urchin_density)


col_order <- c("trial_id", "round", "sp", "abundance", "size_class", "size", "weight", "kelp_in", "kelp_out", "time_ran")
pf <- pf[, col_order]
rf <- rf[, col_order] #reordering the columns so that the red and purple dataframes match 

df <- bind_rows(pf, rf) %>%
  group_by(trial_id, sp, abundance, kelp_in, kelp_out, round, time_ran, size_class) %>%
  rename(urchin_density = abundance, test_diameter = size, mass = weight) %>%
  mutate(size_class = case_when(
    size_class == "1" ~ "10-19", 
    size_class == "2" ~ "20-29", 
    size_class == "3" ~ "30-39", 
    size_class == "4" ~ "40-49", 
    size_class == "5" ~ "50-59", 
    size_class == "6" ~ "60-70", 
    size_class == "10.0-29" ~ "10-29", 
    size_class == "30-49" ~ "30-49", 
    size_class == "50-69" ~ "50-69", 
    size_class == "70-89" ~ "70-89", 
    size_class == "90-109" ~ "90-109")) %>%
  ungroup() %>%
  select(-round)


write.csv(df, file = "data/size_data/combined_sizedata.csv", row.names = F, quote = F)
