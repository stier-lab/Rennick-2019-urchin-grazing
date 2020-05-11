# Rennick-2019-urchin-grazing
Kelp consumption trials in May and summer 2019

*Predicting Herbivory Rate of Red and Purple Sea Urchins on Giant Kelp in the Santa Barbara Channel as a function of density size and available detrital supply*

This repository houses code and data related to a manuscript (XXXX) that explores how the density and size of purple urchins influences their herbivory rate and collective herbivory pressure within coastal kelp forests. Using both observational data gathered from 40 transects in the Santa Barbara Channel and experimental data collected from lab-based herbivory trials, we 1) produced 4 model predictions to represent the relationships between purple and red urchin size and density, and herbivory rate, 2) applied these model preictions to the observational data collected by the LTER in order to predict herbivory pressure through time and space, within the Santa Barbara Channel, and 3) relate those model estimates of herbivory pressure to observed trends in kelp biomass and investiage how available detrital supply may be influencing the liklihood of an ecocystem transition from a kelp-domianted state to an urchin- dominated state as the result of herbivory pressure.

This repo is maintained by Mae Rennick (GitHub: [@maerennick](https://github.com/maerennick)) at the University of California, Santa Barbara in the Department of Ecology, Evolution, & Marine Biology.

# Code

file name | description 
---|-----------
density_analysis.R | This script generates figure 2 which fits sigmoidal and linear regression models to the experimental denisty data (data/density_experiment/derived/urchin_density_data_cleaned.csv)in order to predict how the density of red and purple urchins will affect herbivory rate. All of the models produced in this script were tested and compared using AIC to determine the most appropriate model for this dataset. 
size_analysis.R | This script generates figure 3. It utilizes experimental size data (Purple_size_data.csv , red_size_data_2.csv , and red_size_data_round2.csv) to find a best fit model to predict how the size of red and purple urchins will affect herbivory rate. All models were compared using AIC and anova.
histograms.R | This script uses observational data of purple and red urchin density and size obtained from the Santa Barbara LTER (LTE_Urchin_All_Years_20190611.csv) in order to determine size distribution and frequency of red and purple urhcins across five coastal sites within the Santa Barbara Channel. This script additionally includes calculations to determine mean and range of urchin density and size, and density comparisons between red and purple urchins. 
observational_analysis.R | This script generates a model that fits the relationship between urchin density and kelp biomass based on observational data recorded by the Santa Barbara LTER (Annual_All_Species_Biomass_at_transect.csv). This data was then used to track the realtionship through time and space through a generated mixed effects model. 
spatio-temporal_analysis.R | This script maps predicted herbivory pressure derived from the size analayis (size_analysis.R) and denisty analysis (density_analysis.R) of red and purple urchins thorugh time across 11 sites in the Santa Barbara Channel by applying it to the observational data reported by the Santa Barbara LTER for red and purple urchin size and density. Additionally this script utilizes observational data of available detrital supply to track the relationship between available detrital supply, predicted herbivory pressure, and live kelp biomass. 
spatio_temporal_maps.R | This script maps predicted herbivory pressure across nine coastal sites in the Santa Barbara Channel. It visually represents strength of predicted herbivory pressure determined by the size and density analysis of red and purple urchins paired with the observational size and density data from the LTER, across time and space. 


# LTER Data ?
*/data/LTER*

Data files exceed 100 MB, but can be downloaded from the Santa Barbara Coastal Long Term Ecological Research webpage at: https://sbclter.msi.ucsb.edu


 
