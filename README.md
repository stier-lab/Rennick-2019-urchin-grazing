# Rennick-2019-urchin-grazing
Kelp consumption trials in May and summer 2019

*Predicting Herbivory Rate of Red and Purple Sea Urchins on Giant Kelp in the Santa Barbara Channel as a function of density and available detrital supply*

This repository houses code and data related to a manuscript (XXXX) that explores how the density  of purple urchins influences their herbivory rate and collective herbivory pressure within coastal kelp forests. Using both observational data gathered from 40 transects in the Santa Barbara Channel and experimental data collected from lab-based herbivory trials, we 1) produced 4 model predictions to represent the relationships between purple and red urchin density, and herbivory rate, 2) applied these model preictions to the observational data collected by the LTER in order to predict herbivory pressure through time and space, within the Santa Barbara Channel, and 3) relate those model estimates of herbivory pressure to observed trends in kelp biomass and investiage how available detrital supply may be influencing the liklihood of an ecocystem transition from a kelp-domianted state to an urchin- dominated state as the result of herbivory pressure.

This repo is maintained by Mae Rennick (GitHub: [@maerennick](https://github.com/maerennick)) at the University of California, Santa Barbara in the Department of Ecology, Evolution, & Marine Biology.

# Code

file name | description 
---|-----------
density_analysis.R | This script generates figure 2 which fits sigmoidal and linear regression models to the experimental denisty data (data/density_experiment/derived/urchin_density_data_cleaned.csv)in order to predict how the density of red and purple urchins will affect herbivory rate. All of the models produced in this script were tested and compared using AIC to determine the most appropriate model for this dataset. 
biomasshistograms.R | This script generates figure 1 uses observational data of purple and red urchin density obtained from the Santa Barbara LTER (LTE_Urchin_All_Years_20190611.csv) in order to determine the frequency of red and purple urhcins across 11 coastal sites within the Santa Barbara Channel. This script additionally includes calculations to determine mean and range of urchin density, and density comparisons between red and purple urchins. 
ggmapping.R | This scrip generates figure 3 (panel a) by mapping the California coast line and ploting predicted herbivory pressure (interaction strength between kelp biomass and urhcin biomass) derived from the denisty analysis (density_analysis.R) of red and purple urchins thorugh time across 11 sites in the Santa Barbara Channel by applying it to the observational data reported by the Santa Barbara LTER for red and purple urchin density.
spatio-temporal.R | This script generates figure 3 (pannels b and c) 4 and 5. Figure 3 was constructed by ploting predicted consumption derived from the denisty analysis (density_analysis.R) of red and purple urchins thorough time across 11 sites in the Santa Barbara Channel by applying it to the observational data reported by the Santa Barbara LTER for red and purple urchin density (Annual_All_Species_Biomass_at_transect.csv). This data was then used to track the realtionship through time and space through a generated mixed effects model. Additionally this script utilizes observational data of kelp biomass to track living kelp biomass thorough time across 11 sites in the Santa Barbara Channel which was extracted from observational data reported by the Santa Barbara LTER. To construct figure 4, we then plotted the LTER estimates of kelp biomass and combined urchin biomass at each site, and then incoorperated  available detrital supply in panel b to track the relationship between available detrital supply, predicted herbivory pressure, and live kelp biomass. Lastly, we used this data to compare changes in kelp biomass at all 11 sites when detrital supply was higher and lower than predictied consumption rate to generate figure 5. 
spatio_temporal_maps.R | This script maps predicted herbivory pressure across nine coastal sites in the Santa Barbara Channel. It visually represents strength of predicted herbivory pressure determined by the size and density analysis of red and purple urchins paired with the observational size and density data from the LTER, across time and space to generate supplementary figure A1.


# LTER Data 2020
*/data/LTER*

Santa Barbara Coastal LTER and D. Reed. 2020. SBC LTER: Reef: Long-term experiment: biomass of kelp forest species, ongoing since 2008 ver 6. Environmental Data Initiative. 

Data files exceed 100 MB, but can be downloaded from the Santa Barbara Coastal Long Term Ecological Research webpage at: https://doi.org/10.6073/pasta/47db4ee01f516b0a47b7c585fd552645 

# LTER Data 2021
*/data/LTER*

Santa Barbara Coastal LTER, D. Reed, and R. Miller. 2021. SBC LTER: Reef: Annual time series of biomass for kelp forest species, ongoing since 2000 ver 10. Environmental Data Initiative.

Data files exceed 100 MB, but can be downloaded from the Santa Barbara Coastal Long Term Ecological Research webpage at:https://doi.org/10.6073/pasta/f1cf070648d7654ada052835afb2cfe9 (Accessed 2021-03-16).




 
