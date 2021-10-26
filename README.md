# Detrial supply suppresses deforestation to maintain healthy kelp forest ecosystems

**Mae Rennick\*, Bartholomew DiFiore\*, Joseph Curtis, Daniel Reed, and Adrian Stier**  
\* Equal contributors

This repository houses code and data related to [Rennick et al. 20XX](insert link to published manuscript here)

## Abstract
The biomass of primary producers is ultimately determined by the rate at which they grow and how quickly they die or are eaten. In some instances, consumers disturb communities, dramatically reducing primary producer biomass and driving communities towards alternative ecosystem states. Yet the mechanisms behind consumer-driven disturbances remain unclear. Here, we tested how changes in the biomass and foraging habits of an important group of marine herbivores— sea urchins— triggers disturbance to a foundation species— giant kelp (Macrocystis pyrifera)— on temperate rocky reefs. Specifically, we demonstrated that urchin foraging rates increased linearly with urchin biomass in mesocosm foraging trials. We applied these results to 21 years of spatially explicit observational data on kelp community dynamics to test the long-standing hypothesis that deforestation of kelp communities occurs when the consumptive capacity of urchin populations exceeds the production of kelp detritus. We found a 50-fold reduction in giant kelp biomass when the consumptive capacity of sea urchins exceeded detrital production. Furthermore, the annual change in kelp biomass was dissociated from urchin biomass when detrital kelp was plentiful but as detrital supply became more limited, the annual change in kelp became increasingly negative. Our results suggest that the balance between detrital production and the consumptive capacity of sea urchins is critical to predicting when and where consumers are likely to force state transitions in kelp forest ecosystems.

This repo is maintained by Bart DiFiore (GitHub: [@bartdifiore](https://github.com/bartdifiore)) at the University of California, Santa Barbara in the Department of Ecology, Evolution, & Marine Biology.

# Code

file name | description 
---|-----------
density_analysis.R | This script generates figure 2 which fits sigmoidal and linear regression models to the experimental denisty data (data/density_experiment/derived/urchin_density_data_cleaned.csv)in order to predict how the density of red and purple urchins will affect herbivory rate. All of the models produced in this script were tested and compared using AIC to determine the most appropriate model for this dataset. 
biomasshistograms.R | This script generates figure 1 uses observational data of purple and red urchin density obtained from the Santa Barbara LTER (LTE_Urchin_All_Years_20190611.csv) in order to determine the frequency of red and purple urhcins across 11 coastal sites within the Santa Barbara Channel. This script additionally includes calculations to determine mean and range of urchin density, and density comparisons between red and purple urchins. 
ggmapping.R | This scrip generates figure 3 (panel a) by mapping the California coast line and ploting predicted herbivory pressure (interaction strength between kelp biomass and urhcin biomass) derived from the denisty analysis (density_analysis.R) of red and purple urchins thorugh time across 11 sites in the Santa Barbara Channel by applying it to the observational data reported by the Santa Barbara LTER for red and purple urchin density.
spatio-temporal.R | This script generates figure 3 (pannels b and c) 4 and 5. Figure 3 was constructed by ploting predicted consumption derived from the denisty analysis (density_analysis.R) of red and purple urchins thorough time across 11 sites in the Santa Barbara Channel by applying it to the observational data reported by the Santa Barbara LTER for red and purple urchin density (Annual_All_Species_Biomass_at_transect.csv). This data was then used to track the realtionship through time and space through a generated mixed effects model. Additionally this script utilizes observational data of kelp biomass to track living kelp biomass thorough time across 11 sites in the Santa Barbara Channel which was extracted from observational data reported by the Santa Barbara LTER. To construct figure 4, we then plotted the LTER estimates of kelp biomass and combined urchin biomass at each site, and then incoorperated  available detrital supply in panel b to track the relationship between available detrital supply, predicted herbivory pressure, and live kelp biomass. Lastly, we used this data to compare changes in kelp biomass at all 11 sites when detrital supply was higher and lower than predictied consumption rate to generate figure 5. 
spatio_temporal_maps.R | This script maps predicted herbivory pressure across nine coastal sites in the Santa Barbara Channel. It visually represents strength of predicted herbivory pressure determined by the size and density analysis of red and purple urchins paired with the observational size and density data from the LTER, across time and space to generate supplementary figure A1.

# Add links to data here





 
