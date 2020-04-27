# Rennick-2019-urchin-grazing
Kelp consumption trials in May and summer 2019

*Predicting Herbivory Rate of Red and Purple Sea Urchins on Giant Kelp in the Santa Barbara Channel as a function of density size and available detrital supply*

This repository houses code and data related to a manuscript (XXXX) that explores how the density and size of purple urchins influences their herbivory rate and collective herbivory pressure within coastal kelp forests. Using both observational data gathered from 40 transects in the Santa Barbara Channel and experimental data collected from lab-based herbivory trials, we 1) produced 4 model predictions to represent the relationships between purple and red urchin size and density, and herbivory rate, 2) applied these model preictions to the observational data collected by the LTER in order to predict herbivory pressure through time and space, within the Santa Barbara Channel, and 3) relate those model estimates of herbivory pressure to observed trends in kelp biomass and investiage how available detrital supply may be influencing the liklihood of an ecocystem transition from a kelp-domianted state to an urchin- dominated state as the result of herbivory pressure.

This repo is maintained by Mae Rennick (GitHub: [@maerennick](https://github.com/maerennick)) at the University of California, Santa Barbara in the Department of Ecology, Evolution, & Marine Biology.

#Code

file name | description 
---|-----------

density_analysis.R | This script utilizes best fit model to predict how the density of red and purple urchins will affect herbivory rate. Models were tested and compared using AIC to determine the most appropriate model for this dataset. This script was used to create figure 2. 
size_analysis.R | This script utilizes experimental data collected from kelp consumption trials conducted at the University of Santa Barbara California during the summer of 2019, to find a best fit model to predict how the size of red and purple urchins will affect herbivory rate.
