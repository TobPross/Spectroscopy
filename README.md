# Spectroscopy
This repository contains the Spectroscopy script that is designed to
- read and process spectral files
- create a prediction model via partial least squares regression
- visualise the prediction model
- run outlier detection based on F-statistics
- apply the prediction model to new spectral data

In order to find the optimal prediction model, the script uses a brute force approach by running multiple plsr models with randomly choosen parameters on multiple cpu cores.
The script was tested with ASD FieldSpec4 and Bruker MPA spectral files.

## Dependencies
This script was created and tested with:
- R version 4.2.2 (2022-10-31 ucrt) -- "Innocent and Trusting"
- plantspec 1.0
- foreach 1.52
- doParallel 1.0.17
- dplyr 1.1.0
- tidyr 1.3
- ggplot 3.4.1
