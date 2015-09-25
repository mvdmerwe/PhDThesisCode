# Thesis supplementary material

Code and scripts to accompany my thesis titled 'An optimisation approach for assigning resources to defensive tasks during wildfires'.

## Model formulations

Filename | Description
--------|---------------------------------------------------
AP.cmpl | Asset protection model formulation in CMPL format.
COPTW.cmpl | Formulation for the cooperative orienteering problem with time windows in CMPL format.
Rerouting1.cmpl | Rerouting model formulation with secondary objectives and constraints presented in \S 5.3.1.
Rerouting2.cmpl | Rerouting model formulation with secondary objectives and constraints presented in  \S 5.3.2.
Rerouting3.cmpl | Rerouting model formulation  with secondary objectives and constraints presented in \S 5.3.3.
SP.m | Stochastic model formulation in Maltab's m-file format.
TOPTW.cmpl | Traditional formulation for the team orienteering problem with time windows in CMPL format.
TOPTWnew.cmpl | New formulation for the team orienteering problem with time windows presented in Chapter 3.   

## Data files

Filename | Description
--------|---------------------------------------------------
Rerouting1.csv | Pre-disruption data used for the rerouting model demonstration in Chapter~5.\\ \rowcol
Rerouting2.csv | Post-disruption data used for the rerouting model demonstration in Chapter~5.\\ 
Stochastic.csv | Data used for the stochastic model demonstration in Chapter 6.\\ \rowcol
VehicleProperties.csv | The vehicle capabilities and starting position data file used in Chapters~5~and~6.\\ 

## Scripts

Filename | Description
--------|---------------------------------------------------
FetchDriveTime.php | Takes a list of longitude latitude coordinates in csv format as input. The output is drive time and distance matrices which are constructed using the Google Distance Matrix Service.
PlotSolutionGoogleMap.m | Matlab script for creating maps with the routes of vehicles overlayed. For example see Figure 3.2(b). 
