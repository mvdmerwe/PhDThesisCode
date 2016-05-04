
# An optimisation approach for assigning resources to defensive tasks during wildfires

## Abstract
All over the world, wildfires have a big economic, social and environmental impact. It is expected that climate change will result in more frequent, large, catastrophic wildfires. Responding to these large wildfires is a difficult task with high stakes. Incident management teams (IMTs) managing the response to large, escaped wildfires operate in high-pressure environments where they must make complex, time-critical decisions under fast moving, changing conditions.
  
Past research on providing decision support to IMTs focused on modelling initial attack, fire line construction, pre-incident deployment and longer-term planning. However, on days of extreme fire weather, when large fires are burning in hot, dry and windy conditions, fire suppression may be both ineffective and unsafe. The aim of this thesis is to address the problem of assigning resources to alternative tasks besides direct fire suppression.

A description of the wildfire resource assignment problem is presented. A mixed-integer programming model is formulated to capture features that are unique to the problem of protecting assets during wildfires. The formulated model generalises the team orienteering problem with time windows, allowing for mixed vehicle types, interchangeable and complementary vehicle capabilities, and travel times which are determined by vehicle specific speed and road network information. The protection requirements of locations are defined in terms of vehicle capabilities. 

Two approaches are presented to deal with the dynamic nature of wildfire planning: a dynamic rerouting approach and a two-stage stochastic programming approach. The rerouteing approach is appropriate when disruptions are unexpected. The aim is to reassign vehicle in a manner that minimises changes to current vehicle assignment. The stochastic approach uses likelihood estimates for fire spread scenarios. Initial vehicle assignments are made in the first stage with the opportunity for adjustments in the second stage based on observed fire-weather outcomes.

The proposed approaches resulted in a set of complementary models for wildfire resource assignment. They can, among other, account for mixed vehicle capabilities, handle unexpected changes and incorporate fire spread scenario likelihoods. The models are computationally feasible and have the potential to provide real-time decision support to IMTs.

## Thesis supplementary material

Code and scripts to accompany my thesis. 

### Model formulations

Filename | Description
--------|---------------------------------------------------
AP.cmpl | Asset protection model formulation in CMPL format.
COPTW.cmpl | Formulation for the cooperative orienteering problem with time windows in CMPL format.
Rerouting1.cmpl | Rerouting model formulation with secondary objectives and constraints presented in Section 5.3.1.
Rerouting2.cmpl | Rerouting model formulation with secondary objectives and constraints presented in  Section 5.3.2.
Rerouting3.cmpl | Rerouting model formulation  with secondary objectives and constraints presented in Section 5.3.3.
SP.m | Stochastic model formulation in Maltab's m-file format.
TOPTW.cmpl | Traditional formulation for the team orienteering problem with time windows in CMPL format.
TOPTWnew.cmpl | New formulation for the team orienteering problem with time windows presented in Chapter 3.   

### Data files

Filename | Description
--------|---------------------------------------------------
Rerouting1.csv | Pre-disruption data used for the rerouting model demonstration in Chapter 5.
Rerouting2.csv | Post-disruption data used for the rerouting model demonstration in Chapter 5.
Stochastic.csv | Data used for the stochastic model demonstration in Chapter 6.
VehicleProperties.csv | The vehicle capabilities and starting position data file used in Chapters 5 and 6. 

### Algorithms
GA_COPTW.m | Genetic algorithm for the cooperative orienteering problem

### Scripts

Filename | Description
--------|---------------------------------------------------
FetchDriveTime.php | Takes a list of longitude latitude coordinates in csv format as input. The output is drive time and distance matrices which are constructed using the Google Distance Matrix Service.
PlotSolutionGoogleMap.m | Matlab script for creating maps with the routes of vehicles overlayed. For example see Figure 3.2(b). 
