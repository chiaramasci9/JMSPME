# JMSPME
This repository contains the tools to fit Joint Multinomial Semi-Parametric Mixed-Effects (JMSPME) models. JMSPME models assume the random-effects relative to each response category to be discretely distributed on an a-priori unknown number of support points. This modelling induces a natural clustering of the highest-level units. 

The repository contains three R files:
- JMSPME.R implements the R function to fit the model. In the initial setting of the parameters, the user choose the type of mixed-effects structure to consider (only random intercept, only random slope, random intercept and slope)

- auxiliary.R implements an auxiliary function for the represetation a posteriori of the identified mass points of discrete random effects

- JMSPME_simulation_study implements a simulation study for a 3-categories multinomial responses and 3 different data generating processes.
