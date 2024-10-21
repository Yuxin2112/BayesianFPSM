# Bayesian Flexible Parametric Relative Survival Model, with an application on Queensland cancer data 

This project aims to develop a new Bayesian spatial flexible parametric survival
model (Bayesian spatial FPSM) to obtain small area estimates of Loss of Life Expectancy
(LLE) and Crude probability of death (Cr), to validate the model via a simulation study, and
to further validate it against a real-world application.

## Key contributors

Author: Yuxin Huang

Supervisors:
Helen Thompson
Susanna Cramb
Jessica Cameron
Peter Badde
Kerrie Mengersen

## Input data

The raw simulated cancer data and population data, stored in path `Data\Simulated cancer data (raw)` 
The path `Data\Simulated cancer data (clean)` contained the clean data wihch ready to used.


## Analysis

- Data preperation: 'popmort_poisson.do'

- Fitting flexible parametric model: 'stpm2_winbugs.do'

- Generating MCMC chain(s) : 'mcmc_winbugs.R'

## Repository content

- `Data`: list of data used in the analysis
- `Analysis/R`: collection of `R` scripts containing the analysis code
- `Analysis/STATA`: collection of `STATA` scripts containing the analysis code
- `Reports/`: collection of code used to generate the outputs
- `Outputs/`: collection of outputs
