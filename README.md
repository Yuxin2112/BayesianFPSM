# Bayesian Flexible Parametric Relative Survival Model, with an application on Queensland cancer data 

This project aims to develop a new Bayesian spatial flexible parametric survival
model (Bayesian spatial FPSM) to obtain small area estimates of Loss of Life Expectancy
(LLE) and Crude probability of death (Cr), to validate the model via a simulation study, and
to further validate it against a real-world application.

## Key contributors

Author: Yuxin Huang

Supervisors:
A/Prof. Helen Thompson,
A/Prof. Susanna Cramb,
Dr. Jessica Cameron,
Prof. Peter Badde,
Distinguished Prof. Kerrie Mengersen

## Input data

The raw simulated cancer data and population data, stored in path `Data\Simulated cancer data (raw)` ,
The path `Data\Simulated cancer data (clean)` contained the clean data which was used in analysis.\

The folder `Data\Geographical data` contains SLA information, which used in spatial modelling, and map creation.


## Analysis

### Step 1: Data preperation
'STATA/popmort_poisson.do' -- fitting possion model to population data, to generate population mortality and survival rate\
'STATA/adjm.do' -- Create adjacency matrix 


- Fitting flexible parametric model: 'stpm2_winbugs.do'

- Generating MCMC chain(s) : 'mcmc_winbugs.R'

## Repository content

- `Data`: list of data used in the analysis
- `Analysis/R`: collection of `R` scripts containing the analysis code
- `Analysis/STATA`: collection of `STATA` scripts containing the analysis code
- `Reports/`: collection of code used to generate the outputs
- `Outputs/`: collection of outputs
