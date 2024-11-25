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

## Repository content

- `Data`: list of data used in the analysis
- `Analysis/R`: collection of `R` scripts containing the analysis code
- `Analysis/STATA`: collection of `STATA` scripts containing the analysis code
- `Reports/`: collection of code used to generate the outputs
- `Outputs/`: collection of outputs

## Data

### raw data

Folder: `Data/Simulated cancer data (raw)` containing two raw datasets:
- cancer.txt -  contains simulated data to represent 
              a hypothetical gender-specific cancer

Data columns: cancer.txt (15051 records)

|id	|	Unique identifier for each case.|
|sex|	All values are 2 (one gender).|
|year|		Year of diagnosis. Ranges from 1997 to 2007.|
|site10group|	All values are 1 (hypothetical cancer).|
|dx_date|		Date of diagnosis (DDMMMYYYY).|
|dth_date|	Date of death (DDMMMYYYY).|
|age|		Age at diagnosis (years).|
|bagegroup|	Broad age group at diagnosis (1=0-49 years, 2=50-69 years, 3=70-89 years).Missing means age was 90+ years.|
|grid	|	Represents the geographical area of residence at diagnosis (values 1 to 478). |
|fu_date	|	Date of censoring (31Dec2007) or date of death (DDMMMYYYY).|
|exit	|	Same as fu_date (DDMMMYYYY).|


- pop_dths.txt - contains simulated data to represent hypothetical 
                 population mortality files

Data columns: pop_dths.txt (99902 records)

year		Ranges from 1997 to 2007.

sex		All values are 2 (one gender).

agegroup	Values from 1 to 19 representing 5-year age groups (0-4,5-9...,90+).

grid		Represents the geographical area of residence at diagnosis (values 1 to 478). 

count		Number of deaths (all causes).

pop		Population size.

### clean data
Folder `Data\Simulated cancer data (clean)` contained the clean data which was used in analysis.

Folder `Data\Geographical data` contains SLA information, which used in spatial modelling, and map creation.


## Analysis

### Step 1: Data preperation
- `STATA/popmort_poisson.do`:  fitting possion model to population data, to generate population mortality and survival rate (`popmort_work.dta`)

- `STATA/adjm.do`: Create adjacency matrix  (`adj_matrix.dta`)

### Step 2: Fitting model
(model file: `Analysis/Model_YH_gamma0.5.txt`)
-  `stpm2_winbugs.do`: runs Bayesian analysis in WinBugs from STATA

-  `mcmc_winbugs.R`: runs Bayesian analysis in WinBugs from R

-  `mcmc_diagnosis.R`: perform Geweke test on MCMC chain, and generate trace plots and density plots

- `lle_est.do` - obtains estimates for the Loss in Life Expectancy 

- `cr_est.do` - obtains estimates for the crude probability of death  

### Step 3: Generate outputs
- `results_mapping.R`: create maps of estimation


