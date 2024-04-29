clear all
capture log close
local todaydate: di %tdCYND date(c(current_date),"DMY")
log using "C:\Users\YuxinHuang\Work\Stata\Log\master_`todaydate'", replace text
set more off


/******************************************************************************************
	   Program: master.do
	   
	   Project: Examine spatial variation in cumulative incidence function due to cancer and the loss in life expectancy.  
	   
	   Purpose: Create a master file with all steps 
	           
	   
	   Created: 2012-12-20 by Yuliya Leontyeva  
	   Updated: 2021-4-18 by Yuxin Huang
	   
	   Software: STATA 17
   *****************************************************************************************/  
// Set the timer of the programm
timer on 1   

// Set the working directory

// Working directory for Yuliya
*cd "C:\Users\leontyev\Work\Study\Data" 
*cd "C:\Users\YuliyaLeontyeva\OneDrive - Cancer Council Queensland\Desktop\Today"
*sysdir set PLUS "C:\Users\YuliyaLeontyeva\ado\plus"


// working directory for Yuxin - CCQ
cd "C:\Users\YuxinHuang\Work\Stata\do file"
sysdir set PLUS "C:\Users\YuxinHuang\ado\plus"
sysdir set PERSONAL "C:\Users\YuxinHuang\ado\plus"

// working directory for Yuxin - QUT


// Step 1: Create general population mortality files. The number of death and the total population is calculated by borrowing information from the neighbouring areas, the population mortality rates were obtained using flexible poisson model

do popmort_Poisson.do

// Step 2: Prepare cancer cohort file 

do cancer.do

// Step 3: Sensitivity analysis to choose an appropriate flexible parametric relative survival model. Obtain initial values for the model's parameters

do sens_anal.do 

// Step 3 Obtaining estimates with WunBugs 

do stpm2_winbugs.do

// Step 4 Evaluation of the Bayesian model

*do Bayes_assess.do

// Step 5 Obtaining estimates 

do lle_est_final.do

do cif_est.do

// Step 6 Create visalization

*do StataMapping_QLD.do

// Turn of the timer and list it:

timer off 1
timer list 
log close	