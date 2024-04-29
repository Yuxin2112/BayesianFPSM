clear all
capture log close
local todaydate: di %tdCYND date(c(current_date),"DMY")
log using "C:\Users\YuxinHuang\Work\Stata\Log\cancer_`todaydate'", replace text
set more off


/******************************************************************************************
	   Program: cancer.do
	   
	   Project: Examine spatial variation in cumulative incidence function due to cancer and the loss in life expectancy.  
	   
	   Purpose: Prepare cancer cohort data set
	   
	   Data:   simulated raw data from Github: Rawdata
	           
	   Created: 2012-12-02 by Yuliya Leontyeva  
	   
	   Software: STATA 17
   *****************************************************************************************/  
// Set the timer of the programm
timer on 1   

// Set the working directory


cd "C:\Users\YuxinHuang\Work\Stata\do file"

sysdir set PLUS "C:\Users\YuxinHuang\ado\plus"
sysdir set PERSONAL "C:\Users\YuxinHuang\ado\plus"

clear all
set more off
set scheme sj

//Load cancer data: 

use "RAW\cancer.dta", clear

// control that there are no any duplicates
*duplicates report id
*sort id
 

// Convert strings to date variables:

gen exit1=date(exit,"DMY")
gen dx_date1=date(dx_date,"DMY") // diagnosis time
gen dth_date1 = date(dth_date,"DMY") // date of death 

drop if bagegroup == . 
// 394 observations are deleted 

keep if age >= 15
keep if age < 89

// Check whether there are observations, where death time == time of diagnosis:

gen ind = 1 if dx_date1 == exit1 & death == 1   
tab ind

// There are 10 observations 

// Add 0.5 days to the death date for those who has date of death the same as date of diagnosis

replace exit1 = exit1 + 0.5 if ind == 1


// keep only necessary for the analysis variables:

keep id sex year dx_date1 exit1 death age bagegroup grid  

// stset the data set 
stset exit1, f(death) origin(time dx_date1) id(id) scale(365.24)

keep if _st == 1

// 361 observations have date of death before date of diagnosis 

// Calculate restricted cubic splines for ln(t):

*gen lnt = log(_t)

*rcsgen lnt, gen(rcs) dgen(drcs) df(3)


// Attained age
gen _age = min(floor(age + _t), 89) 

// Attained year 
gen _year = min(floor(year + _t), 2007)

// Merge with population file using sex, bagegroup, year & grid and keep only those who are matched :
merge m:1 sex _age _year grid using "C:\Users\YuxinHuang\Work\Stata\Data\WORK\popmort.dta", nogen keep(match)

save cancer_work.dta,replace