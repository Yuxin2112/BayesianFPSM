clear all
capture log close
local todaydate: di %tdCYND date(c(current_date),"DMY")
log using "Z:\Australia\Australia\Study\Log\winbugs_`todaydate'", replace text
set more off

/******************************************************************************************
	   Program: stpm2_winbugs.do
	   
	   Project: Examine spatial variation in cumulative incidence function due to cancer and the loss in life expectancy.  
	   
	   Purpose: Obtain estimates with WinBugs
	           
	   
	   Created 21/11/2022 by Yuliya Leontyeva
       Updated 09/12/2022 by Yuliya Leontyeva 
               12/12/2022 by Yuliya Leontyeva
			   13/01/2023 by Yuliya Leontyeva
			   22/02/2023 by Yuliya Leontyeva
	   
	   Software: STATA 17
   *****************************************************************************************/  
// Set the timer of the programm
timer on 1   
clear all 
local path "C:\Users\n11117761\Work"

// Preparation of the data for WinBugs 


// rename default frame:

frame rename default work 

//Load prepared cancer data: 
clear frames  
use `path'\Data\CleanData\cancer_work.dta, clear
sort grid 
gen ind = _n // This variable will be used to link different data sets

// We have to include extra variables from adjusted matrix data set:

frame create temp 
frame temp {
	*use `path'\Data\RawData\adj_matrix.dta, clear 
	use `path'\Data\RawData\adj_matrix_our.dta, clear
	sort grid adj_grid
	gen ind = _n
}

// First, add a variable adj_grid to the existing data set.
// adj_grid variable contains ID numbers of the adjacent areas
// for each area 

// Link two data frames with id:

frlink 1:1 ind, frame(temp)
frget adj_grid, from(temp)

drop temp // drop a variable connecting two frames 

// create a vector of length 478 (the total number of areas) giving # of neighbours for each area:
frame temp {
	sort grid
	by grid: gen num = _N
	by grid: keep if _n == 1
	replace ind = grid 
}

frlink 1:1 ind, frame(temp)

frget num, from(temp)
drop temp 

// Save the data in WinBugs fromat as a list
// WinBugs does not allow variabels with underscore, so we tell WinBugs to rename it

// Update 09/12/2022 Based on the sensitivity analysis, the chosen model has 5 df for baseline and splines of continuous age with 2 df:

// create splines for age 

rcsgen age, df(2) gen(cage) orthog 

// Fit the model to obtain initial values 
stpm2 cage1 cage2, df(5) scale(hazard) bhazard(rate)

// save estimates in locals:
// I am not sure yet whether to use them :) We try without them 
/*
local cage1 `= [xb][cage1]'
local cage2 `= [xb][cage2]'

local rcs1 `= [xb][_rcs1]'
local rcs2 `= [xb][_rcs2]'
local rcs3 `= [xb][_rcs3]'
local rcs4 `= [xb][_rcs4]'
local rcs5 `= [xb][_rcs5]'
*/

cd "C:\Users\n11117761\Work\Data\WinBugs"

local N = _N

qui wbslist (var grid _d, name(grid d) format(%3.0f)) (var _rcs1 _rcs2 _rcs3 _rcs4 _rcs5 _d_rcs1 _d_rcs2 _d_rcs3 _d_rcs4 _d_rcs5 cage1 cage2 rate, name(rcs1 rcs2 rcs3 rcs4 rcs5 drcs1 drcs2 drcs3 drcs4 drcs5 cage1 cage2 rate) format(%4.3f)) (var _t,name(t) format(%4.3f)) (vector num if num !=.,format(%2.0f)) (vector adj_grid if adj_grid != . , name(adj) format(%2.0f)) (sumNumNeigh = 2724) (Nsla = 478) (datarows = `N') using Data_our.txt, replace

 
// We can take initial values from stpm2 model // Create initial values for chain 1:

// Based on Susanna's code: 
*wbslist (tauu = 1, tauv = 1, sigmau = 1, sigmav = 1, gamma = c(-1.3,1,0,0,0,0), beta = c(0.3,-0.1)) using init1.txt,  replace 

// Jessica suggested to use initial values for three chains:


frame create temp2
set seed 1368053
frame temp2 {
forvalues i = 1/3 {

forvalues j = 1/2 {
	
    local beta`j' = rnormal(0,0.00001)
}

forvalues k = 1/6 {
local gamma`k' = rnormal(0,0.00001)
}

local u = rnormal(0,1)

local v = rnormal(0,1)


local sigmau = runiform(0.01,20)
local tauu = `sigmau'^-2

local sigmav = runiform(0.01,20)
local tauv = `sigmav'^-2

wbslist (tauu = `tauu', tauv = `tauv', sigmau = `sigmau', sigmav = `sigmav', gamma = c(`gamma1', `gamma2', `gamma3', `gamma4', `gamma5', `gamma6'), beta = c(`beta1',`beta2')) using init`i'.txt,  replace 

}
}

/*
// Create the script to run in WinBugs 
wbsscript using script.txt, ///
    model("C:\Users\leontyev\Work\Study\Data\WinBugs\Model.txt") /// 
    data("C:\Users\leontyev\Work\Study\Data\WinBugs\Data.txt")  ///
    inits("C:\Users\leontyev\Work\Study\Data\WinBugs\init1.txt") ///
	log("C:\Users\leontyev\Work\Study\Data\WinBugs\log.txt")   ///
    coda("C:\Users\leontyev\Work\Study\Data\WinBugs\out") ///
	set(tauu tauv gamma beta) ///
    burn(5000) update(5000) path("`c(pwd)'") noquit replace 
*/
	// Susanna uses 150000 burn, should we use the same?
	
wbsscript using script.txt, ///
    model("Model.txt") /// 
    data("Data.txt")  ///
	initsfile(init1.txt + init2.txt + init3.txt)
	log("log.txt")   ///
    coda("out") ///
	set(tauu tauv gamma beta) ///
    burnin(10000) update(5000) thin(10) path("`c(pwd)'") noquit replace 
	
		
	// Run WinBugs within Stata 
/*
	wbsrun using `c(pwd)'/Study/Data/WinBugs/script.txt , ///
             background /// run WinBugs in the background, while Stata is in control
			 time /// displys the time of running the command
*/

wbsrun using script.txt 

// Read data from coda files created by WinBugs into Stata 			 
wbscoda using WinBugs/out, clear 

// save the data set in .dta file:

save "winbugs_est.dta", replace



// Turn of the timer and list it:

timer off 1
timer list 
log close	




