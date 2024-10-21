clear all
capture log close
local todaydate: di %tdCYND date(c(current_date),"DMY")
log using "C:\Users\n11117761\Work\Log\popmort_Poisson`todaydate'", replace text
set more off


/******************************************************************************************
	   Program: popmort.do
	   
	   Project: Examine spatial variation in cumulative incidence function due to cancer and the loss in life expectancy.  
	   
	   Purpose: Create expected mortality rates by borrowing information from the neighbouring areas with flexible Poisson model
	   
	   Data:    raw mortality file: "C:\Users\leontyev\Work\Study\Data\RawData\pop_dths.txt"
	           
	   
	   Created: 2022-05-01 by Susanna Cramb  
	   Updated: 2022-12-19 by Yuliya Leontyeva
	            2023-01-13 by Yuliya Leontyeva 
				2023-01-23 by Yuliya Leontyeva
				2023-02-09 by Yuliya Leontyeva
				
				2023-02-25 by Yuxin Huang
				2023-08-14 by Yuxin Huang
	   
	   Software: STATA 17
   *****************************************************************************************/  
// Set the timer of the programm
timer on 1 

// Download a package for calculating general population mortality rates 
cd "C:\Users\n11117761\Work\Data"  // QUT
sysdir set PLUS "C:\Users\n11117761\ado\plus"
sysdir set PERSONAL "C:\Users\n11117761\ado\plus"

*net from http://www.stata.com/stb/stb59/  
*net install ssa14 

// Define necessary locals:

// number of areas in whole population mortality file 
clear all 
local n_grid 478

local path "C:\Users\n11117761\Work\Data\RawData"


// Import population mortality file:
use "`path'\pop_dths", clear 



// In cancer cohort, the population is older 15 years, therefore we can drop agegroups == 1,2,3

drop if agegroup < 4

// We still need to create raw numbers based on aggregated 5-year age and 5-year calendar period, becaue we cannot fit Poisson model when expose (# of people at risk) = 0

// Aggregate data by 5-year calendar period (similar to Susanna's but could be discussed):

gen yeargroup = 1 if inrange(year,1997,2002)
replace yeargroup = 2 if inrange(year,2003,2007)

sort yeargroup grid sex agegroup

collapse (sum) pop count, by(grid yeargroup sex agegroup)
tempfile temp
save `temp'


// Read the adjacency matrix 

	*import delimited "`path'\adj_matrix.txt", clear 
	use "`path'\adj_matrix", clear 
	
	forvalues val = 1/`n_grid' {
	preserve
    qui keep if grid ==`val' // leave only numbers of neighbouring areas for a specific area  
    keep adj_grid
    rename adj_grid grid
    sort grid
	// merge with cancer cohort 
    qui merge 1:m grid using `temp'
	qui keep if _merge==3 | grid==`val' // we keep all observations for a given region and all its neighbours
	 
	 // calculate the total number of deaths and the total population in these areas
	 collapse (sum) pop count, by(sex yeargroup agegroup)
     gen grid=`val'
	 
	 tempfile gridcount`val'
	 qui save `gridcount`val''
	 
	 restore
	}
	 
	 
	 // merge all temporary files:
	
	use `gridcount1', clear
     
forvalues val = 2/`n_grid' {
    append using `gridcount`val''
    }

// sort the data
sort sex grid yeargroup agegroup  



// explore the range of # of deaths and # of population:

codebook count pop 	

// create observed annual mortality rates: 
gen obs_rate= count / pop


// After having updated Susanna's adj_mat we obtained that count == 0 only for agegroup == 12, yeargroup == 2 and grid == 410

// save this data set, otherwise it take a lot of time to rerun the code:

save "C:\Users\n11117761\Work\Data\Temp\grid_agg", replace 	


////////////////////////////////////////////////////////////////////////////////

use "C:\Users\n11117761\Work\Data\Temp\grid_agg", clear  	

local n_grid 478

// create midage in each age category:
gen midage=agegroup*5 - 3 

codebook midage

// midage from 17-92

replace midage=89 if midage>90

codebook midage

// Create splines of midage with knots at (18 25 50 75 87 89) // can be discussed

rcsgen midage, gen(agercs) knots(18 25 50 75 87 89) // It means that we have 4 df, because two are boundary notes

 
// Fit Poisson model with splines of midage stratified by grid & calendar period

forvalues grid = 1 / `n_grid' {
	forvalues yeargr = 1/2 {
	preserve
	

qui poisson count agercs* if grid == `grid' & yeargroup == `yeargr', exp(pop)
 
// Create the general population mortality rates for each age from 15-89 and calendar year 1997-2007:

clear
range age 15 89 75

if `yeargr' == 1 {
	range year 1997 2002 6
}
else {
	range year 2003 2007 5
}

fillin age year 
drop if age == . | year == .
gen sex = 2
gen grid = `grid'

rcsgen age, gen(agercs) knots(18 25 50 75 87 89)

predict rate, ir

drop agercs*

tempfile rate`grid'`yeargr'
save `rate`grid'`yeargr''

restore
}
}

// Merge all the grids:

use `rate11', clear
append using `rate12'

forvalues grid = 2 / `n_grid' {
	forvalues yeargr = 1/2 {
	append using `rate`grid'`yeargr''
}
}

drop _fillin
sort grid year age 

gen prob = exp(-rate) 


save "C:\Users\n11117761\Work\Data\CleanData\popmort", replace 	

timer off 1
timer list 
log close	