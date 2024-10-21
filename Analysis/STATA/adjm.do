clear all
capture log close
local todaydate: di %tdCYND date(c(current_date),"DMY")
log using "C:\Users\n11117761\Work\Log\adjm_`todaydate'", replace text
set more off


/******************************************************************************************
	   Program: adjm.do
	   
	   Project: Examine spatial variation in the loss in life expectancy.  
	   
	   Purpose: Create adjacency matrix 
	   
	   Data:    SA2_2021_AUST_SHP_GDA2020  - the contemporary data 
	            1259030002_sla06aaust_shape - SLA areas used by Susanna
	   
	   Created: 2023-01-19 by Yuliya Leontyeva  
	   Updated: 2023-02-07 by Yuliya Leontyeva
	   
	   Software: STATA 17
   *****************************************************************************************/  
// Set the timer of the programm
timer on 1   

// Change working directory:

*cd "Z:\Australia\Australia\Study\Data\RawData"
cd "C:\Users\n11117761\Work\Data\MapInfo file"  // working directory for Yuxin Huang
sysdir set PLUS "C:\Users\n11117761\ado\plus"
sysdir set PERSONAL "C:\Users\n11117761\ado\plus"  

// The shape zip file can be downloaded from Australian Bureau of Statistics for Statistical Areas Level 2:

// https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/access-and-downloads/digital-boundary-files


// Convert to Stata format shapefiles:
/*
unzipfile SA2_2021_AUST_SHP_GDA2020.zip, replace 
spshape2dta SA2_2021_AUST_GDA2020, replace saving(SA2_2021)
use SA2_2021, clear 
*/

// Use Susanna's code:

unzipfile 1259030002_sla06aaust_shape.zip, replace 
spshape2dta SLA06aAUST, replace saving(SLA06)
use SLA06, clear

// In this shape file _CX and _CY are longitude and latitude values correspondingly
// I checke the coordinates on Wikipedia 
// It means that we have to tell Stata about that

spset, modify coordsys(latlong, kilometers)



// Create a data set only for Queensland:

* keep if STE_NAME21 == "Queensland" // for 2021 data
// 548 SA2 areas in Queensland

keep if STATE_CODE == "3"
// 478 observations
codebook


/* This is for 2021 data 
// keep only variables that we need (I do not which we need, I just assume)

keep _ID _CX _CY SA2_CODE21 SA2_NAME21 CHG_FLAG21 CHG_LBL21 GCC_CODE21 GCC_NAME21 STE_CODE21 STE_NAME21 AREASQKM21

// From Jessica's code:

// Remove Chistmas, Cocos and Norfolk Island

drop if SA2_NAME21 == "Christmas Island" | SA2_NAME21 == "Cocos Island" | SA2_NAME21 == "Norfolk Island"
// My understanding that they are not in Queensland area 
drop if SA2_NAME21 == "Lord Howe Island" // also from outside Queensland

drop if _CX == . | _CY == .
// 2 observations are deleted 

*/

// Thus, we have the same number of SLA as in Susanna's code. 
// It means that NO islands were excluded (islands, low population or missing info in _CX or _CY)

// I have checked a simulated data with general population, there are no SLAs where annual average population is less than 5. The smallest number is 11.


// we have to create grid variabel, the same as in cancer cohort to be able to match then:
gen grid = _n

// verify that a new variable grid really does uniquely identify the observations
bysort grid: assert _N !=1

save "C:\Users\n11117761\Work\Data\MapInfo file\SLA06_qln", replace 
spcompress // so that shape file contains the geographical area only for Queensland

// Tell sp to use the common ID variable:

spset grid, modify replace 
save "C:\Users\n11117761\Work\Data\MapInfo file\SLA06_qln", replace 

// Create an adjacency matrix:

spmatrix create contiguity W // our matrix contains 18 islands

spmatrix summarize W, gen(adj_grid)

save "C:\Users\n11117761\Work\Data\MapInfo file\SLA06_qln", replace 

// Yuxin save the file as SLA06_QLD.dta 

// Create a variable with these areas, i.e. where adj_grid == 0 
list grid if adj_grid == 0 

spmatrix export W using contig.txt


// Define the nearest areas (as Yuxin has done it):

ssc install geodist

// create a data set only with areas without neighbours:

frame copy default temp
frame temp{
	keep if adj_grid == 0
	save SLA06_qln_noneighb, replace
}

// create a data set only with those who have neighbours:

frame default{
	keep if adj_grid != 0
	rename * *_good
	rename _* _*_good
	save SLA06_qln_neighb, replace
}

*use Z:\Australia\Australia\Study\Data\RawData\SLA06_qln_noneighb, clear

frame change temp 
cross using "SLA06_qln_neighb.dta" // merge each SLA without neighbourg with all other areas
geodist _CY _CX _CY_good _CX_good,gen(dist) sphere

// In the data _CX is a longitude and _CY is lattitude 
// sort by areas without neighbours and distance

sort grid dist


// Keep the nearest neighbourgh:

by grid: keep if _n == 1

// Ok, now for each area we have found the nearest neighbourgh


// save this data 
save Z:\Australia\Australia\Study\Data\RawData\SLA06_qln_new_neighb, replace

