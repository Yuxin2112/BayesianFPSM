clear all 
capture log close
local todaydate: di %tdCYND date(c(current_date),"DMY")
*log using "Z:\Australia\Australia\Study\Log\adj3_`todaydate'", replace text
log using "C:\Users\n11117761\Work\Log\adj3_`todaydate'", replace text
set more off

/******************************************************************************************
	   Program: contig_matrix.do
	   
	   Project: Examine spatial variation in cumulative incidence function due to cancer and the loss in life expectancy.  
	   
	   Purpose: Update weights in weigning matrix by taking into account new neigbours
	            forSLA without neigbours. Create a data set with SLA and all 
				neigbours for each SLA (this data set will be used in WinBugs)
	           
	   
	   Created: 2023-03-08 by Yuliya Leontyeva  
	   Updated: 
	   Obtained: 2023-3-09 by Yuxin Huang
	   
	   Software: STATA 17
   *****************************************************************************************/  

  // Create a data set with 19 SLA (listed by Yuxin) and assign to each SLA all the others as neigbours
  // Most of these areas are those without neigbours except 116 and 442. But there is no any problem with convergency for 116 and 442. 
  
  // According to Susanna's matrix 116 has a neighb 407; 442 - 448 
  // According to our matrix 116 - 394; 442 - 374 & 452 ( abit confused why two)
  
  
cd "C:\Users\n11117761\Work\Data"  // working directory for Yuxin Huang
sysdir set PLUS "C:\Users\n11117761\ado\plus"
sysdir set PERSONAL "C:\Users\n11117761\ado\plus"  
  
set obs 19
generate var1 = 4 in 1
replace var1 = 18 in 2
replace var1 = 229 in 3
replace var1 = 237 in 4
replace var1 = 291 in 5
replace var1 = 297 in 6
replace var1 = 317 in 7
replace var1 = 340 in 8
replace var1 = 354 in 9
replace var1 = 359 in 10
replace var1 = 374 in 11
replace var1 = 390 in 12
replace var1 = 413 in 13
replace var1 = 414 in 14
replace var1 = 416 in 15
replace var1 = 459 in 16
replace var1 = 463 in 17
replace var1 = 477 in 18
replace var1 = 478 in 19

gen var2 = var1
fillin var1 var2

drop if _fillin == 0
drop _fillin
rename var1 grid
rename var2 adj_grid



. sort grid

*. save Z:\Australia\Australia\Study\Data\RawData\temp3_grid, replace

  
. save temp3_grid, replace
   
   
  
// Make changes in the weightening matrix 

// Import original weigning matrix, created before in adjm.do

spmatrix import C using "C:\Users\n11117761\Work\Data\MapInfo file\contig.txt", replace 

spmatrix matafromsp W id = C

// use the data set with grid, which do not have neighbourgs and their nearest neighbours:

use temp3_grid, clear 

// here grid is SLA without neighbourgs and adj_grid is the nearest neigbourgh
// Change their weight from 0 to a standardised 0.15 

mata
data = st_data(.,("grid","adj_grid")) 

for (k=1; k<=rows(data); k++) {

for (i=1; i<=rows(W); i++) {
for (j=1; j<=cols(W); j++) {
if (i == data[k,1] & j == data[k,2]) W[i,j] = .154200334037381
if (i == data[k,2] & j == data[k,1]) W[i,j] = .154200334037381
 }
 }
}

*st_matrix("W", W) // save as a matrix


end

// Dimension of W matrix is r1 - r478, c1-c478

*mat list W



mata

adj = J(3000, 2, .)

r=1

for (i=1; i<=rows(W); i++) {
  for (j=1; j<=cols(W); j++) {
   if (W[i,j] != 0 ) { 
	 adj[r,1] = i 
	 adj[r,2] = j
	 
	r++
    }
	
}
}
st_matrix("AdjM", adj)

end 

clear 

// save matrix as Stata data set:
svmat AdjM, name(grid)

// you store back into matrix C
spmatrix spfrommata C = W id, replace

// Now export:

*spmatrix export C using contig_new.txt, replace 

// You can save in Stata format as well 
rename grid1 grid 
rename grid2 adj_grid 

drop if grid == .

local n = _N

// add also an extra neigbourg to Brisbie island 


insobs 4
replace grid = 189 if _n == `n' +1
replace grid = 192 if _n == `n' +2
replace grid = 116 if _n ==  `n' +3
replace grid = 442 if _n == `n' +4

replace adj_grid = 192 if _n == `n' +1
replace adj_grid = 189 if _n == `n' +2
replace adj_grid = 407 if _n == `n' +3
replace adj_grid = 448 if _n == `n' +4 

save "C:\Users\n11117761\Work\Data\RawData\adj_matrix3", replace 

/*
// Create a data set with the numbers of neigbours for each SLA:

bysort grid: gen num = _N

bysort grid: keep if _n ==1

keep grid num

// Save the matrix:

save "Z:\Australia\Australia\Study\Data\RawData\adj", replace 

*/

log close 