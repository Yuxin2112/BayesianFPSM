clear all
capture log close
local todaydate: di %tdCYND date(c(current_date),"DMY")
log using "C:\Users\n11117761\Work\Log\lle_est_`todaydate'", replace text
set more off


/******************************************************************************************
	   Program: lle_est.do
	   
	   Project: Examine spatial variation in cumulative incidence function due to cancer and the loss in life expectancy.  
	   
	   Purpose: Obtain the estimates of the LLE using the 2,5: 50 & 97,5 quantiles bayesian estimates
	   Data:    bc cohort: "C:\Users\leontyev\Work\Study\Data\cancer_work.dta"
	            estimates: "C:\Users\leontyev\Work\Study\Data\quantiles_est.dta"
	   
	   Created: 2022-12-15 by Yuliya Leontyeva 
	   Updated: 2023-03-06 by Yuliya Leontyeva
	   
	   updated: 2023-03-13 by Yuxin Huang
	   updated: 2023-3-17 by Yuxin Huang
	   
	   Software: STATA 17
   *****************************************************************************************/  
// Set the timer of the programm
timer on 1   


cd "C:\Users\n11117761\Work\Data"  // QUT
sysdir set PLUS "C:\Users\n11117761\ado\plus"
sysdir set PERSONAL "C:\Users\n11117761\ado\plus"

local path "C:\Users\n11117761\Work\Data"

clear frames 
clear mata


// Necessary mata programm:		
mata:

    mata set matastrict off
	
	void gq(string scalar weightsname, string scalar nodesname)
{
	n =  strtoreal(st_local("nodes"))
	i = range(1,n,1)'
	i1 = range(1,n-1,1)'
	muzero = 2
	a = J(1,n,0)
	b = i1:/sqrt(4 :* i1:^2 :- 1)	
	A= diag(a)
	for(j=1;j<=n-1;j++){
		A[j,j+1] = b[j]
		A[j+1,j] = b[j]
	}	
	symeigensystem(A,vec,nodes)
	weights = (vec[1,]:^2:*muzero)'
	weights = weights[order(nodes',1)]
	nodes = nodes'[order(nodes',1)']
	st_matrix(weightsname,weights)
	st_matrix(nodesname,nodes)
}
end


// Step 1: // Our estimates from the Bayesian analysis are saved in the winbugs_est.dta data set. In Stata after running WinBugs the data contains posterior distributions of each parameter in the model. Save those estimates in a frame called estimates:

frame rename default estimates
frame estimates {
	*use "quantiles_est.dta", clear
    use "WinBugs\further_analysis\model2\quantiles_est_agegeroup.dta", clear  
	gen ind = _n
	save temp_quantile_est.dta,replace
}

// rename to be accepted by Stata



// Create another frame with cancer cohort:


frame create data 
frame change data
frame data: use `path'\CleanData\cancer.dta, clear

// Merge with popmort file

// Attained age
gen _age = min(floor(age + _t), 89) 

// Attained year 
gen _year = min(floor(year + _t), 2007)

// Merge with population file using sex, bagegroup, year & grid and keep only those who are matched :
merge m:1 sex _age _year grid using "`path'\CleanData\popmort.dta", nogen keep(match)

// create categorical variable: agegroupaca

gen agegroupaca1 = 1 if inrange(age,15,54)
replace agegroupaca1 = 0 if inrange(age,55,89)

gen agegroupaca2 = 1 if inrange(age,55,64)
replace agegroupaca2 = 0 if inrange(age,15,54)
replace agegroupaca2 = 0 if inrange(age,65,89)

gen agegroupaca3 = 1 if inrange(age,65,74)
replace agegroupaca3 = 0 if inrange(age,15,64)
replace agegroupaca3 = 0 if inrange(age,75,89)

gen agegroupaca4 = 1 if inrange(age,75,89)
replace agegroupaca4 = 0 if inrange(age,15,74)
// Create splines of ln(_t) and save knots and Rmatrix to use afterwords:

gen lnt = log(_t)
rcsgen lnt, gen(rcs) df(5) orthog 

global ltknots = r(knots)
matrix ltR = r(R)

drop rcs*


// Create splines of age and save knots and Rmatrix to use afterwords:

rcsgen age, df(2) gen(cage) orthog 
global aknots = r(knots)
matrix aR = r(R)

// create dummy variable: agegroupaca for fit model2



// Copy the existing data set (so that not to destroy) for further manupulations

frame copy data temp 

frame change temp

// keep only necessary variables: 

keep sex year age grid cage* agegroupaca*
sort grid sex year age 
gen ind = _n


// Define all necessary local variables:

local t0 0           // the beginning of the follow-up period
local tmax 70        // 10-year restricted mean survival 
local nodes 30       // # of  intervals follow-up time is divided, usually 30 is enough  
local meanexp exp    // avarege survival for cancer patients if they did not have cancer
local meanobs obs    // average observed survival for cancer patients 
local survprob prob  // name of a covariate in the population life tables 
local maxage 89      // max age in the population life tables
local maxyear 2007   // max year in the population life tables
local grid 478       // # of regions


tempvar S_star_`t0'
gen `S_star_`t0'' =  1  // expected survival at the beginning of follow-up

local using "`path'\CleanData\popmort.dta" // population life tables 

// Loop over each year until tmax 
// Merge on interval specific expected survival at every year
// Calculate cumulative expected survival

    
	local b = `t0' + 1
	
	forvalues i=`b'/`tmax' {		
		
		tempvar S_star_`i'
		
		local j=`i'-1
	
		 gen _age = floor(min(age + `j',`maxage')) 
		 gen _year= floor(min(year  + `j', `maxyear'))  
		
		sort sex grid _year _age

		qui merge m:1 sex grid _year _age using `using', nogen keep(1 3) keepusing(`survprob')

		gen `S_star_`i''=`S_star_`j'' * `survprob'  // cumulative survival 
		
		drop _age _year `survprob' 			
		
		}
	
			
	******************************************************************
	* Find the nodes and the weights for the numerical integration   *
	* using the guassquad command and then calculate the time points *
	******************************************************************

	
    tempname weightsmat nodesmat
	mata gq("`weightsmat'","`nodesmat'") // find optimal abscissas and weights for the adaptive Legendre-Gaussian integration for a specified # of nodes 
	
	// The integral limits must be converted to [-1:1] to apply the Gaussian integration, i.e. we have to calculate:
	// ((b-a)/2)*x_i + (b - a)/2

	forvalues i=1/`nodes' {			// loop over all nodes to create the time points
		local t`i'= (`tmax'-`t0')*0.5*el(`nodesmat',`i',1)+(`tmax'+`t0')*0.5
		tempvar tvar`i'
		 gen `tvar`i''=`t`i''
	}
	
	****************************************************************************
	* Calculate cumulative expected survival at every time point of interest,  *
	* multiply with the weights, and do the integration (summation)        *
	****************************************************************************
	*set trace on 
		forvalues i=1/`nodes' {
		tempvar S`i' S_W_`i'
		local floort=floor(`t`i'')
		local ceilt=ceil(`t`i'')
		local dist=`t`i''-`floort'
	    gen `S`i''= `S_star_`floort''-(`S_star_`floort''-`S_star_`ceilt'')*`dist' 
		gen `S_W_`i''= `S`i''*el(`weightsmat',`i',1) 
		drop `S`i'' 
	}
	
	/*Sum up to get the mean expected survival*/
	local SW_list `S_W_1'
	forvalues i=2/`nodes' {
		local SW_list `SW_list' `S_W_`i''
	}
	
	 egen `meanexp' = rowtotal(`SW_list') 
	 replace `meanexp' = `meanexp'*(`tmax'-`t0')*0.5 
	
	
		
// Obtain average of observed survival:
// Loop over each time point and save R(t)S*(t) (the integrand) for each time point       	

// loop over all regions (478)

sort ind 

// merge with the data set with Bayesian estimates 

merge 1:1 ind using "temp_quantile_est.dta", nogen  

*set trace on 

/*
forvalues j = 1/`grid' {
		
	tempvar lnt1
	
	/*if grid == `j'*/ gen `lnt1' = log(`tvar1')
	
	rcsgen `lnt1' if grid == `j', gen(rcs) knots($ltknots) rmatrix(ltR)
	drop `lnt1'
	
	// loop over all observations in the quantiles data set and create local variables:

        forvalues ind = 1/2 {
				local b`ind'_2 = q2[`ind']
				*di `b`ind'_2'
				local b`ind'_50 = q50[`ind']
				*di `b`ind'_50'
				local b`ind'_97 = q97[`ind']
				*di `b`ind'_97'
			}
				
					   
	    forvalues ind = 3/8 {
			local s = `ind' - 2 
				local gamma`s'_2 = q2[`ind']
				*di `gamma`s'_2'
				local gamma`s'_50 = q50[`ind']
				*di `gamma`s'_50'
				local gamma`s'_97 = q97[`ind']		
				*di `gamma`s'_97'
		} 
		
		
		
		// obtain estimates of random effects:
		        local k = 8 +`j'
				local u`j'_2 = q2[`k']
				local u`j'_50 = q50[`k']
				local u`j'_97 = q97[`k']				
			
		        
				local l = 486 + `j'
				local v`j'_2 = q2[`l']
				local v`j'_50 = q50[`l']
				local v`j'_97 = q97[`l']				
		 	
		
	// obtain log cumulative hazard as a linear predictor at all quantilies at first time point 
	*tempvar lnH1_2 lnH1_50 lnH1_97
		
	foreach q of numlist 2 50 97 {
	tempvar lnH1_`q'
	
	*set trace on 
	gen double `lnH1_`q'' = `gamma1_`q'' + `gamma2_`q''*rcs1 + `gamma3_`q''*rcs2 + `gamma4_`q''*rcs3 + `gamma5_`q''*rcs4 + `gamma6_`q''*rcs5 + `b1_`q''*cage1 + `b2_`q''*cage2 + `u`j'_`q'' + `v`j'_`q''  
	
}
	drop rcs*
	
	// obtain relative survival at first time point at different quantilies:
	
	tempvar rs1_2 rs1_50 rs1_97 predictstat_2 predictstat_50 predictstat_97
	
	 gen double `rs1_2' = exp(-exp(`lnH1_2'))
	 gen double `rs1_50' = exp(-exp(`lnH1_50'))
	 gen double `rs1_97' = exp(-exp(`lnH1_97'))
	 
	 
	 gen double `predictstat_2' = `rs1_2'*`S_W_1'
	 gen double `predictstat_50' = `rs1_50'*`S_W_1'
	 gen double `predictstat_97' = `rs1_97'*`S_W_1'
	
	 drop `lnH1_2' `lnH1_50' `lnH1_97' `rs1_2' `rs1_50' `rs1_97'  
	
		
	forvalues i=2/`nodes' {
		
		// also a loop of three quantiles
		
		foreach q of numlist 2 50 97 {
			
		 tempvar lnt`i'
		 gen `lnt`i'' = log(`tvar`i'')
		
		rcsgen `lnt`i'' if grid == `j', gen(rcs) knots($ltknots) rmatrix(ltR) 
		drop `lnt`i''
		
		
		tempvar lnH`i'_`q' 
		
		gen double `lnH`i'_`q'' = `gamma1_`q'' + `gamma2_`q''*rcs1 + `gamma3_`q''*rcs2 + `gamma4_`q''*rcs3 + `gamma5_`q''*rcs4 + `gamma6_`q''*rcs5 + `b1_`q''*cage1 + `b2_`q''*cage2 + `u`j'_`q'' + `v`j'_`q''  
	    
		drop rcs*
		
		tempvar rs`i'_`q'
		gen double `rs`i'_`q'' = exp(-exp(`lnH`i'_`q''))
	    
			
		local predictstat_`q' `predictstat_`q'' + (`rs`i'_`q''*`S_W_`i'')
				
		drop `lnH`i'_`q'' 
		 
		}
	   
	   }
	   
	   	   
	   // Calculate mean observed survival for each quantile
	   
 foreach q of numlist 2 50 97 {	
 gen double `meanobs'_`q' = 0.5*`tmax'*(`predictstat_`q'') // approximate integral 
 
 forvalues i=2/`nodes' {
 	drop `rs`i'_`q''
 }
 }

 */
 
forvalues j = 1/`grid' {
		
	tempvar lnt1
	
	/*if grid == `j'*/ gen `lnt1' = log(`tvar1')
	
	rcsgen `lnt1' if grid == `j', gen(rcs) knots($ltknots) rmatrix(ltR)
	drop `lnt1'
	
	// loop over all observations in the quantiles data set and create local variables:

        forvalues ind = 1/6 {
				local gamma`ind'_2 = q2[`ind']
				*di `gamma`ind'_2'
				local gamma`ind'_50 = q50[`ind']
				*di `gamma`ind'_50'
				local gamma`ind'_97 = q97[`ind']		
				*di `gamma`ind'_97'
			}
				
					   
	    forvalues ind = 7/9 {
			local s = `ind' - 6 
				local b`s'_2 = q2[`ind']
				*di `b`s'_2'
				local b`s'_50 = q50[`ind']
				*di `b`s'_50'
				local b`s'_97 = q97[`ind']
				*di `b`s'_97'
		} 
		
		
		
		// obtain estimates of random effects:
		        local k = 9 +`j'
				local u`j'_2 = q2[`k']
				local u`j'_50 = q50[`k']
				local u`j'_97 = q97[`k']				
			
		        
				local l = 487 + `j'
				local v`j'_2 = q2[`l']
				local v`j'_50 = q50[`l']
				local v`j'_97 = q97[`l']				
		 	
		
	// obtain log cumulative hazard as a linear predictor at all quantilies at first time point 
	*tempvar lnH1_2 lnH1_50 lnH1_97
		
	foreach q of numlist 2 50 97 {
	tempvar lnH1_`q'
	
	*set trace on 
	gen double `lnH1_`q'' = `gamma1_`q'' + `gamma2_`q''*rcs1 + `gamma3_`q''*rcs2 + `gamma4_`q''*rcs3 + `gamma5_`q''*rcs4 + `gamma6_`q''*rcs5 + `b1_`q''*agegroupaca2 + `b2_`q''*agegroupaca3 + `b3_`q''*agegroupaca4 + `u`j'_`q'' + `v`j'_`q''  
	
}
	drop rcs*
	
	// obtain relative survival at first time point at different quantilies:
	
	tempvar rs1_2 rs1_50 rs1_97 predictstat_2 predictstat_50 predictstat_97
	
	 gen double `rs1_2' = exp(-exp(`lnH1_2'))
	 gen double `rs1_50' = exp(-exp(`lnH1_50'))
	 gen double `rs1_97' = exp(-exp(`lnH1_97'))
	 
	 
	 gen double `predictstat_2' = `rs1_2'*`S_W_1'
	 gen double `predictstat_50' = `rs1_50'*`S_W_1'
	 gen double `predictstat_97' = `rs1_97'*`S_W_1'
	
	 drop `lnH1_2' `lnH1_50' `lnH1_97' `rs1_2' `rs1_50' `rs1_97'  
	
		
	forvalues i=2/`nodes' {
		
		// also a loop of three quantiles
		
		foreach q of numlist 2 50 97 {
			
		 tempvar lnt`i'
		 gen `lnt`i'' = log(`tvar`i'')
		
		rcsgen `lnt`i'' if grid == `j', gen(rcs) knots($ltknots) rmatrix(ltR) 
		drop `lnt`i''
		
		
		tempvar lnH`i'_`q' 
		
		gen double `lnH`i'_`q'' = `gamma1_`q'' + `gamma2_`q''*rcs1 + `gamma3_`q''*rcs2 + `gamma4_`q''*rcs3 + `gamma5_`q''*rcs4 + `gamma6_`q''*rcs5 + `b1_`q''*agegroupaca2 + `b2_`q''*agegroupaca3 + `b3_`q''*agegroupaca4 + `u`j'_`q'' + `v`j'_`q''  
	    
		drop rcs*
		
		tempvar rs`i'_`q'
		gen double `rs`i'_`q'' = exp(-exp(`lnH`i'_`q''))
	    
			
		local predictstat_`q' `predictstat_`q'' + (`rs`i'_`q''*`S_W_`i'')
				
		drop `lnH`i'_`q'' 
		 
		}
	   
	   }
	   
	   	   
	   // Calculate mean observed survival for each quantile
	   
 foreach q of numlist 2 50 97 {	
 gen double `meanobs'_`q' = 0.5*`tmax'*(`predictstat_`q'') // approximate integral 
 
 forvalues i=2/`nodes' {
 	drop `rs`i'_`q''
 }
 }
 
 

preserve
 
 // save results for each region separately 
keep if grid == `j' 
tempfile temp_`j'
save `temp_`j''

restore 
drop `meanobs'_2 `meanobs'_50 `meanobs'_97 
		
	} 		

// join together all regions 
use `temp_1', clear
	
forvalues j = 2/`grid' {	
append using `temp_`j''  
}


gen lle_2 = exp - obs_2 
gen lle_50 = exp - obs_50
gen lle_97 = exp - obs_97


drop _*

// save data set:

*save lle_est.dta, replace
save lle_est_agegroup.dta, replace
	
// Turn of the timer and list it:

timer off 1
timer list 
log close	

