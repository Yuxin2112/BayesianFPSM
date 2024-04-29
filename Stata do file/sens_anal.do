clear all
capture log close
local todaydate: di %tdCYND date(c(current_date),"DMY")
log using "C:\Users\leontyev\Work\Study\Log\sens_anal_`todaydate'", replace text
set more off


/******************************************************************************************
	   Program: sens_anal.do
	   
	   Project: Examine spatial variation in cumulative incidence function due to cancer and the loss in life expectancy.  
	   
	   Purpose: Choose a "better" flexible parametric relative survival model 
	           
	   
	   Created: 2012-12-20 by Yuliya Leontyeva  
	   Updated: 2023-4-19 by Yuxin Huang
	   
	   Software: STATA 17
   *****************************************************************************************/  
// Set the timer of the programm
timer on 1   

// Set the working directory

*cd "C:\Users\leontyev\Work\Study\Data" 
*cd "C:\Users\YuliyaLeontyeva\OneDrive - Cancer Council Queensland\Desktop\Today"
*sysdir set PLUS "C:\Users\YuliyaLeontyeva\ado\plus"

// working directory for Yuxin - CCQ
cd "C:\Users\YuxinHuang\Work"
sysdir set PLUS "C:\Users\YuxinHuang\ado\plus"
sysdir set PERSONAL "C:\Users\YuxinHuang\ado\plus"

// Use a prepaired data set:

use cancer_work, clear 

// Yuxin code:

// Create dummy variables:

tab bagegroup, gen(bagegroup)


forvalues i = 2/5 {
	*xi:mfp, select(0.05): 
	stpm2 bagegroup2 bagegroup3, df(`i') scale(hazard) bhazard(rate)
	*di "df=" `i', "AIC = " %12.3f e(AIC) "BIC = " %12.3f e(BIC)
	estimates store cage_`i'
	
	
	forvalues j = 2/3 {
		rcsgen age, gen(sage) df(`j')
		*xi:mfp, select(0.05): 
		stpm2 sage*, df(`i') scale(hazard) bhazard(rate)
		drop sage*
		estimates store sage_`i'_`j'
	}
}


collect: estimates table *,  b(%7.2f) stats(dfbase ll aic bic) drop(_d*) title(Sensitivity analysis)


collect dims 
collect layout (colname) (rowname) (cmdset)
collect preview

// Export as a Latex file 
collect export "table", as(tex) replace


// Turn of the timer and list it:

timer off 1
timer list 
log close	
