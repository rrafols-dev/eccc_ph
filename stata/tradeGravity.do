********************************************************************************
*Step 1, merging datasets
********************************************************************************
*left hand side variables:
use "$a\migr_distance_province", clear
	merge 1:1 o_id d_id using "$a\migr_xchars_province.dta"
	assert _merge==3
	drop _merge
	tempfile x
	save `x'

*rhs:
use "$a\trade_provincepairs.dta", clear
	merge m:1 o_id d_id using `x'
	count if o_id==d_id & _merge==1 // all
	foreach var of varlist sameregion sameisland diff_lat diff_long distance {
		replace `var'=0 if _merge==1
	}
	drop _merge

*Create vars
gen logdist=log(distance)
replace logdist=0 if logdist==.
label var logdist "Log distance"
gen ihsdist=asinh(distance)
label var ihsdist "Distance"

*generate fixed effects:
egen ofe = group (o_id)
egen dfe = group (d_id)
egen oyfe = group (year o_id)
egen dyfe = group (year d_id)	
gen homebias=cond(o_id==d_id,1,0)
label var homebias "Hometown bias"
drop if year==2005



********************************************************************************
* Step 2: PPML -Appendix Table / calibration goods gravity
********************************************************************************
set more off
est clear
cd "$output"
		
foreach out in value_php {
		eststo b1: ppmlhdfe `out' ihsdist sameisland sameregion homebias diff_long diff_lat, absorb(oyfe dyfe) cluster(o_id d_id) keepsingleton
		estadd local ofe Y 
		estadd local dfe Y
		estadd local ofeonly N 
		estadd local dfeonly N	
		eststo b2: ppmlhdfe `out' ihsdist sameisland sameregion homebias diff_long diff_lat if year==2000, absorb(ofe dfe) cluster(o_id d_id) keepsingleton
		estadd local ofe N 
		estadd local dfe N
		estadd local ofeonly Y 
		estadd local dfeonly Y			
		eststo b3: ppmlhdfe `out' ihsdist sameisland sameregion homebias diff_long diff_lat if year==2010, absorb(ofe dfe) cluster(o_id d_id) keepsingleton
		estadd local ofe N 
		estadd local dfe N
		estadd local ofeonly Y 
		estadd local dfeonly Y
		esttab b*  using "tradeihsgravity_`out'.tex", label replace noconstant drop("_cons") ///
				stat(ofe dfe ofeonly dfeonly N chi2 r2_p, fmt(%9.0f %9.0f  %9.0f   %9.3f %9.3f )  ///
				labels("Origin x Year FE" "Dest. x Year FE" "Origin FE" "Dest. FE" "Obs." "Wald \$\chi^{2}\$" "Pseudo \$R^{2}\$")) ///
				star(* 0.10 ** 0.05 *** 0.01) b(%9.3f) t(%9.3f) se  compress booktabs 

}	

