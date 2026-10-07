
********************************************************************************
*Distribution of Migration Costs
********************************************************************************
use "$a\migr_gravity_ROW.dta", clear
keep o_id d_id distance
duplicates drop _all, force
gen ihsdist=asinh(distance)
gen betas=1.171 // <--input coefficients here
gen betau=1.315 // <--
gen dist_s=betas*ihsdist
gen dist_u=betau*ihsdist
cd "$output"
scatter dist_s dist_u
twoway (kdensity dist_s, color(pink) lpattern(dash)) ///
       (kdensity dist_u, lcolor(edkblue) lwidth(medium) fcolor(none)), ///
		legend(order(1 "Skilled" 2 "Low-Skilled") ring(4) position(6) row(1)) xtitle("Migration Costs") ///
	   xlabel(0(2)10, nogrid) ylabel(, nogrid) 
graph export migcost_ROW_ihs.png, replace 
  
*rearrange data to do execute test
preserve 
	rename dist_s dist
	gen skilled=1
	tempfile z
	save `z'
restore
	rename dist_u dist
	gen skilled=0
	append using `z'
count 
bys skilled: sum dist
ksmirnov dist, by(skilled)	  
/*

Two-sample Kolmogorovâ€“Smirnov test for equality of distribution functions

Smaller group             D     p-value  
---------------------------------------
0                    0.0000       1.000
1                   -0.4424       0.000
Combined K-S         0.4424       0.000

Note: Ties exist in combined dataset;
      there are 1949168 unique values out of 5126370 observations.



*/
