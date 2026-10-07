use "$a\migr_gravity_ROW.dta", clear
gen homebias=cond(o_id==d_id,1,0)

gen migsh=wa/o_wpop
gen migshsk=wa_skilled/o_wpop_skilled
gen migshun=wa_unskilled/o_wpop_unskilled
gen migshsk2=wa_skilled/o_wpop
gen migshun2=wa_unskilled/o_wpop

sum wa wa_skilled wa_unskilled  migsh migshsk migshun  migshsk2 migshun2
*other functional forms of depvar:
gen logdist=log(distance)
gen logwas=log(wa_skilled)
gen logwau=log(wa_unskilled)
gen logwa=log(wa)
gen logmigsh=log(wa/o_wpop)
gen logmigshsk=log(wa_skilled/o_wpop_skilled)
gen logmigshun=log(wa_unskilled/o_wpop_unskilled)
gen logpopo=log(o_wpop)
replace logdist=0 if logdist==.
gen ihsdist=asinh(distance)
egen ofe = group (o_id)
egen dfe = group (d_id)
egen oyfe = group (year o_id)
egen dyfe = group (year d_id)


est clear
cd "$output"
		
local i=1
foreach out in migsh  migshsk migshun { //
		eststo b`i': ppmlhdfe `out' ihsdist sameisland sameprov homebias diff_long diff_lat, absorb(feo_`out'=oyfe fed_`out'=dyfe) cluster(o_id d_id) keepsingleton
		estadd local ofe Y 
		estadd local dfe Y
		
		esttab b`i'  using "outmaincov2_`i'_ihs.tex", label replace noconstant drop("_cons") ///
				stat(ofe dfe N chi2 r2_p, fmt(%9.0f %9.0f  %9.0f   %9.3f %9.3f )  ///
				labels("Origin x Year FE" "Dest. x Year FE" "Obs." "Wald \$\chi^{2}\$" "Pseudo \$R^{2}\$")) ///
				star(* 0.10 ** 0.05 *** 0.01) b(%9.3f) t(%9.3f) se  compress booktabs 

local i=`i'+1
}	
*when all regressions are completed, compile to one tex:
esttab b*  using "outmaincov2_migr_ihs.tex", label replace noconstant drop("_cons") ///
	stat(ofe dfe N chi2 r2_p, fmt(%9.0f %9.0f  %9.0f   %9.3f %9.3f )  ///
	labels("Origin x Year FE" "Dest. x Year FE" "Obs." "Wald \$\chi^{2}\$" "Pseudo \$R^{2}\$")) ///
	star(* 0.10 ** 0.05 *** 0.01) b(%9.3f) t(%9.3f) se  compress booktabs mtitle("All" "Skilled" "Unskilled")
*and delete the yearly version
forval i=1/3 { 
	erase "outmaincov2_`i'_ihs.tex"
}



*Referee response/robustness:
capture noi drop feo* fed*
drop if d_prov=="88" | o_prov=="88" // foreign origin and destinations
local i=1
foreach out in migsh  migshsk migshun { //
		eststo b`i': ppmlhdfe `out' ihsdist sameisland sameprov homebias diff_long diff_lat, absorb(feo_`out'=oyfe fed_`out'=dyfe) cluster(o_id d_id) keepsingleton
		estadd local ofe Y 
		estadd local dfe Y
		
		esttab b`i'  using "outmaincov2_`i'_noForeign_ihs.tex", label replace noconstant drop("_cons") ///
				stat(ofe dfe N chi2 r2_p, fmt(%9.0f %9.0f  %9.0f   %9.3f %9.3f )  ///
				labels("Origin x Year FE" "Dest. x Year FE" "Obs." "Wald \$\chi^{2}\$" "Pseudo \$R^{2}\$")) ///
				star(* 0.10 ** 0.05 *** 0.01) b(%9.3f) t(%9.3f) se  compress booktabs 

local i=`i'+1
}	
*when all regressions are completed, compile to one tex:
esttab b*  using "outmaincov2_migr_noForeign_ihs.tex", label replace noconstant drop("_cons") ///
	stat(ofe dfe N chi2 r2_p, fmt(%9.0f %9.0f  %9.0f   %9.3f %9.3f )  ///
	labels("Origin x Year FE" "Dest. x Year FE" "Obs." "Wald \$\chi^{2}\$" "Pseudo \$R^{2}\$")) ///
	star(* 0.10 ** 0.05 *** 0.01) b(%9.3f) t(%9.3f) se  compress booktabs mtitle("All" "Skilled" "Unskilled")
*and delete the yearly version
forval i=1/3 { 
	erase "outmaincov2_`i'_noForeign_ihs.tex"
}





