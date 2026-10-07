use "$a\cf_results_migFull.dta", clear

forval i=5/5{
gen ineq_sim`i'=Ws_sim`i'/Wu_sim`i'
}

lookfor sim1
foreach var in `r(varlist)'{
	label var `var' "no SLR yet"
}

lookfor sim2
foreach var in `r(varlist)'{
	label var  `var' "SLR extreme rcp85"
}

lookfor sim3
foreach var in `r(varlist)'{
	label var `var' "SLR median rcp85"
}

lookfor sim4
foreach var in `r(varlist)'{
	label var  `var' "SLR extreme rcp45"
}

lookfor sim5
foreach var in `r(varlist)'{
	label var `var' "SLR median rcp45"
}

lookfor sim6
foreach var in `r(varlist)'{
	label var `var' "SLR uniform 1.4 meter"
}

lookfor sim7
foreach var in `r(varlist)'{
	label var `var' "SLR uniform 2.0 meter"
}

lookfor sim8
foreach var in `r(varlist)'{
	label var `var' "SLR uniform 3.0 meter"
}

keep scenario  *sim5 
replace scenario="baseline" if scenario=="SLRwTemp"
drop if strpos(scenario, "newCity") | strpos(scenario, "coastProtect")

foreach var of varlist Y_sim5-ineq_sim5{
	replace `var'=100*(`var'-1)
}
format Y_sim5-ineq_sim5  %9.2f


order scenario W_sim5 Ws_sim5 Wu_sim5 Y ineq_sim5, first
* Then transpose
sxpose, clear firstnames force
destring *, force replace

gen varname=""
replace varname = "Welfare, aggregate" in 1
replace varname = "Welfare, skilled" in 2
replace varname = "Welfare, unskilled" in 3
replace varname = "Output, aggregate" in 4
replace varname = "Inequality" in 5


order varname baseline  SLRwTemp_reducedMig SLRwTemp_noMigCost  SLRwTemp_prodOnly SLRwTemp_amenOnly SLRonly 
format baseline  SLRwTemp_reducedMig SLRwTemp_noMigCost  SLRwTemp_prodOnly SLRwTemp_amenOnly SLRonly   %9.2f
listtex varname baseline  SLRwTemp_prodOnly SLRwTemp_amenOnly SLRonly  SLRwTemp_reducedMig SLRwTemp_noMigCost using "$output/maincf.tex", end("\\") replace	
