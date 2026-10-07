use "$a\cf_results_migFull.dta", clear

forval i=5/5{
gen ineq_sim`i'=Ws_sim`i'/Wu_sim`i'
}

lookfor sim5
foreach var in `r(varlist)'{
	label var `var' "SLR median rcp45"
}

replace scenario="baseline" if scenario=="SLRwTemp"
keep scenario  *sim5 
keep if scenario=="baseline" | strpos(scenario, "newCity") | strpos(scenario, "coastProtect")


foreach var of varlist Y_sim5-ineq_sim5{
	replace `var'=100*(`var'-1)
}
format Y_sim5-ineq_sim5  %9.2f

order scenario W_sim5 Ws_sim5 Wu_sim5 Y ineq_sim5

* Then transpose
sxpose, clear firstnames force
destring *, force replace

gen varname=""
replace varname = "Welfare, aggregate" in 1
replace varname = "Welfare, skilled" in 2
replace varname = "Welfare, unskilled" in 3
replace varname = "Output, aggregate" in 4
replace varname = "Inequality" in 5

order varname baseline SLRwTemp_coastProtect SLRwTemp_newCity_boost20 
format baseline SLRwTemp_coastProtect SLRwTemp_newCity_boost20 %9.2f
listtex varname baseline SLRwTemp_coastProtect SLRwTemp_newCity_boost20 using "$output/adaption.tex", end("\\") replace	


