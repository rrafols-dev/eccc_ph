* ..extract cf results from 3 folders
local folders "upperbound lowerbound 2050"

local firstLoop=1
foreach f of local folders {

	di as text "===================="
	di as text "Processing folder: `f'"
	di as text "===================="

	* set path for this folder
	global csv_temp "${user}\outputs_csv\mig_full/`f'"

	*===================================*
	* BASELINE OBJECTS
	*===================================*
	insheet using "${user}\outputs_csv\mig_full\u_resultspy_new_2010.csv", clear delimit(",") 
	drop v1
	ds
	foreach v of varlist `r(varlist)' {
		local lbl : variable label `v'
		if ("`lbl'"!="") {
			local newname = strtoname("`lbl'")
			rename `v' `newname'
		}
	}
	tempfile x
	save `x'

	import delimited "${user}\outputs_csv\mig_full\resultspy_new_2010.csv", varnames(1) clear
	ds
	foreach v of varlist `r(varlist)' {
		local lbl : variable label `v'
		if ("`lbl'"!="") {
			local newname = strtoname("`lbl'")
			rename `v' `newname'
		}
	}
	merge 1:1 id using `x'
	drop _merge v1

	foreach var of varlist ws-normbarTu {
		format `var' %12.4f
	}

	sort id
	keep id Ls Lu Y u_Lso u_Luo
	rename u_Lso baseline_Ls
	rename u_Luo baseline_Lu
	rename Y baseline_Y

	tempfile baseline_weights
	gen counter=_n
	save `baseline_weights'


	*===================================*
	* SCENARIO FILES (in paper: slr scenario 5 + ave rcp)
	*===================================*	
	insheet using "$csv_temp/cf_SLRwTemp_sim5.csv", clear delimit(",") 

	ds
	foreach v of varlist `r(varlist)' {
		local lbl : variable label `v'
		if ("`lbl'"!="") {
			local newname = strtoname("`lbl'")
			rename `v' `newname'
		}
	}

	keep Ws_hat Wu_hat Ls_hat Lu_hat Y_hat Y_hat2
	rename * *_sim5
	gen counter=_n
	merge 1:1 counter using `baseline_weights'
	assert _merge==3
	drop _merge

	tempfile baseline_weights
	save `baseline_weights', replace
	


	*===================================*
	* NEW ALLOCATIONS
	*===================================*
	foreach n in "s" "u"{
		gen L`n'_new_sim5=baseline_L`n'*L`n'_hat_sim5
		drop L`n'_hat_sim5
	}


	gen totalL=baseline_Ls+baseline_Lu


	preserve
		collapse (rawsum) baseline_Ls baseline_Lu Ls_* Lu_* totalL 
		sum
	restore


	*===================================*
	* WELFARE
	*===================================*
	foreach n in "s" "u"{
		replace W`n'_hat_sim5= baseline_L`n' * W`n'_hat_sim5
	}
	gen W_sim5 = Ws_hat_sim5 + Wu_hat_sim5
	replace Y_hat_sim5 = Y_hat_sim5*( baseline_Ls + baseline_Lu)


	collapse (rawsum) Ws_* Wu_* Ls_new_* Lu_new_*  W_sim* totalL Y_hat_* baseline_Ls baseline_Lu

		replace W_sim5 = W_sim5/totalL
		gen Y_sim5 = Y_hat_sim5/totalL
		gen Ws_sim5 =Ws_hat_sim5 / baseline_Ls
		gen Wu_sim5 =Wu_hat_sim5 /baseline_Lu

		drop Ws_hat_sim5  Wu_hat_sim5  Lu_new_sim5 Ls_new_sim5 Y_hat_sim5
	
	drop totalL
	order Y_sim* W_sim* Ws_sim* Wu_sim*

	gen scenario="SLRwTemp"
	gen folder="`f'"

	if `firstLoop'==1 {
		save "$a\cf_robustness_migFull.dta", replace
	}
	else{
		order folder scenario, first
		append using "$a\cf_robustness_migFull.dta"
		save "$a\cf_robustness_migFull.dta", replace
	}

	local firstLoop=`firstLoop'+1
}


*******************
* Generate table 
*******************
use "$a\cf_results_migFull.dta", clear
keep if scenario=="SLRwTemp"
replace scenario="baseline" 
append using "$a\cf_robustness_migFull.dta"
gen ineq_sim5=Ws_sim5/Wu_sim5

replace scenario=folder if folder!=""
keep scenario *sim5 
foreach var of varlist Y_sim5-ineq_sim5{
	replace `var'=100*(`var'-1)
}


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


order varname baseline upperbound lowerbound _var2
format baseline upperbound lowerbound _var2   %9.1f
listtex varname baseline upperbound lowerbound _var2  using "$output/robustnesscf.tex", end("\\") replace	
