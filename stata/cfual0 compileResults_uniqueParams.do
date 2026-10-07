global csv 		"${user}\outputs_csv\mig_full"
global scenario "SLRonly SLRwTemp SLRwTemp_prodOnly SLRwTemp_amenOnly SLRwTemp_noMigCost SLRwTemp_reducedMig SLRwTemp_coastProtect SLRwTemp_newCity_boost20"


local firstLoop=1
foreach cfscenario of global scenario{
	*inverted objects:
	insheet using "$csv/u_resultspy_new_2010.csv", clear delimit(",") 
	drop v1
	ds
	foreach v of varlist `r(varlist)' {
		local lbl : variable label `v'
		if "`lbl'" != "" {
			local newname = strtoname("`lbl'")
			rename `v' `newname'
		}
	}
	tempfile x
	save `x'
	import delimited "$csv/resultspy_new_2010.csv", varnames(1) clear
		ds
		foreach v of varlist `r(varlist)' {
			local lbl : variable label `v'
			if "`lbl'" != "" {
				local newname = strtoname("`lbl'")
				rename `v' `newname'
			}
		}
		merge 1:1 id using `x'
		drop _merge v1
	foreach var of varlist ws-normbarTu{
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

	
	forval i=5/5 {
	insheet using "$csv/cf_`cfscenario'_sim`i'.csv", clear delimit(",") 
	ds
	foreach v of varlist `r(varlist)' {
		local lbl : variable label `v'
		if "`lbl'" != "" {
			local newname = strtoname("`lbl'")
			rename `v' `newname'
		}
	}
	keep Ws_hat Wu_hat Ls_hat Lu_hat Y_hat Y_hat2
	rename * *_sim`i'
	gen counter=_n
	merge 1:1 counter using `baseline_weights'
	assert _merge==3
	drop _merge
	tempfile baseline_weights
	save `baseline_weights', replace
	}
	
	*calculate  the new allocations of Ls/Lu from hat-ratios: 
	forval i=5/5 {
		foreach n in "s" "u"{
			gen L`n'_new_sim`i'=baseline_L`n'*L`n'_hat_sim`i'
			drop L`n'_hat_sim`i'
		}
	}
	gen totalL=baseline_Ls+baseline_Lu
	corr Ls* baseline_Ls 
	corr Lu* baseline_Lu 
	sleep 1000
	
	*Checks (make sure they have the same distribution splits):
	preserve
		collapse (rawsum) baseline_Ls baseline_Lu Ls_* Lu_* totalL 
		sum // raw-totals the same
	restore
	
	forval i=5/5 {
		foreach n in "s" "u"{
			replace W`n'_hat_sim`i'= baseline_L`n' * W`n'_hat_sim`i' 
		}
		gen W_sim`i' = Ws_hat_sim`i' + Wu_hat_sim`i' // sum across all locations 
		replace Y_hat_sim`i' = Y_hat_sim`i'*( baseline_Ls + baseline_Lu) //new dist Ls_new_sim`i'+ Lu_new_sim`i'
	}
	collapse (rawsum) Ws_* Wu_* Ls_new_* Lu_new_*  W_sim* totalL Y_hat_* baseline_Ls baseline_Lu
	
	forval i=5/5{
		replace W_sim`i'=W_sim`i'/totalL // rescale back to per-capita
		gen Y_sim`i'=Y_hat_sim`i'/totalL // rescale back to per-capita
		gen Ws_sim`i'=Ws_hat_sim`i'/ baseline_Ls // rescale back to per-capita
		gen Wu_sim`i'=Wu_hat_sim`i'/baseline_Lu // rescale back to per-capita
		drop Ws_hat_sim`i'  Wu_hat_sim`i'  Lu_new_sim`i' Ls_new_sim`i' Y_hat_sim`i'
	}
	drop totalL
	order Y_sim* W_sim* Ws_sim* Wu_sim*
	gen scenario="`cfscenario'"
	
	if `firstLoop'==1 {
		save "$a\cf_results_migFull.dta", replace
	}
	else{
		order scenario, first
		append using "$a\cf_results_migFull.dta"
		save "$a\cf_results_migFull.dta", replace
	}
	
	local firstLoop=`firstLoop'+1
}



