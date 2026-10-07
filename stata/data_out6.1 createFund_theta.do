global csvP		"${user}\inputs_csv\"



********************************************************************************
* Create hats
********************************************************************************	
use "$a\locFundamentals.dta", clear	
egen region=group(h_reg)
egen prov=group(h_prov)
egen mun=group(h_mun h_prov)
foreach var of varlist barTs barTu barAs barAu {
	replace `var'=log(`var')
}
ds *_mean
foreach v of varlist `r(varlist)' {
    summarize `v'
	replace `v'=`r(mean)' if `v'==.
}


pwcorr barTs barTu barAs barAu `r(varlist)'

*Create hats
foreach var of varlist barTs barTu  {
	reg `var' elev_mean tri_mean soilWatercont_mean
	predict `var'_hat
}

foreach var of varlist barAs barAu  {
	reg `var' soilPH_mean
	predict `var'_hat
}


keep barTs_hat barTu_hat barAs_hat barAu_hat id_2 year
reshape wide *_hat, i(id_2) j(year)
sort id_2
export delimited "$csvP\exogObjects_resid.csv", replace
