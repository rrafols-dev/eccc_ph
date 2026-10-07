
********************************************************************************
*part 1: create migration matrices from cf 2000 -> 2010 (2010 origin -> 2010 dest)
********************************************************************************

import delimited "${user}\inputs_csv\Lod_sk_2010.csv", varnames(1) clear 
drop wa_skilled8887
drop if o_id==8887
gen order_id=_n
order order_id, first
replace order_id=order_id-1
egen wa_skilled = rowtotal(wa_skilled1- wa_skilled9993)

*generate stayers from the diagonal elements: 
gen rstayers_skilled=.
		forval i=0/1600{
			qui sum o_id if order_id==`i'
			qui replace rstayers_skilled=wa_skilled`r(mean)' if order_id==`i'
		}
drop wa_skilled1- wa_skilled9993
tempfile z_skilled
save `z_skilled', replace


import delimited "${user}\inputs_csv\Lod_un_2010.csv", varnames(1) clear 
drop wa_unskilled8887
drop if o_id==8887
gen order_id=_n
order order_id, first
replace order_id=order_id-1
egen wa_unskilled = rowtotal(wa_unskilled1- wa_unskilled9993)
*generate stayers from the diagonal elements: 
gen rstayers_unskilled=.
		forval i=0/1599 {
			qui sum o_id if order_id==`i'
			qui replace rstayers_unskilled=wa_unskilled`r(mean)' if order_id==`i'
		}
drop wa_unskilled1- wa_unskilled9993
tempfile z_unskilled
save `z_unskilled', replace


foreach name in "skilled" "unskilled" {
	di "`name'"
	import delimited "${user}/outputs_csv/mig_full/validateTHETA_migration_`name'.csv", varnames(1) clear
	drop origin
	egen rowtotal=rowtotal(dest_0-dest_1599)
	sum rowtotal // range of min-max should be around 1
	gen order_id=_n
	replace order_id=order_id-1
	merge 1:1 order_id using `z_`name''
	forval i=0/1599 {
		qui replace dest_`i'=round(dest_`i'*wa_`name')
	}

	order order_id o_id, first
	forval i=0/1599 {
		qui sum o_id if order_id==`i'
		rename dest_`i' dest_`r(mean)'ZZ
	}
	rename *ZZ *
	drop rowtotal wa_`name' _merge
	
	*generate stayers:
		gen stayers_`name'=.
		forval i=0/1599 {
			qui sum o_id if order_id==`i'
			qui replace stayers_`name'=dest_`r(mean)' if order_id==`i'
			qui replace dest_`r(mean)'=. if order_id==`i'
		}	
	corr rstayers_`name' stayers_`name'
	save "$p\validate_mig_`name'.dta", replace

}


********************************************************************************
*part 2: bring back foreign pair
********************************************************************************
**FOR SKILLED:
*row:
import delimited "${user}\inputs_csv\Lod_sk_2010.csv", varnames(1) clear 
	keep  if o_id==8887
	gen order_id=1600
	order order_id, first
	order wa_skilled8887, last
	rename wa_skilled* dest_*
	gen stayers_skilled=dest_8887
	gen rstayers_skilled=dest_8887
	replace dest_8887=. if o_id==8887
	append using  "$p\validate_mig_skilled.dta"
	drop dest_8887
	sort order_id
	save "$p\validate_mig_skilled.dta", replace
	
*column:
import delimited "${user}\inputs_csv\Lod_sk_2010.csv", varnames(1) clear 
	keep o_id wa_skilled8887
	rename wa_skilled8887 dest_8887
	merge 1:1 o_id using  "$p\validate_mig_skilled.dta"
	order dest_8887, after(dest_9993)
	drop _merge	
	save "$p\validate_mig_skilled.dta", replace

**FOR UNSKILLED:
*row:
import delimited "${user}\inputs_csv\Lod_un_2010.csv", varnames(1) clear 
	keep  if o_id==8887
	gen order_id=1600
	order order_id, first
	order wa_unskilled8887, last
	rename wa_unskilled* dest_*
	gen stayers_unskilled=dest_8887
	gen rstayers_unskilled=dest_8887
	replace dest_8887=. if o_id==8887
	append using  "$p\validate_mig_unskilled.dta"
	drop dest_8887
	sort order_id
	save "$p\validate_mig_unskilled.dta", replace
	
*column:
import delimited "${user}\inputs_csv\Lod_un_2010.csv", varnames(1) clear 
	keep o_id wa_unskilled8887
	rename wa_unskilled8887 dest_8887
	merge 1:1 o_id using  "$p\validate_mig_unskilled.dta"
	order dest_8887, after(dest_9993)
	drop _merge	
	save "$p\validate_mig_unskilled.dta", replace
	
********************************************************************************
*part 3: generate other flows
********************************************************************************	
*because RF results include foreign, must do here too:
foreach name in "skilled" "unskilled" {
	use "$p\validate_mig_`name'.dta", clear
	
		egen o_wpop_`name'=rowtotal(dest_1- dest_8887)
		qui replace o_wpop_`name'=o_wpop_`name'+stayers_`name'
		egen outflow_wa_`name'=rowtotal(dest_1- dest_8887)

		gen inflow_wa_`name'=.
		forval i=0/1600 {
			qui sum o_id if order_id==`i'
			egen temp=total(dest_`r(mean)')
			qui replace inflow_wa_`name'=temp if order_id==`i'
			drop temp
		}
		drop dest_1- dest_9993 order_id
		gen net_wa_`name' = outflow_wa_`name' - inflow_wa_`name' 
	save "$p\validate_mig_`name'.dta", replace
}



********************************************************************************
*part 4: combine
********************************************************************************


use  "$p\validate_mig_skilled.dta", clear
	merge 1:1 o_id using "$p\validate_mig_unskilled.dta", nogen
*create rates:
gen outmigsk_mod=(outflow_wa_skilled/o_wpop_skilled)
gen outmigun_mod=(outflow_wa_unskilled/o_wpop_unskilled)
gen inmig2sk_mod=(inflow_wa_skilled/o_wpop_skilled)
gen inmig2un_mod=(inflow_wa_unskilled/o_wpop_unskilled)
gen netmig2sk_mod=(net_wa_skilled/o_wpop_skilled)
gen netmig2un_mod=(net_wa_unskilled/o_wpop_unskilled)

*inmig/netmig rates as good as data qlty: trim rates above 1 (due to 0 stayers)
foreach var of varlist outmigsk_mod-netmig2un_mod { 
	qui replace `var'=. if `var'>1 //100% migrate share is a data qlty issue
}
gen year=2010
order o_id year, first
rename o_id d_id
drop dest_8887

rename *_wa* *_pred*
rename *_wpop* *_pred*

drop if d_id==8887   // foreign row
save "$p\validateTHETA_processing.dta", replace

*clean up:
erase "$p\validate_mig_skilled.dta" 
erase "$p\validate_mig_unskilled.dta"
	
	
