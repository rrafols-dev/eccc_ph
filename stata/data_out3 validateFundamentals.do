

********************************************************************************
*Step 1: Are recovered objects correlated with observed data? --Exogenous--
********************************************************************************	
use "$a\locFundamentals.dta", clear	
egen region=group(h_reg)
egen prov=group(h_prov)
egen mun=group(h_mun h_prov)
foreach var of varlist barTs barTu barAs barAu {
	replace `var'=log(`var')
}

mata: mata clear
local i = 1
lookfor _mean
foreach var of varlist  `r(varlist)' lat lon coastal {
	if "`var'"!="coastal" {
		replace `var'=log(`var')		
	}
	reg `var'  barAs i.year  i.prov, nocons vce(cluster prov)
		capture noi outreg, keep(barAs)  rtitle("`: var label `var''") stats(b se) ///
				 store(row`i')  starlevels(10 5 1) starloc(1) 
		capture noi	outreg, replay(simple) append(row`i') ctitles("", Skilled) ///
				store(simple) note("")
	reg `var'  barAu   i.year i.prov , nocons   vce(cluster prov)
		capture noi outreg, keep(barAu)  rtitle("`: var label `var''") stats(b se) ///
				 store(row`i')  starlevels(10 5 1) starloc(1) 
		capture noi	outreg, replay(simple3) append(row`i') ctitles("", Low-skilled ) ///
				store(simple3) note("")

	reg `var'  barTs i.year i.prov, nocons  vce(cluster prov)
		capture noi outreg, keep(barTs)  rtitle("`: var label `var''") stats(b se) ///
				 store(row`i')  starlevels(10 5 1) starloc(1) 
		capture noi	outreg, replay(simple2) append(row`i') ctitles("", Skilled) ///
				store(simple2) note("")
	reg `var'  barTu i.year i.prov , nocons  vce(cluster prov)
		capture noi outreg, keep(barTu)  rtitle("`: var label `var''") stats(b se) ///
				 store(row`i')  starlevels(10 5 1) starloc(1) 
		capture noi	outreg, replay(simple4) append(row`i') ctitles("", Low-skilled ) ///
				store(simple4) note("")				
	local ++i
	
}
	outreg, replay(simple) merge(simple3) store(simple)	
	outreg, replay(simple2) merge(simple4) store(simple2)		
	outreg, replay(simple) merge(simple2) store(simple)	
	
	outreg using "$output/temp_exogamenities_migFull.tex", ///
		replay(simple) tex nocenter note("") fragment plain replace 


********************************************************************************
*Step 2: Are recovered objects correlated with observed data? --Endogenous--
********************************************************************************
*Productivities (should validate against TFP-correlates!)
use "$a\locFundamentals.dta", clear
merge 1:1 id_2 year using "$a\LFS_imputedMeans_mun.dta", keep(1 3) nogen
merge m:1 h_prov year using "$a\tradeXports_prov.dta", keep(1 3) nogen
	*Davao Occidental (86) is part of Davao del Sur (24)
	sum value_php if h_prov=="24" & year==2000
		replace value_php=`r(mean)' if h_prov=="86" & year==2000 & value_php==.
	sum value_php if h_prov=="24" & year==2010
		replace value_php=`r(mean)' if h_prov=="86" & year==2010 & value_php==.

egen region=group(h_reg)
egen prov=group(h_prov)
egen mun=group(h_mun h_prov)

gen shcollege = p_16to64_t/p_16to64
label var value_php "Log exports"
label var shcollege "Share with college education"
label var manuf "Worker shares in manufacturing"
label var services "Worker shares in services"
label var agri "Worker shares in agriculture"
label var selfemployed "Share of own-account workers"
label var unpaidworker "Share of unpaid workers"

foreach var of varlist Ts Tu {
	replace `var'=log(`var')
}

mata: mata clear
local i = 1
foreach var of varlist shcollege unpaidworker selfemployed services manuf agri {
	replace `var'=log(`var')
	qui reg `var'  Ts i.year i.prov , nocons
		capture noi outreg, keep(Ts)  rtitle("`: var label `var''") stats(b se) ///
				 store(row`i')  starlevels(10 5 1) starloc(1) 
		capture noi	outreg, replay(simple) append(row`i') ctitles("", Skilled) ///
				store(simple) note("")
	qui reg `var'  Tu i.year i.prov, nocons
		capture noi outreg, keep(Tu)  rtitle("`: var label `var''") stats(b se) ///
				 store(row`i')  starlevels(10 5 1) starloc(1) 
		capture noi	outreg, replay(simple3) append(row`i') ctitles("", Low-skilled ) ///
				store(simple3) note("")

	local ++i
	
}

	outreg, replay(simple) merge(simple3) store(simple)	
	outreg using "$output/temp_productivities_migFull.tex", ///	
		replay(simple) tex nocenter note("") fragment plain replace 
		
********************************************************************************
*Amenities
use "$a\locFundamentals.dta", clear
egen region=group(h_reg)
egen prov=group(h_prov)
egen mun=group(h_mun h_prov)

foreach var of varlist As Au   {
	replace `var'=log(`var')
}

replace electricity=landline_cable if year==2010
replace n_stores = n_stores + 1


label var weaktenure "Share of residents with insecure housing tenure"
label var slum2 "Share of households in informal settlements"
label var healthctr "Share of villages with a health center"
label var  college "Share of villages with universities"
label var hospital "Share of villages with a hospital"
label var electricity "Share of villages with electricity/ICT"
label var planned "Share of villages with street or grid pattern"
label var n_stores "No. wholesale or department stores per 1000 residents"
label var hwayaccess "Share of villages with highway access"
label var  worship "Share of villages with a church or mosque"
mata: mata clear
local i = 1
foreach var of varlist planned electricity hwayaccess  postofc college hospital  worship weaktenure slum2 n_stores {
	replace `var'=log(`var')
	qui reg `var'  As i.year i.prov , nocons  vce(cluster prov)
		capture noi outreg, keep(As)  rtitle("`: var label `var''") stats(b se) ///
				 store(row`i')  starlevels(10 5 1) starloc(1) 
		capture noi	outreg, replay(simple) append(row`i') ctitles("", Skilled) ///
				store(simple) note("")
	qui reg `var'  Au i.year i.prov, nocons  vce(cluster prov)
		capture noi outreg, keep(Au)  rtitle("`: var label `var''") stats(b se) ///
				 store(row`i')  starlevels(10 5 1) starloc(1) 
		capture noi	outreg, replay(simple3) append(row`i') ctitles("", Low-skilled ) ///
				store(simple3) note("")
	local ++i
	
}
	outreg, replay(simple) merge(simple3) store(simple)	
	outreg using "$output/temp_amenities_migFull.tex", ///	
		replay(simple) tex nocenter note("") fragment plain replace 
