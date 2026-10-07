********************************************************************************
*Step 1: merge map with model composites  (python output)
********************************************************************************
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
	tempfile x
	save `x', replace
insheet using "$csv_main\UAs.csv", clear delimit(",")
merge m:1 id using `x'
assert _merge==3
drop _merge  id
foreach var of varlist ws-normbarTu{
	format `var' %12.4f
}
sort id_2


rename id_2 id
merge m:1 id using "$m\Raw\munmap.dta"
	assert _merge!=1
	drop _merge
	
cd "$output"
*Test equilibrium objects:
scatter wages_pred wages,  saving(z1.gph, replace)
scatter wageu_pred wageu,  saving(z2.gph, replace)
scatter Ls_pred Ls,  saving(z3.gph, replace)
scatter Lu_pred Lu,  saving(z4.gph, replace)	
graph combine z1.gph z1.gph z3.gph z4.gph , row(2) title("Predicted vs Observed") ///
graphregion(color(white) lcolor(white) style(none) lwidth(none) lpattern(blank)) ///
plotregion(fcolor(white) lcolor(white)) title(, fcolor(white) lcolor(white))


*welfare
spmap Wu using "$m\Raw\munmapc", id(id) fcolor(RdBu)  ///
 clnumber(8)  ndocolor(none) graphregion(color(white) lcolor(white) lwidth(none)) ///
 ocolor(none none none none none none none none)  ///
 legend(on col(1) pos(2) symysize(4) symxsize(4) size(small) colgap(20) bmargin(1))
graph export "welfare_unskilled.png", replace  


spmap Ws using "$m\Raw\munmapc", id(id) fcolor(RdBu) ///
 clnumber(8)  ndocolor(none) graphregion(color(white) lcolor(white) lwidth(none)) ///
 ocolor(none none none none none none none none) ///
 legend(on col(1) pos(2) symysize(4) symxsize(4) size(small) colgap(20) bmargin(1))
graph export "welfare_skilled.png", replace  

*exog amenities (fundamental productivity)
format barAu barAs %12.2f
spmap barAu using "$m\Raw\munmapc", id(id) fcolor(RdBu) ///
 clnumber(8)  ndocolor(none) graphregion(color(white) lcolor(white) lwidth(none)) ///
 ocolor(none none none none none none none none) ///
 legend(on col(1) pos(2) symysize(4) symxsize(4) size(small) colgap(20) bmargin(1))
graph export "exogamen_unskilled.png", replace 

spmap  barAs using "$m\Raw\munmapc", id(id) fcolor(RdBu)  ///
 clnumber(8)  ndocolor(none) graphregion(color(white) lcolor(white) lwidth(none)) ///
 ocolor(none none none none none none none none) ///
 legend(on col(1) pos(2) symysize(4) symxsize(4) size(small) colgap(20) bmargin(1))
graph export "exogamen_skilled.png", replace  


*exog productivity
format barTu barTs %12.2f
spmap Tu using "$m\Raw\munmapc", id(id) fcolor(RdBu)  ///
 clnumber(8)  ndocolor(none) graphregion(color(white) lcolor(white) lwidth(none)) ///
 ocolor(none none none none none none none none)  ///
 legend(on col(1) pos(2) symysize(4) symxsize(4) size(small) colgap(20) bmargin(1))
 graph export "exogprod_unskilled.png", replace 
format barTs barTs %12.2f
spmap Ts using "$m\Raw\munmapc", id(id) fcolor(RdBu)  ///
 clnumber(8)  ndocolor(none) graphregion(color(white) lcolor(white) lwidth(none)) ///
 ocolor(none none none none none none none none)  ///
 legend(on col(1) pos(2) symysize(4) symxsize(4) size(small) colgap(20) bmargin(1))
 graph export "exogprod_skilled.png", replace  
 