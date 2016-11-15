clear
set memory 100m
/*
infile FirmSizeCapital FirmSizeProduction FirmSizeSales ProductionCapacity Demand Productivity using "./data2.txt", clear
generate float lnFirmSizeProduction=ln(FirmSizeProduction) if FirmSizeProduction > 0
generate float lnFirmSizeSales=ln(FirmSizeSales) if FirmSizeSales > 0
generate float lnFirmSizeCapital = ln(FirmSizeCapital) if FirmSizeProduction > 0
save "./firmsize.dta", replace
*/
/*
histogram FirmSizeProduction if FirmSizeProduction > 0, bin(60) normal freq
graph save "./firmsizeproduction.gph",  replace 
/* histogram lnFirmSizeProduction if FirmSizeProduction > 0, bin(60) normal freq
graph save "./lnfirmsizeproduction.gph",  replace */
/*histogram FirmSizeSales if FirmSizeSales > 0, bin(60) normal freq
graph save "./firmsizesales.gph",  replace */
/* histogram lnFirmSizeSales if FirmSizeSales > 0, bin(60) normal freq
graph save "./lnfirmsizesales.gph",  replace */ 
histogram FirmSizeCapital if FirmSizeCapital > 0, bin(60) normal freq
graph save "./firmsizecaptital.gph",  replace 
/* histogram lnFirmSizeCapital if FirmSizeCapital > 0, bin(60) normal freq
graph save "./lnfirmsizecaptital.gph",  replace */

clear
set memory 20m
*/
infile time Yc gdpDeflator firmNum1 firmNum2 Utility UtilityPL UtilityOH WageBill SalaryBill EmploymentProdL EmploymentOHL Restarts RandD RandD1 RandD2 avgProfit avgProfitAboveAvgRandD avgProfitBelowAvgRandD avgRecentProfitType1 avgRecentProfitType2 using "./data1.txt"
/* egen float time = fill(1,2,) */
tset time
tssmooth ma Ycsmoothed1 = Yc, window(12,1,0) /* (100,1 0) */ /*  12 weeks or 100 days = 1 quarter */
gen Ycsmoothed = Ycsmoothed1 if _n > 12
tssmooth ma RandDsmoothed1 = RandD, window(12,1,0) /* (100,1 0) */
gen RandDsmoothed = RandDsmoothed1 if _n > 12
tssmooth ma RandD1smoothed1 = RandD1, window(12,1,0) /* (100,1 0) */
gen RandD1smoothed = RandD1smoothed1 if _n > 12
tssmooth ma RandD2smoothed1 = RandD2, window(12,1,0) /* (100,1 0) */
gen RandD2smoothed = RandD2smoothed1 if _n > 12

tssmooth ma Utilitysmoothed1 = UtilityPL, window(12,1,0) /* (100,1 0) */
gen UtilityPLsmoothed = Utilitysmoothed1 if _n > 12
tssmooth ma Utilitysmoothed2 = UtilityOH, window(12,1,0) /* (100,1 0) */
gen UtilityOHsmoothed = Utilitysmoothed2 if _n > 12


save "./output.dta", replace
twoway (line Yc time if time>0, clwidth(vthin)), xlabel(#5, alternate) ylabel(#3, alternate)  /* yscale(log) */
graph save "./yc.gph",  replace 
/*twoway (line nomGDP time if time>0, clwidth(vthin)), xlabel(#5, alternate) ylabel(#3, alternate)  /* yscale(log) */
graph save "./nomGDP.gph",  replace 
twoway (line realGDP time if time>0, clwidth(vthin)), xlabel(#5, alternate) ylabel(#3, alternate)  yscale(log)
graph save "./realGDP.gph",  replace 
twoway (line gdpDeflator time if time>0, clwidth(vthin)), xlabel(#5, alternate) ylabel(#3, alternate)  /* yscale(log) */
graph save "./gdpDeflator.gph",  replace 
twoway (line WageBill time if time>0, clwidth(vthin)), xlabel(#5, alternate) ylabel(#3, alternate)  /* yscale(log) */
graph save "./wageBill.gph",  replace 
*/

twoway (line Ycsmoothed time if time>0, clwidth(vthin)), xlabel(#5, alternate) ylabel(#3, alternate) /* yscale(log) */
graph save "./ycsmoothed.gph",  replace
/* 
twoway (line Yc time if time>0, clwidth(vthin)) (line Ycsmoothed time if time>0, clwidth(vthin)), ylabel(#3, alternate) 
graph save "./yc-ycsmoothed.gph", replace
twoway (line Utility time if time>0, clwidth(vthin)) , xlabel(#5, alternate) ylabel(#3, alternate)
graph save "./utility.gph", replace
*/
twoway (line UtilityPLsmoothed time if time>0, clwidth(vthin)) , xlabel(#5, alternate) ylabel(#3, alternate)
graph save "./utilityPLsmoothed.gph", replace
twoway (line UtilityOHsmoothed time if time>0, clwidth(vthin)) , xlabel(#5, alternate) ylabel(#3, alternate)
graph save "./utilityOHsmoothed.gph", replace

gen U2overU1 = UtilityOH/UtilityPL
twoway (line U2overU1 time if time>0, clwidth(vthin)) 
//xlabel (#5, alternate) ylabel(#3, alternate))
graph save "./U2overU1.gph", replace 

twoway (line Utilitysmoothed1 time if time>0, clwidth(vthin)) (line Utilitysmoothed2 time if time>0, clwidth(vthin))
graph save "./U1andU2smoothed.gph", replace

twoway (line UtilityPL time if time>0, clwidth(vthin)) (line UtilityOH time if time>0, clwidth(vthin))
graph save "./U1andU2.gph", replace

/*twoway (line Abar time if time>0, clwidth(vthin)), xlabel(#5, alternate) ylabel(#3, alternate)
graph save "./thetabar-abar.gph",  replace
twoway (line UtilityOH time if time>0, clwidth(vthin)) , xlabel(#5, alternate) ylabel(#3, alternate)
graph save "./utilityOH.gph", replace
twoway (line UtilityPL time if time>0, clwidth(vthin)) , xlabel(#5, alternate) ylabel(#3, alternate)
graph save "./utilityPL.gph", replace

graph combine "ycsmoothed" "utilityPLsmoothed" "utilityOHsmoothed"
graph save "./dynamics3.gph",  replace
*/

twoway (line firmNum1 time if time>0, clwidth(vthin)) (line firmNum2 time if time>0, clwidth(vthin))
//, xlabel(#5, alternate) ylabel(#3, alternate)
graph save "./n1-n2.gph", replace

twoway(line avgProfit time if time>0, clwidth(vthin) xlabel(#5, alternate) ylabel(#3, alternate))
twoway (line avgProfitAboveAvgRandD time, clwidth(vthin)) (line avgProfitBelowAvgRandD time, clwidth(vthin)) (line avgProfit time, clwidth(vthin))
graph save "./profitsRandD.gph", replace

twoway (line avgRecentProfitType1 time, clwidth(vthin)) (line avgRecentProfitType2 time, clwidth(vthin)) (line avgProfit time, clwidth(vthin))
graph save "./profitsByMarket.gph", replace

twoway (line EmploymentProdL time, clwidth(vthin)) (line EmploymentOHL time, clwidth(vthin))
graph save "./employment.gph", replace

twoway (line RandD1 time if time>0, clwidth(vthin))(line RandD2 time if time>0, clwidth(vthin)) 
//, xlabel(#5, alternate) ylabel(#3, alternate) , xlabel(#5, alternate) ylabel(#3, alternate)
graph save "./RandDByMarket.gph", replace

twoway (line Restarts time if time>0, clwidth(vthin) xlabel(#5, alternate) ylabel(#3, alternate)) 
graph save "./restarts.gph",  replace
twoway (line RandD time if time>0, clwidth(vthin) xlabel(#5, alternate) ylabel(#3, alternate)) 
graph save "./randd.gph",  replace
twoway (line RandDsmoothed time if time>0, clwidth(vthin) xlabel(#5, alternate) ylabel(#3, alternate)) 
graph save "./randdsmoothed.gph",  replace


graph combine "ycsmoothed" "U1andU2" "U2overU1" "n1-n2" "profitsByMarket" "RandDByMarket"
graph save "./dynamics3.gph",  replace


/*graph combine  "yc" "ycsmoothed" "nomGDP" "realGDP" "gdpDeflator" "utility" "thetabar-abar" "restarts" */
/*graph combine  "yc-ycsmoothed" "utility" "thetabar-abar" "restarts */ 
/*graph combine "ycsmoothed" "gdpDeflator" "utility" "firmsizesales" */
/*graph combine "ycsmoothed" "restarts" "utility" "firmsizeproduction"*/
/*graph combine "ycsmoothed" "randdsmoothed" "utility" "firmsizeproduction"*/ 
/*graph combine "ycsmoothed" "randdsmoothed" "utility" "wageBill" */
/*
graph combine "ycsmoothed" "utilitysmoothed" 
graph save "./dynamics1.gph",  replace
graph combine "ycsmoothed" "utility" 
graph save "./dynamics2.gph",  replace
*/


save "./output.dta", replace





