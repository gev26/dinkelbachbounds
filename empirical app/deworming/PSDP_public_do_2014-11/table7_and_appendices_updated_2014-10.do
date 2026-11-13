/*MIGUEL AND KREMER (2004)*/
/*DATA ANALYSIS DO FILES: TABLE 7 AND APPENDICES - UPDATED WITH FINAL DATA AND CODE*/
/*VERSION: OCTOBER 2014*/

clear
set more off
version 12.1

* Incorporate namelist information
	use "$da/namelist" 
	keep if visit==981 
	keep pupid sex elg98 elg99 sch98v1 wgrp* sap* totpar98 std Isem* 

* Create standard indicators based upon std98v1
	gen std_fs = std if (std>-1 & std<9) 
	replace std_fs = -1 if (std==55) 
	tab std_fs, gen(Istd) 
	summ Istd* 
	drop Istd10 std_fs 

	sort pupid 
	save "$dt/f1", replace 
	clear 

* Incorporate compliance data
	use "$da/comply" 
	keep pupid any98 any99 
	sort pupid 
	merge pupid using "$dt/f1" 
	tab _merge 
	
* Keep those with compliance data
	drop if _merge==1 
	drop _merge 
	sort pupid 
	save "$dt/f1", replace 
	clear 

* Incorporate school-level information
	use "$da/schoolvar" 
	keep *sch* *pop* mk96_s z_inf98 distlake 
* Normalize the test score
	replace mk96_s = mk96_s/0.8318*0.4357 

	rename schid sch98v1 
	sort sch98v1 
	save "$dt/f2", replace 
	use "$dt/f1" 
	sort sch98v1 
	merge sch98v1 using "$dt/f2" 
	tab _merge  
	drop _merge 

	sort pupid 
	save "$dt/f1", replace 

* Incorporate parasitological test data
	use "$da/wormed" 
	keep pupid numics98 numinf98 hw98_ics al98_who sm98_who tt98_ics ///
		numics99 numinf99 hw99_ics al99_who sm99_who tt99_ics ///
		any_hw99 any_sm99 any_al99 any_tt99 ///
		al99 any_ics98 any_98 any_hw98 any_al98 any_sm98 any_tt98 ///
		any_ics99 any_geo99_updated 
	sort pupid 
	merge pupid using "$dt/f1" 
	tab _merge 

* Keep those with parasitological test data
	drop if _merge==1 
	drop _merge 
	sort pupid 
	save "$dt/f1", replace 

* Add new controls, for G1 and G2 together
	gen sc12_3km_updated = sch1_3km_updated + sch2_3km_updated 
	gen sc12_36k_updated = sch1_36k_updated + sch2_36k_updated 
	gen po12_3km_updated = pop1_3km_updated + pop2_3km_updated 
	gen po12_36k_updated = pop1_36k_updated + pop2_36k_updated 
	gen wgrp12 = max(wgrp1, wgrp2) 

* Treatment group - population density interactions
	gen Ipop1_3_updated   = wgrp1 * pop1_3km_updated 
	gen IpopT_3_updated   = wgrp1 * popT_3km_updated 
	gen Ipop1_36_updated  = wgrp1 * pop1_36k_updated 
	gen IpopT_36_updated  = wgrp1 * popT_36k_updated 
	gen IIpop1_3_updated  = wgrp2 * pop1_3km_updated 
	gen IIpopT_3_updated  = wgrp2 * popT_3km_updated 
	gen IIpop1_36_updated = wgrp2 * pop1_36k_updated 
	gen IIpopT_36_updated = wgrp2 * popT_36k_updated 

* Incorporate ratios of Group 1 / Total Pupils
	gen ratio_03_updated = 0  
	replace ratio_03_updated = pop1_3km_updated/popT_3km_updated if popT_3km_updated~=0 
	gen ratio_36_updated = 0 
	replace ratio_36_updated = pop1_36k_updated/popT_36k_updated if popT_36k_updated~=0 

	drop if sch98v1>300 
	save "$dt/f1", replace 

* Weight each school by its total initial namelist population
	use "$dt/f1" 
	collapse (count) nsch=pupid (count) npar=numics99, by(sch98v1) 
	sort sch98v1 
	save "$dt/f2", replace 
	use "$dt/f1", clear 
	sort sch98v1 
	merge sch98v1 using "$dt/f2" 
	tab _merge  
	drop _merge 
	gen indiv_weight = nsch/npar 
	save "$dt/f1", replace 

* Generate selection into treatment indicator variable
	gen select=0 
	replace select=1 if (wgrp==1 & any98==1) | (wgrp==2 & any99==1) 
	replace select=. if (any98==. & any99==.) 
	gen Iwgrp1_select = wgrp1*select 

* Generate eligibility interaction
	gen Iwgrp1_elg = wgrp1*elg98 

* Generate difference in infection rates between 1998, 1999
	gen ics_d = any_ics99 - any_ics98 
	gen sm_d  = sm99_who - sm98_who 
	summ *_d, detail 

* Define sample of infection micro-regressions
	keep if (sch98v1~=. & (any_ics99~=. | any_ics98~=.) & Istd1~=. & pop1_3km_updated~=. & wgrp~=. & distlake~=.) 
	save "$dt/f2", replace 
	drop if (any_ics99==. | sm99_who==. | any_geo99_updated==.) 
	save "$dt/f1", replace 

* Summarize average local densities from pupil - needed for calculations on page 187
* Assumption 1: can extend regression results for group 1&2 students in standards 3-7 
*	who have worm and pupq data to the entire student population
* Assumption 2: can extend results of spillovers onto group 2 students to group 3 students

	* Average number of G1 individuals within 0-3 and 3-6 km as reported in published paper;
		collapse (mean) pop1_3km_original pop1_3km_updated pop1_36k_original pop1_36k_updated [aw=indiv_weight] if (wgrp==1 | wgrp==2) 
		summ 

	* Average number of G1 individuals within 0-3 and 3-6 km for comparison schools;
		use "$dt/f1", clear 
		collapse (mean) pop1_3km_original pop1_3km_updated pop1_36k_original pop1_36k_updated [aw=indiv_weight] if wgrp==2 
		summ 

	* Average number of G1 individuals within 0-3 and 3-6 km for treatment schools;
		use "$dt/f1", clear 
		collapse (mean) pop1_3km_original pop1_3km_updated pop1_36k_original pop1_36k_updated [aw=indiv_weight] if wgrp==1 
		summ 

	* Number of pupils in group 1 schools
		use "$da/schoolvar", clear 
		collapse (sum) pup_pop if wgrp==1 
		list 

	* Number of pupils in group 2&3 schools
		use "$da/schoolvar", clear 
		collapse (sum) pup_pop if wgrp==2 | wgrp==3 
		list 

* Define regression controls 
	use "$dt/f1", clear 
	global x_base = "sap* Istd4-Istd9 mk96_s" 

**** TABLE 7
	local i=1
	foreach var in any_ics99 sm99_who any_geo99_updated { 
		display "******** TABLE 7, (`i')********" 
		dprobit `var' wgrp1 pop1_3km_updated popT_3km_updated pop1_36k_updated popT_36k_updated $x_base [pw=indiv_weight] if (wgrp==1 | wgrp==2), robust cluster(sch98v1) 
		local i=`i'+1
		display "******** TABLE 7, (`i')********" 
		dprobit `var' wgrp1 pop1_3km_updated popT_3km_updated pop1_36k_updated popT_36k_updated select Iwgrp1_select $x_base [pw=indiv_weight] if (wgrp==1 | wgrp==2), robust cluster(sch98v1) 
		local i=`i'+1
		display "******** TABLE 7, (`i')********" 
		dprobit `var' wgrp1 pop1_3km_updated popT_3km_updated pop1_36k_updated popT_36k_updated Ipop1_3_updated Ipop1_36_updated $x_base [pw=indiv_weight] if (wgrp==1 | wgrp==2), robust cluster(sch98v1) 
		local i=`i'+1
		} 

**** APPENDIX TABLE A2
	use "$dt/f2", clear 
	collapse any98 any99 *po* $x_base indiv_weight wgrp*, by(sch98v1) 
	* APPENDIX TABLE A2, (1)
	regress any98 pop1_3km_updated pop1_36k_updated popT_3km_updated popT_36k_updated sap* Istd4-Istd9 mk96_s if wgrp==1 [pw=indiv_weight], robust 
 	* APPENDIX TABLE A2, (2)
	regress any99 po12_3km_updated po12_36k_updated popT_3km_updated popT_36k_updated sap* Istd4-Istd9 mk96_s if wgrp==1 | wgrp==2 [pw=indiv_weight], robust 
	clear 

**** APPENDIX TABLE A3
	use "$dt/f1", clear 
	* APPENDIX TABLE A3, (1)
	dprobit any_ics99 wgrp1 pop1_3km_updated popT_3km_updated pop1_36k_updated popT_36k_updated $x_base [pw=indiv_weight] if (wgrp==1 | wgrp==2), robust cluster(sch98v1) 
	*Note that column (2) requires school-level GPS data 
	* APPENDIX TABLE A3, (3)
	dprobit any_ics99 wgrp1 ratio_03_updated popT_3km_updated ratio_36_updated popT_36k_updated $x_base [pw=indiv_weight] if (wgrp==1 | wgrp==2), robust cluster(sch98v1) 
	* APPENDIX TABLE A3, (4)
	dprobit any_ics99 any_ics98 pop1_3km_updated popT_3km_updated pop1_36k_updated popT_36k_updated $x_base [pw=indiv_weight] if (wgrp==1), robust cluster(sch98v1) 
	* APPENDIX TABLE A3, (5)
	dprobit sm99_who wgrp1 pop1_3km_updated popT_3km_updated pop1_36k_updated popT_36k_updated $x_base [pw=indiv_weight] if (wgrp==1 | wgrp==2), robust cluster(sch98v1) 
	*Note that column (6) requires school-level GPS data 
	* APPENDIX TABLE A3, (7)
	dprobit sm99_who wgrp1 ratio_03_updated popT_3km_updated ratio_36_updated popT_36k_updated $x_base [pw=indiv_weight] if (wgrp==1 | wgrp==2), robust cluster(sch98v1) 
	* APPENDIX TABLE A3, (8)
	dprobit sm99_who sm98_who pop1_3km_updated popT_3km_updated pop1_36k_updated popT_36k_updated $x_base [pw=indiv_weight] if (wgrp==1), robust cluster(sch98v1) 

**** APPENDIX TABLE A4
	* APPENDIX TABLE A4, (1)
	dprobit any_ics99 wgrp1 pop1_3km_updated popT_3km_updated pop1_36k_updated popT_36k_updated select Iwgrp1_select $x_base [pw=indiv_weight] if (wgrp==1 | wgrp==2), robust cluster(sch98v1) 
	* APPENDIX TABLE A4, (2)
	regress any_ics99 wgrp1 pop1_3km_updated popT_3km_updated pop1_36k_updated popT_36k_updated select Iwgrp1_select $x_base (wgrp1 pop1_3km_updated popT_3km_updated pop1_36k_updated popT_36k_updated elg98 Iwgrp1_elg $x_base) [aw=indiv_weight] if (wgrp==1 | wgrp==2), robust cluster(sch98v1) 
	clear 

******* UPDATED P.187-188 CALCULATIONS
	use "$dt/f1", clear 
	global x_base = "sap* Istd4-Istd9 mk96_s" 

	* Average cross-school externality effect
	* Regression 1
		dprobit any_ics99 wgrp1 pop1_3km_updated popT_3km_updated pop1_36k_updated popT_36k_updated $x_base [pw=indiv_weight] if (wgrp==1 | wgrp==2), robust cluster(sch98v1) 
		probit  any_ics99 wgrp1 pop1_3km_updated popT_3km_updated pop1_36k_updated popT_36k_updated $x_base [pw=indiv_weight] if (wgrp==1 | wgrp==2), robust cluster(sch98v1) 
		matrix B=e(b) 
		matrix V=e(V) 	
	* Calculation of average cross-school externality effect (p.187)
		di (448*0.21 + 1108*0.05)/1000 
	* Calculation of average cross-school externality effect standard error (p.187)
		gen C = 1 
		* Take avg RHS values
		collapse (mean) any_ics99 wgrp1 pop1_3km_updated popT_3km_updated pop1_36k_updated popT_36k_updated $x_base C [aw=indiv_weight] if wgrp12==1  
		summ 
		* First for treatment schools
		replace wgrp1 = 0 /*looking only at externalities, not direct treatment*/
		mkmat wgrp1 - C, mat(T1) 
		* Then for comparison schools
		replace wgrp1 = 0 
		replace pop1_3km_updated = 0 
		replace pop1_36k_updated = 0 
		mkmat wgrp1 - C, mat(T0) 
		* Calculation
		matrix t1 = T1*B' 
		matrix t0 = T0*B' 
		matrix tau = normprob(t1[1,1]) - normprob(t0[1,1]) 
		matrix tau2 = t1[1,1] - t0[1,1] 
		matrix tau_v = (T1-T0)*V*(T1-T0)' 
		matrix tau_se2 = sqrt(tau_v[1,1]) 
		matrix tau_se = tau_se2[1,1]*tau[1,1]/tau2[1,1] 
		matrix list tau 
		matrix list tau_se 

	* Table VII, Regression 3 for treatment and comparison school externality effects
		use "$dt/f1", clear 
		dprobit any_ics99 wgrp1 pop1_3km_updated popT_3km_updated pop1_36k_updated popT_36k_updated Ipop1_3_updated Ipop1_36_updated $x_base [pw=indiv_weight] if (wgrp==1 | wgrp==2), robust cluster(sch98v1) 
		probit  any_ics99 wgrp1 pop1_3km_updated popT_3km_updated pop1_36k_updated popT_36k_updated Ipop1_3_updated Ipop1_36_updated $x_base [pw=indiv_weight] if (wgrp==1 | wgrp==2), robust cluster(sch98v1) 
		matrix V=e(V) 
		matrix B=e(b) 	

	* Avg cross school externality for control schools (p.187)
		di (496*0.10 + 1027*0.08)/1000 

	* Avg cross school externality for treatment schools (p.187)
		di (398*0.10 + 1192*0.08 + 398*0.27 + 1192*0.01)/1000 

	* Standard error for "true reduction" of 9 pp from paper (p.188)
		gen C = 1 
		* Take avg RHS values
		collapse (mean) any_ics99 wgrp1 pop1_3km_updated popT_3km_updated pop1_36k_updated popT_36k_updated Ipop1_3_updated Ipop1_36_updated $x_base C [aw=indiv_weight] if wgrp==1 
		summ 
		* First for treatment schools
		replace wgrp1 = 1 
		replace pop1_3km_updated = 496  /*setting to comparison school mean to give avg cross-school elasticity for comparison schools*/
		replace pop1_36k_updated = 1027 
		mkmat wgrp1 - C, mat(T1) 
		* Then for comparison schools
		replace wgrp1 = 0 
		replace pop1_3km_updated = 0 
		replace pop1_36k_updated = 0 
		replace Ipop1_3_updated = 0 
		replace Ipop1_36_updated = 0 
		mkmat wgrp1 - C, mat(T0) 
		* Calculation
		matrix t1 = T1*B' 
		matrix t0 = T0*B' 
		matrix tau = normprob(t1[1,1]) - normprob(t0[1,1]) 
		matrix tau2 = t1[1,1] - t0[1,1] 
		matrix tau_v = (T1-T0)*V*(T1-T0)' 
		matrix tau_se2 = sqrt(tau_v[1,1]) 
		matrix tau_se = tau_se2[1,1]*tau[1,1]/tau2[1,1] 
		matrix list tau 
		matrix list tau_se 
	