/*MIGUEL AND KREMER (2004)*/
/*DATA ANALYSIS DO FILES: TABLE 9 AND APPENDIX 4*/
/*VERSION: OCTOBER 2014*/

clear
set more off
version 12.1

* Incorporate data from Namelist
	use "$da/namelist" 
	keep pupid dupid visit std 
	reshape wide std, i(pupid dupid) j(visit) 
	sort pupid dupid 
	save "$dt/f1", replace 
	clear 
	use "$da/namelist" 
	keep if visit==981 
	keep pupid dupid elg98 elg99 sch98v1 wgrp* totpar98 totpar99 sap* yrbirth 
	sort pupid dupid 
	merge pupid dupid using "$dt/f1" 
	tab _merge  
	drop _merge 
	save "$dt/f1", replace 
	clear 

* Incorporate exam data
	use "$da/test" 
	keep pupid dpexam ics98 ics99 tksep98t std98t3 stddp98 schdp98 
	drop if tksep98==0 
	rename std98t3 stdtst 
	replace stdtst=stddp98 if stdtst==. 
	drop stddp98 
	replace dpexam=0 if (tksep98==1 & dpexam==1) 
	replace dpexam=0 if (dpexam==1 & ics98==.) 
	drop if ics98==. & ics99==. 
	sort pupid 
	save "$dt/f2", replace 
	use "$dt/f1", clear 
	sort pupid 
	merge pupid using "$dt/f2" 
	tab _merge  
	drop _merge 
	sort pupid 
	replace sch98v1 = schdp98 if (sch98v1==. & dpexam==1) 
	drop if sch98v1==. & dpexam~=1 
	save "$dt/f1", replace 

* Incorporate an indicator for whether pupil responded to 1998 Pupil Questionnaire
	use "$da/pupq" 
	keep pupid sex_98_9 
	gen Ipupque = . 
	replace Ipupque = 1 if sex~=. 
	drop sex_98_9 
	sort pupid 
	save "$dt/f2", replace 
	clear 
	use "$dt/f1" 
	merge pupid using "$dt/f2" 
	tab _merge  
	drop _merge 
	sort pupid 
	save "$dt/f1", replace 

* Create measure of standard in 1998 (include drop-outs during the year among those not taking exams)
	gen     std = std981 
	replace std = std982 if (std==.) 
	replace std = stdtst if (std==.) 
	gen std98 = std if ((std>0 & std<9) | std==55) 
	save "$dt/f1", replace 

* Create an indicator for standard based on 1998 test taken
	tab stdtst, gen(Istdtst) 
	drop Istdtst1 
	save "$dt/f1", replace 

* Incorporate mock score and average infection rate
	use "$da/schoolvar" 
	keep schid wgrp mk96_s z9899_78 z9899_34 z9899_56 
	* Normalize mock score
	replace mk96_s = mk96_s*(0.4357)/(0.8318) 
	rename schid sch98v1 
	sort sch98v1 
	save "$dt/f9", replace 
	clear 
	use "$dt/f1" 
	sort sch98v1 
	merge sch98v1 using "$dt/f9" 
	tab _merge  
	drop _merge 
	sort pupid 
	gen     p1 = z9899_34 
	replace p1 = z9899_56 if (stdtst==5 | stdtst==6 | std98==5 | std98==6) 
	replace p1 = z9899_78 if (stdtst==7 | stdtst==8 | std98==7 | std98==8) 

* Create an indicator for Standards 1&2
	gen     Istd2 = 0 
	replace Istd2 = 1 if (std98==2 | std98==1) 

* Create treatment indicators
	gen     t198 = 1 if wgrp==1 
	replace t198 = 0 if (wgrp==2 | wgrp==3) 
	replace t198 = . if wgrp==. 
	
	gen     t199 = 1 if wgrp==2 
	replace t199 = 0 if (wgrp==1 | wgrp==3) 
	replace t199 = . if wgrp==. 

	gen     t298 = 0 
	replace t298 = . if wgrp==. 

	gen     t299 = 1 if wgrp==1 
	replace t299 = 0 if (wgrp==2 | wgrp==3) 
	replace t299 = . if wgrp==. 
	save "$dt/f1", replace 

* Drop duplicates and reshape to panel dataset
	sort pupid 
	duplicates drop pupid, force 
	compress 
	reshape long ics t1 t2 elg, i(pupid) j(year) 

* Create a measure of total participation by year
	gen     totpart = . 
	replace totpart = totpar98 if year==98 
	replace totpart = totpar99 if year==99 

* Create SAP-year interactions
	gen     Y98sap1 = 0 
	replace Y98sap1 = sap1 if year==98 
	gen     Y98sap2 = 0 
	replace Y98sap2 = sap2 if year==98 
	gen     Y98sap3 = 0 
	replace Y98sap3 = sap3 if year==98 
	gen     Y98sap4 = 0 
	replace Y98sap4 = sap4 if year==98 
	save "$dt/f2", replace 
	clear 

* Incorporate local pupid population density measures
	use "$da/schoolvar" 
	keep schid pop*3km* pop*36k* 
	rename schid sch98v1 
	sort sch98v1 
	save "$dt/f9", replace 
	clear 
	use "$dt/f2" 
	sort sch98v1 
	merge sch98v1 using "$dt/f9" 
	tab _merge  
	drop _merge 

	gen     popTR_3km_updated = pop1_3km_updated 
	replace popTR_3km_updated = pop1_3km_updated + pop2_3km_updated if year==99 
	gen     popTR_36k_updated = pop1_36k_updated 
	replace popTR_36k_updated = pop1_36k_updated + pop2_36k_updated if year==99 

	gen     popTR_3km_original = pop1_3km_original 
	replace popTR_3km_original = pop1_3km_original + pop2_3km_original if year==99 
	gen     popTR_36k_original = pop1_36k_original 
	replace popTR_36k_original = pop1_36k_original + pop2_36k_original if year==99 

	keep pupid ics t* e* p* sap* Y98sap* sch98v1 mk96_s Istd* Ipupque year dpexam totpart pop* 
	save "$dt/table10", replace 

**** TABLE 10, COLUMN 1
	regress ics totpart elg mk96_s p1 sap* Y98sap* Istdtst* Istd2 if dpexam==0 & t1~=., robust cluster(sch98v1) 
	sum ics if totpart!=. & elg!=. & mk96_s!=. & p1!=. & sap1!=. & sap2!=. & sap3!=. & sap4!=. & Istd2!=.  ///
		& Y98sap1!=. & Y98sap2!=. & Y98sap3!=. & Y98sap4!=. & Istdtst2!=. & Istdtst3!=. ///
		& Istdtst4!=. &Istdtst5!=. & Istdtst6!=. & Istdtst7!=. & dpexam==0 & t1~=. 

**** TABLE 10, COLUMN 2: UPDATED
	summ popTR*updated if dpexam==0 & totpart~=. 
	regress ics t1 t2 elg 	mk96_s p1 sap* Y98sap* Istdtst* Istd2 ///
		popTR_3km_updated popT_3km_updated popTR_36k_updated popT_36k_updated ///
		if (dpexam==0 & totpart~=.), robust cluster(sch98v1) 
	* Average cross-school externality: UPDATED
		lincom 629*popTR_3km_updated + 1711*popTR_36k_updated 
	sum ics if totpart!=. & elg!=. & mk96_s!=. & p1!=. & sap1!=. & sap2!=. & sap3!=. & sap4!=. & Istd2!=. ///
		& Y98sap1!=. & Y98sap2!=. & Y98sap3!=. & Y98sap4!=. & Istdtst2!=. & Istdtst3!=. ///
		& Istdtst4!=. &Istdtst5!=. & Istdtst6!=. & Istdtst7!=. & dpexam==0 & t1~=. 
	
**** TABLE 10, COLUMN 2: ORIGINAL
	summ popTR*original if dpexam==0 & totpart~=. 
	regress ics t1 t2 elg 	mk96_s p1 sap* Y98sap* Istdtst* Istd2 ///
		popTR_3km_original popT_3km_original popTR_36k_original popT_36k_original ///
		if (dpexam==0 & totpart~=.), robust cluster(sch98v1) 
	* Average cross-school externality: ORIGINAL
		lincom 636*popTR_3km_original + 1190*popTR_36k_original 
		*matrix list e(V) 

	sum ics if totpart!=. & elg!=. & mk96_s!=. & p1!=. & sap1!=. & sap2!=. & sap3!=. & sap4!=. & Istd2!=. ///
		& Y98sap1!=. & Y98sap2!=. & Y98sap3!=. & Y98sap4!=. & Istdtst2!=. & Istdtst3!=. ///
		& Istdtst4!=. &Istdtst5!=. & Istdtst6!=. & Istdtst7!=. & dpexam==0 & t1~=. 

**** TABLE 10, COLUMN 3: UPDATED
	regress ics t1 t2 elg mk96_s p1 sap* Y98sap* Istdtst* Istd2 popTR_3km_updated popT_3km_updated popTR_36k_updated popT_36k_updated ///
		if (dpexam==0 & Ipupque~=. & totpart~=.), robust cluster(sch98v1) 
	sum ics if (dpexam==0 & Ipupque~=. & totpart~=. & elg!=. & mk96_s!=. & p1!=. & sap1!=. & sap2!=. & sap3!=. & sap4!=. ///
		& Istd2!=. & Y98sap1!=. & Y98sap2!=. & Y98sap3!=. & Y98sap4!=. & Istdtst2!=. ///
		& Istdtst3!=. & Istdtst4!=. &Istdtst5!=. & Istdtst6!=. & Istdtst7!=.) 

**** TABLE 10, COLUMN 3: ORIGINAL
	regress ics t1 t2 elg mk96_s p1 sap* Y98sap* Istdtst* Istd2 ///
		popTR_3km_original popT_3km_original popTR_36k_original popT_36k_original ///
		if (dpexam==0 & Ipupque~=. & totpart~=.), robust cluster(sch98v1) 
	sum ics if (dpexam==0 & Ipupque~=. & totpart~=. & elg!=. & mk96_s!=. & p1!=. & sap1!=. & sap2!=. & sap3!=. & sap4!=. ///
		& Istd2!=. & Y98sap1!=. & Y98sap2!=. & Y98sap3!=. & Y98sap4!=. & Istdtst2!=. ///
		& Istdtst3!=. & Istdtst4!=. &Istdtst5!=. & Istdtst6!=. & Istdtst7!=.) 
