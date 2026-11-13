/*MIGUEL AND KREMER (2004)*/
/*DATA ANALYSIS DO FILES: TABLE 1*/
/*VERSION: OCTOBER 2014*/

clear
set more off
version 12.1

* Start with Namelist data
	use "$da/namelist" 

* Each school is a distinct data point, weighted by number of pupils
	keep if visit==981 
	collapse sex elg98 stdgap yrbirth wgrp* (count) np=pupid, by (sch98v1) 

**** TABLE 1: PANEL A
	bys wgrp: summ sex elg98 stdgap yrbirth [aw=np] 
	foreach var in sex elg98 stdgap yrbirth { 
		regress `var' wgrp1 wgrp2 [aw=np] 
		} 

* Incorporate data from Pupil Questionnaire
	keep sch98v1 wgrp* 
	rename sch98v1 schid 
	sort schid 
	save "$dt/f1", replace 
	clear 
	use "$da/pupq" 

* Only keep pupils with 1998 data
	drop if pupdate_98_1=="" &  schid_98_2==. 

* Incorporate treatment group variable
	gen schid = schid_98_2 
	sort schid 
	merge schid using "$dt/f1" 
	tab _merge  
	drop _merge 

* Create measure of pre-program school attendance based on # days absent in previous four weeks
	gen preatt_98 = (20-absdays_98_6)/20 

* Create indicator for "Household Has Livestock"
	gen     Ilivestock_98 = (cows_98_23==1 | goats_98_24==1 | sheep_98_25==1 | pigs_98_26==1) 
	replace Ilivestock_98 = . if (cows_98_23==. | goats_98_24==. | sheep_98_25==. | pigs_98_26==.) 

* Create indicator for "Child Sick Often"
	gen     Isoften_98 = (fallsick_98_37==3) 
	replace Isoften_98 = . if fallsick_98_37==. 

* Create indicator for "Child Clean"
	gen     Iclean_98 = (clean_98_15==1) 
	replace Iclean_98 = . if clean_98_15==. 

* Each school is a distinct data point, weighted by number of pupils
collapse preatt_98 havelatr_98_33 Ilivestock_98 waz_98 bloodst_98_58 Isoften_98 malaria_98_48 Iclean_98 ///
	wgrp* (count) np = pupid, by(schid) 

**** TABLE 1: PANEL B
	bys wgrp: summ preatt_98 havelatr_98_33 Ilivestock_98 waz_98 bloodst_98_58 Isoften_98 malaria_98_48 Iclean_98 [aw=np] 
	foreach var in preatt_98 havelatr_98_33 Ilivestock_98 waz_98 bloodst_98_58 Isoften_98 malaria_98_48 Iclean_98 { 
		regress `var' wgrp1 wgrp2 [aw=np] 
		} 
	clear 

* Use School data
	use "$da/schoolvar" 

* Create worm group indicators
	gen wgrp1 = (wgrp==1) 
	gen wgrp2 = (wgrp==2) 
	gen wgrp3 = (wgrp==3) 

* Normalize 1996 mock tests to be in units of individual std dev, equivalent to 1998, 1999
	replace mk96_s = mk96_s*(0.4357)/(0.8318) 

**** TABLE 1: PANEL C
	bys wgrp: summ mk96_s distlake pup_pop latr_pup z_inf98 pop1_3km_updated pop1_36k_updated popT_3km_updated popT_36k_updated  
	foreach var in mk96_s distlake pup_pop latr_pup z_inf98 pop1_3km_updated pop1_36k_updated popT_3km_updated popT_36k_updated { 
			regress `var' wgrp1 wgrp2 
		} 
