/*MIGUEL AND KREMER (2004)*/
/*DATA ANALYSIS DO FILES: TABLE 5*/
/*VERSION: OCTOBER 2014*/

clear
set more off
version 12.1

* Incorporate data on eligibility and parasitological exams
	use "$da/namelist" 
	keep if visit==981 
	keep pupid sch98v1 wgrp elg98 elg99 
	sort pupid 
	save "$dt/f1", replace 

	use "$da/wormed" 
	sort pupid 
	merge pupid using "$dt/f1" 
	tab _merge 
	keep if _merge==3 
	drop _merge 

* Generate low HB indicator
	gen hb100 = (hb<100) 
	replace hb100=. if hb==. 

* Restrict sample appropriately
	keep if any_ics99~=. 

**** TABLE 5, PANEL A 
	summ any_ics98 any_ics99 hw99_ics al99_who sm99_who tt99_ics if wgrp==1
	summ any_ics98 any_ics99 hw99_ics al99_who sm99_who tt99_ics if wgrp==2
	foreach var in any_ics99 hw99_ics al99_who sm99_who tt99_ics { 
		regress `var' wgrp1 if (wgrp==1 | wgrp==2), robust cluster(sch98v1) 
		} 

**** TABLE 5, PANEL B: HB AND ANEMIA ONLY
	summ hb hb10 if wgrp==1
	summ hb hb10 if wgrp==2
	foreach var in hb hb10 { 
		regress `var' wgrp1 if (wgrp==1 | wgrp==2), robust cluster(sch98v1) 
		} 

* Incorporate data on treatment group and from pupil questionnaire
	use "$da/namelist", clear 
	keep if visit==981 
	keep pupid sch98v1 wgrp* 
	sort pupid 
	save "$dt/f1", replace 

	use "$da/pupq" 
	sort pupid 
	merge pupid using "$dt/f1" 
	tab _merge  
	drop _merge 

* Indicator for "child sick often"
	gen     Isoften_99 = (soften_99_39==3) 
	replace Isoften_99 = . if soften_99_39==. 

* Indicator for "child clean"
	gen     Iclean_99 = (clean_99_13==1) 
	replace Iclean_99 = . if clean_99_13==. 

* Indicator for "child wears shoes"
	gen     Ishoes_99 = (shoes_99_10==1 | shoes_99_10==2) 
	replace Ishoes_99 = . if shoes_99_10==. 

	keep wgrp* sch98v1 slastwk_99_40 Isoften_99 haz99 waz_99 Iclean_99 Ishoes_99 dayswat_99_36 
	
**** TABLE 5, PANELS B AND C
	summ slastwk_99_40 Isoften_99 haz99 waz_99 Iclean_99 Ishoes_99 dayswat_99_36 if wgrp==1
	summ slastwk_99_40 Isoften_99 haz99 waz_99 Iclean_99 Ishoes_99 dayswat_99_36 if wgrp==2 | wgrp==3
	foreach var in slastwk_99_40 Isoften_99 haz99 waz_99 Iclean_99 Ishoes_99 dayswat_99_36 { 
		regress `var' wgrp1, robust cluster(sch98v1) 
		} 
