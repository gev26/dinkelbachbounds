/*MIGUEL AND KREMER (2004)*/
/*DATA ANALYSIS DO FILES: TABLE 6*/
/*VERSION: OCTOBER 2014*/

clear
set more off
version 12.1

* Incorporate namelist information
	use "$da/namelist" 
	keep if visit==981 
	keep pupid elg98 elg99 sch98v1 wgrp* totpar98 stdgap std 
	sort pupid 
	save "$dt/f1", replace 
	clear 

*Incorporate compliance information	
	use "$da/comply" 
	keep pupid any98 any99 
	sort pupid 
	merge pupid using "$dt/f1" 
	tab _merge 
	drop if _merge==1 
	drop _merge 
	sort pupid 
	save "$dt/f1", replace 
	clear 

*Incorporate infection information	
	use "$da/wormed" 
	keep pupid numics98 numics99 hw99_ics al99_who sm99_who tt99_ics any_ics98 any_ics99 
	sort pupid 
	merge pupid using "$dt/f1" 
	tab _merge 
	drop if _merge==1 
	drop _merge 
	sort pupid 
	save "$dt/f1", replace 
	clear 

*Incorporate Pupil Questionnaire information	
	use "$da/pupq" 
	keep pupid havelatr_98_33 waz_98 malaria_98_48 clean_98_15 
	sort pupid 
	merge pupid using "$dt/f1" 
	tab _merge 
	drop if _merge==1 
	drop _merge 
	sort pupid 

* Construct indicator for "child clean"
	gen     Iclean_98 = (clean_98_15==1) 
	replace Iclean_98 = . if clean_98_15==. 

* Restrict to those with non-missing eligibility data
	keep if elg98~=. 
	save "$dt/f1", replace 

* Generate indicator for in sample both 1998, 1999
	gen     sample99 = . 
	replace sample99 = 0 if (any_ics98~= . & any_ics99==.) 
	replace sample99 = 1 if (any_ics98~= . & any_ics99~=.) 

* Only consider groups 1,2 children
	keep if wgrp==1 | wgrp==2 


**** TABLE 6, PANEL A
	
	* Any moderate-heavy infection, Group 1, Treated 1998
		summ any_ics98 if (wgrp==1 & any98==1) & sample99==1 & elg98==1
	
	* Any moderate-heavy infection, Group 1, Untreated 1998
		summ any_ics98 if (wgrp==1 & any98==0) & sample99==1 & elg98==1
	
	* Proportion of 1998 parasitological, Group 1, Treated 1998
		summ sample99 if (wgrp==1 & any98==1) & any_ics98~=. & elg98==1
		
	* Proportion of 1998 parasitological, Group 1, Untreated 1998
		summ sample99 if (wgrp==1 & any98==0) & any_ics98~=.  & elg98==1

	foreach var in havelatr_98_33 stdgap waz_98 malaria_98_48 Iclean_98 { 
		summ `var' if (wgrp==1 & any98==1) & any_ics99~=. & elg98==1
		summ `var' if (wgrp==1 & any98==0) & any_ics99~=. & elg98==1 
		summ `var' if (wgrp==2 & any99==1) & any_ics99~=. & elg98==1
		summ `var' if (wgrp==2 & any99==0) & any_ics99~=. & elg98==1 

	* G1 treated 1998 vs G2 treated 1999
		regress `var' wgrp1 if ((wgrp==1 & any98==1) | (wgrp==2 & any99==1)) & any_ics99~=. & elg98==1, robust cluster(sch98v1) 
		
	* G1 untreated 1998 vs G2 untreated 1999
		regress `var' wgrp1 if ((wgrp==1 & any98==0) | (wgrp==2 & any99==0)) & any_ics99~=. & elg98==1, robust cluster(sch98v1) 
		}


**** TABLE 6, PANEL B

	*Girls <13 years, and all boys
		foreach var in any_ics99 hw99_ics al99_who sm99_who tt99_ics { 
			summ `var' if (wgrp==1 & any98==1) & any_ics99~=. & elg98==1
			summ `var' if (wgrp==1 & any98==0) & any_ics99~=. & elg98==1 
			summ `var' if (wgrp==2 & any99==1) & any_ics99~=. & elg98==1 
			summ `var' if (wgrp==2 & any99==0) & any_ics99~=. & elg98==1 
			regress `var' wgrp1 if ((wgrp==1 & any98==1) | (wgrp==2 & any99==1)) & any_ics99~=. & elg98==1, robust cluster(sch98v1) 
			regress `var' wgrp1 if ((wgrp==1 & any98==0) | (wgrp==2 & any99==0)) & any_ics99~=. & elg98==1, robust cluster(sch98v1)
		} 

	*Girls >=13 years
		summ any_ics98 if (wgrp==1 & any98==1) & sample99==1 & elg98==0
		summ any_ics98 if (wgrp==1 & any98==0) & sample99==1 & elg98==0
		summ any_ics99 if (wgrp==1 & any98==1) & any_ics99~=. & elg98==0
		summ any_ics99 if (wgrp==1 & any98==0) & any_ics99~=. & elg98==0
		summ any_ics99 if (wgrp==2 & any99==1) & any_ics99~=. & elg98==0
		summ any_ics99 if (wgrp==2 & any99==0) & any_ics99~=. & elg98==0
		regress any_ics99 wgrp1 if ((wgrp==1 & any98==1) | (wgrp==2 & any99==1)) & any_ics99~=. & elg98==0, robust cluster(sch98v1) 
		regress any_ics99 wgrp1 if ((wgrp==1 & any98==0) | (wgrp==2 & any99==0)) & any_ics99~=. & elg98==0, robust cluster(sch98v1)


* Incorporate school participation information
	use "$da/namelist", clear 
	sort pupid 
	keep pupid visit std98v1 prs wgrp* sch98v1 elg98 
	keep if (visit>981 & visit<993) 
	save "$dt/temp", replace 

*Incorporate compliance information	
	use "$da/comply", clear 
	keep pupid any98 any99 
	sort pupid 
	merge pupid using "$dt/temp" 
	tab _merge 
	drop if _merge==1 
	drop _merge 

*Exclude pupils	not in grades 1-7
	drop if std<1 | std>7 
	keep if elg98==1 
	save "$dt/f2", replace 

	
**** TABLE 6, PANEL C	

	keep if wgrp==1 & any98==1 
	collapse prs (count) np=pupid, by(sch98v1) 
	sum prs [aw=np] 
	clear 

	use "$dt/f2" 
	keep if wgrp==1 & any98==0 
	collapse prs (count) np=pupid, by(sch98v1) 
	sum prs [aw=np] 
	clear 

	use "$dt/f2" 
	keep if wgrp==2 & any99==1 
	collapse prs (count) np=pupid, by(sch98v1) 
	sum prs [aw=np] 
	clear 

	use "$dt/f2" 
	keep if wgrp==2 & any99==0 
	collapse prs (count) np=pupid, by(sch98v1) 
	sum prs [aw=np] 
	clear 

	use "$dt/f2" 
	keep if ((wgrp==1 & any98==1) | (wgrp==2 & any99==1)) 
	collapse prs wgrp1 (count) np=pupid, by(sch98v1) 
	regress prs wgrp1 [aw=np], robust 
	clear 

	use "$dt/f2" 
	keep if ((wgrp==1 & any98==0) | (wgrp==2 & any99==0)) 
	collapse prs wgrp1 (count) np=pupid, by(sch98v1) 
	regress prs wgrp1 [aw=np], robust 
