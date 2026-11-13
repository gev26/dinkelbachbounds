/*MIGUEL AND KREMER (2004)*/
/*DATA ANALYSIS DO FILES: TABLE 2*/
/*VERSION: OCTOBER 2014*/

clear
set more off
version 12.1

* Incorporate child gender and age information with parasitological exam data
	use "$da/namelist" 
	collapse sex yrbirth sch98v1, by (pupid) 
	sort pupid 
	save "$dt/f1", replace 
	use "$da/wormed" 
	drop sch98v1 
	sort pupid 
	merge pupid using "$dt/f1" 
	tab _merge  
	drop _merge 
	save "$dt/f2", replace 
	clear 

* Incorporate distance to lake and treatment group information
	use "$da/schoolvar" 
	keep schid distlake wgrp 
	rename schid sch98v1 
	sort sch98v1 
	save "$dt/f1", replace 
	clear 
	use "$dt/f2" 
	drop wgrp1 
	sort sch98v1 
	merge sch98v1 using "$dt/f1" 
	tab _merge  
	drop _merge 
	drop if hw98==. 

* Change units for average infection intensity variables from 100 milligrams to grams
	replace hw98 = hw98*10 
	replace al98 = al98*10 
	replace sm98 = sm98*10 
	replace tt98 = tt98*10 

**** TABLE 1, Column (1): Prevalence of infection
	* Hookwork
		summ any_hw98 
	* Roundworm
		summ any_al98 	
	* Schistosomiasis
		summ any_sm98 
	* Schistosomiasis, schools < 5km
		summ any_sm98 if distlake < 5 
	* Whipworm
		summ any_tt98 
	* At least one infection
		summ any_98 
		* Born since 1985
			summ any_98 if (yrbirth>=1985 & yrbirth~=.) 
		* Born before 1985
			summ any_98 if yrbirth<1985 
		* Female
			summ any_98 if sex==0
		* Male
			summ any_98 if sex==1
	* At least two & three infections
		gen atleast2=0 if numinf98!=.
		replace atleast2=1 if numinf98>=2 & numinf98!=.
		tab atleast2
		gen atleast3=0 if numinf98!=.
		replace atleast3=1 if numinf98>=3 & numinf98!=.
		tab atleast3

**** TABLE 1, Column (2): Prevalence of moderate-heavy infection
	* Hookworm
		summ hw98_ics 
	* Roundworm
		summ al98_who 
	* Schistosomiasis
		summ sm98_who 
	* Schistosomiasis, schools < 5km
		summ sm98_who if distlake < 5 
	* Whipworm
		summ tt98_ics 
	* At least one infection
		summ any_ics98 
		* Born since 1985
			summ any_ics98 if (yrbirth>=1985 & yrbirth~=.) 
		* Born before 1985
			summ any_ics98 if yrbirth<1985 
		* Female 
			summ any_ics98 if sex==0
		* Male
			summ any_ics98 if sex==1
	* At least two & three infections
		gen atleast2i=0 if numics98!=.
		replace atleast2i=1 if numics98>=2 & numics98!=.
		tab atleast2i
		gen atleast3i=0 if numics98!=.
		replace atleast3i=1 if numics98>=3 & numics98!=.
		tab atleast3i

**** TABLE 1, Column (3): Average worm load
	* Hookworm
		summ hw98 
	* Roundworm
		summ al98 
	* Schistosomiasis
		summ sm98 
	* Schistosomiasis, schools < 5km
		summ sm98 if distlake < 5 
	* Whipworm
		summ tt98
