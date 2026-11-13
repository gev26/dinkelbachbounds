/*MIGUEL AND KREMER (2004)*/
/*DATA ANALYSIS DO FILES: TABLE 3*/
/*VERSION: OCTOBER 2014*/

clear
set more off
version 12.1

* Merge Compliance data with Namelist data
	use "$da/namelist" 
	rename schid sch 
	drop Tmonths 
	reshape wide date sch std obs prs Isem*, i(pupid dupid) j(visit) 
	sort pupid 
	save "$dt/f1", replace 
	use "$da/comply" 
	sort pupid 
	merge pupid using "$dt/f1" 
	aorder
	
/* 	Schools mass treated with PZQ 1998: 158, 206, 207, 215, 220, 268
	Schools mass treated with PZQ 1999: 158, 206, 207, 215, 220, 268 
		AND 116, 182, 183, 203, 208, 221, 262, 263, 269, 284	*/ 
* 	Schools missing 1999 treatment information: 133, 274

**** TABLE 3, TOP PANEL
	
	* Drop individudals with missing group information
		drop if wgrp==. 

	* Girls < 13
		* Any medical treatment in 1998
			bys wgrp: tab any98 if elg98 == 1
		* Albendazole
			bys wgrp: tab a981  if elg98 == 1 
		* Praziquantel
			tab p98 if elg98 == 1 &  wgrp1==1 & psch98==1
			tab p98 if elg98 == 1 & (wgrp2==1) 
			tab p98 if elg98 == 1 & (wgrp3==1) 
		* Albendazole
			bys wgrp: tab a982  if elg98 == 1 

	* Girls >= 13
		* Any medical treatment in 1998
			bys wgrp: tab any98 if elg98==0 
		* Albendazole
			bys wgrp: tab a981 if elg98==0 
		* Praziquantel
			tab p98 if elg98==0 & wgrp1==1 & psch98==1
			tab p98 if elg98==0 & (wgrp2==1) 
			tab p98 if elg98==0 & (wgrp3==1) 
		* Albendazole
			bys wgrp: tab a982 if elg98==0 


**** TABLE 3, MIDDLE PANEL

	* Drop pupils in standard 8 in 1998
		drop if std981==8 

	* Girls < 13
		* Any medical treatment in 1999
			bys wgrp: tab any99 if elg99 == 1 
		* Albendazole
			bys wgrp: tab a991  if elg99 == 1 
		* Praziquantel
			bys wgrp: tab p99   if elg99 == 1 & psch99==1
		* Albendazole
			bys wgrp: tab a992  if elg99 == 1 

	* Girls >= 13
		* Any medical treatment in 1999
			bys wgrp: tab any99 if elg99 == 0 
		* Albendazole
			bys wgrp: tab a991  if elg99 == 0 
		* Praziquantel
			bys wgrp: tab p99   if elg99 == 0 & psch99==1
		* Albendazole
			bys wgrp: tab a992  if elg99 == 0


**** TABLE 3, BOTTOM PANEL

	* Drop pupils not present during 1999 school visits
		drop if (totprs99==0 | totprs99==.) 

	* Girls < 13
		* Any medical treatment in 1999
			bys wgrp: tab any99 if elg99 == 1 
		* Albendazole
			bys wgrp: tab a991  if elg99 == 1 
		* Praziquantel
			bys wgrp: tab p99   if elg99 == 1 & psch99==1 
		* Albendazole
			bys wgrp: tab a992  if elg99 == 1 

	* Girls >= 13
		* Any medical treatment in 1999
			bys wgrp: tab any99 if elg99 == 0 
		* Albendazole
			bys wgrp: tab a991  if elg99 == 0 
		* Praziquantel
			bys wgrp: tab p99   if elg99 == 0 & psch99==1 
		* Albendazole
			bys wgrp: tab a992  if elg99 == 0 
