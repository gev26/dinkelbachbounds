/*MIGUEL AND KREMER (2004)*/
/*DATA ANALYSIS DO FILES: TABLE 9 AND APPENDIX 4 - ORIGINAL VERSION*/
/*VERSION: OCTOBER 2014*/

clear
set more off
version 12.1
#delimit ;

* Incorporate Namelist data;
	use "$da\namelist_old";

* Create further school assistance controls; 
	gen Y98 = 0;
	replace Y98 = 1 if (visit>980 & visit<990);
	gen Y98sap1 = sap1*Y98;
	gen Y98sap2 = sap2*Y98;
	gen Y98sap3 = sap3*Y98;
	gen Y98sap4 = sap4*Y98;
	save "$dt\f9", replace;

* Incorporate school and zonal variables;
	use "$da\schoolvar";
	keep schid mk96_s z9899* pop1_3km* pop2_3km* popT_3km* pop1_36k* pop2_36k* popT_36k* distlake;
	rename schid sch98v1;
	sort sch98v1;
	save "$dt\f1", replace;
	clear;
	use "$dt\f9";
	sort sch98v1;
	merge sch98v1 using "$dt\f1";
	sum _merge mk96_s; drop _merge;

* Normalize and adjust mock scores to individual units;
	replace mk96_s = mk96*0.4357/0.8318;

* Generate year measure;
	gen yr = .;
	replace yr = 1 if (visit>=981 & visit<993); 
	replace yr = 2 if (visit>992 & visit<999);

* Create treatment indicators;
	* First year of treatment;
	gen t1 = 0;
	replace t1 = 1 if (wgrp==1 & visit>981 & visit<993);
	replace t1 = 1 if (wgrp==2 & visit>992 & visit<999);
	replace t1 = . if wgrp==.;
	* Second year of treatment;
	gen t2 = 0;
	replace t2 = 1 if (wgrp==1 & visit>992 & visit<999);
	replace t2 = .  if wgrp==.;

* Other indicators;
	gen t1e = elg*t1;
	gen t2e = elg*t2;

* Create standard-specific measure of zonal infection rate;
	gen p1 = z9899_34;
	replace p1 = z9899_56 if (std98v1==5 | std98v1==6);
	replace p1 = z9899_78 if (std98v1==7 | std98v1==8);
	drop z9899*;

* Create standard indicators, based on 1998 visit 1 standard;
	gen std_fs = std98v1 if (std98v1>-1 & std98v1<9);
	replace std_fs = -1 if (std98v1==55);
	tab std_fs, gen(Istd);
	summ Istd*;
	drop Istd10 std_fs;
	save "$dt\f9", replace;
	clear;

* Incorporate compliance data;
	use "$da\comply_old";
	keep pupid any98 any99;
	sort pupid;
	save "$dt\f2", replace;
	clear;
	use "$dt\f9";
	sort pupid;
	merge pupid using "$dt\f2";
	tab _merge; drop _merge;

* THREE PERIODS: pre-treatment (981), Year 1 (982-992), Year 2 (993-998);
	* First year of treatment indicator;
	gen treat_y1 = 0 if visit == 981;
	replace treat_y1 = any98 if (visit>981 & visit<993);
	replace treat_y1 = any99 if (any98==0 & visit>992 & visit<999);
	replace treat_y1 = 0 if (any98==1 & visit>992 & visit<999);
	* Second year of treatment indicator;
	* If treated in first and second years;
	gen treat_y2 = 0 if visit<993;
	replace treat_y2 = 1 if (any98==1 & any99==1 & visit>992 & visit<999);
	* Pupils who were treated in 1998 and not in 1999 are counted as not treated;
	replace treat_y2 = 0 if (any98==1 & any99==0 & visit>992 & visit<999);
	replace treat_y2 = 0 if (any98==0 & visit>992 & visit<999);
	compress;
	save "$dt\f9", replace;

* Collapse data by pupil and year, where YEAR1 = 982-992, YEAR2 = 983-998;
	sort pupid yr;
	collapse (mean) sch98v1 prs t1 t2 elg p1 mk96_s 
		Y98sap1 Y98sap2 Y98sap3 Y98sap4 sap1 sap2 sap3 sap4 
		Istd1 Istd2 Istd3 Istd4 Istd5 Istd6 Istd7 Istd8 Istd9 
		Isem1 Isem2 Isem3 pop1_3km* pop2_3km* pop1_36k* pop2_36k* popT_36k* popT_3km* 
		any98 any99 wgrp (sum) obs 
		if (t1~=. & elg~=. & sch98v1~=. &  mk96_s~=. &  p1~=. & Istd2~=. & pop1_3km_original~=.), by(pupid yr);
	keep e* t* p* sap* Y98sap* sch98v1 prs* mk96_s Istd* Isem* pop* sch* obs yr;

* Create measure of population in the vicinity in the year of treatment;
	gen pop_3km_original = pop1_3km_original;
	gen pop_36k_original = pop1_3km_original; 
		/*Note that this line above has a coding error: should be gen pop_36k_original = pop1_36k_original. 
		Leaving in coding error in order to replicate original resuls*/;
	replace pop_3km_original = pop_3km_original + pop2_3km_original if yr==2;
	replace pop_36k_original = pop_36k_original + pop2_36k_original if yr==2;

* Create an indicator for whether the school received treatment;
	gen t_any = 0;
	replace t_any=1 if (t1==1 | t2==1);
	replace t_any=. if t1==. | t2==.;
	save "$dt\table9a", replace;

**** TABLE 9, COLUMN 1;	
	sum prs [aw=obs] if  (t1~=. & elg~=. & sch98v1~=. & mk96_s~=. & p1~=. & Istd2~=. & pop1_3km_original~=.);
	regress prs t_any elg p1 mk96_s Y98sap* sap* Istd* Isem* [aw=obs] if (t1~=. & elg~=. & sch98v1~=. & mk96_s~=. & p1~=. & Istd2~=. 
		& pop1_3km_original~=.), robust cluster(sch98v1);

**** TABLE 9, COLUMN 2;	
regress prs t1 t2 elg p1 mk96_s Y98sap* sap* Istd* Isem* [aw=obs] if (t1~=. & elg~=. & sch98v1~=. & mk96_s~=. & p1~=. & Istd2~=. 
	& pop1_3km_original~=.), robust cluster(sch98v1);

**** TABLE 9, COLUMN 3;	
regress prs t1 t2 elg p1 mk96_s Y98sap* sap* Istd* Isem* pop_3km_original popT_3km_original pop_36k_original popT_36k_original [aw=obs] 
	if (t1~=. & elg~=. & sch98v1~=. & mk96_s~=. & p1~=. & Istd2~=. & pop1_3km_original~=.), robust cluster(sch98v1);

* Locals for population density for school attendance sample;
	summ pop_3km_original [aw=obs] if (t1~=. & elg~=. & sch98v1~=. & mk96_s~=. & p1~=. & Istd2~=. & pop1_3km_original~=.);
	local num_3km = round(r(mean));
	summ pop_36k_original [aw=obs] if (t1~=. & elg~=. & sch98v1~=. & mk96_s~=. & p1~=. & Istd2~=. & pop1_3km_original~=.);
	local num_36k = round(r(mean));

**** AVERAGE SPILLOVER EFFECT CALCULATIONS;
	use "$dt\table9a", clear;
	regress prs t_any elg p1 mk96_s Y98sap* sap* Istd* Isem* pop_3km_original popT_3km_original pop_36k_original popT_36k_original [aw=obs] 
		if (t1~=. & elg~=. & sch98v1~=. & mk96_s~=. & p1~=. & Istd2~=. & pop1_3km_original~=.), robust cluster(sch98v1);
	regress prs t_any elg p1 mk96_s Y98sap* sap* Istd* Isem* pop_3km_original popT_3km_original pop_36k_original popT_36k_original [aw=obs] if (t1~=. & elg~=. & sch98v1~=. & mk96_s~=. & p1~=. & Istd2~=. & pop1_3km_original~=.), robust cluster(sch98v1); 
	matrix V = e(V); 
	matrix B = e(b) ;
	gen C = 1; 
	* Take average RHS values;
	collapse (mean) t_any elg p1 mk96_s Y98sap* sap* Istd* Isem* pop_3km_original popT_3km_original pop_36k_original popT_36k_original C [aw=obs] if (t1~=. & elg~=. & sch98v1~=. & mk96_s~=. & p1~=. & Istd2~=. & pop1_3km_original~=. & prs~=.); 
	summ; 
	* First for treatment school;
	replace t_any= 0 ;
	mkmat t_any-C, mat(T1) ;
	* Then for comparison schools;
	replace t_any= 0 ;
	replace pop_3km_original = 0 ;
	replace pop_36k_original = 0 ;
	mkmat t_any-C, mat(T0) ;
	* Calculation;
	matrix t1 = T1*B' ;
	matrix t0 = T0*B' ;
	matrix tau = t1[1,1] - t0[1,1] ;
	matrix tau_v = (T1-T0)*V*(T1-T0)'; 
	matrix tau_se = sqrt(tau_v[1,1]) ;
	matrix list tau;  /*Average spillover effect*/
	matrix list tau_se;  /*Std dev*/
	clear ;

**** AVERAGE SPILLOVER EFFECT AND DIRECT EFFECT;
	use "$dt\table9a", clear;
	regress prs t_any elg p1 mk96_s Y98sap* sap* Istd* Isem* pop_3km_original popT_3km_original pop_36k_original popT_36k_original [aw=obs] 
		if (t1~=. & elg~=. & sch98v1~=. & mk96_s~=. & p1~=. & Istd2~=. & pop1_3km_original~=.), robust cluster(sch98v1);
	matrix V = e(V);
	matrix B = e(b);
	gen C = 1;
	* Take average RHS values;
	collapse (mean) t_any elg p1 mk96_s Y98sap* sap* Istd* Isem* pop_3km_original popT_3km_original pop_36k_original popT_36k_original C 
		[aw=obs] if (t1~=. & elg~=. & sch98v1~=. & mk96_s~=. & p1~=. & Istd2~=. & pop1_3km_original~=. & prs~=.);
	* First for treatment schools;
	replace t_any= 1;
	mkmat t_any-C, mat(T1);
	* Then for comparison schools;	
	replace t_any= 0;
	replace pop_3km_original = 0;
	replace pop_36k_original = 0;
	mkmat t_any-C, mat(T0);
	* Calculation; 
	matrix t1 = T1*B';
	matrix t0 = T0*B';
	matrix tau = t1[1,1] - t0[1,1];
	matrix tau_v = (T1-T0)*V*(T1-T0)';
	matrix tau_se = sqrt(tau_v[1,1]);
	matrix list tau;
	matrix list tau_se;
	clear;

* Incorporate data on number of infections in 1999;
	use "$da\wormed";
	keep pupid numics99;
	gen ics99 = numics99;
	replace ics99 = 1 if (numics>1 & numics<5);
	sort pupid;
	save "$dt\f4", replace;
	clear;
	use "$dt\f9";
	sort pupid;
	merge pupid using "$dt\f4";
	tab _merge; drop _merge;
	sort pupid;
	save "$dt\f9", replace;

* Only consider year 1 of treatment;
	keep if visit<993 & visit>981;
	sort pupid yr;
	collapse (mean) sch98v1 prs t1 t2 elg p1 mk96_s t1e t2e 
		Y98sap1 Y98sap2 Y98sap3 Y98sap4 sap1 sap2 sap3 sap4 
		Istd1 Istd2 Istd3 Istd4 Istd5 Istd6 Istd7 Istd8 Istd9 
		pop1_3km_original pop2_3km_original pop1_36k_original pop2_36k_original popT_36k_original popT_3km_original 
		ics99 any98 any99 wgrp distlake
		(sum) obs
		if (t1~=. & elg~=. & sch98v1~=. &  mk96_s~=. &  p1~=. & Istd2~=. & pop1_3km_original~=.), by(pupid yr);

* Generate selection into treatment indicator variable, and interaction with first year treatment;
	gen select=.;
	replace select=1 if (wgrp==1 & any98==1) | (wgrp==2 & any99==1);
	replace select=0 if (wgrp==1 & any98==0) | (wgrp==2 & any99==0);
	gen t1_select = t1*select;
	save "$dt\table9b", replace;

**** TABLE 9, COLUMN 4;
	sum prs [aw=obs] if (t1~=. & elg~=. & sch98v1~=. & mk96_s~=. & p1~=. & Istd2~=. & pop1_3km_original~=. & select~=. & (wgrp==1 | wgrp==2));
	regress prs t1 elg p1 mk96_s Y98sap* Istd* [aw=obs] if (t1~=. & elg~=. & sch98v1~=. & mk96_s~=. & p1~=. & Istd2~=. & pop1_3km_original~=. 
		& select~=. & (wgrp==1 | wgrp==2)), robust cluster(sch98v1);

**** TABLE 9, COLUMN 5;
regress prs t1 select t1_select
         p1 mk96_s Y98sap* Istd* pop1_3km_original popT_3km_original pop1_36k_original popT_36k_original 
         [aw=obs]
         if (t1~=. & elg~=. & sch98v1~=. & mk96_s~=. & p1~=. & Istd2~=. & pop1_3km_original~=. & select~=.  & (wgrp==1 | wgrp==2)),
         robust cluster(sch98v1);

**** TABLE 9, COLUMN 6;
sum prs [aw=obs] if (t1~=. & elg~=. & sch98v1~=. & mk96_s~=. & p1~=. & Istd2~=. & ics99~=. & pop1_3km_original~=.);
regress prs ics99 elg mk96_s p1 Y98sap* Istd* popT_36k_original popT_3km_original [aw=obs]
	if (t1~=. & elg~=. & sch98v1~=. & mk96_s~=. & p1~=. & Istd2~=. & ics99~=. & pop1_3km_original~=.),
	robust cluster(sch98v1);

**** TABLE 9, COLUMN 7 - SCHOOL AVERAGES;
collapse prs ics99 t1 pop1_3km_original pop1_36k_original popT_36k_original popT_3km_original elg mk96_s p1 Y98sap* Istd* wgrp (sum) obs
	if (t1~=. & elg~=. & sch98v1~=. & mk96_s~=. & p1~=. & Istd2~=. & ics99~=. & pop1_3km_original~=.),
	by (sch98v1);
regress prs ics99 elg mk96_s p1 Y98sap* Istd* popT_36k_original popT_3km_original
	(t1 pop1_3km_original pop1_36k_original popT_36k_original popT_3km_original elg mk96_s p1 Y98sap* Istd*)
	[aw=obs] if (wgrp==1 | wgrp==2), robust;
clear;
use "$dt\table9b";

**** APPENDIX TABLE A4, COLUMN 3;
regress prs t1 select t1_select
         p1 mk96_s Y98sap* Istd* pop1_3km_original popT_3km_original pop1_36k_original popT_36k_original 
         [aw=obs]
         if (t1~=. & elg~=. & sch98v1~=. & mk96_s~=. & p1~=. & Istd2~=. & pop1_3km_original~=. & select~=.  & (wgrp==1 | wgrp==2)),
         robust cluster(sch98v1);

**** APPENDIX TABLE A4, COLUMN 4;
regress prs t1 select t1_select
         p1 mk96_s Y98sap* Istd* pop1_3km_original popT_3km_original pop1_36k_original popT_36k_original 
         (t1 elg t1e
         p1 mk96_s Y98sap* Istd* pop1_3km_original popT_3km_original pop1_36k_original popT_36k_original)
         [aw=obs]
         if (t1~=. & elg~=. & sch98v1~=. & mk96_s~=. & p1~=. & Istd2~=. & pop1_3km_original~=. & select~=.
         & (wgrp==1 | wgrp==2)),
         robust cluster(sch98v1);
