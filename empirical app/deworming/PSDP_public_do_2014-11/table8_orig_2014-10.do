/*MIGUEL AND KREMER (2004)*/
/*DATA ANALYSIS DO FILES: TABLE 8 - ORIGINAL VERSION*/
/*VERSION: OCTOBER 2014*/

clear
set more off
version 12.1
#delimit ;

*Prepare Namelist data;
	use "$da\namelist_old";
	sort pupid visit;
	keep if visit==981;
	keep pupid std;
	rename std std98v1;
	sort pupid;
	save "$dt\f1", replace;

	use "$da\namelist_old";
	sort pupid;
	merge pupid using "$dt\f1";
	tab _merge; drop _merge;
	sort pupid;
	keep pupid visit std std98v1 prs wgrp* sch98v1 elg sex yrbirth;
	save "$dt\f1", replace;
	clear;

**** TABLE 8;

	* PANEL A, ROW 1;
	use "$dt\f1";
	keep if (visit>981 & visit<993 & elg==1);
	collapse prs wgrp* (count) np = pupid, by(sch98v1);
	summ prs [aw=np];
	bys wgrp: summ prs [aw=np];
	regress prs wgrp1 [aw=np];
	regress prs wgrp1 wgrp2 [aw=np];
	clear;

	* PANEL A, ROW 2;
	use "$dt\f1";
	keep if (visit>981 & visit<993 & elg==0);
	collapse prs wgrp* (count) np = pupid, by(sch98v1);
	summ prs [aw=np];
	bys wgrp: summ prs [aw=np];
	regress prs wgrp1 [aw=np];
	regress prs wgrp1 wgrp2 [aw=np];
	clear;

	* PANEL A, ROW 3;
	use "$dt\f1";
	keep if (visit>981 & visit<993 & elg==1 & (std98v1==0 | std98v1==1 | std98v1==2));
	collapse prs wgrp* (count) np = pupid, by(sch98v1);
	summ prs [aw=np];
	bys wgrp: summ prs [aw=np];
	regress prs wgrp1 [aw=np];
	regress prs wgrp1 wgrp2 [aw=np];
	clear;

	* PANEL A, ROW 4;
	use "$dt\f1";
	keep if (visit>981 & visit<993 & elg==1 & (std98v1==3 | std98v1==4 | std98v1==5));
	collapse prs wgrp* (count) np = pupid, by(sch98v1);
	summ prs [aw=np];
	bys wgrp: summ prs [aw=np];
	regress prs wgrp1 [aw=np];
	regress prs wgrp1 wgrp2 [aw=np];
	clear;

	* PANEL A, ROW 5;
	use "$dt\f1";
	keep if (visit>981 & visit<993 & elg==1 & (std98v1==6 | std98v1==7 | std98v1==8));
	collapse prs wgrp* (count) np = pupid, by(sch98v1);
	summ prs [aw=np];
	bys wgrp: summ prs [aw=np];
	regress prs wgrp1 [aw=np];
	regress prs wgrp1 wgrp2 [aw=np];
	clear;

	* PANEL A, ROW 6;
	use "$dt\f1";
	keep if (visit>981 & visit<993 & elg==1 & (std98v1==55));
	collapse prs wgrp* (count) np = pupid, by(sch98v1);
	summ prs [aw=np];
	bys wgrp: summ prs [aw=np];
	regress prs wgrp1 [aw=np];
	regress prs wgrp1 wgrp2 [aw=np];
	clear;

	* PANEL A, ROW 7;
	use "$dt\f1";
	keep if (visit>981 & visit<993 & elg~=. & sex==0);
	collapse prs wgrp* (count) np = pupid, by(sch98v1);
	summ prs [aw=np];
	bys wgrp: summ prs [aw=np];
	regress prs wgrp1 [aw=np];
	regress prs wgrp1 wgrp2 [aw=np];
	clear;

	* PANEL A, ROW 8;
	use "$dt\f1";
	keep if (visit>981 & visit<993 & elg~=. & sex==1);
	collapse prs wgrp* (count) np = pupid, by(sch98v1);
	summ prs [aw=np];
	bys wgrp: summ prs [aw=np];
	regress prs wgrp1 [aw=np];
	regress prs wgrp1 wgrp2 [aw=np];
	clear;


	* PANEL B, ROW 1;
	use "$dt\f1";
	keep if (visit>992 & visit<999 & elg==1);
	collapse prs wgrp* (count) np = pupid, by(sch98v1);
	summ prs [aw=np];
	bys wgrp: summ prs [aw=np];
	regress prs wgrp1 wgrp2 [aw=np];
	clear;

	* PANEL B, ROW 2;
	use "$dt\f1";
	keep if (visit>992 & visit<999 & elg==0);
	collapse prs wgrp* (count) np = pupid, by(sch98v1);
	summ prs [aw=np];
	bys wgrp: summ prs [aw=np];
	regress prs wgrp1 wgrp2 [aw=np];
	clear;

	* PANEL B, ROW 3;
	use "$dt\f1";
	keep if (visit>992 & visit<999 & elg==1 & (std98v1==0 | std98v1==1 | std98v1==2));
	collapse prs wgrp* (count) np = pupid, by(sch98v1);
	summ prs [aw=np];
	bys wgrp: summ prs [aw=np];
	regress prs wgrp1 wgrp2 [aw=np];
	clear;

	* PANEL B, ROW 4;
	use "$dt\f1";
	keep if (visit>992 & visit<999 & elg==1 & (std98v1==3 | std98v1==4 | std98v1==5));
	collapse prs wgrp* (count) np = pupid, by(sch98v1);
	summ prs [aw=np];
	bys wgrp: summ prs [aw=np];
	regress prs wgrp1 wgrp2 [aw=np];
	clear;

	* PANEL B, ROW 5;
	use "$dt\f1";
	keep if (visit>992 & visit<999 & elg==1 & (std98v1==6 | std98v1==7 | std98v1==8));
	collapse prs wgrp* (count) np = pupid, by(sch98v1);
	summ prs [aw=np];
	bys wgrp: summ prs [aw=np];
	regress prs wgrp1 wgrp2 [aw=np];
	clear;

	* PANEL B, ROW 6;
	use "$dt\f1";
	keep if (visit>992 & visit<999 & elg==1 & (std98v1==55));
	collapse prs wgrp* (count) np = pupid, by(sch98v1);
	summ prs [aw=np];
	bys wgrp: summ prs [aw=np];
	regress prs wgrp1 wgrp2 [aw=np];
	clear;

	* PANEL B, ROW 7;
	use "$dt\f1";
	keep if (visit>992 & visit<999 & elg~=. & sex==0);
	collapse prs wgrp* (count) np = pupid, by(sch98v1);
	summ prs [aw=np];
	bys wgrp: summ prs [aw=np];
	regress prs wgrp1 wgrp2 [aw=np];
	clear;

	* PANEL B, ROW 8;
	use "$dt\f1";
	keep if (visit>992 & visit<999 & elg~=. & sex==1);
	collapse prs wgrp* (count) np = pupid, by(sch98v1);
	summ prs [aw=np];
	bys wgrp: summ prs [aw=np];
	regress prs wgrp1 wgrp2 [aw=np];
	clear;
