/*MIGUEL AND KREMER (2004)*/
/*DATA ANALYSIS DO FILES: TABLE 4*/
/*VERSION: OCTOBER 2014*/

clear
set more off
#delimit ;
version 12.1;

* Incorporate information on pupils;
use "$da\namelist";
rename schid sch;
drop Tmonths;
reshape wide date sch std obs prs Isem*, i(pupid dupid) j(visit);
keep pupid sch98* sch99*;
save "$dt\f2", replace;
clear;

* Incorporate information on treatment group of each school;
use "$da\schoolvar";
keep schid wgrp;
gen W1=(wgrp==1);
gen W2=(wgrp==2);
gen W3=(wgrp==3);
drop wgrp;
save "$dt\f1", replace;
clear;

* Assign treatment group to school attended each visit;
foreach num in 98 99 {;
	use "$dt\f1";
	rename schid sch`num'1;
	sort sch`num'1;
	save "$dt\f1", replace;
	clear;
	local i=1;
	while (`i'<9)
		{;
		use "$dt\f2";
		sort sch`num'`i';
		merge sch`num'`i' using "$dt\f1";
		tab _merge; drop _merge;
		rename W1 W1`num'`i';
		rename W2 W2`num'`i';
		rename W3 W3`num'`i';
		sort sch`num'`i';
		save "$dt\f2", replace;
		clear;
		local j=`i'+1;
		use "$dt\f1";
		rename sch`num'`i' sch`num'`j';
		sort sch`num'`j';
		save "$dt\f1", replace;
		clear;
		local i=`j';
		};
	use "$dt\f1";
	rename sch`num'9 schid;
	save "$dt\f1", replace;
	};

use "$dt\f2";

* Generate indicators for group 1, 2 and 3 schools by year;
gen G198 = max(W1982, W1983, W1984, W1985, W1986, W1987, W1988);
gen G298 = max(W2982, W2983, W2984, W2985, W2986, W2987, W2988);
gen G398 = max(W3982, W3983, W3984, W3985, W3986, W3987, W3988);
gen G199 = max(W1991, W1992, W1993, W1994, W1995, W1996, W1997, W1998);
gen G299 = max(W2991, W2992, W2993, W2994, W2995, W2996, W2997, W2998);
gen G399 = max(W3991, W3992, W3993, W3994, W3995, W3996, W3997, W3998);


** GROUP 1 OUT-TRANSFERS;

* Transfers of group 1 pupils to other group 1 schools during 1998;
gen trs11_98 = 0;
replace trs11_98 = 1 if
	((W1981==1 & W1982==1 & sch981~=sch982) |
	(W1981==1 & W1983==1 & sch981~=sch983) |
	(W1981==1 & W1984==1 & sch981~=sch984) |
	(W1981==1 & W1985==1 & sch981~=sch985) |
	(W1981==1 & W1986==1 & sch981~=sch986) |
	(W1981==1 & W1987==1 & sch981~=sch987) |
	(W1981==1 & W1988==1 & sch981~=sch988) );
replace trs11_98 = . if (W1981==. | W1981==0 | G198==.);

* Transfers of group 1 pupils to group 2 schools during 1998;
gen trs12_98 = .;
replace trs12_98 = 1 if (W1981==1 & G298==1);
replace trs12_98 = 0 if (W1981==1 & G298==0);

gen trs13_98 = .;
replace trs13_98 = 1 if (W1981==1 & G398==1);
replace trs13_98 = 0 if (W1981==1 & G398==0);

* Transfers of group 1 pupils to other group 1 schools during 1999;
gen trs11_99 = 0;
replace trs11_99 = 1 if
	((W1981==1 & W1991==1 & sch981~=sch991) |
	(W1981==1 & W1992==1 & sch981~=sch992) |
	(W1981==1 & W1993==1 & sch981~=sch993) |
	(W1981==1 & W1994==1 & sch981~=sch994) |
	(W1981==1 & W1995==1 & sch981~=sch995) |
	(W1981==1 & W1996==1 & sch981~=sch996) |
	(W1981==1 & W1997==1 & sch981~=sch997) |
	(W1981==1 & W1998==1 & sch981~=sch998) );
replace trs11_99 = . if (W1981==. | W1981==0 | G199==.);

gen trs12_99 = .;
replace trs12_99 = 1 if (W1981==1 & G299==1);
replace trs12_99 = 0 if (W1981==1 & G299==0);

gen trs13_99 = .;
replace trs13_99 = 1 if (W1981==1 & G399==1);
replace trs13_99 = 0 if (W1981==1 & G399==0);


** GROUP 2 OUT-TRANSFERS;

* Transfers of group 2 pupils to other group 2 schools during 1998;
gen trs22_98 = 0;
replace trs22_98 = 1 if
	((W2981==1 & W2982==1 & sch981~=sch982) |
	(W2981==1 & W2983==1 & sch981~=sch983) |
	(W2981==1 & W2984==1 & sch981~=sch984) |
	(W2981==1 & W2985==1 & sch981~=sch985) |
	(W2981==1 & W2986==1 & sch981~=sch986) |
	(W2981==1 & W2987==1 & sch981~=sch987) |
	(W2981==1 & W2988==1 & sch981~=sch988) );
replace trs22_98 = . if (W2981==. | W2981==0 | G298==.);

* Transfers of group 2 pupils to group 1 schools during 1998;
gen trs21_98 = .;
replace trs21_98 = 1 if (W2981==1 & G198==1);
replace trs21_98 = 0 if (W2981==1 & G198==0);

gen trs23_98 = .;
replace trs23_98 = 1 if (W2981==1 & G398==1);
replace trs23_98 = 0 if (W2981==1 & G398==0);

* Transfers of group 2 pupils to other group 2 schools during 1999;
gen trs22_99 = 0;
replace trs22_99 = 1 if 
	((W2981==1 & W2991==1 & sch981~=sch991) |
	(W2981==1 & W2992==1 & sch981~=sch992) |
	(W2981==1 & W2993==1 & sch981~=sch993) |
	(W2981==1 & W2994==1 & sch981~=sch994) |
	(W2981==1 & W2995==1 & sch981~=sch995) |
	(W2981==1 & W2996==1 & sch981~=sch996) |
	(W2981==1 & W2997==1 & sch981~=sch997) |
	(W2981==1 & W2998==1 & sch981~=sch998) );
replace trs22_99 = . if (W2981==. | W2981==0 | G299==.);

gen trs21_99 = .;
replace trs21_99 = 1 if (W2981==1 & G199==1);
replace trs21_99 = 0 if (W2981==1 & G199==0);

gen trs23_99 = .;
replace trs23_99 = 1 if (W2981==1 & G399==1);
replace trs23_99 = 0 if (W2981==1 & G399==0);


** GROUP 3 OUT-TRANSFERS;

* Transfers of group 3 pupils to other group 3 schools during 1998;
gen trs33_98 = 0;
replace trs33_98 = 1 if
	((W3981==1 & W3982==1 & sch981~=sch982) |
	(W3981==1 & W3983==1 & sch981~=sch983) |
	(W3981==1 & W3984==1 & sch981~=sch984) |
	(W3981==1 & W3985==1 & sch981~=sch985) |
	(W3981==1 & W3986==1 & sch981~=sch986) |
	(W3981==1 & W3987==1 & sch981~=sch987) |
	(W3981==1 & W3988==1 & sch981~=sch988) );
replace trs33_98 = . if (W3981==. | W3981==0 | G398==.);

* Transfers of group 3 pupils to group 1 schools during 1998;
gen trs31_98 = .;
replace trs31_98 = 1 if (W3981==1 & G198==1);
replace trs31_98 = 0 if (W3981==1 & G198==0);

gen trs32_98 = .;
replace trs32_98 = 1 if (W3981==1 & G298==1);
replace trs32_98 = 0 if (W3981==1 & G298==0);

* Transfers of group 3 pupils to other group 3 schools during 1999;
gen trs33_99 = 0;
replace trs33_99 = 1 if
	((W3981==1 & W3991==1 & sch981~=sch991) |
	(W3981==1 & W3992==1 & sch981~=sch992) |
	(W3981==1 & W3993==1 & sch981~=sch993) |
	(W3981==1 & W3994==1 & sch981~=sch994) |
	(W3981==1 & W3995==1 & sch981~=sch995) |
	(W3981==1 & W3996==1 & sch981~=sch996) |
	(W3981==1 & W3997==1 & sch981~=sch997) |
	(W3981==1 & W3998==1 & sch981~=sch998) );
replace trs33_99 = . if (W3981==. | W3981==0 | G399==.);

gen trs31_99 = .;
replace trs31_99 = 1 if (W3981==1 & G199==1);
replace trs31_99 = 0 if (W3981==1 & G199==0);

gen trs32_99 = .;
replace trs32_99 = 1 if (W3981==1 & G299==1);
replace trs32_99 = 0 if (W3981==1 & G299==0);


** TABLE 4: TOTAL TRANSFERS;
collapse trs*;
foreach num in
	98 99
{;
gen tot1_`num' = trs11_`num' + trs21_`num' + trs31_`num';
gen tot2_`num' = trs12_`num' + trs22_`num' + trs32_`num';
gen tot3_`num' = trs13_`num' + trs23_`num' + trs33_`num';
};
summ trs* tot1* tot2* tot3*;
