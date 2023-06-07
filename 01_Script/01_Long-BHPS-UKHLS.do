/*****************************************************************************************
* MERGING INDIVIDUAL FILES FROM HARMONISED BHPS AND UKHLS IN LONG FORMAT                 *
* To match individual level files from the harmonised BHPS and Understanding Society     *
* in long format, you need to remove the wave prefixes in the two sets of files and      *
* generate a wave identifier that works across both sets of files. The pidp will         *
* work as the unique cross-wave identifier across both sets of files. This code only     *
* keeps individuals who took part in BHPS and drops those who joined as part of          *
* Understanding Society.                                                                 *
*****************************************************************************************/

// change current file location
cd "C:\work\Forschung\Climate Change_Replication\02_Data"

// assign global macro to refer to Understanding Society data
global ukhls "C:\work\Forschung\Climate Change_Replication\00_Input_data\UKDA-6931-stata\stata\stata13_se"

// assign global macro to refer to geo identifiers
global geo "C:\work\Forschung\Climate Change_Replication\00_Input_data\UKDA-6670-stata\stata\stata13"

// assign global macros for the lists of waves
global BHPSwaves "a b c d e f g h i j k l m n o p q r"
global UKHLSwaves_bh "b c d e f g h i j" // since BHPS respondents did not take 
									 // part in Wave 1, begin at Wave 2
									 // - update this to include 
									 // new waves as they are released
global UKHLSwaves "a b c d e f g h i j" // use all UKHLS waves
global UKHLSno 10	// number of waves of UKHLS data	



*-------------------------------------------------------------
*------------------- Individual level ------------------------			
*-------------------------------------------------------------				 

// loop through the waves of bhps
foreach w of global BHPSwaves {

	// find the wave number
	local waveno=strpos("abcdefghijklmnopqrstuvwxyz","`w'")

	// open the individual file for that wave
	// use pidp b`w'_age_dv b`w'_paygu_dv using "$ukhls\bhps_w`waveno'\b`w'_indresp_protext", clear
	use "$ukhls\bhps_w`waveno'\b`w'_indresp_protect", clear
	
	// remove the wave prefix
	rename b`w'_* *
	
	// keep relevant variables
	
	* identifiers
	global vars_id pidp hidp pno memorig 
	
	* meta data 
	global var_met buno buno_dv sampst ivfio istrtdatd istrtdatm istrtdaty intdatd_dv intdatm_dv intdaty_dv hhmem indin01_lw lrwtuk1 indinus_xw indinub_xw xrwtuk1
	
	* locality
	global vars_loc gor_dv mvever mvmnth mvyr plnew movest origadd adcts addrmov_dv movdir lkmove xpmove xpmvmnth xpmvyr plnowm plnowy4 urban_dv
	
	* demographics
	// global vars_dem age_dv ukborn bornuk_dv plbornc_cc sex sex_dv mastat marstat jbstat race racel_dv nchild_dv ethn_dv
	global vars_dem age age_dv sex plbornc yr2uk4 mastat marstat_dv race racel_dv depchl_dv fnspno mnspno mnspid fnspid
	
	* socec
	// global vars_socec hiqual_dv jbsoc90_cc jbsoc00_cc jbnssec_dv jbnssec8_dv jbnssec5_dv jbnssec3_dv fimnnet_dv fimnlabnet_dv a_fimnlabgrs_dv paygu_dv
	global vars_socec jbstat jbsemp jbhrs jbot jbttwt basnsa basrate basrest paynu_dv paygu_dv fimnlabgrs_dv fimngrs_dv jbnssec_dv jbnssec8_dv jbnssec5_dv finnow save fiyrinvinc_dv fiyrdia
	
	* educ
	global vars_educ edasp edtype feend hiqual_dv nhiqual_dv qfhigh_dv scend school
	
	* health and satisfaction
	//global vars_he sf12mcs_dv sf12pcs_dv health scghq1_dv scghq2_dv
	global vars_he sf1 scsf1 sf12mcs_dv sf12pcs_dv health scghq1_dv scghq2_dv sclfsat1 sclfsat2 sclfsato lfsat1 lfsat2 lfsato scwemwb scwemwbb scghqa scghqb scghqc scghqd scghqe scghqf scghqg scghqh scghqi scghqj scghqk scghql slp_qual hrs_slph hrs_slpm med_slp scslp_qual scmed_slp scsf2a scsf2b scsf3a scsf3b scsf4a scsf4b scsf5 scsf6a scsf6b scsf6c scsf7
	
	// variables efficacy and trust
	global vars_eff scwemwb scwemwbb scghqb scghqf scghqd scghqh se1 se2 se3 se4 se5 se6 se7 se8 se9 se10 riskb scriska scriskb sctrust ivcoop iv4 volun volfreq chargv charfreq
	
	* individual and family background
	// global vars_bkg scend_dv j1soc00_cc maid macob maedqf masoc90_cc masoc00_cc masoc10_cc paid pacob paedqf pasoc90_cc pasoc00_cc pasoc10_cc
	global vars_bkg manssec_dv panssec_dv manssec8_dv panssec8_dv paju maju pacob macob maid macob maedqf masoc90_cc masoc00_cc masoc10_cc paid pacob paedqf pasoc90_cc pasoc00_cc pasoc10_cc
	
	* household
	global vars_hh fihhmn nchild_dv husits huboss hoh howlng
	
	* attitudes
	global vars_att scopfama scopfamb scopfamd scopfamf scopfamh oprlg1 oprlg2 oprlg3 vote1 vote2 vote3 vote4 vote5 vote6 vote7 vote8 vote3_all poleff1 poleff2 poleff3 poleff4 oppola oppolb oppolc oppold
	
	* environment variables
	global vars_env1 scenv_bccc scenv_bcon scenv_brit scenv_canc scenv_ccls scenv_cfit scenv_chwo scenv_crex scenv_crlf scenv_dstr scenv_exag scenv_fitl scenv_ftst scenv_futr scenv_grn scenv_meds scenv_noot scenv_nowo scenv_pmep scenv_pmre scenv_tlat orgm3 orga3
	global vars_env2 openv1 openv2 openv3 openv4 openva openvb openvc
	
	* Climate change variables
	global vars_clim opcca opccb opccc opccd opcce opccf scopecl200 scopecl30
	
	* environmental behaviour variables
	global vars_envb1 envhabit1 envhabit2 envhabit3 envhabit4 envhabit5 envhabit6 envhabit7 envhabit8 envhabit9 envhabit10 envhabit11
	global vars_envb2 grnlfa grnlfb grnlfc grnlfd grnlfe grnlff grnlfg grnlfh trcarfq trbikefq caruse
	
	* Neighbourhood
	global vars_nb simarea nbrcoh1 nbrcoh2 nbrcoh3 nbrcoh4 nbrcoh_dv nbrcohdk_dv nbrsnci_dv crdark crburg crcar crdrnk crgraf crmugg crrace crteen crvand lknbrd llknbrd locchd locsera locserap locseras locserb locserc locserd locsere opngbha opngbhb opngbhc opngbhd opngbhe opngbhf opngbhg opngbhh scopngbha scopngbhb scopngbhc scopngbhd scopngbhe scopngbhf scopngbhg scopngbhh
	
	* Accomodation
	global vars_acc lfsat3 lkmovy netuse netpuse
	
	* Commuting
	global vars_com jsttwtb jbttwt jbpl twkdiff1 twkdiff2 twkdiff3 twkdiff4 twkdiff5 twkdiff6 twkdiff7 twkdiff8 twkdiff9
	
	* Distances
	global vars_dis distmov distmov_dv jsworkdis workdis mafar pafar chfar mlivedistf mlivedist
	
	* Language
	global vars_lan iv6d englang natidb engspk engtel engform readdif formdif teldif spkdif
	
	* keep variables (continue if not available)
	// keep $vars_id $vars_dem $vars_loc vars_educ $vars_socec $vars_he $vars_bkg $vars_hh $vars_att $var_met
	global varlist1 $vars_id
	global varlist $vars_dem $vars_loc $vars_educ $vars_socec $vars_he $vars_eff $vars_bkg $vars_hh $vars_att $var_met $vars_env1 $vars_env2 $vars_clim $vars_envb1 $vars_envb2 $vars_nb $vars_acc $vars_com $vars_dis $vars_lan
	foreach v of global varlist {
		capture confirm var `v'							//Var exists?
		if !_rc {
			global varlist1 ${varlist1} `v'	//only existing vars in varlist1
		}
	}		
	keep $varlist1
	
	// generate a variable which records the wave number
	gen wave=`waveno'
	
	// gen year
	gen year=1990+`waveno'
	
		
	**** Add dependet child variable from indall data for bhps
	merge 1:1 pidp using "$ukhls\bhps_w`waveno'\b`w'_indall_protect", keepusing(b`w'_depchl_dv)
	drop if _merge == 2
	drop _merge
	rename b`w'_depchl_dv depchl_dv
	
	
	// save the file for future use
	save tmp_b`w'_indresp, replace
}

// loop through the relevant waves of Understanding Society
foreach w of global UKHLSwaves {

	// find the wave number
	local waveno=strpos("abcdefghijklmnopqrstuvwxyz","`w'")
	
	// open the individual level file for that wave
	// use pidp pid `w'_age_dv `w'_paygu_dv using "$ukhls/ukhls_w`waveno'/`w'_indresp_protect", clear
	use "$ukhls/ukhls_w`waveno'/`w'_indresp_protect", clear
	
	/*
	// keep the individual if they have a pid - ie were part of BHPS
	// individuals have pid==-8 (inapplicable) if they were not part of BHPS
	keep if pid>0
	
	// drop the pid variable
	drop pid
	*/
	
	// remove the wave prefix
	rename `w'_* *
	
	// keep relevant variables
	* keep variables (continue if not available)
	// keep $vars_id $vars_loc $vars_dem $vars_socec $vars_he $vars_bkg $vars_hh
	global varlist1 $vars_id
	global varlist $vars_dem $vars_loc $vars_educ $vars_socec $vars_he $vars_eff $vars_bkg $vars_hh $vars_att $var_met $vars_env1 $vars_env2 $vars_clim $vars_envb1 $vars_envb2 $vars_nb $vars_acc $vars_com $vars_dis $vars_lan
	foreach v of global varlist {
		capture confirm var `v'							//Var exists?
		if !_rc {
			global varlist1 ${varlist1} `v'	//only existing vars in varlist1
		}
	}		
	keep $varlist1

	// generate a variable which records the wave number + 17 
	// - treating wave 2 ukhls as wave 19 of bhps --> TR: changed to 18!
	gen wave=`waveno'+18
	
	// gen year
	gen year=1990+`waveno'+18
	
	// save the file for future use
	save tmp_`w'_indresp, replace
}

// loop through the waves of bhps
foreach w of global BHPSwaves {
	
	// first time through the loop
	if "`w'"=="a" {
	
		// reopen the first file created
		use tmp_ba_indresp, clear
		
	// following times through the loop	
	} 
	else {	
		
		// append each file in turn
		append using tmp_b`w'_indresp
	}
}

// loop through the waves of ukhls from Wave 1
foreach w of global UKHLSwaves {
	
	// append each file in turn
	append using tmp_`w'_indresp
}

// create labels for the wave variable
// loop through the waves of bhps
foreach n of numlist 1/18 {

	// add a label for each wave number in turn
	lab def wave `n' "BHPS Wave `n'", modify
}

// loop through the waves of ukhls 
// (using the global macro UKHLSno to define the last wave)
foreach n of numlist 1/$UKHLSno {
	
	// calculate which label value this label will apply to
	local waveref=`n'+18
	
	// add a label for each wave in turn
	lab def wave `waveref' "UKHLS Wave `n'", modify
}

// apply the label to the wave variable
lab val wave wave

// check how many observations are available from each wave
tab wave

// order
order $vars_id $varlist
order pidp hidp pno year wave memorig

// Sort by id year
sort pidp year



// Use xwavedat to update person-constant variables
global var_con birthm birthy ukborn bornuk_dv plbornc scend_dv feend_dv school_dv racel_dv ethn_dv generation yr2uk4 evermar_dv anychild_dv paju maju pacob macob maid macob maedqf masoc90_cc masoc00_cc masoc10_cc paid pacob paedqf pasoc90_cc pasoc00_cc pasoc10_cc psnenub_xd

/* // not neccessary with update replace option
foreach v of global var_con {
	capture confirm var `v'							//Var exists?
	if !_rc {
		drop `v'	//drop current var
	}
}	
*/

merge m:1 pidp using "$ukhls\ukhls_wx\xwavedat_protect", keepusing(pidp $var_con)  update replace
drop if _merge<=2
drop _merge


// order
order $vars_id $varlist
order pidp hidp pno year wave memorig gor_dv urban_dv sex age age_dv birthm birthy bornuk_dv plbornc yr2uk4 racel_dv ethn_dv generation mastat marstat_dv evermar_dv anychild_dv hiqual_dv nhiqual_dv qfhigh_dv scend_dv feend_dv school_dv plnew origadd movdir mvmnth mvyr plnowm plnowy4 lkmove xpmove xpmvmnth xpmvyr 

// Sort by id year
sort pidp year

// Code negatives as missing values
mvdecode _all, mv(-21/-10=.a\-9=.\-8/-7=.b\-2/-1=.)



// Harmonise identical variables from BHPS and UKHLS

* Cross-secitonal Weights
gen ukpopweight = indinub_xw
replace ukpopweight = xrwtuk1 if wave <= 18
replace ukpopweight = indinus_xw if ukpopweight ==. & wave == 19


* Longitunial Weights
gen longweight01 = indin01_lw
replace longweight01 = lrwtuk1 if wave <= 18

drop indinub_xw xrwtuk1 indin01_lw lrwtuk1

* Marital status
recode mastat (7 = 1) (8 = 4) (9 = 5) (10 = 3)
replace marstat_dv = mastat if wave <= 18
drop mastat mastat_dv

/*
* Address change
replace origadd = plnew if wave <= 18
drop plnew
*/

* Climate change effect
replace scopecl200 = opccf if wave <= 18
drop opccf
replace scopecl30 = opcce if wave <= 18
drop opcce

* Environmental behaviour
replace envhabit1 = grnlfa if wave <= 18
drop grnlfa
replace envhabit2 = grnlfb if wave <= 18
drop grnlfb
replace envhabit3 = grnlfc if wave <= 18
drop grnlfc
replace envhabit4 = grnlfd if wave <= 18
drop grnlfd
replace envhabit5 = grnlfe if wave <= 18
drop grnlfe
* grnlff missing in UKHLS
replace envhabit6 = grnlfg if wave <= 18
drop grnlfg
replace envhabit7 = grnlfh if wave <= 18
drop grnlfh

* label hiqual
label val hiqual_dv b_nhiqual_dv


// save the file containing all waves
save all_indresp, replace

// erase each temporary file using loops
foreach w of global BHPSwaves {
	erase tmp_b`w'_indresp.dta
}
foreach w of global UKHLSwaves {
	erase tmp_`w'_indresp.dta
}



*------------------------------------------------------------
*------------------- Household level ------------------------
*------------------------------------------------------------

// loop through the waves of bhps
foreach w of global BHPSwaves {

	// find the wave number
	local waveno=strpos("abcdefghijklmnopqrstuvwxyz","`w'")

	// open the individual file for that wave
	use "$ukhls\bhps_w`waveno'\b`w'_hhresp_protect", clear
	
	// remove the wave prefix
	rename b`w'_* *
	capture confirm var origadd						//Var exists?
	if !_rc {
		rename origadd hhorigadd	//rename (same var in indresp
	}
	
	// keep relevant variables
	
	* identifiers
	global vars_hid hidp
	
	* demographics 
	global var_hdem hhsize hhtype hhtype_dv agechy_dv nkids_dv nch02_dv nch34_dv nch511_dv nch1215_dv nemp_dv ncouple_dv fihhmngrs_dv carown carval region
	
	* locaility 
	global var_hloc hhmove hsivlw hhorigadd crburg crcar crdrnk crgraf crmugg crrace crrubsh crteen crvand

	* accomodation
	global var_haccom tenure_dv hsownd hsownd_bh hsval hscost mglife mgold mgnew hsroom hsprbg hsprbh hsprbi hsprbj hsprbp hsprbq rent rent_dv rentgrs_dv rent1 rent2 rent3 rent4 rent5 rent6 rent7 xphsdb
	global var_haccom2 mgynot mgynot_bh mgextra houscost1_dv houscost2_dv 
	
	* environmental
	global var_henv grimyn noisyn
	
	* Electricity and gas
	global var_elic fuelhave1 fuelhave2 fuelhave3 xpgasy gaspay xpelecy elecpay xpduely fuelduel duelpay heatch heatyp hheat xpoily
	
	* reference person
	global var_head hrpid hrpno
	
	* keep variables (continue if not available)
	global varlist1 $vars_hid
	global varlist $var_hdem $var_hloc $var_haccom $var_haccom2 $var_henv $var_elic $var_head 
	foreach v of global varlist {
		capture confirm var `v'							//Var exists?
		if !_rc {
			global varlist1 ${varlist1} `v'	//only existing vars in varlist1
		}
	}		
	keep $varlist1
	
	// generate a variable which records the wave number
	gen wave=`waveno'
	
	// gen year
	gen year=1990+`waveno'	
	
	// save the file for future use
	save tmp_b`w'_hhresp, replace
}

// loop through the relevant waves of Understanding Society
foreach w of global UKHLSwaves {

	// find the wave number
	local waveno=strpos("abcdefghijklmnopqrstuvwxyz","`w'")
	
	// open the individual level file for that wave
	use "$ukhls/ukhls_w`waveno'/`w'_hhresp_protect", clear
	
	// remove the wave prefix
	rename `w'_* *
	capture confirm var origadd						//Var exists?
	if !_rc {
		rename origadd hhorigadd	//rename (same var in indresp
	}
	
	
	// keep relevant variables
	* keep variables (continue if not available)
	global varlist1 $vars_hid
	global varlist $var_hdem $var_hloc $var_haccom $var_haccom2 $var_henv $var_elic $var_head 
	foreach v of global varlist {
		capture confirm var `v'							//Var exists?
		if !_rc {
			global varlist1 ${varlist1} `v'	//only existing vars in varlist1
		}
	}		
	keep $varlist1

	// generate a variable which records the wave number + 17 
	// - treating wave 2 ukhls as wave 19 of bhps --> TR: changed to 18!
	gen wave=`waveno'+18
	
	// gen year
	gen year=1990+`waveno'+18
	
	// save the file for future use
	save tmp_`w'_hhresp, replace
}

// loop through the waves of bhps
foreach w of global BHPSwaves {
	
	// first time through the loop
	if "`w'"=="a" {
	
		// reopen the first file created
		use tmp_ba_hhresp, clear
		
	// following times through the loop	
	} 
	else {	
		
		// append each file in turn
		append using tmp_b`w'_hhresp
	}
}

// loop through the waves of ukhls from Wave 1
foreach w of global UKHLSwaves {
	
	// append each file in turn
	append using tmp_`w'_hhresp
}

// create labels for the wave variable
// loop through the waves of bhps
foreach n of numlist 1/18 {

	// add a label for each wave number in turn
	lab def wave `n' "BHPS Wave `n'", modify
}

// loop through the waves of ukhls 
// (using the global macro UKHLSno to define the last wave)
foreach n of numlist 1/$UKHLSno {
	
	// calculate which label value this label will apply to
	local waveref=`n'+18
	
	// add a label for each wave in turn
	lab def wave `waveref' "UKHLS Wave `n'", modify
}

// apply the label to the wave variable
lab val wave wave

// check how many observations are available from each wave
tab wave

// order
order $vars_hid $varlist
order hidp year wave

// Sort by id year
sort hidp year


// Code negatives as missing values
mvdecode _all, mv(-21/-10=.a\-9=.\-8/-7=.b\-2/-1=.)


// save the file containing all waves
save all_hhresp, replace

// erase each temporary file using loops
foreach w of global BHPSwaves {
	erase tmp_b`w'_hhresp.dta
}
foreach w of global UKHLSwaves {
	erase tmp_`w'_hhresp.dta
}






*--------------------------------------------------------------------
*------------------- Prepare geo identifiers ------------------------
*--------------------------------------------------------------------

// loop through the waves of bhps
foreach w of global BHPSwaves {

	// find the wave number
	local waveno=strpos("abcdefghijklmnopqrstuvwxyz","`w'")

	// open the individual file for that wave
	use "$geo\bhps\b`w'_lsoa01_protect", clear
	
	// remove the wave prefix
	rename b`w'_* *
	
	// Drop hip
	keep hidp lsoa01
	
	// generate a variable which records the wave number
	gen wave=`waveno'
	
	// gen year
	gen year=1990+`waveno'
	
	// save the file for future use
	save tmp_b`w'_lsoa, replace
}

// loop through the relevant waves of Understanding Society
foreach w of global UKHLSwaves {

	// find the wave number
	local waveno=strpos("abcdefghijklmnopqrstuvwxyz","`w'")
	
	// open the individual level file for that wave
	use "$geo/ukhls/`w'_lsoa01_protect", clear
	
	// remove the wave prefix
	rename `w'_* *
	
	// Drop hip
	keep hidp lsoa01

	// generate a variable which records the wave number + 17 
	// - treating wave 2 ukhls as wave 19 of bhps --> TR: changed to 18!
	gen wave=`waveno'+18
	
	// gen year
	gen year=1990+`waveno'+18
	
	// save the file for future use
	save tmp_`w'_lsoa, replace
}

// loop through the waves of bhps
foreach w of global BHPSwaves {
	
	// first time through the loop
	if "`w'"=="a" {
	
		// reopen the first file created
		use tmp_ba_lsoa, clear
		
	// following times through the loop	
	} 
	else {	
		
		// append each file in turn
		append using tmp_b`w'_lsoa
	}
}

// loop through the waves of ukhls from Wave 1
foreach w of global UKHLSwaves {
	
	// append each file in turn
	append using tmp_`w'_lsoa
}

// create labels for the wave variable
// loop through the waves of bhps
foreach n of numlist 1/18 {

	// add a label for each wave number in turn
	lab def wave `n' "BHPS Wave `n'", modify
}

// loop through the waves of ukhls 
// (using the global macro UKHLSno to define the last wave)
foreach n of numlist 1/$UKHLSno {
	
	// calculate which label value this label will apply to
	local waveref=`n'+18
	
	// add a label for each wave in turn
	lab def wave `waveref' "UKHLS Wave `n'", modify
}

// apply the label to the wave variable
lab val wave wave

// check how many observations are available from each wave
tab wave

// Sort by id year
sort hidp year


// save the file containing all waves
save all_lsoa, replace

// erase each temporary file using loops
foreach w of global BHPSwaves {
	erase tmp_b`w'_lsoa.dta
}
foreach w of global UKHLSwaves {
	erase tmp_`w'_lsoa.dta
}



*-------------------------------------------------------
*------------------- Merge data ------------------------
*-------------------------------------------------------

// Load ind
use all_indresp, clear

// merge hh data
merge m:1 hidp year wave using all_hhresp
drop if _merge==2
drop _merge

// merge hh data
merge m:1 hidp year wave using all_lsoa
drop if _merge==2
drop _merge

/*
// Add lastnm mother and father identifier for bhps from youth questionnaire
merge m:1 pidp using all_youth_lastnm
drop if _merge==2
drop _merge

* replace with NA for UKHLS obs
replace fnspid_bh = . if wave > 18
replace mnspid_bh = . if wave > 18
*/

// Add egoalt parents pnos for bhps
merge 1:1 pidp wave using bhps_egoalt_childpno
drop if _merge==2
drop _merge

* Replace NA with zero (not in hh, similar to UKHLS coding)
replace mnspno_bh = 0 if mnspno_bh == . & wave <= 18 
replace fnspno_bh = 0 if fnspno_bh == . & wave <= 18 


// Sort by id year
sort hidp year


// Harmonise hh head identifier
gen head = .
replace head = 0 if hrpid != pidp & hrpid != .
replace head = 1 if hrpid == pidp

replace head = 2 - hoh if wave <= 18

drop hrpid hoh

// Harmonise benefit unit number
replace buno_dv = buno if wave <= 18

drop buno

// Order
order pidp-wave head

// Sort
sort pidp wave 

// save the file containing all waves
save all_bhpsukhls_stata, replace 

// save the file containing all waves
saveold all_bhpsukhls, replace nolabel version(12) 

