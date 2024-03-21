********************************************************************************
***-------------------------HIV-risk behaviours study------------------------***
***--------------------------Hybrid model regression-------------------------***
********************************************************************************

********************************************************************************
* Which protective factors are linked with a reduction in HIV-risk behaviours
* among adolescents? 

********************************************************************************
* Within-person: examines variation between individuals
* Between-group: differences between groups (e.g. exposed vs. unexposed)
* Hybrid model: combines within- and between-person estimates

********************************************************************************

***-------------------Preparing variables for hybrid models------------------***

* Mean age of males at wave 1: 13.28397
* Mean age of females at wave 1: 14.35995

* Recode age variables
gen Age_i = Age_c
recode Age_i (0/19=0)(20/25=1)

bysort ID: egen Age_min = min(Age_c) 
gen c_Age_c = Age_min - 13.28397 if Girl_i == 0
replace c_Age_c = Age_min - 14.35995 if Girl_i == 1

gen time = Age_c - Age_min

* Recode outliers
recode GrTot_c (7/28=6)

* Generate within person centred variables

*m: between-person coefficients (average level of explanatory variable)
*d: within-person coefficients (deviation from one's average), eliminated time-invariant confounding

foreach var of varlist ScEnrl_i Food_i Food_c GrTot_c PosP_c MonP_c ComP_c /*Age_c Age_i Age_cat*/ { 
bysort ID: egen m123_`var'=mean(`var') 
gen d123_`var'=`var'-m123_`var' 
bysort ID: egen m23_`var'=mean(`var') if wave == 2 | wave == 3
gen d23_`var'=`var'-m23_`var' if wave == 2 | wave == 3
bysort ID: egen m12_`var'=mean(`var') if wave == 1 | wave == 2
gen d12_`var'=`var'-m12_`var' if wave == 1 | wave == 2
}

* Defining panel data (panelvar & timevar)
xtset ID wave

* Label outcome variables
label variable SxUnpr_i "Condomless sex"
label variable SxPrt2_i "Multiple sexual partners"
label variable SxDisp_i "Age-disparate sex"
label variable TSx6m_i "Transactional sex"
label variable SxSubs_i "Sex on substances{sup:a}"

label variable d123_ScEnrl_i "W_School_enrolment"
label variable m123_ScEnrl_i "B_School_enrolment"
label variable d123_Food_c "W_Food_security"
label variable m123_Food_c "B_Food_security"
label variable d123_GrTot_c "W_Number of grants"
label variable m123_GrTot_c "B_Number of grants"
label variable d123_PosP_c "W_Positive parenting"
label variable m123_PosP_c "B_Positive parenting"
label variable d123_MonP_c "W_Parental supervision"
label variable m123_MonP_c "B_Parental supervision"
label variable d123_ComP_c "W_Parental communication"
label variable m123_ComP_c "B_Parental communication"
label variable HIV_i "HIV status"

***---------------------Estimate ICC for study outcomes----------------------***

* After double checking, the Stata ICC command automatically calculates this
* value correctly for binary models: https://www.stata.com/manuals/meestaticc.pdf
* "The intraclass correlation for this model is where γ = σ21 for a mixed-effects
* linear regression, γ = 1 for a mixed-effects probit and ordered probit 
* regression, γ = π2/3 for a mixed-effects logistic and ordered logistic regression"

foreach var in SxPrt2_i TSx6m_i SxDisp_i SxUnpr_i SxSubs_i{
* Females
xtmelogit `var' || ID: if Girl == 1, or refineopts(iterate(0)) variance
estat icc
mat icc = r(icc2)
mat iccci =  r(ci2)
mat Gicc`var' = (icc, iccci)
* Males
xtmelogit `var' || ID: if Girl == 0, or refineopts(iterate(0)) variance
estat icc
mat icc = r(icc2)
mat iccci =  r(ci2)
mat Bicc`var' = (icc, iccci)
}
* Generate table of ICCs
* Females
matrix GICCAll = (GiccSxPrt2_i \ GiccTSx6m_i \ GiccSxDisp_i \ GiccSxUnpr_i \ GiccSxSubs_i)
matrix list GICCAll
* Males
matrix BICCAll = (BiccSxPrt2_i \ BiccTSx6m_i \ BiccSxDisp_i \ BiccSxUnpr_i \ BiccSxSubs_i)
matrix list BICCAll

***----------------------------Univariable regression------------------------***

foreach var in SxPrt2_i TSx6m_i SxDisp_i SxUnpr_i{
foreach pred in ScEnrl_i Food_c GrTot_c PosP_c MonP_c ComP_c{
* Females
logit `var' d123_`pred' m123_`pred' if Girl_i == 1, or
mat e = r(table)
mat b = e[1, 1..2]'
mat lb = e[5, 1..2]'
mat ub = e[6, 1..2]'
mat p = e[4, 1..2]'
mat G`var'`pred' = (b, lb, ub, p)
* Males
logit `var' d123_`pred' m123_`pred' if Girl_i == 0, or
mat e = r(table)
mat b = e[1, 1..2]'
mat lb = e[5, 1..2]'
mat ub = e[6, 1..2]'
mat p = e[4, 1..2]'
mat B`var'`pred' = (b, lb, ub, p)
}
}
foreach pred in ScEnrl_i Food_c GrTot_c PosP_c MonP_c ComP_c{
* Females
logit SxSubs_i d23_`pred' m23_`pred' if Girl_i == 1, or
mat e = r(table)
mat b = e[1, 1..2]'
mat lb = e[5, 1..2]'
mat ub = e[6, 1..2]'
mat p = e[4, 1..2]'
mat GSxSubs_i`pred' = (b, lb, ub, p)
* Males
logit SxSubs_i d23_`pred' m23_`pred' if Girl_i == 0, or
mat e = r(table)
mat b = e[1, 1..2]'
mat lb = e[5, 1..2]'
mat ub = e[6, 1..2]'
mat p = e[4, 1..2]'
mat BSxSubs_i`pred' = (b, lb, ub, p)
}
* Generate table of odds ratios for outcomes seperately
* Females
matrix GUNIUNPR = (GSxUnpr_iScEnrl_i\ GSxUnpr_iFood_c\ GSxUnpr_iGrTot_c\ GSxUnpr_iPosP_c\ GSxUnpr_iMonP_c\ GSxUnpr_iComP_c)
matrix GUNISXPRT2 = (GSxPrt2_iScEnrl_i\ GSxPrt2_iFood_c\ GSxPrt2_iGrTot_c\ GSxPrt2_iPosP_c\ GSxPrt2_iMonP_c\ GSxPrt2_iComP_c)
matrix GUNIDISP = (GSxDisp_iScEnrl_i\ GSxDisp_iFood_c\ GSxDisp_iGrTot_c\ GSxDisp_iPosP_c\ GSxDisp_iMonP_c\ GSxDisp_iComP_c)
matrix GUNITSX6M = (GTSx6m_iScEnrl_i\ GTSx6m_iFood_c\ GTSx6m_iGrTot_c\ GTSx6m_iPosP_c\ GTSx6m_iMonP_c\ GTSx6m_iComP_c)
matrix GUNISXSUBS = (GSxSubs_iScEnrl_i\ GSxSubs_iFood_c\ GSxSubs_iGrTot_c\ GSxSubs_iPosP_c\ GSxSubs_iMonP_c\ GSxSubs_iComP_c)
* Males
matrix BUNIUNPR = (BSxUnpr_iScEnrl_i\ BSxUnpr_iFood_c\ BSxUnpr_iGrTot_c\ BSxUnpr_iPosP_c\ BSxUnpr_iMonP_c\ GSxUnpr_iComP_c)
matrix BUNISXPRT2 = (BSxPrt2_iScEnrl_i\ BSxPrt2_iFood_c\ BSxPrt2_iGrTot_c\ BSxPrt2_iPosP_c\ BSxPrt2_iMonP_c\ BSxPrt2_iComP_c)
matrix BUNIDISP = (BSxDisp_iScEnrl_i\ BSxDisp_iFood_c\ BSxDisp_iGrTot_c\ BSxDisp_iPosP_c\ BSxDisp_iMonP_c\ BSxDisp_iComP_c)
matrix BUNITSX6M = (BTSx6m_iScEnrl_i\ BTSx6m_iFood_c\ BTSx6m_iGrTot_c\ BTSx6m_iPosP_c\ BTSx6m_iMonP_c\ BTSx6m_iComP_c)
matrix BUNISXSUBS = (BSxSubs_iScEnrl_i\ BSxSubs_iFood_c\ BSxSubs_iGrTot_c\ BSxSubs_iPosP_c\ BSxSubs_iMonP_c\ BSxSubs_iComP_c)
* Generate table of odds ratios for all outcomes together
* Females
matrix GUNIAll = (GUNISXPRT2, GUNITSX6M, GUNIDISP, GUNIUNPR, GUNISXSUBS)
matrix list GUNIAll
* Males
matrix BUNIAll = (BUNISXPRT2, BUNITSX6M, BUNIDISP, BUNIUNPR, BUNISXSUBS)
matrix list BUNIAll

***--------------------------Multivariable regression------------------------***

* Outcomes at all three timepoints
foreach var in SxPrt2_i TSx6m_i SxDisp_i SxUnpr_i{
* Females
xtmelogit `var' d123_GrTot_c m123_GrTot_c d123_PosP_c m123_PosP_c d123_MonP_c m123_MonP_c d123_ComP_c m123_ComP_c d123_ScEnrl_i m123_ScEnrl_i d123_Food_c m123_Food_c HIV_i Rural_i House_i HHSiz_c MOrph_i POrph_i Age_c wave if Girl_i == 1 || ID:, or variance
est store model
predict predicted_`var'1, mu
roctab `var' predicted_`var'1
mat G`var'roc = r(area)
est replay model, or
mat e = r(table)
mat b = e[1, 1...]'
mat lb = e[5, 1...]'
mat ub = e[6, 1...]'
mat p = e[4, 1...]'
mat G`var' = (b, lb, ub, p)
test d123_GrTot_c=m123_GrTot_c
mat t1 = (r(chi2), r(p))
test d123_PosP_c=m123_PosP_c
mat t2 = (r(chi2), r(p))
test d123_MonP_c=m123_MonP_c
mat t3 = (r(chi2), r(p))
test d123_ComP_c=m123_ComP_c
mat t4 = (r(chi2), r(p))
test d123_ScEnrl_i=m123_ScEnrl_i
mat t5 = (r(chi2), r(p))
test d123_Food_c=m123_Food_c
mat t6 = (r(chi2), r(p))
mat G`var'test = (t1\t2\t3\t4\t5\t6)
* Males
xtmelogit `var' d123_GrTot_c m123_GrTot_c d123_PosP_c m123_PosP_c d123_MonP_c m123_MonP_c d123_ComP_c m123_ComP_c d123_ScEnrl_i m123_ScEnrl_i d123_Food_c m123_Food_c HIV_i Rural_i House_i HHSiz_c MOrph_i POrph_i Age_c wave if Girl_i == 0 || ID:, or variance
est store model
predict predicted_`var'0, mu
roctab `var' predicted_`var'0
mat B`var'roc = r(area)
est replay model, or
mat e = r(table)
mat b = e[1, 1...]'
mat lb = e[5, 1...]'
mat ub = e[6, 1...]'
mat p = e[4, 1...]'
mat B`var' = (b, lb, ub, p)
test d123_GrTot_c=m123_GrTot_c
mat t1 = (r(chi2), r(p))
test d123_PosP_c=m123_PosP_c
mat t2 = (r(chi2), r(p))
test d123_MonP_c=m123_MonP_c
mat t3 = (r(chi2), r(p))
test d123_ComP_c=m123_ComP_c
mat t4 = (r(chi2), r(p))
test d123_ScEnrl_i=m123_ScEnrl_i
mat t5 = (r(chi2), r(p))
test d123_Food_c=m123_Food_c
mat t6 = (r(chi2), r(p))
mat B`var'test = (t1\t2\t3\t4\t5\t6)
}

* Sex after substance use at T2 and T3
* Females
xtmelogit SxSubs_i d23_GrTot_c m23_GrTot_c d23_PosP_c m23_PosP_c d23_MonP_c m23_MonP_c d23_ComP_c m23_ComP_c d123_ScEnrl_i m123_ScEnrl_i d123_Food_c m123_Food_c HIV_i Rural_i House_i HHSiz_c MOrph_i POrph_i Age_c wave if Girl_i == 1 || ID:, or variance
est store model
predict predicted_SxSubs_i1, mu
roctab SxSubs_i predicted_SxSubs_i1
mat GSxSubs_iroc = r(area)
est replay model, or
mat e = r(table)
mat b = e[1, 1...]'
mat lb = e[5, 1...]'
mat ub = e[6, 1...]'
mat p = e[4, 1...]'
mat GSxSubs_i = (b, lb, ub, p)
test d23_GrTot_c=m23_GrTot_c
mat t1 = (r(chi2), r(p))
test d23_PosP_c=m23_PosP_c
mat t2 = (r(chi2), r(p))
test d23_MonP_c=m23_MonP_c
mat t3 = (r(chi2), r(p))
test d23_ComP_c=m23_ComP_c
mat t4 = (r(chi2), r(p))
test d123_ScEnrl_i=m123_ScEnrl_i
mat t5 = (r(chi2), r(p))
test d123_Food_c=m123_Food_c
mat t6 = (r(chi2), r(p))
mat GSxSubs_itest = (t1\t2\t3\t4\t5\t6)
* Males
xtmelogit SxSubs_i d23_GrTot_c m23_GrTot_c d23_PosP_c m23_PosP_c d23_MonP_c m23_MonP_c d23_ComP_c m23_ComP_c d123_ScEnrl_i m123_ScEnrl_i d123_Food_c m123_Food_c HIV_i Rural_i House_i HHSiz_c MOrph_i POrph_i Age_c wave if Girl_i == 0 || ID:, or variance
est store model
predict predicted_SxSubs_i0, mu
roctab SxSubs_i predicted_SxSubs_i0
mat BSxSubs_iroc = r(area)
est replay model, or
mat e = r(table)
mat b = e[1, 1...]'
mat lb = e[5, 1...]'
mat ub = e[6, 1...]'
mat p = e[4, 1...]'
mat BSxSubs_i = (b, lb, ub, p)
test d23_GrTot_c=m23_GrTot_c
mat t1 = (r(chi2), r(p))
test d23_PosP_c=m23_PosP_c
mat t2 = (r(chi2), r(p))
test d23_MonP_c=m23_MonP_c
mat t3 = (r(chi2), r(p))
test d23_ComP_c=m23_ComP_c
mat t4 = (r(chi2), r(p))
test d123_ScEnrl_i=m123_ScEnrl_i
mat t5 = (r(chi2), r(p))
test d123_Food_c=m123_Food_c
mat t6 = (r(chi2), r(p))
mat BSxSubs_itest = (t1\t2\t3\t4\t5\t6)

* Generate table of odds ratios for all outcomes together
* Females
matrix GMULTI = (GSxPrt2_i, GTSx6m_i, GSxDisp_i, GSxUnpr_i, GSxSubs_i)
mat list GMULTI 
* Males
matrix BMULTI = (BSxPrt2_i, BTSx6m_i, BSxDisp_i, BSxUnpr_i, BSxSubs_i)
mat list BMULTI 

* Summarise tests of equality for deviation and mean coefficients
* Females
matrix GMULTI_HHTest = (GSxPrt2_itest, GTSx6m_itest, GSxDisp_itest, GSxUnpr_itest, GSxSubs_itest)
matrix list GMULTI_HHTest
* Males
matrix BMULTI_HHTest = (BSxPrt2_itest, BTSx6m_itest, BSxDisp_itest, BSxUnpr_itest, BSxSubs_itest)
matrix list BMULTI_HHTest

* Summarise ROC tests
* Females
matrix GMULTIAllRoc= (GSxPrt2_iroc, GTSx6m_iroc, GSxDisp_iroc, GSxUnpr_iroc, GSxSubs_iroc)
matrix list GMULTIAllRoc
* Males
matrix BMULTIAllRoc = (BSxPrt2_iroc, BTSx6m_iroc, BSxDisp_iroc, BSxUnpr_iroc, BSxSubs_iroc)
matrix list BMULTIAllRoc

/// Model 2

* Females
foreach var in SxPrt2_i TSx6m_i SxDisp_i SxUnpr_i SxSubs_i{
if "`var'" == "SxPrt2_i" {
xtmelogit `var' GrTot_c PosP_c d123_MonP_c m123_MonP_c ComP_c ScEnrl_i Food_c HIV_i Rural_i House_i HHSiz_c MOrph_i POrph_i Age_c wave if Girl_i == 1 || ID: , or variance
}
else{
if "`var'" == "SxUnpr_i" {	
xtmelogit `var' GrTot_c PosP_c d123_MonP_c m123_MonP_c ComP_c d123_ScEnrl_i m123_ScEnrl_i d123_Food_c m123_Food_c HIV_i Rural_i House_i HHSiz_c MOrph_i POrph_i Age_c wave if Girl_i == 1 || ID: , or variance
}
else{
if "`var'" == "SxSubs_i" {
xtmelogit `var' GrTot_c PosP_c d23_MonP_c m23_MonP_c ComP_c ScEnrl_i Food_c HIV_i Rural_i House_i HHSiz_c MOrph_i POrph_i Age_c wave if Girl_i == 1 || ID: , or variance
}
else{
if "`var'" == "TSx6m_i" | "`var'" == "SxDisp_i" {
xtmelogit `var' GrTot_c PosP_c MonP_c ComP_c ScEnrl_i Food_c HIV_i Rural_i House_i HHSiz_c MOrph_i POrph_i Age_c wave if Girl_i == 1 || ID: , or variance
}
}
}
}
est store GModel`var'
predict predicted_`var'11, mu
roctab `var' predicted_`var'11
mat G`var'roc1 = r(area)
est replay GModel`var', or
mat e = r(table)
mat e = r(table)
mat b = e[1, 1...]'
mat se = e[2, 1...]'
mat lb = e[5, 1...]'
mat ub = e[6, 1...]'
mat p = e[4, 1...]'
mat G`var' = (b, lb, ub, p)
}

* Males
foreach var in SxPrt2_i TSx6m_i SxDisp_i SxUnpr_i SxSubs_i{
if "`var'" == "SxPrt2_i" | "`var'" == "TSx6m_i" | "`var'" == "SxDisp_i" | "`var'" == "SxSubs_i"{
xtmelogit `var' GrTot_c PosP_c MonP_c ComP_c ScEnrl_i Food_c HIV_i Rural_i House_i HHSiz_c MOrph_i POrph_i Age_c wave if Girl_i == 0 || ID:, or variance
}
else{
if "`var'" == "SxUnpr_i" {	
xtmelogit `var' GrTot_c PosP_c MonP_c ComP_c d123_ScEnrl_i m123_ScEnrl_i Food_c HIV_i Rural_i House_i HHSiz_c MOrph_i POrph_i Age_c wave if Girl_i == 0 || ID: , or variance
}
}
est store Bmodel`var'
predict predicted_`var'00, mu
roctab `var' predicted_`var'00
mat B`var'roc1 = r(area)
mat B= r(area)
est replay Bmodel`var', or
mat e = r(table)
mat b = e[1, 1...]'
mat lb = e[5, 1...]'
mat ub = e[6, 1...]'
mat p = e[4, 1...]'
mat B`var' = (b, lb, ub, p)
}

* Generate table of odds ratios for all outcomes together
* Females
mat list GSxPrt2_i
mat list GTSx6m_i
mat list GSxDisp_i
mat list GSxUnpr_i
mat list GSxSubs_i
* Males
mat list BSxPrt2_i
mat list BTSx6m_i
mat list BSxDisp_i
mat list BSxUnpr_i
mat list BSxSubs_i

* Summarise ROC tests
* Females
matrix GMULTIAllRoc1= (GSxPrt2_iroc1, GTSx6m_iroc1, GSxDisp_iroc1, GSxUnpr_iroc1, GSxSubs_iroc1)
matrix list GMULTIAllRoc1
* Males
matrix BMULTIAllRoc1= (BSxPrt2_iroc1, BTSx6m_iroc1, BSxDisp_iroc1, BSxUnpr_iroc1, BSxSubs_iroc1)
matrix list BMULTIAllRoc1

***------------------------------------End-----------------------------------***