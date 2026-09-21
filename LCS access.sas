/**********************************************************************************************************************
SAS program for geographic access to lung cancer screening and early-stage lung cancer diagnosis
Data sources:  SEER database submitted in November 2021 (follow-up through December 31, 2020) and linked Medicaid Analytic eXtracts Personal Summary files
Exposure: Geographic access to lung cancer screening centers
 - Lung Cancer Screening Centers of Excellence (N=1947) was identified using data from the GO2 Foundation for Lung Cancer and the American College of 
Radiology (2018). 
 - Geographic access to LCS centers was quantified using a refined two-step floating catchment area method (Shrestha P et al. Geographic Access to 
Cancer Care and Treatment and Outcomes of Early-Stage Non–Small Cell Lung Cancer. JAMA Network Open 2025;8(3):e251061-e251061. Lian M et al. Geographic Access to 
Cancer Care and Breast Cancer Treatment in Low-Income Women. Med Care 2025;63(9):694-702. 
 - Catchment areas were defined using Euclidean distance thresholds of 10, 30, 60 miles for metropolitan, urban, and rural areas, respectively.
 - Census tract populations from the 2010 US Census were used to derive population-weighted centroids.
 - Tract-level access scores were aggregated to the county level.
 - Geographic access scores were categorized into quartiles based on their distribution in the study population, with Quartile 1 representing the 
lowest level of geographic access.
Outcome variable: Early-stage lung cancer diagnosis
 - Localized disease identified using the SEER Combined Summary Stage (2004+).
Covariates: age, sex, race/ethnicity, health insurance, non-metropolitan residence, county-level socioeconomic deprivation, county-level vehicle 
ownership, year of diagnosis, and histologic subtype.
 - Health insurance at diagnosis was determined using SEER primary payer information.
 - Medicaid enrollment before or at diagnosis was ascertained from the Medicaid Analytic eXtracts Personal Summary files. Individuals with Medicaid 
coverage limited to breast/cervical cancer treatment, pregnancy-related services, or family planning services were classified as not having Medicaid coverage.
 - Non-metropolitan residence was defined using the Rural-Urban Continuum Codes.
 - County-level socioeconomic deprivation was measured using a composite deprivation index derived from 21 variables from the 2008-2012 American 
Community Surveys (Lian M. et al. Statistical Assessment of Neighborhood Socioeconomic Deprivation Environment in Spatial Epidemiologic Studies. 
Open J Stat 2016;6(3):436-442.) Counties were categorized into quintiles based on the national distribution of the index, with higher quintiles 
indicating greater socioeconomic deprivation.
 - County-level vehicle ownership was obtained from the 2015-2019 American Community Surveys and categorized as high or low based on the median 
percentage of households without a vehicle.
 - County-level public transit availability (2016-2018) was obtained from the National Neighborhood Data Archive. Metropolitan counties were 
classified as having limited transportation resources if they had both low vehicle ownership and low public transit availability, defined as below 
the median number of public transit stops per square mile. 
Statistical models:
 - Modified Poisson regression with robust variance estimation and county-level clustering was used to estimate multivariable-adjusted prevalence 
ratios.
 - Linear trends were assessed by modeling the continuous geographic access score. 
***********************************************************************************************************************/

/******read in SEER data (seer_in) using the algorithm provided by SEER******/
data seer1;set seer_in;
 if year_of_diagnosis in ('2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019');
 if sequence_number in ('00','01');
 if primary_site in ('C340','C341','C342','C343','C348','C349');
run;
data seer2;set seer1;
age=Agerecodewithsingleages_and_100*1;
/*race*/
race1=RACE_ETHNICITY;
race2=RACE_RECODE_W_B_AI_API;
race3=NHIA_DERIVED_HISP_ORIGIN;
race4=RaceandoriginrecodeNHWNHBNHAIAN;
/*SEER registries*/
reg=seer_registry*1;   
reg1=SEERregistrywithCAandGAaswholes;
/*time of diagnosis: dxm (month), dxy (year)*/
dxm=MONTH_OF_DIAGNOSIS_RECODE;
dxy=YEAR_OF_DIAGNOSIS*1;
dxdate1=(dxy-2006)*12;
if dxm=. then dxm=6;
dxdate=dxdate1+dxm;
dxdate2=dxy*100+dxm;
/*Histology*/
histo3=HISTOLOGIC_TYPE_ICD_O_3;  
/*grade: grade1 (through 2017), grade2 (clinical 2018+), grade3 (pathological 2018+)*/
grade1=grade_thru_2017;   
grade2=grade_clinical_2018;
grade3=grade_pathological_2018;  
/*reporting source*/
report=TYPE_OF_REPORTING_SOURCE;  
/*tumor size*/
size1=CSTUMOR_SIZE_EXT_EVAL_2004_2015;
size2=CS_TUMOR_SIZE_2004_2015;
size3=tumor_size_summary_2016;
/*lymph nodes*/
node1=cs_lymph_nodes_2004_2015;
nodeex=REGIONAL_NODES_EXAMINED_1988;
/*cancer stage: stage1(2004-2015, 6th ed),stage2(2010-2015, 7th ed),stage3(2004-2017),stage4(1973-2015),stage5(1998-2017),stage6(2004+),stage7(2016-2017,7th ed),stage8(2018+)*/
stage1=DERIVEDAJCCSTAGEGROUP6THED20042*1;
stage2=DERIVEDAJCCSTAGEGROUP7THED20102*1;
stage3=SEERCOMBINEDSUMMARYSTAGE2000200*1;
stage4=SEER_HISTORIC_STAGE_A_1973_2015*1;
stage5=SUMMARY_STAGE_2000_1998_2017*1;
stage6=combined_summary_stage_2004;
stage7=derivedseercmbstg_grp_2016_2017;
stage8=derivedEOD2018_stage_group_2018;
tnmt1=Derived_AJCC_T_6th_ed_2004_2015;
tnmn1=Derived_AJCC_N_6th_ed_2004_2015;
tnmm1=Derived_AJCC_M_6th_ed_2004_2015;
tnmt2=Derived_AJCC_T_7th_ed_2010_2015;
tnmn2=Derived_AJCC_N_7th_ed_2010_2015;
tnmm2=Derived_AJCC_M_7th_ed_2010_2015;
tnmt7=DerivedSEERCombined_T_2016_2017;
tnmn7=DerivedSEERCombined_N_2016_2017;
tnmm7=DerivedSEERCombined_M_2016_2017;
tnmt8=Derived_EOD_2018_T_2018;
tnmn8=Derived_EOD_2018_N_2018;
tnmm8=Derived_EOD_2018_M_2018;
tnmother1=SeparateTumorNodulesIpsilateral;  /*Separate Tumor Nodules Ipsilateral Lung Recode (2010+)*/
tnmother2=Pleural_Effusion_Recode_2010;  /*Pleural Effusion Recode (2010+)*/
tnmother3=CS_tumor_size_2004_2015*1; /*CS Tumor size (2004-2015)*/
tnmother4=CSsitespecificfactor120042017va; /*CS Site-specific factor 1(2004-2017 varying by schema)*/
tnmother5=CSsitespecificfactor220042017va; /*CS Site-specific factor 2(2004-2017 varying by schema)*/ 
tnmother6=Tumor_Size_Summary_2016; /*Tumor Size Summary (2016+)*/
/*treatment*/
treatm=MONTH_THERAPY_STARTED;   /*month of first therapy*/
treaty=YEAR_THERAPY_STARTED;
treatf=DATE_THERAPY_STARTED_FLAG;
txdate1=(treaty-2006)*12;treatm=treatm*1;
txdate=txdate1+treatm;
treato=OTHER_CANCER_DIRECTED_THERAPY;  /*other cancer-directed therapy*/
chemo=CHEMOTHERAPY_RECODE_YES_NO_UNK*1;  /*chemo*/
rad1=radiation; /*radiation therapy*/
surg1=RX_SUMM_SURG_PRIM_SITE_1998;
treatseq1=RX_SUMM_SURG_RAD_SEQ*1;
treatseq2=RX_SUMM_SYSTEMIC_SURG_SEQ*1;
behavior=BEHAVIOR_RECODE_FOR_ANALYSIS;  /*behavior code*/
/*county FIPs;early FIPs of Hawaii were changed*/
cty=state*1000+county;
if CTY=15911 then CTY=15001;
if CTY=15912 then CTY=15003;
if CTY=15913 then CTY=15005;
if CTY=15914 then CTY=15007;
if CTY=15915 then CTY=15009;
if CTY=15900 then CTY=15009;
cty1=county_at_dx_geocode_1990;
cty2=county_at_dx_geocode_2000;
cty3=county_at_dx_geocode_2010;
/*outcome: event1(cause specific death),event2(other cause death),event3(vital status),time1(survival months)*/
event1=SEERCAUSESPECIFICDEATHCLASSIFIC;  /*cause specific death*/
event2=SEEROTHERCAUSEOFDEATHCLASSIFICA;  /*other cause death*/
deathm=SEER_DATEOFDEATH_MONTH;  /*death month*/ 
deathy=SEER_DATEOFDEATH_YEAR;  /*death year*/
event3=VITALSTATUSRECODESTUDYCUTOFFUSE;  /*vital status*/
time1=survival_months*1;
followm=MONTH_OF_LAST_FOLLOW_UP_RECODE;
followy=YEAR_OF_LAST_FOLLOW_UP_RECODE;
/*rural*/
rural03=RURAL_URBAN_CONTINUUM_CODE_2003;
rural13=RURAL_URBAN_CONTINUUM_CODE_2013;
/*insurance*/
insurance1=INSURANCE_RECODE_2007;
insurance2=PRIMARY_PAYER_AT_DX*1;
/*laterality*/
later=laterality;

aa=1;
keep patient_id reg reg1 age mari race1-race4 dxm dxy dxdate histo3 grade1-grade3 report size1-size3 node1 stage1-stage8 treatm treaty treatf txdate 
     treato chemo rad1 surg1 treatseq1 treatseq2 node1 nodeex behavior cty cty1-cty3 event1-event3 deathm deathy time1 followm followy
	 rural03 rural13 insurance1 insurance2 later aa dxdate2 sex Medicaid_Flag tnmt1 tnmt2 tnmt7 tnmt8 tnmn1 tnmn2 tnmn7 tnmn8 tnmm1 tnmm2 tnmm7 tnmm8 
     tnmother1-tnmother6
	 ;
run;proc sort;by patient_id;run;

/******read in Medicaid enrollment files (taf_base2014-taf_base2020, maxdata2005-maxdata2015) using the algorithm provided by SEER and 
identify eligible Medicaid enrollment******/
data seer_id;set seer2;keep patient_id aa dxy dxm;run;proc sort;by patient_id;run;

data taf1;set taf_base2014-taf_base2020;
year=RFRNC_YR*1;state=STATE_CD;
array a1{12} MDCD_ENRLMT_DAYS_01-MDCD_ENRLMT_DAYS_12;  /*Medicaid enrollment days*/
array a2{12} enr1-enr12;
array a3{12} DUAL_ELGBL_CD_01-DUAL_ELGBL_CD_12;        /*Dual eligible code*/
array a4{12} MC_PLAN_TYPE_CD_01-MC_PLAN_TYPE_CD_12;    /*managed care plan type*/
array a5{12} dual1-dual12;
array a6{12} dual21-dual32;
array a7{12} RSTRCTD_BNFTS_CD_01-RSTRCTD_BNFTS_CD_12;  /*scope of medicaid benefits*/
array a8{12} medcad1-medcad12;
array a9{12} ELGBLTY_GRP_CD_01-ELGBLTY_GRP_CD_12;      /*eligibility group code*/
array a10{12} medcad21-medcad32;
array a11{12} meddate1-meddate12;
array a12{12} mccdate1-mccdate12;
array a13{12} ffsdate1-ffsdate12;

do i=1 to 12;
  if a1{i}>=28 then a2{i}=1;else a2{i}=0;   /*Medicaid enrollment based on Medicaid enrollment days*/

  if a3{i} in ('02','04','08') then a5{i}=1;  /*full duals*/
  if a3{i} in ('01','03','05','06') then a5{i}=2;  /*partial duals*/
  if a3{i} in ('00') then a5{i}=3;  /*not Medicare*/
  if a3{i} in ('09') then a5{i}=4; /*Medicare not Medicaid*/

  if a4{i} in ('01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','60','70','80') then a6{i}=1; /*Medicaid through MCO*/
  if a4{i} in ('01','02','03','04','05','06','07','15','16','17','18','19','20','60','70','80') then a6{i}=2; /*Medicaid through MCO - exclude mental health, substance use disorders, and dental*/

  if a7{i} in ('0') then a8{i}=0;  /*not Medicaid*/
  if a7{i} in ('1') then a8{i}=1;  /*full scope of Medicaid*/
  if a7{i} in ('2','3','5','E') then a8{i}=2;  /*Medicaid but restricted benefits*/
  if a7{i} in ('7') then a8{i}=3;  /*Medicaid under an alternative package of benchmark-equivalent coverage as enacted by the Deficit Reduction Act of 2005. adult medicaid expansion eligibls not medically frail who have limited/restricted ABP*/  
  if a7{i} in ('B') then a8{i}=4;  /*Medicaid using Health Opportunity Account*/
  if a7{i} in ('D') then a8{i}=5;  /*Medicaid under "Money Follows the Person" rebalancing demonstration*/
  if a7{i} in ('C') then a8{i}=6;  /*CHIP dental*/
  if a7{i} in ('4','6','A','F') then a8{i}=7; /*restricted benefits not related to cancer care*/
  
  if a9{i} in ('11','12','13','14','15','16','17','18','19','20','21','22','24','37','39','40','41','47','48','49','51','52','59','60') then a10{i}=1;  /*disabled/blind*/
  if a9{i} in ('23','25','26','43') then a10{i}=2;  /*Medicaid+Medicare*/
  if a9{i} in ('01','02','03','04','27','28','32','33','36','38','42','44','46','56','67','71','72','73','74','75') then a10{i}=3;  /*low-income*/
  if a9{i} in ('05','06','07','08','09','29','30','31','34','35','36','45','50','53','54','55','61','62','63','64','65','66''68','69','70','76') then a10{i}=4;  /*other categories unrelated to cancer care*/

  if a2{i}=1 and a10{i}^=4 and a8{i} in (1,2,3,4,5,.) then a11{i}=i;   /*Medicaid enrollment related to cancer care*/
  if a11{i}>. and a6{i}=2 then a12{i}=i;   /*Medicaid managed care enrollment related to cancer care*/
  if a11{i}>. and a12{i}=. then a13{i}=i;   /*Medicaid FFS related to cancer care*/
end;
elig_mons=0;
do j=1 to 12;
 if a11{j}>0 then elig_mons=elig_mons+1;
end;
if elig_mons=12 then do;start1=year*100+1;end1=year*100+12;end;
if 0<elig_mons<12 then do;start1=year*100+coalesce(of meddate1-meddate12);end1=year*100+coalesce(of meddate12-meddate1);end; /*select the first non-missing month as start month and the last non-missing value as end month*/

keep patient_id year dual1-dual12 dual21-dual32 medcad1-medcad12 medcad21-medcad32 meddate1-meddate12 mccdate1-mccdate12 ffsdate1-ffsdate12 state elig_mons start1 end1; 
run;proc sort;by patient_id year;run;

data max;set maxdata2005-maxdata2015;
year=MAX_YR_DT*1;state=STATE_CD;
array a1{12} EL_DAYS_EL_CNT_1-EL_DAYS_EL_CNT_12;
array a2{12} enr1-enr12;
array a3{12} EL_MDCR_DUAL_MO_1-EL_MDCR_DUAL_MO_12;
array a4{12} EL_PHP_TYPE_1_1-EL_PHP_TYPE_1_12;
array a5{12} dual1-dual12;
array a6{12} dual21-dual32;
array a7{12} EL_RSTRCT_BNFT_FLG_1-EL_RSTRCT_BNFT_FLG_12;  /*scope of medicaid benefits*/
array a8{12} medcad1-medcad12;
array a9{12} MAX_ELG_CD_MO_1-MAX_ELG_CD_MO_12;      /*eligibility group code*/ 
array a10{12} medcad21-medcad32;
array a11{12} meddate1-meddate12;
array a12{12} mccdate1-mccdate12;
array a13{12} ffsdate1-ffsdate12;

do i=1 to 12;
  if a1{i}>=28 then a2{i}=1;else a2{i}=0;   /*Medicaid enrollment based on Medicaid enrollment days*/

  if a3{i} in ('01','03','05','06','07','51','53','55','56','57') then a5{i}=2;  /*partial duals*/
  if a3{i} in ('02','04','08','52','54','58') then a5{i}=1;  /*full duals*/
  if a3{i} in ('00') then a5{i}=3;  /*not Medicare*/
  if a3{i} in ('09','59') then a5{i}=5;  /*duals-type unknown*/
  if a3{i} in ('50') then a5{i}=4;  /*Medicare not Medicaid*/

  if a4{i} in ('01','02','03','04','05','06','07','08') then a6{i}=1; /*Medicaid through MCO*/
  if a4{i} in ('01','05','06','07','08') then a6{i}=2; /*Medicaid through MCO - exclude mental health, substance use disorders, and dental*/

  if a9{i} in ('12','22','32','42','52') then a10{i}=1;  /*disabled/blind*/
  if a9{i} in ('11','21','31','41','51') then a10{i}=2;  /*aged*/
  if a9{i} in ('15','17','25','35','45','55') then a10{i}=3;  /*low-income*/
  if a9{i}='00' then a10{i}=0;  /*not Medicaid*/
  if a9{i} in ('99','ZZ') then a10{i}=9; /*unknown eligibility*/
  if a9{i} in ('14','16','24','34','3A','44','48','54') then a10{i}=4; /*other categories unrelated to cancer care*/

  if a7{i} in ('0') then a8{i}=0;  /*not Medicaid*/
  if a7{i} in ('1') then a8{i}=1;  /*full scope of Medicaid*/
  if a7{i} in ('2','3','5','X','Y','Z') then a8{i}=2;  /*Medicaid but restricted benefits*/
  if a7{i} in ('7') then a8{i}=3;  /*Medicaid under an alternative package of benchmark-equivalent coverage as enacted by the Deficit Reduction Act of 2005*/  
  if a7{i} in ('B') then a8{i}=4;  /*Medicaid using Health Opportunity Account*/
  if a7{i} in ('8') then a8{i}=5;  /*Medicaid under "Money Follows the Person" rebalancing demonstration*/
  if a7{i} in ('C') then a8{i}=6;  /*CHIP dental*/
  if a7{i} in ('W') then a8{i}=7; /*medicaid premium assistance only*/
  if a7{i} in ('4','6','A') then a8{i}=8;  /*restricted benefits not related to cancer care*/
  if a7{i} in ('9') then a8{i}=9;  /*benefit restriction unknown*/

  if a2{i}=1 and a10{i} in (1,2,3,9,.) and a8{i} in (1,2,3,4,5,7,9,.) then a11{i}=i; /*Medicaid enrollment related to cancer care*/
  if a11{i}>. and a6{i}=2 then a12{i}=i;   /*Medicaid managed care enrollment related to cancer care*/
  if a11{i}>. and a12{i}=. then a13{i}=i;   /*Medicaid FFS related to cancer care*/
end;
elig_mons=0;
do j=1 to 12;
 if a11{j}>0 then elig_mons=elig_mons+1;
end;
if elig_mons=12 then do;start1=year*100+1;end1=year*100+12;end;
if 0<elig_mons<12 then do;start1=year*100+coalesce(of meddate1-meddate12);end1=year*100+coalesce(of meddate12-meddate1);end; /*select the first non-missing month as start month and the last non-missing value as end month*/

keep patient_id year dual1-dual12 dual21-dual32 medcad1-medcad12 medcad21-medcad32 meddate1-meddate12 mccdate1-mccdate12 ffsdate1-ffsdate12 state elig_mons start1 end1; 
run;proc sort;by patient_id year;run;

data taf1a;set taf1 max;run;proc sort;by patient_id year;run;

data taf2;merge seer_id taf1a;by patient_id;if aa=1;
diff=dxy-year;
if -1<=diff<=2;

/*dual status and medicaid status at cancer diagnosis*/
array a1{12} dual1-dual12;
array a2{12} dual21-dual32;
array a3{12} meddate1-meddate12;
array a4{12} medcad1-medcad12;
array a5{12} medcad21-medcad32;
array a6{12} mccdate1-mccdate12;
array a7{12} ffsdate1-ffsdate12;

if dxy=year then do;
 do i=1 to 12;
  if dxm=i and a1{i} in (1) then dual13=1;  /*full duals*/
  if dxm=i and a1{i} in (2) then dual13=2;  /*partial duals*/
  if dxm=i and a1{i} in (4) then dual13=4;  /*Medicare not Medicaid*/
  if dxm=i and a1{i} in (3) then dual13=5;  /*Not Medicare*/

  if dxm=i and a3{i}>0 then medicaid1=1;   /*medicaid - used to determine enrollment timeline*/
  if dxm=i and a6{i}>0 then medicaid4=1;   /*Medicaid managed care related to cancer care*/
  if dxm=i and a7{i}>0 then medicaid5=1;   /*Medicaid FFS related to cancer care*/

  if dxm=i and a4{i}=1 then medicaid2=1; /*full scope of Medicaid*/
  if dxm=i and a4{i}=2 then medicaid2=2; /*restricted medicaid benefit*/
  if dxm=i and a4{i}=3 then medicaid2=3;/*Medicaid under an alternative package of benchmark-equivalent coverage as enacted by the Deficit Reduction Act of 2005*/  
  if dxm=i and a4{i}=4 then medicaid2=4;  /*Medicaid using Health Opportunity Account*/
  if dxm=i and a4{i}=5 then medicaid2=5;  /*Medicaid under "Money Follows the Person" rebalancing demonstration*/
  if dxm=i and a4{i}=7 then medicaid2=7; /*medicaid premium assistance only*/
  if dxm=i and a4{i}=9 then medicaid2=9;  /*benefit restriction unknown*/

  if dxm=i and a5{i}=1 then medicaid3=1;  /*disabled/blind*/
  if dxm=i and a5{i}=2 then medicaid3=2;  /*aged*/
  if dxm=i and a5{i}=3 then medicaid3=3;  /*low-income*/
  if dxm=i and a5{i}=9 then medicaid3=9;  /*unknown*/
 end;
end;
keep patient_id medicaid1-medicaid5 dual13 year meddate1-meddate12 mccdate1-mccdate12 ffsdate1-ffsdate12 dxy dxm diff;
run; proc sort;by patient_id year;run;

data taf3;set taf2;
array a1{36} meddate1-meddate12 mccdate1-mccdate12 ffsdate1-ffsdate12 ;
array a2{36} meddate13-meddate24 mccdate13-mccdate24 ffsdate13-ffsdate24;
do i=1 to 36; if a1{i}>0 then a2{i}=(year-2005)*12+a1{i}; end;
keep patient_id year dxy dxm diff meddate13-meddate24 mccdate13-mccdate24 ffsdate13-ffsdate24;run;proc sort;by patient_id;run;

data taf4 taf4a;set taf3;by patient_id;
 if not (first.patient_id and last.patient_id) then output taf4;else output taf4a;   /*taf4a includes cases with only one observation*/
run;
data taf5;set taf4;
 array a1{36} meddate13-meddate24 mccdate13-mccdate24 ffsdate13-ffsdate24;
 array a2{36} meddate_at1-meddate_at12 mccdate_at1-mccdate_at12 ffsdate_at1-ffsdate_at12;
 array a3{36} meddate_bf1-meddate_bf12 mccdate_bf1-mccdate_bf12 ffsdate_bf1-ffsdate_bf12;;
 array a4{36} meddate_af1-meddate_af12 mccdate_af1-mccdate_af12 ffsdate_af1-ffsdate_af12;;

 if diff=0 then do;
  do i=1 to 36;
   a2{i}=a1{i};
  end;
 end;
 if diff=-1 then do;
  do i=1 to 36;
   a4{i}=a1{i};
  end;
 end;
 if diff in (1,2) then do;
  do i=1 to 36;
   a3{i}=a1{i};
  end;
 end;
drop i meddate13-meddate24 mccdate13-mccdate24 ffsdate13-ffsdate24; 
 run;proc sort;by patient_id;run;
data taf6;set taf5;if diff=0;
 keep patient_id meddate_at1-meddate_at12 mccdate_at1-mccdate_at12 ffsdate_at1-ffsdate_at12 dxy dxm;
run;proc sort;by patient_id;run;
data taf7;set taf5;if diff=-1;dxyaf=dxy;dxmaf=dxm;
 keep patient_id meddate_af1-meddate_af12 mccdate_af1-mccdate_af12 ffsdate_af1-ffsdate_af12 dxyaf dxmaf;
run;proc sort;by patient_id;run;
data taf8a;set taf5;if diff=2;dxybf1=dxy;dxmbf1=dxm;
 keep patient_id meddate_bf1-meddate_bf12 mccdate_bf1-mccdate_bf12 ffsdate_bf1-ffsdate_bf12 dxybf1 dxmbf1;
run;proc sort;by patient_id;run;
data taf8b;set taf5;if diff=1;dxybf2=dxy;dxmbf2=dxm; 
 array a1{36} meddate_bf1-meddate_bf12 mccdate_bf1-mccdate_bf12 ffsdate_bf1-ffsdate_bf12; 
 array a2{36} meddate_bf13-meddate_bf24 mccdate_bf13-mccdate_bf24 ffsdate_bf13-ffsdate_bf24;
 do i=1 to 36;a2{i}=a1{i};end;
 keep patient_id meddate_bf13-meddate_bf24 mccdate_bf13-mccdate_bf24 ffsdate_bf13-ffsdate_bf24 dxybf2 dxmbf2;
run;proc sort;by patient_id;run;
data taf9;merge taf6 taf7 taf8a taf8b;by patient_id;run;
data taf10 taf11;set taf9;by patient_id;
 if first.patient_id and last.patient_id then output taf10;else output taf11;  /*taf11 includes the cases with more than one records*/
 run;

 /*define Medicaid enrollment patterns*/
data taf12;set taf10;
 array a1{12} meddate_at1-meddate_at12;
 array a2{48} meddate_bf1-meddate_bf24 meddate_at1-meddate_at12 meddate_af1-meddate_af12;
 array a3{48} meddate_bf25-meddate_bf48 meddate_at13-meddate_at24 meddate_af13-meddate_af24;

if dxy=. and dxybf1>. then dxy=dxybf1;if dxy=. and dxybf2>. then dxy=dxybf2;if dxy=. and dxyaf>. then dxy=dxyaf;
if dxm=. and dxmbf1>. then dxm=dxmbf1;if dxm=. and dxmbf2>. then dxm=dxmbf2;if dxm=. and dxmaf>. then dxm=dxmaf;

do n=1 to 48; if a2{n}>0 then a3{n}=1;if a2{n} in (0,.) then a3{n}=0;end;
cc=0;dd=0;ee=0;ff=0;
do i=1 to 12;
 if i=dxm then do;if a1{i}>0 then meddx=1;if a1{i} in (0,.) then meddx=0;  /*medicaid enrollment at diagnosis*/
  do l=(dxm+23) to (dxm+12) by -1;   /*count months of Medicaid enrollment in 12 months before cancer diagnosis*/
   cc=sum(cc,a3{l});
  end;
  do m=(dxm+25) to (dxm+36);    /*count months of Medicaid enrollment in 12 months after diagnosis*/
   dd=sum(dd,a3{m});
  end;
  do h=(dxm+25) to (dxm+30);    /*count months of Medicaid enrollment in 6 months after diagnosis*/
   ee=sum(ee,a3{h});
  end;
  do k=(dxm+25) to (dxm+26);    /*count months of Medicaid enrollment in 2 months after diagnosis*/
   ff=sum(ff,a3{k});
  end;
 end;
end;
if cc>10 then pattern1=1;if 0<cc<=10 then pattern1=0;if cc=0 then pattern1=2;label pattern1='0=discontinuous enroll within 1 year before diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if dd>10 then pattern2=1;if 0<dd<=10 then pattern2=0;if dd=0 then pattern2=2;label pattern2='0=discontinuous enroll withn 1 year after diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if ee=6 then pattern3=1;if 0<ee<6 then pattern3=0;if ee=0 then pattern3=2;label pattern3='0=discontinuous enroll within 6 months after diagnosis,1=continuous enroll (=6 mons),2=no enroll';
if ff=2 then pattern4=1;if 0<ff<2 then pattern4=0;if ff=0 then pattern4=2;label pattern4='0=discontinuous enroll within 2 months after diagnosis,1=continuous enroll (=2 mons),2=no enroll';
if meddx=1 then do;
 if pattern1=1 and pattern2=1 then pattern5=1;if pattern1=1 and pattern2=0 then pattern5=2;if pattern1=0 and pattern2=1 then pattern5=3; 
 if pattern1=0 and pattern2=0 then pattern5=4;if pattern1=2 and pattern2=1 then pattern5=5;if pattern1=2 and pattern2=0 then pattern5=6;
 if pattern1=1 and pattern2=2 then pattern5=7;if pattern1=0 and pattern2=2 then pattern5=8;if pattern1=2 and pattern2=2 then pattern5=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern2=1 then pattern5=9;if pattern1=1 and pattern2=0 then pattern5=10;if pattern1=0 and pattern2=1 then pattern5=11; 
 if pattern1=0 and pattern2=0 then pattern5=12;if pattern1=2 and pattern2=1 then pattern5=13;if pattern1=2 and pattern2=0 then pattern5=14;
 if pattern1=1 and pattern2=2 then pattern5=15;if pattern1=0 and pattern2=2 then pattern5=16;if pattern1=2 and pattern2=2 then pattern5=17;
end;
label pattern5='1=continuous enroll 1yr before,at and 1yr after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern3=1 then pattern6=1;if pattern1=1 and pattern3=0 then pattern6=2;if pattern1=0 and pattern3=1 then pattern6=3; 
 if pattern1=0 and pattern3=0 then pattern6=4;if pattern1=2 and pattern3=1 then pattern6=5;if pattern1=2 and pattern3=0 then pattern6=6;
 if pattern1=1 and pattern3=2 then pattern6=7;if pattern1=0 and pattern3=2 then pattern6=8;if pattern1=2 and pattern3=2 then pattern6=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern3=1 then pattern6=9;if pattern1=1 and pattern3=0 then pattern6=10;if pattern1=0 and pattern3=1 then pattern6=11; 
 if pattern1=0 and pattern3=0 then pattern6=12;if pattern1=2 and pattern3=1 then pattern6=13;if pattern1=2 and pattern3=0 then pattern6=14;
 if pattern1=1 and pattern3=2 then pattern6=15;if pattern1=0 and pattern3=2 then pattern6=16;if pattern1=2 and pattern3=2 then pattern6=17;
end;
label pattern6='1=continuous enroll 1yr before,at and 6months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern4=1 then pattern7=1;if pattern1=1 and pattern4=0 then pattern7=2;if pattern1=0 and pattern4=1 then pattern7=3; 
 if pattern1=0 and pattern4=0 then pattern7=4;if pattern1=2 and pattern4=1 then pattern7=5;if pattern1=2 and pattern4=0 then pattern7=6;
 if pattern1=1 and pattern4=2 then pattern7=7;if pattern1=0 and pattern4=2 then pattern7=8;if pattern1=2 and pattern4=2 then pattern7=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern4=1 then pattern7=9;if pattern1=1 and pattern4=0 then pattern7=10;if pattern1=0 and pattern4=1 then pattern7=11; 
 if pattern1=0 and pattern4=0 then pattern7=12;if pattern1=2 and pattern4=1 then pattern7=13;if pattern1=2 and pattern4=0 then pattern7=14;
 if pattern1=1 and pattern4=2 then pattern7=15;if pattern1=0 and pattern4=2 then pattern7=16;if pattern1=2 and pattern4=2 then pattern7=17;
end;
label pattern7='1=continuous enroll 1yr before,at and 2months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 

gg=0;
do o=(dxm+11) to dxm by -1;   /*count months of Medicaid enrollment in 13-24 months before cancer diagnosis*/
 gg=sum(gg,a3{o});
end;

bb=0;
if pattern5 in (1,2,7,9,10,15) then do;
 do j=1 to 12;
  if j=dxm then do;
   if gg>10 and a3{j}>0 then medenr_mo=a2{j};
   if gg>10 and a3{j}=0 then medenr_mo=a2{j+1};
   if 0<gg<=10 then do;
    if a3{j+11}=0 and a3{j+12}>0 then medenr_mo=a2{j+12};
	if a3{j+11}=0 and a3{j+12}=0 and a3{j+13}>0 then medenr_mo=a2{j+13};
	if a3{j+11}>0 then do;
     do v=(dxm+11) to dxm by -1 until(a3{v}=0);
	  if a2{v}>0 then bb=a2{v};
     end;
	medenr_mo=bb;
    end;
   end;
   if gg=0 and a3{j+13}>0 then medenr_mo=a2{j+13};
   if gg=0 and a3{j+13}=0 then medenr_mo=a2{j+14};
  end;
 end;
end;
if pattern5 in (3,4,8) then do;
 do p=1 to 12;
  if p=dxm then do;
   if a3{p+23}=0 then medenr_mo=a1{p};
   if a3{p+23}>0 then do;
    do q=(dxm+23) to (dxm+12) by -1 until(a3{q}=0);
     if a2{q}>0 then bb=a2{q};
    end;
    medenr_mo=bb;
   end;
  end;
 end;
end;
if pattern5 in (5,6) then do;
 do r=1 to 12;
  if r=dxm then medenr_mo=a1{r};
 end;
end;
if pattern5=11 then do;
 do c=1 to 12;
  if c=dxm then do;
   if a3{c+23}>0 then do;
    do d=(dxm+23) to (dxm+12) by -1 until(a3{d}=0);
     if a2{d}>0 then bb=a2{d};
    end;
    medenr_mo=bb;
   end;
   if a3{c+23}=0 and a3{c+25}>0 then medenr_mo=a2{c+25};
   if a3{c+23}=0 and a3{c+25}=0 then medenr_mo=a2{c+26};
  end;
 end; 
end;
if pattern5=12 then do;
 do e=1 to 12;
  if e=dxm then do;
   if a3{e+23}>0 then do;
    do f=(dxm+23) to (dxm+12) by -1 until(a3{f}=0);
     if a2{f}>0 then bb=a2{f};
    end;
   end;
   if a3{e+23}=0 then do;
    do g=(dxm+25) to (dxm+36) until(a3{g}>0);
     if a2{g}>0 then bb=a2{g};
    end;
   end;
   medenr_mo=bb;
  end;
 end;
end;
if pattern5 in (13,14) then do; 
 do s=1 to 12;
  if s=dxm then do;
   if a3{s+25}=0 then do;
    do t=(dxm+25) to (dxm+36) until(a3{t}>0);
     if a2{t}>0 then bb=a2{t};
    end;
    medenr_mo=bb;
   end;
   if a3{s+25}>0 then medenr_mo=a2{s+25};
  end;
 end;
end;
if pattern5=18 then do;
 do u=1 to 12;
  if u=dxm then medenr_mo=a1{u};
 end;
end;
drop c d e f g h i j k l m n o p q r s t u v bb cc dd ee ff gg;
run;
data taf12a;set taf12;keep patient_id medenr_mo pattern5-pattern7;run;proc sort;by patient_id;run;

proc freq data=taf11 noprint;tables patient_id/out=taf11a;run;
data taf11b;set taf11a;if count=2;b=1;keep patient_id b;run;proc sort;by patient_id;run;
proc sort data=taf11;by patient_id;run;
data taf13;merge taf11 taf11b;by patient_id;if b=1;drop b;run;
data taf13a;set taf13;if mod(_n_,2)=1;
 array a1{48} meddate_bf1-meddate_bf24 meddate_at1-meddate_at12 meddate_af1-meddate_af12;
 array a2{48} meddate1_bf1-meddate1_bf24 meddate1_at1-meddate1_at12 meddate1_af1-meddate1_af12;
 do i=1 to 48;a2{i}=a1{i};end;
 drop meddate_bf1-meddate_bf24 meddate_at1-meddate_at12 meddate_af1-meddate_af12 i;
run;proc sort;by patient_id;run;
data taf13b;set taf13;if mod(_n_,2)=0;
 array a1{48} meddate_bf1-meddate_bf24 meddate_at1-meddate_at12 meddate_af1-meddate_af12;
 array a2{48} meddate2_bf1-meddate2_bf24 meddate2_at1-meddate2_at12 meddate2_af1-meddate2_af12;
 do i=1 to 48;a2{i}=a1{i};end;
 keep patient_id meddate2_bf1-meddate2_bf24 meddate2_at1-meddate2_at12 meddate2_af1-meddate2_af12;
run;proc sort;by patient_id;run;
data taf13c;merge taf13a taf13b;by patient_id;
 array a1{48} meddate1_bf1-meddate1_bf24 meddate1_at1-meddate1_at12 meddate1_af1-meddate1_af12;
 array a2{48} meddate2_bf1-meddate2_bf24 meddate2_at1-meddate2_at12 meddate2_af1-meddate2_af12;
 array a3{48} meddate_bf1-meddate_bf24 meddate_at1-meddate_at12 meddate_af1-meddate_af12;
 do i=1 to 48;
  if a1{i}=.|a2{i}=. then a3{i}=sum(a1{i},a2{i});
  if a1{i}>=0 and a2{i}>=0 then a3{i}=min(a1{i},a2{i});
 end;
 drop i meddate1_bf1-meddate1_bf24 meddate1_at1-meddate1_at12 meddate1_af1-meddate1_af12 meddate2_bf1-meddate2_bf24 meddate2_at1-meddate2_at12 meddate2_af1-meddate2_af12;
run;
data taf13d;set taf13c;
 array a1{12} meddate_at1-meddate_at12;
 array a2{48} meddate_bf1-meddate_bf24 meddate_at1-meddate_at12 meddate_af1-meddate_af12;
 array a3{48} meddate_bf25-meddate_bf48 meddate_at21-meddate_at32 meddate_af21-meddate_af32;

if dxy=. and dxybf1>. then dxy=dxybf1;if dxy=. and dxybf2>. then dxy=dxybf2;if dxy=. and dxyaf>. then dxy=dxyaf;
if dxm=. and dxmbf1>. then dxm=dxmbf1;if dxm=. and dxmbf2>. then dxm=dxmbf2;if dxm=. and dxmaf>. then dxm=dxmaf;

do n=1 to 48; if a2{n}>0 then a3{n}=1;if a2{n} in (0,.) then a3{n}=0;end;
cc=0;dd=0;ee=0;ff=0;
do i=1 to 12;
 if i=dxm then do;if a1{i}>0 then meddx=1;if a1{i} in (0,.) then meddx=0;  /*medicaid enrollment at diagnosis*/
  do l=(dxm+23) to (dxm+12) by -1;   /*count months of Medicaid enrollment in 12 months before cancer diagnosis*/
   cc=sum(cc,a3{l});
  end;
  do m=(dxm+25) to (dxm+36);    /*count months of Medicaid enrollment in 12 months after diagnosis*/
   dd=sum(dd,a3{m});
  end;
  do h=(dxm+25) to (dxm+30);    /*count months of Medicaid enrollment in 6 months after diagnosis*/
   ee=sum(ee,a3{h});
  end;
  do k=(dxm+25) to (dxm+26);    /*count months of Medicaid enrollment in 2 months after diagnosis*/
   ff=sum(ff,a3{k});
  end;
 end;
end;
if cc>10 then pattern1=1;if 0<cc<=10 then pattern1=0;if cc=0 then pattern1=2;label pattern1='0=discontinuous enroll within 1 year before diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if dd>10 then pattern2=1;if 0<dd<=10 then pattern2=0;if dd=0 then pattern2=2;label pattern2='0=discontinuous enroll withn 1 year after diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if ee=6 then pattern3=1;if 0<ee<6 then pattern3=0;if ee=0 then pattern3=2;label pattern3='0=discontinuous enroll within 6 months after diagnosis,1=continuous enroll (=6 mons),2=no enroll';
if ff=2 then pattern4=1;if 0<ff<2 then pattern4=0;if ff=0 then pattern4=2;label pattern4='0=discontinuous enroll within 2 months after diagnosis,1=continuous enroll (=2 mons),2=no enroll';
if meddx=1 then do;
 if pattern1=1 and pattern2=1 then pattern5=1;if pattern1=1 and pattern2=0 then pattern5=2;if pattern1=0 and pattern2=1 then pattern5=3; 
 if pattern1=0 and pattern2=0 then pattern5=4;if pattern1=2 and pattern2=1 then pattern5=5;if pattern1=2 and pattern2=0 then pattern5=6;
 if pattern1=1 and pattern2=2 then pattern5=7;if pattern1=0 and pattern2=2 then pattern5=8;if pattern1=2 and pattern2=2 then pattern5=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern2=1 then pattern5=9;if pattern1=1 and pattern2=0 then pattern5=10;if pattern1=0 and pattern2=1 then pattern5=11; 
 if pattern1=0 and pattern2=0 then pattern5=12;if pattern1=2 and pattern2=1 then pattern5=13;if pattern1=2 and pattern2=0 then pattern5=14;
 if pattern1=1 and pattern2=2 then pattern5=15;if pattern1=0 and pattern2=2 then pattern5=16;if pattern1=2 and pattern2=2 then pattern5=17;
end;
label pattern5='1=continuous enroll 1yr before,at and 1yr after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern3=1 then pattern6=1;if pattern1=1 and pattern3=0 then pattern6=2;if pattern1=0 and pattern3=1 then pattern6=3; 
 if pattern1=0 and pattern3=0 then pattern6=4;if pattern1=2 and pattern3=1 then pattern6=5;if pattern1=2 and pattern3=0 then pattern6=6;
 if pattern1=1 and pattern3=2 then pattern6=7;if pattern1=0 and pattern3=2 then pattern6=8;if pattern1=2 and pattern3=2 then pattern6=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern3=1 then pattern6=9;if pattern1=1 and pattern3=0 then pattern6=10;if pattern1=0 and pattern3=1 then pattern6=11; 
 if pattern1=0 and pattern3=0 then pattern6=12;if pattern1=2 and pattern3=1 then pattern6=13;if pattern1=2 and pattern3=0 then pattern6=14;
 if pattern1=1 and pattern3=2 then pattern6=15;if pattern1=0 and pattern3=2 then pattern6=16;if pattern1=2 and pattern3=2 then pattern6=17;
end;
label pattern6='1=continuous enroll 1yr before,at and 6months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern4=1 then pattern7=1;if pattern1=1 and pattern4=0 then pattern7=2;if pattern1=0 and pattern4=1 then pattern7=3; 
 if pattern1=0 and pattern4=0 then pattern7=4;if pattern1=2 and pattern4=1 then pattern7=5;if pattern1=2 and pattern4=0 then pattern7=6;
 if pattern1=1 and pattern4=2 then pattern7=7;if pattern1=0 and pattern4=2 then pattern7=8;if pattern1=2 and pattern4=2 then pattern7=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern4=1 then pattern7=9;if pattern1=1 and pattern4=0 then pattern7=10;if pattern1=0 and pattern4=1 then pattern7=11; 
 if pattern1=0 and pattern4=0 then pattern7=12;if pattern1=2 and pattern4=1 then pattern7=13;if pattern1=2 and pattern4=0 then pattern7=14;
 if pattern1=1 and pattern4=2 then pattern7=15;if pattern1=0 and pattern4=2 then pattern7=16;if pattern1=2 and pattern4=2 then pattern7=17;
end;
label pattern7='1=continuous enroll 1yr before,at and 2months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 

gg=0;
do o=(dxm+11) to dxm by -1;   /*count months of Medicaid enrollment in 13-24 months before cancer diagnosis*/
 gg=sum(gg,a3{o});
end;

bb=0;
if pattern5 in (1,2,7,9,10,15) then do;
 do j=1 to 12;
  if j=dxm then do;
   if gg>10 and a3{j}>0 then medenr_mo=a2{j};
   if gg>10 and a3{j}=0 then medenr_mo=a2{j+1};
   if 0<gg<=10 then do;
    if a3{j+11}=0 and a3{j+12}>0 then medenr_mo=a2{j+12};
	if a3{j+11}=0 and a3{j+12}=0 and a3{j+13}>0 then medenr_mo=a2{j+13};
	if a3{j+11}>0 then do;
     do v=(dxm+11) to dxm by -1 until(a3{v}=0);
	  if a2{v}>0 then bb=a2{v};
     end;
	medenr_mo=bb;
    end;
   end;
   if gg=0 and a3{j+13}>0 then medenr_mo=a2{j+13};
   if gg=0 and a3{j+13}=0 then medenr_mo=a2{j+14};
  end;
 end;
end;
if pattern5 in (3,4,8) then do;
 do p=1 to 12;
  if p=dxm then do;
   if a3{p+23}=0 then medenr_mo=a1{p};
   if a3{p+23}>0 then do;
    do q=(dxm+23) to (dxm+12) by -1 until(a3{q}=0);
     if a2{q}>0 then bb=a2{q};
    end;
    medenr_mo=bb;
   end;
  end;
 end;
end;
if pattern5 in (5,6) then do;
 do r=1 to 12;
  if r=dxm then medenr_mo=a1{r};
 end;
end;
if pattern5=11 then do;
 do c=1 to 12;
  if c=dxm then do;
   if a3{c+23}>0 then do;
    do d=(dxm+23) to (dxm+12) by -1 until(a3{d}=0);
     if a2{d}>0 then bb=a2{d};
    end;
    medenr_mo=bb;
   end;
   if a3{c+23}=0 and a3{c+25}>0 then medenr_mo=a2{c+25};
   if a3{c+23}=0 and a3{c+25}=0 then medenr_mo=a2{c+26};
  end;
 end; 
end;
if pattern5=12 then do;
 do e=1 to 12;
  if e=dxm then do;
   if a3{e+23}>0 then do;
    do f=(dxm+23) to (dxm+12) by -1 until(a3{f}=0);
     if a2{f}>0 then bb=a2{f};
    end;
   end;
   if a3{e+23}=0 then do;
    do g=(dxm+25) to (dxm+36) until(a3{g}>0);
     if a2{g}>0 then bb=a2{g};
    end;
   end;
   medenr_mo=bb;
  end;
 end;
end;
if pattern5 in (13,14) then do; 
 do s=1 to 12;
  if s=dxm then do;
   if a3{s+25}=0 then do;
    do t=(dxm+25) to (dxm+36) until(a3{t}>0);
     if a2{t}>0 then bb=a2{t};
    end;
    medenr_mo=bb;
   end;
   if a3{s+25}>0 then medenr_mo=a2{s+25};
  end;
 end;
end;
if pattern5=18 then do;
 do u=1 to 12;
  if u=dxm then medenr_mo=a1{u};
 end;
end;
drop c d e f g h i j k l m n o p q r s t u v bb cc dd ee ff gg;
run;
data taf13e;set taf13d;keep patient_id medenr_mo pattern5-pattern7;run;proc sort;by patient_id;run;

data taf11c;set taf11a;if count=3;b=1;keep patient_id b;run;proc sort;by patient_id;run;
data taf14;merge taf11 taf11c;by patient_id;if b=1;drop b;run;
data taf14a;set taf14;if mod(_n_,3)=1;
 array a1{48} meddate_bf1-meddate_bf24 meddate_at1-meddate_at12 meddate_af1-meddate_af12;
 array a2{48} meddate1_bf1-meddate1_bf24 meddate1_at1-meddate1_at12 meddate1_af1-meddate1_af12;
 do i=1 to 48;a2{i}=a1{i};end;
 drop meddate_bf1-meddate_bf24 meddate_at1-meddate_at12 meddate_af1-meddate_af12 i;
run;proc sort;by patient_id;run;
data taf14b;set taf14;if mod(_n_,3)=2;
 array a1{48} meddate_bf1-meddate_bf24 meddate_at1-meddate_at12 meddate_af1-meddate_af12;
 array a2{48} meddate2_bf1-meddate2_bf24 meddate2_at1-meddate2_at12 meddate2_af1-meddate2_af12;
 do i=1 to 48;a2{i}=a1{i};end;
 keep patient_id meddate2_bf1-meddate2_bf24 meddate2_at1-meddate2_at12 meddate2_af1-meddate2_af12;
run;proc sort;by patient_id;run;
data taf14c;set taf14;if mod(_n_,3)=0;
 array a1{48} meddate_bf1-meddate_bf24 meddate_at1-meddate_at12 meddate_af1-meddate_af12;
 array a2{48} meddate3_bf1-meddate3_bf24 meddate3_at1-meddate3_at12 meddate3_af1-meddate3_af12;
 do i=1 to 48;a2{i}=a1{i};end;
 keep patient_id meddate3_bf1-meddate3_bf24 meddate3_at1-meddate3_at12 meddate3_af1-meddate3_af12;
run;proc sort;by patient_id;run;
data taf14d;merge taf14a taf14b taf14c;by patient_id;
 array a1{48} meddate1_bf1-meddate1_bf24 meddate1_at1-meddate1_at12 meddate1_af1-meddate1_af12;
 array a2{48} meddate2_bf1-meddate2_bf24 meddate2_at1-meddate2_at12 meddate2_af1-meddate2_af12;
 array a3{48} meddate3_bf1-meddate3_bf24 meddate3_at1-meddate3_at12 meddate3_af1-meddate3_af12;
 array a4{48} meddate_bf1-meddate_bf24 meddate_at1-meddate_at12 meddate_af1-meddate_af12;
 do i=1 to 48;
  if (a1{i}=. and a2{i}=. and a3{i}>=0)|(a1{i}=. and a2{i}>=0 and a3{i}=.)|(a1{i}>=0 and a2{i}=. and a3{i}=.) then a4{i}=sum(a1{i},a2{i},a3{i});
  if (a1{i}>=0 and a2{i}>=0)|(a1{i}>=0 and a3{i}>=0)|(a2{i}>=0 and a3{i}>=0) then a4{i}=min(a1{i},a2{i},a3{i});
 end;
 drop i meddate1_bf1-meddate1_bf24 meddate1_at1-meddate1_at12 meddate1_af1-meddate1_af12 
        meddate2_bf1-meddate2_bf24 meddate2_at1-meddate2_at12 meddate2_af1-meddate2_af12
        meddate3_bf1-meddate3_bf24 meddate3_at1-meddate3_at12 meddate3_af1-meddate3_af12;
run;
data taf14e;set taf14d;
 array a1{12} meddate_at1-meddate_at12;
 array a2{48} meddate_bf1-meddate_bf24 meddate_at1-meddate_at12 meddate_af1-meddate_af12;
 array a3{48} meddate_bf25-meddate_bf48 meddate_at21-meddate_at32 meddate_af21-meddate_af32;

if dxy=. and dxybf1>. then dxy=dxybf1;if dxy=. and dxybf2>. then dxy=dxybf2;if dxy=. and dxyaf>. then dxy=dxyaf;
if dxm=. and dxmbf1>. then dxm=dxmbf1;if dxm=. and dxmbf2>. then dxm=dxmbf2;if dxm=. and dxmaf>. then dxm=dxmaf;

do n=1 to 48; if a2{n}>0 then a3{n}=1;if a2{n} in (0,.) then a3{n}=0;end;
cc=0;dd=0;ee=0;ff=0;
do i=1 to 12;
 if i=dxm then do;if a1{i}>0 then meddx=1;if a1{i} in (0,.) then meddx=0;  /*medicaid enrollment at diagnosis*/
  do l=(dxm+23) to (dxm+12) by -1;   /*count months of Medicaid enrollment in 12 months before cancer diagnosis*/
   cc=sum(cc,a3{l});
  end;
  do m=(dxm+25) to (dxm+36);    /*count months of Medicaid enrollment in 12 months after diagnosis*/
   dd=sum(dd,a3{m});
  end;
  do h=(dxm+25) to (dxm+30);    /*count months of Medicaid enrollment in 6 months after diagnosis*/
   ee=sum(ee,a3{h});
  end;
  do k=(dxm+25) to (dxm+26);    /*count months of Medicaid enrollment in 2 months after diagnosis*/
   ff=sum(ff,a3{k});
  end;
 end;
end;
if cc>10 then pattern1=1;if 0<cc<=10 then pattern1=0;if cc=0 then pattern1=2;label pattern1='0=discontinuous enroll within 1 year before diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if dd>10 then pattern2=1;if 0<dd<=10 then pattern2=0;if dd=0 then pattern2=2;label pattern2='0=discontinuous enroll withn 1 year after diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if ee=6 then pattern3=1;if 0<ee<6 then pattern3=0;if ee=0 then pattern3=2;label pattern3='0=discontinuous enroll within 6 months after diagnosis,1=continuous enroll (=6 mons),2=no enroll';
if ff=2 then pattern4=1;if 0<ff<2 then pattern4=0;if ff=0 then pattern4=2;label pattern4='0=discontinuous enroll within 2 months after diagnosis,1=continuous enroll (=2 mons),2=no enroll';
if meddx=1 then do;
 if pattern1=1 and pattern2=1 then pattern5=1;if pattern1=1 and pattern2=0 then pattern5=2;if pattern1=0 and pattern2=1 then pattern5=3; 
 if pattern1=0 and pattern2=0 then pattern5=4;if pattern1=2 and pattern2=1 then pattern5=5;if pattern1=2 and pattern2=0 then pattern5=6;
 if pattern1=1 and pattern2=2 then pattern5=7;if pattern1=0 and pattern2=2 then pattern5=8;if pattern1=2 and pattern2=2 then pattern5=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern2=1 then pattern5=9;if pattern1=1 and pattern2=0 then pattern5=10;if pattern1=0 and pattern2=1 then pattern5=11; 
 if pattern1=0 and pattern2=0 then pattern5=12;if pattern1=2 and pattern2=1 then pattern5=13;if pattern1=2 and pattern2=0 then pattern5=14;
 if pattern1=1 and pattern2=2 then pattern5=15;if pattern1=0 and pattern2=2 then pattern5=16;if pattern1=2 and pattern2=2 then pattern5=17;
end;
label pattern5='1=continuous enroll 1yr before,at and 1yr after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern3=1 then pattern6=1;if pattern1=1 and pattern3=0 then pattern6=2;if pattern1=0 and pattern3=1 then pattern6=3; 
 if pattern1=0 and pattern3=0 then pattern6=4;if pattern1=2 and pattern3=1 then pattern6=5;if pattern1=2 and pattern3=0 then pattern6=6;
 if pattern1=1 and pattern3=2 then pattern6=7;if pattern1=0 and pattern3=2 then pattern6=8;if pattern1=2 and pattern3=2 then pattern6=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern3=1 then pattern6=9;if pattern1=1 and pattern3=0 then pattern6=10;if pattern1=0 and pattern3=1 then pattern6=11; 
 if pattern1=0 and pattern3=0 then pattern6=12;if pattern1=2 and pattern3=1 then pattern6=13;if pattern1=2 and pattern3=0 then pattern6=14;
 if pattern1=1 and pattern3=2 then pattern6=15;if pattern1=0 and pattern3=2 then pattern6=16;if pattern1=2 and pattern3=2 then pattern6=17;
end;
label pattern6='1=continuous enroll 1yr before,at and 6months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern4=1 then pattern7=1;if pattern1=1 and pattern4=0 then pattern7=2;if pattern1=0 and pattern4=1 then pattern7=3; 
 if pattern1=0 and pattern4=0 then pattern7=4;if pattern1=2 and pattern4=1 then pattern7=5;if pattern1=2 and pattern4=0 then pattern7=6;
 if pattern1=1 and pattern4=2 then pattern7=7;if pattern1=0 and pattern4=2 then pattern7=8;if pattern1=2 and pattern4=2 then pattern7=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern4=1 then pattern7=9;if pattern1=1 and pattern4=0 then pattern7=10;if pattern1=0 and pattern4=1 then pattern7=11; 
 if pattern1=0 and pattern4=0 then pattern7=12;if pattern1=2 and pattern4=1 then pattern7=13;if pattern1=2 and pattern4=0 then pattern7=14;
 if pattern1=1 and pattern4=2 then pattern7=15;if pattern1=0 and pattern4=2 then pattern7=16;if pattern1=2 and pattern4=2 then pattern7=17;
end;
label pattern7='1=continuous enroll 1yr before,at and 2months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 

gg=0;
do o=(dxm+11) to dxm by -1;   /*count months of Medicaid enrollment in 13-24 months before cancer diagnosis*/
 gg=sum(gg,a3{o});
end;

bb=0;
if pattern5 in (1,2,7,9,10,15) then do;
 do j=1 to 12;
  if j=dxm then do;
   if gg>10 and a3{j}>0 then medenr_mo=a2{j};
   if gg>10 and a3{j}=0 then medenr_mo=a2{j+1};
   if 0<gg<=10 then do;
    if a3{j+11}=0 and a3{j+12}>0 then medenr_mo=a2{j+12};
	if a3{j+11}=0 and a3{j+12}=0 and a3{j+13}>0 then medenr_mo=a2{j+13};
	if a3{j+11}>0 then do;
     do v=(dxm+11) to dxm by -1 until(a3{v}=0);
	  if a2{v}>0 then bb=a2{v};
     end;
	medenr_mo=bb;
    end;
   end;
   if gg=0 and a3{j+13}>0 then medenr_mo=a2{j+13};
   if gg=0 and a3{j+13}=0 then medenr_mo=a2{j+14};
  end;
 end;
end;
if pattern5 in (3,4,8) then do;
 do p=1 to 12;
  if p=dxm then do;
   if a3{p+23}=0 then medenr_mo=a1{p};
   if a3{p+23}>0 then do;
    do q=(dxm+23) to (dxm+12) by -1 until(a3{q}=0);
     if a2{q}>0 then bb=a2{q};
    end;
    medenr_mo=bb;
   end;
  end;
 end;
end;
if pattern5 in (5,6) then do;
 do r=1 to 12;
  if r=dxm then medenr_mo=a1{r};
 end;
end;
if pattern5=11 then do;
 do c=1 to 12;
  if c=dxm then do;
   if a3{c+23}>0 then do;
    do d=(dxm+23) to (dxm+12) by -1 until(a3{d}=0);
     if a2{d}>0 then bb=a2{d};
    end;
    medenr_mo=bb;
   end;
   if a3{c+23}=0 and a3{c+25}>0 then medenr_mo=a2{c+25};
   if a3{c+23}=0 and a3{c+25}=0 then medenr_mo=a2{c+26};
  end;
 end; 
end;
if pattern5=12 then do;
 do e=1 to 12;
  if e=dxm then do;
   if a3{e+23}>0 then do;
    do f=(dxm+23) to (dxm+12) by -1 until(a3{f}=0);
     if a2{f}>0 then bb=a2{f};
    end;
   end;
   if a3{e+23}=0 then do;
    do g=(dxm+25) to (dxm+36) until(a3{g}>0);
     if a2{g}>0 then bb=a2{g};
    end;
   end;
   medenr_mo=bb;
  end;
 end;
end;
if pattern5 in (13,14) then do; 
 do s=1 to 12;
  if s=dxm then do;
   if a3{s+25}=0 then do;
    do t=(dxm+25) to (dxm+36) until(a3{t}>0);
     if a2{t}>0 then bb=a2{t};
    end;
    medenr_mo=bb;
   end;
   if a3{s+25}>0 then medenr_mo=a2{s+25};
  end;
 end;
end;
if pattern5=18 then do;
 do u=1 to 12;
  if u=dxm then medenr_mo=a1{u};
 end;
end;
drop c d e f g h i j k l m n o p q r s t u v bb cc dd ee ff gg;
run;
data taf14f;set taf14e;keep patient_id medenr_mo pattern5-pattern7;run;proc sort;by patient_id;run;
data taf5a;set taf4a;
 array a1{12} meddate13-meddate24;
 array a2{12} at1-at12;
 array a3{12} bf1-bf12;
 array a4{12} bf13-bf24;
 array a5{12} af1-af12;
 array a6{48} bf1-bf24 at1-at12 af1-af12;
 array a7{48} ff1-ff48;
cc=0;dd=0;ee=0;ff=0;
do i=1 to 12;
 if diff=0 then a2{i}=a1{i};
 if diff=1 then a4{i}=a1{i};
 if diff=2 then a3{i}=a1{i};
 if diff=-1 then a5{i}=a1{i};
end;
do n=1 to 48;
 if a6{n}>0 then a7{n}=1;if a6{n} in (0,.) then a7{n}=0;
end;
do h=1 to 12;
 if h=dxm and a2{h}>0 then meddx=1;
end;
do k=(dxm+23) to (dxm+12) by -1;  /*count months of Medicaid enrollment in 12 months before cancer diagnosis*/
 cc=sum(cc,a7{k});
end;
do l=(dxm+25) to (dxm+36);  /*count months of Medicaid enrollment in 12 months after diagnosis*/
 dd=sum(dd,a7{l});
end;
do m=(dxm+25) to (dxm+30); /*count months of Medicaid enrollment in 6 months after diagnosis*/
 ee=sum(ee,a7{m});
end;
do o=(dxm+25) to (dxm+26); /*count months of Medicaid enrollment in 2 months after diagnosis*/
 ff=sum(ff,a7{o});
end;

if cc>10 then pattern1=1;if 0<cc<=10 then pattern1=0;if cc=0 then pattern1=2;label pattern1='0=discontinuous enroll within 1 year before diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if dd>10 then pattern2=1;if 0<dd<=10 then pattern2=0;if dd=0 then pattern2=2;label pattern2='0=discontinuous enroll withn 1 year after diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if ee=6 then pattern3=1;if 0<ee<6 then pattern3=0;if ee=0 then pattern3=2;label pattern3='0=discontinuous enroll within 6 months after diagnosis,1=continuous enroll (=6 mons),2=no enroll';
if ff=2 then pattern4=1;if 0<ff<2 then pattern4=0;if ff=0 then pattern4=2;label pattern4='0=discontinuous enroll within 2 months after diagnosis,1=continuous enroll (=2 mons),2=no enroll';
if meddx=1 then do;
 if pattern1=1 and pattern2=1 then pattern5=1;if pattern1=1 and pattern2=0 then pattern5=2;if pattern1=0 and pattern2=1 then pattern5=3; 
 if pattern1=0 and pattern2=0 then pattern5=4;if pattern1=2 and pattern2=1 then pattern5=5;if pattern1=2 and pattern2=0 then pattern5=6;
 if pattern1=1 and pattern2=2 then pattern5=7;if pattern1=0 and pattern2=2 then pattern5=8;if pattern1=2 and pattern2=2 then pattern5=18;
end;
if meddx in (0,.) then do;
 if pattern1=1 and pattern2=1 then pattern5=9;if pattern1=1 and pattern2=0 then pattern5=10;if pattern1=0 and pattern2=1 then pattern5=11; 
 if pattern1=0 and pattern2=0 then pattern5=12;if pattern1=2 and pattern2=1 then pattern5=13;if pattern1=2 and pattern2=0 then pattern5=14;
 if pattern1=1 and pattern2=2 then pattern5=15;if pattern1=0 and pattern2=2 then pattern5=16;if pattern1=2 and pattern2=2 then pattern5=17;
end;
label pattern5='1=continuous enroll 1yr before,at and 1yr after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern3=1 then pattern6=1;if pattern1=1 and pattern3=0 then pattern6=2;if pattern1=0 and pattern3=1 then pattern6=3; 
 if pattern1=0 and pattern3=0 then pattern6=4;if pattern1=2 and pattern3=1 then pattern6=5;if pattern1=2 and pattern3=0 then pattern6=6;
 if pattern1=1 and pattern3=2 then pattern6=7;if pattern1=0 and pattern3=2 then pattern6=8;if pattern1=2 and pattern3=2 then pattern6=18;
end;
if meddx in (0,.) then do;
 if pattern1=1 and pattern3=1 then pattern6=9;if pattern1=1 and pattern3=0 then pattern6=10;if pattern1=0 and pattern3=1 then pattern6=11; 
 if pattern1=0 and pattern3=0 then pattern6=12;if pattern1=2 and pattern3=1 then pattern6=13;if pattern1=2 and pattern3=0 then pattern6=14;
 if pattern1=1 and pattern3=2 then pattern6=15;if pattern1=0 and pattern3=2 then pattern6=16;if pattern1=2 and pattern3=2 then pattern6=17;
end;
label pattern6='1=continuous enroll 1yr before,at and 6months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern4=1 then pattern7=1;if pattern1=1 and pattern4=0 then pattern7=2;if pattern1=0 and pattern4=1 then pattern7=3; 
 if pattern1=0 and pattern4=0 then pattern7=4;if pattern1=2 and pattern4=1 then pattern7=5;if pattern1=2 and pattern4=0 then pattern7=6;
 if pattern1=1 and pattern4=2 then pattern7=7;if pattern1=0 and pattern4=2 then pattern7=8;if pattern1=2 and pattern4=2 then pattern7=18;
end;
if meddx in (0,.) then do;
 if pattern1=1 and pattern4=1 then pattern7=9;if pattern1=1 and pattern4=0 then pattern7=10;if pattern1=0 and pattern4=1 then pattern7=11; 
 if pattern1=0 and pattern4=0 then pattern7=12;if pattern1=2 and pattern4=1 then pattern7=13;if pattern1=2 and pattern4=0 then pattern7=14;
 if pattern1=1 and pattern4=2 then pattern7=15;if pattern1=0 and pattern4=2 then pattern7=16;if pattern1=2 and pattern4=2 then pattern7=17;
end;
label pattern7='1=continuous enroll 1yr before,at and 2months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 

gg=0;
do o=(dxm+11) to dxm by -1;   /*count months of Medicaid enrollment in 13-24 months before cancer diagnosis*/
 gg=sum(gg,a7{o});
end;

bb=0;
if pattern5 in (1,2,7,9,10,15) then do;
 do j=1 to 12;
  if j=dxm then do;
   if gg>10 and a7{j}>0 then medenr_mo=a6{j};
   if gg>10 and a7{j}=0 then medenr_mo=a6{j+1};
   if 0<gg<=10 then do;
    if a7{j+11}=0 and a7{j+12}>0 then medenr_mo=a6{j+12};
	if a7{j+11}=0 and a7{j+12}=0 and a7{j+13}>0 then medenr_mo=a6{j+13};
	if a7{j+11}>0 then do;
     do v=(dxm+11) to dxm by -1 until(a7{v}=0);
	  if a6{v}>0 then bb=a6{v};
     end;
	medenr_mo=bb;
    end;
   end;
   if gg=0 and a7{j+13}>0 then medenr_mo=a6{j+13};
   if gg=0 and a7{j+13}=0 then medenr_mo=a6{j+14};
  end;
 end;
end;
if pattern5 in (3,4,8) then do;
 do p=1 to 12;
  if p=dxm then do;
   if a7{p+23}=0 then medenr_mo=a2{p};
   if a7{p+23}>0 then do;
    do q=(dxm+23) to (dxm+12) by -1 until(a7{q}=0);
     if a6{q}>0 then bb=a6{q};
    end;
    medenr_mo=bb;
   end;
  end;
 end;
end;
if pattern5 in (5,6) then do;
 do r=1 to 12;
  if r=dxm then medenr_mo=a2{r};
 end;
end;
if pattern5=11 then do;
 do c=1 to 12;
  if c=dxm then do;
   if a7{c+23}>0 then do;
    do d=(dxm+23) to (dxm+12) by -1 until(a7{d}=0);
     if a6{d}>0 then bb=a6{d};
    end;
    medenr_mo=bb;
   end;
   if a7{c+23}=0 and a7{c+25}>0 then medenr_mo=a6{c+25};
   if a7{c+23}=0 and a7{c+25}=0 then medenr_mo=a6{c+26};
  end;
 end; 
end;
if pattern5=12 then do;
 do e=1 to 12;
  if e=dxm then do;
   if a7{e+23}>0 then do;
    do f=(dxm+23) to (dxm+12) by -1 until(a7{f}=0);
     if a6{f}>0 then bb=a6{f};
    end;
   end;
   if a7{e+23}=0 then do;
    do g=(dxm+25) to (dxm+36) until(a7{g}>0);
     if a6{g}>0 then bb=a6{g};
    end;
   end;
   medenr_mo=bb;
  end;
 end;
end;
if pattern5 in (13,14) then do; 
 do s=1 to 12;
  if s=dxm then do;
   if a7{s+25}=0 then do;
    do t=(dxm+25) to (dxm+36) until(a7{t}>0);
     if a6{t}>0 then bb=a6{t};
    end;
    medenr_mo=bb;
   end;
   if a7{s+25}>0 then medenr_mo=a6{s+25};
  end;
 end;
end;
if pattern5=18 then do;
 do u=1 to 12;
  if u=dxm then medenr_mo=a2{u};
 end;
end;
drop c d e f g h i j k l m n o p q r s t u v bb cc dd ee ff gg;
run;
data taf5b;set taf5a;keep patient_id medenr_mo pattern5-pattern7;run;proc sort;by patient_id;run;

data taf15;set taf12a taf13e taf14f taf5b;run;proc sort;by patient_id;run;  /*overall medicaid enrollment timeline*/

/*Medicaid managed care enrollment timeline*/
data mcc1;set taf10;
 array a1{12} mccdate_at1-mccdate_at12;
 array a2{48} mccdate_bf1-mccdate_bf24 mccdate_at1-mccdate_at12 mccdate_af1-mccdate_af12;
 array a3{48} mccdate_bf25-mccdate_bf48 mccdate_at13-mccdate_at24 mccdate_af13-mccdate_af24;

if dxy=. and dxybf1>. then dxy=dxybf1;if dxy=. and dxybf2>. then dxy=dxybf2;if dxy=. and dxyaf>. then dxy=dxyaf;
if dxm=. and dxmbf1>. then dxm=dxmbf1;if dxm=. and dxmbf2>. then dxm=dxmbf2;if dxm=. and dxmaf>. then dxm=dxmaf;

do n=1 to 48; if a2{n}>0 then a3{n}=1;if a2{n} in (0,.) then a3{n}=0;end;
cc=0;dd=0;ee=0;ff=0;
do i=1 to 12;
 if i=dxm then do;if a1{i}>0 then meddx=1;if a1{i} in (0,.) then meddx=0;  /*medicaid enrollment at diagnosis*/
  do l=(dxm+23) to (dxm+12) by -1;   /*count months of Medicaid enrollment in 12 months before cancer diagnosis*/
   cc=sum(cc,a3{l});
  end;
  do m=(dxm+25) to (dxm+36);    /*count months of Medicaid enrollment in 12 months after diagnosis*/
   dd=sum(dd,a3{m});
  end;
  do h=(dxm+25) to (dxm+30);    /*count months of Medicaid enrollment in 6 months after diagnosis*/
   ee=sum(ee,a3{h});
  end;
  do k=(dxm+25) to (dxm+26);    /*count months of Medicaid enrollment in 2 months after diagnosis*/
   ff=sum(ff,a3{k});
  end;
 end;
end;
if cc>10 then pattern1=1;if 0<cc<=10 then pattern1=0;if cc=0 then pattern1=2;label pattern1='0=discontinuous enroll within 1 year before diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if dd>10 then pattern2=1;if 0<dd<=10 then pattern2=0;if dd=0 then pattern2=2;label pattern2='0=discontinuous enroll withn 1 year after diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if ee=6 then pattern3=1;if 0<ee<6 then pattern3=0;if ee=0 then pattern3=2;label pattern3='0=discontinuous enroll within 6 months after diagnosis,1=continuous enroll (=6 mons),2=no enroll';
if ff=2 then pattern4=1;if 0<ff<2 then pattern4=0;if ff=0 then pattern4=2;label pattern4='0=discontinuous enroll within 2 months after diagnosis,1=continuous enroll (=2 mons),2=no enroll';
if meddx=1 then do;
 if pattern1=1 and pattern2=1 then pattern5=1;if pattern1=1 and pattern2=0 then pattern5=2;if pattern1=0 and pattern2=1 then pattern5=3; 
 if pattern1=0 and pattern2=0 then pattern5=4;if pattern1=2 and pattern2=1 then pattern5=5;if pattern1=2 and pattern2=0 then pattern5=6;
 if pattern1=1 and pattern2=2 then pattern5=7;if pattern1=0 and pattern2=2 then pattern5=8;if pattern1=2 and pattern2=2 then pattern5=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern2=1 then pattern5=9;if pattern1=1 and pattern2=0 then pattern5=10;if pattern1=0 and pattern2=1 then pattern5=11; 
 if pattern1=0 and pattern2=0 then pattern5=12;if pattern1=2 and pattern2=1 then pattern5=13;if pattern1=2 and pattern2=0 then pattern5=14;
 if pattern1=1 and pattern2=2 then pattern5=15;if pattern1=0 and pattern2=2 then pattern5=16;if pattern1=2 and pattern2=2 then pattern5=17;
end;
label pattern5='1=continuous enroll 1yr before,at and 1yr after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern3=1 then pattern6=1;if pattern1=1 and pattern3=0 then pattern6=2;if pattern1=0 and pattern3=1 then pattern6=3; 
 if pattern1=0 and pattern3=0 then pattern6=4;if pattern1=2 and pattern3=1 then pattern6=5;if pattern1=2 and pattern3=0 then pattern6=6;
 if pattern1=1 and pattern3=2 then pattern6=7;if pattern1=0 and pattern3=2 then pattern6=8;if pattern1=2 and pattern3=2 then pattern6=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern3=1 then pattern6=9;if pattern1=1 and pattern3=0 then pattern6=10;if pattern1=0 and pattern3=1 then pattern6=11; 
 if pattern1=0 and pattern3=0 then pattern6=12;if pattern1=2 and pattern3=1 then pattern6=13;if pattern1=2 and pattern3=0 then pattern6=14;
 if pattern1=1 and pattern3=2 then pattern6=15;if pattern1=0 and pattern3=2 then pattern6=16;if pattern1=2 and pattern3=2 then pattern6=17;
end;
label pattern6='1=continuous enroll 1yr before,at and 6months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern4=1 then pattern7=1;if pattern1=1 and pattern4=0 then pattern7=2;if pattern1=0 and pattern4=1 then pattern7=3; 
 if pattern1=0 and pattern4=0 then pattern7=4;if pattern1=2 and pattern4=1 then pattern7=5;if pattern1=2 and pattern4=0 then pattern7=6;
 if pattern1=1 and pattern4=2 then pattern7=7;if pattern1=0 and pattern4=2 then pattern7=8;if pattern1=2 and pattern4=2 then pattern7=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern4=1 then pattern7=9;if pattern1=1 and pattern4=0 then pattern7=10;if pattern1=0 and pattern4=1 then pattern7=11; 
 if pattern1=0 and pattern4=0 then pattern7=12;if pattern1=2 and pattern4=1 then pattern7=13;if pattern1=2 and pattern4=0 then pattern7=14;
 if pattern1=1 and pattern4=2 then pattern7=15;if pattern1=0 and pattern4=2 then pattern7=16;if pattern1=2 and pattern4=2 then pattern7=17;
end;
label pattern7='1=continuous enroll 1yr before,at and 2months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 

gg=0;
do o=(dxm+11) to dxm by -1;   /*count months of Medicaid enrollment in 13-24 months before cancer diagnosis*/
 gg=sum(gg,a3{o});
end;

bb=0;
if pattern5 in (1,2,7,9,10,15) then do;
 do j=1 to 12;
  if j=dxm then do;
   if gg>10 and a3{j}>0 then medenr_mo=a2{j};
   if gg>10 and a3{j}=0 then medenr_mo=a2{j+1};
   if 0<gg<=10 then do;
    if a3{j+11}=0 and a3{j+12}>0 then medenr_mo=a2{j+12};
	if a3{j+11}=0 and a3{j+12}=0 and a3{j+13}>0 then medenr_mo=a2{j+13};
	if a3{j+11}>0 then do;
     do v=(dxm+11) to dxm by -1 until(a3{v}=0);
	  if a2{v}>0 then bb=a2{v};
     end;
	medenr_mo=bb;
    end;
   end;
   if gg=0 and a3{j+13}>0 then medenr_mo=a2{j+13};
   if gg=0 and a3{j+13}=0 then medenr_mo=a2{j+14};
  end;
 end;
end;
if pattern5 in (3,4,8) then do;
 do p=1 to 12;
  if p=dxm then do;
   if a3{p+23}=0 then medenr_mo=a1{p};
   if a3{p+23}>0 then do;
    do q=(dxm+23) to (dxm+12) by -1 until(a3{q}=0);
     if a2{q}>0 then bb=a2{q};
    end;
    medenr_mo=bb;
   end;
  end;
 end;
end;
if pattern5 in (5,6) then do;
 do r=1 to 12;
  if r=dxm then medenr_mo=a1{r};
 end;
end;
if pattern5=11 then do;
 do c=1 to 12;
  if c=dxm then do;
   if a3{c+23}>0 then do;
    do d=(dxm+23) to (dxm+12) by -1 until(a3{d}=0);
     if a2{d}>0 then bb=a2{d};
    end;
    medenr_mo=bb;
   end;
   if a3{c+23}=0 and a3{c+25}>0 then medenr_mo=a2{c+25};
   if a3{c+23}=0 and a3{c+25}=0 then medenr_mo=a2{c+26};
  end;
 end; 
end;
if pattern5=12 then do;
 do e=1 to 12;
  if e=dxm then do;
   if a3{e+23}>0 then do;
    do f=(dxm+23) to (dxm+12) by -1 until(a3{f}=0);
     if a2{f}>0 then bb=a2{f};
    end;
   end;
   if a3{e+23}=0 then do;
    do g=(dxm+25) to (dxm+36) until(a3{g}>0);
     if a2{g}>0 then bb=a2{g};
    end;
   end;
   medenr_mo=bb;
  end;
 end;
end;
if pattern5 in (13,14) then do; 
 do s=1 to 12;
  if s=dxm then do;
   if a3{s+25}=0 then do;
    do t=(dxm+25) to (dxm+36) until(a3{t}>0);
     if a2{t}>0 then bb=a2{t};
    end;
    medenr_mo=bb;
   end;
   if a3{s+25}>0 then medenr_mo=a2{s+25};
  end;
 end;
end;
if pattern5=18 then do;
 do u=1 to 12;
  if u=dxm then medenr_mo=a1{u};
 end;
end;
drop c d e f g h i j k l m n o p q r s t u v bb cc dd ee ff gg;
run;
data mcc1a;set mcc1;mccenr_mo=medenr_mo;mccpattern5=pattern5;mccpattern6=pattern6;mccpattern7=pattern7;
keep patient_id mccenr_mo mccpattern5-mccpattern7;run;proc sort;by patient_id;run;

data mcc2a;set taf13;if mod(_n_,2)=1;
 array a1{48} mccdate_bf1-mccdate_bf24 mccdate_at1-mccdate_at12 mccdate_af1-mccdate_af12;
 array a2{48} mccdate1_bf1-mccdate1_bf24 mccdate1_at1-mccdate1_at12 mccdate1_af1-mccdate1_af12;
 do i=1 to 48;a2{i}=a1{i};end;
 drop mccdate_bf1-mccdate_bf24 mccdate_at1-mccdate_at12 mccdate_af1-mccdate_af12 i;
run;proc sort;by patient_id;run;
data mcc2b;set taf13;if mod(_n_,2)=0;
 array a1{48} mccdate_bf1-mccdate_bf24 mccdate_at1-mccdate_at12 mccdate_af1-mccdate_af12;
 array a2{48} mccdate2_bf1-mccdate2_bf24 mccdate2_at1-mccdate2_at12 mccdate2_af1-mccdate2_af12;
 do i=1 to 48;a2{i}=a1{i};end;
 keep patient_id mccdate2_bf1-mccdate2_bf24 mccdate2_at1-mccdate2_at12 mccdate2_af1-mccdate2_af12;
run;proc sort;by patient_id;run;
data mcc2c;merge mcc2a mcc2b;by patient_id;
 array a1{48} mccdate1_bf1-mccdate1_bf24 mccdate1_at1-mccdate1_at12 mccdate1_af1-mccdate1_af12;
 array a2{48} mccdate2_bf1-mccdate2_bf24 mccdate2_at1-mccdate2_at12 mccdate2_af1-mccdate2_af12;
 array a3{48} mccdate_bf1-mccdate_bf24 mccdate_at1-mccdate_at12 mccdate_af1-mccdate_af12;
 do i=1 to 48;
  if a1{i}=.|a2{i}=. then a3{i}=sum(a1{i},a2{i});
  if a1{i}>=0 and a2{i}>=0 then a3{i}=min(a1{i},a2{i});
 end;
 drop i mccdate1_bf1-mccdate1_bf24 mccdate1_at1-mccdate1_at12 mccdate1_af1-mccdate1_af12 mccdate2_bf1-mccdate2_bf24 mccdate2_at1-mccdate2_at12 mccdate2_af1-mccdate2_af12;
run;
data mcc2d;set mcc2c;
 array a1{12} mccdate_at1-mccdate_at12;
 array a2{48} mccdate_bf1-mccdate_bf24 mccdate_at1-mccdate_at12 mccdate_af1-mccdate_af12;
 array a3{48} mccdate_bf25-mccdate_bf48 mccdate_at21-mccdate_at32 mccdate_af21-mccdate_af32;

if dxy=. and dxybf1>. then dxy=dxybf1;if dxy=. and dxybf2>. then dxy=dxybf2;if dxy=. and dxyaf>. then dxy=dxyaf;
if dxm=. and dxmbf1>. then dxm=dxmbf1;if dxm=. and dxmbf2>. then dxm=dxmbf2;if dxm=. and dxmaf>. then dxm=dxmaf;

do n=1 to 48; if a2{n}>0 then a3{n}=1;if a2{n} in (0,.) then a3{n}=0;end;
cc=0;dd=0;ee=0;ff=0;
do i=1 to 12;
 if i=dxm then do;if a1{i}>0 then meddx=1;if a1{i} in (0,.) then meddx=0;  /*medicaid enrollment at diagnosis*/
  do l=(dxm+23) to (dxm+12) by -1;   /*count months of Medicaid enrollment in 12 months before cancer diagnosis*/
   cc=sum(cc,a3{l});
  end;
  do m=(dxm+25) to (dxm+36);    /*count months of Medicaid enrollment in 12 months after diagnosis*/
   dd=sum(dd,a3{m});
  end;
  do h=(dxm+25) to (dxm+30);    /*count months of Medicaid enrollment in 6 months after diagnosis*/
   ee=sum(ee,a3{h});
  end;
  do k=(dxm+25) to (dxm+26);    /*count months of Medicaid enrollment in 2 months after diagnosis*/
   ff=sum(ff,a3{k});
  end;
 end;
end;
if cc>10 then pattern1=1;if 0<cc<=10 then pattern1=0;if cc=0 then pattern1=2;label pattern1='0=discontinuous enroll within 1 year before diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if dd>10 then pattern2=1;if 0<dd<=10 then pattern2=0;if dd=0 then pattern2=2;label pattern2='0=discontinuous enroll withn 1 year after diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if ee=6 then pattern3=1;if 0<ee<6 then pattern3=0;if ee=0 then pattern3=2;label pattern3='0=discontinuous enroll within 6 months after diagnosis,1=continuous enroll (=6 mons),2=no enroll';
if ff=2 then pattern4=1;if 0<ff<2 then pattern4=0;if ff=0 then pattern4=2;label pattern4='0=discontinuous enroll within 2 months after diagnosis,1=continuous enroll (=2 mons),2=no enroll';
if meddx=1 then do;
 if pattern1=1 and pattern2=1 then pattern5=1;if pattern1=1 and pattern2=0 then pattern5=2;if pattern1=0 and pattern2=1 then pattern5=3; 
 if pattern1=0 and pattern2=0 then pattern5=4;if pattern1=2 and pattern2=1 then pattern5=5;if pattern1=2 and pattern2=0 then pattern5=6;
 if pattern1=1 and pattern2=2 then pattern5=7;if pattern1=0 and pattern2=2 then pattern5=8;if pattern1=2 and pattern2=2 then pattern5=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern2=1 then pattern5=9;if pattern1=1 and pattern2=0 then pattern5=10;if pattern1=0 and pattern2=1 then pattern5=11; 
 if pattern1=0 and pattern2=0 then pattern5=12;if pattern1=2 and pattern2=1 then pattern5=13;if pattern1=2 and pattern2=0 then pattern5=14;
 if pattern1=1 and pattern2=2 then pattern5=15;if pattern1=0 and pattern2=2 then pattern5=16;if pattern1=2 and pattern2=2 then pattern5=17;
end;
label pattern5='1=continuous enroll 1yr before,at and 1yr after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern3=1 then pattern6=1;if pattern1=1 and pattern3=0 then pattern6=2;if pattern1=0 and pattern3=1 then pattern6=3; 
 if pattern1=0 and pattern3=0 then pattern6=4;if pattern1=2 and pattern3=1 then pattern6=5;if pattern1=2 and pattern3=0 then pattern6=6;
 if pattern1=1 and pattern3=2 then pattern6=7;if pattern1=0 and pattern3=2 then pattern6=8;if pattern1=2 and pattern3=2 then pattern6=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern3=1 then pattern6=9;if pattern1=1 and pattern3=0 then pattern6=10;if pattern1=0 and pattern3=1 then pattern6=11; 
 if pattern1=0 and pattern3=0 then pattern6=12;if pattern1=2 and pattern3=1 then pattern6=13;if pattern1=2 and pattern3=0 then pattern6=14;
 if pattern1=1 and pattern3=2 then pattern6=15;if pattern1=0 and pattern3=2 then pattern6=16;if pattern1=2 and pattern3=2 then pattern6=17;
end;
label pattern6='1=continuous enroll 1yr before,at and 6months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern4=1 then pattern7=1;if pattern1=1 and pattern4=0 then pattern7=2;if pattern1=0 and pattern4=1 then pattern7=3; 
 if pattern1=0 and pattern4=0 then pattern7=4;if pattern1=2 and pattern4=1 then pattern7=5;if pattern1=2 and pattern4=0 then pattern7=6;
 if pattern1=1 and pattern4=2 then pattern7=7;if pattern1=0 and pattern4=2 then pattern7=8;if pattern1=2 and pattern4=2 then pattern7=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern4=1 then pattern7=9;if pattern1=1 and pattern4=0 then pattern7=10;if pattern1=0 and pattern4=1 then pattern7=11; 
 if pattern1=0 and pattern4=0 then pattern7=12;if pattern1=2 and pattern4=1 then pattern7=13;if pattern1=2 and pattern4=0 then pattern7=14;
 if pattern1=1 and pattern4=2 then pattern7=15;if pattern1=0 and pattern4=2 then pattern7=16;if pattern1=2 and pattern4=2 then pattern7=17;
end;
label pattern7='1=continuous enroll 1yr before,at and 2months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 

gg=0;
do o=(dxm+11) to dxm by -1;   /*count months of Medicaid enrollment in 13-24 months before cancer diagnosis*/
 gg=sum(gg,a3{o});
end;

bb=0;
if pattern5 in (1,2,7,9,10,15) then do;
 do j=1 to 12;
  if j=dxm then do;
   if gg>10 and a3{j}>0 then medenr_mo=a2{j};
   if gg>10 and a3{j}=0 then medenr_mo=a2{j+1};
   if 0<gg<=10 then do;
    if a3{j+11}=0 and a3{j+12}>0 then medenr_mo=a2{j+12};
	if a3{j+11}=0 and a3{j+12}=0 and a3{j+13}>0 then medenr_mo=a2{j+13};
	if a3{j+11}>0 then do;
     do v=(dxm+11) to dxm by -1 until(a3{v}=0);
	  if a2{v}>0 then bb=a2{v};
     end;
	medenr_mo=bb;
    end;
   end;
   if gg=0 and a3{j+13}>0 then medenr_mo=a2{j+13};
   if gg=0 and a3{j+13}=0 then medenr_mo=a2{j+14};
  end;
 end;
end;
if pattern5 in (3,4,8) then do;
 do p=1 to 12;
  if p=dxm then do;
   if a3{p+23}=0 then medenr_mo=a1{p};
   if a3{p+23}>0 then do;
    do q=(dxm+23) to (dxm+12) by -1 until(a3{q}=0);
     if a2{q}>0 then bb=a2{q};
    end;
    medenr_mo=bb;
   end;
  end;
 end;
end;
if pattern5 in (5,6) then do;
 do r=1 to 12;
  if r=dxm then medenr_mo=a1{r};
 end;
end;
if pattern5=11 then do;
 do c=1 to 12;
  if c=dxm then do;
   if a3{c+23}>0 then do;
    do d=(dxm+23) to (dxm+12) by -1 until(a3{d}=0);
     if a2{d}>0 then bb=a2{d};
    end;
    medenr_mo=bb;
   end;
   if a3{c+23}=0 and a3{c+25}>0 then medenr_mo=a2{c+25};
   if a3{c+23}=0 and a3{c+25}=0 then medenr_mo=a2{c+26};
  end;
 end; 
end;
if pattern5=12 then do;
 do e=1 to 12;
  if e=dxm then do;
   if a3{e+23}>0 then do;
    do f=(dxm+23) to (dxm+12) by -1 until(a3{f}=0);
     if a2{f}>0 then bb=a2{f};
    end;
   end;
   if a3{e+23}=0 then do;
    do g=(dxm+25) to (dxm+36) until(a3{g}>0);
     if a2{g}>0 then bb=a2{g};
    end;
   end;
   medenr_mo=bb;
  end;
 end;
end;
if pattern5 in (13,14) then do; 
 do s=1 to 12;
  if s=dxm then do;
   if a3{s+25}=0 then do;
    do t=(dxm+25) to (dxm+36) until(a3{t}>0);
     if a2{t}>0 then bb=a2{t};
    end;
    medenr_mo=bb;
   end;
   if a3{s+25}>0 then medenr_mo=a2{s+25};
  end;
 end;
end;
if pattern5=18 then do;
 do u=1 to 12;
  if u=dxm then medenr_mo=a1{u};
 end;
end;
drop c d e f g h i j k l m n o p q r s t u v bb cc dd ee ff gg;
run;
data mcc2e;set mcc2d;mccenr_mo=medenr_mo;mccpattern5=pattern5;mccpattern6=pattern6;mccpattern7=pattern7;
keep patient_id mccenr_mo mccpattern5-mccpattern7;run;proc sort;by patient_id;run;

data mcc3a;set taf14;if mod(_n_,3)=1;
 array a1{48} mccdate_bf1-mccdate_bf24 mccdate_at1-mccdate_at12 mccdate_af1-mccdate_af12;
 array a2{48} mccdate1_bf1-mccdate1_bf24 mccdate1_at1-mccdate1_at12 mccdate1_af1-mccdate1_af12;
 do i=1 to 48;a2{i}=a1{i};end;
 drop mccdate_bf1-mccdate_bf24 mccdate_at1-mccdate_at12 mccdate_af1-mccdate_af12 i;
run;proc sort;by patient_id;run;
data mcc3b;set taf14;if mod(_n_,3)=2;
 array a1{48} mccdate_bf1-mccdate_bf24 mccdate_at1-mccdate_at12 mccdate_af1-mccdate_af12;
 array a2{48} mccdate2_bf1-mccdate2_bf24 mccdate2_at1-mccdate2_at12 mccdate2_af1-mccdate2_af12;
 do i=1 to 48;a2{i}=a1{i};end;
 keep patient_id mccdate2_bf1-mccdate2_bf24 mccdate2_at1-mccdate2_at12 mccdate2_af1-mccdate2_af12;
run;proc sort;by patient_id;run;
data mcc3c;set taf14;if mod(_n_,3)=0;
 array a1{48} mccdate_bf1-mccdate_bf24 mccdate_at1-mccdate_at12 mccdate_af1-mccdate_af12;
 array a2{48} mccdate3_bf1-mccdate3_bf24 mccdate3_at1-mccdate3_at12 mccdate3_af1-mccdate3_af12;
 do i=1 to 48;a2{i}=a1{i};end;
 keep patient_id mccdate3_bf1-mccdate3_bf24 mccdate3_at1-mccdate3_at12 mccdate3_af1-mccdate3_af12;
run;proc sort;by patient_id;run;
data mcc3d;merge mcc3a mcc3b mcc3c;by patient_id;
 array a1{48} mccdate1_bf1-mccdate1_bf24 mccdate1_at1-mccdate1_at12 mccdate1_af1-mccdate1_af12;
 array a2{48} mccdate2_bf1-mccdate2_bf24 mccdate2_at1-mccdate2_at12 mccdate2_af1-mccdate2_af12;
 array a3{48} mccdate3_bf1-mccdate3_bf24 mccdate3_at1-mccdate3_at12 mccdate3_af1-mccdate3_af12;
 array a4{48} mccdate_bf1-mccdate_bf24 mccdate_at1-mccdate_at12 mccdate_af1-mccdate_af12;
 do i=1 to 48;
  if (a1{i}=. and a2{i}=. and a3{i}>=0)|(a1{i}=. and a2{i}>=0 and a3{i}=.)|(a1{i}>=0 and a2{i}=. and a3{i}=.) then a4{i}=sum(a1{i},a2{i},a3{i});
  if (a1{i}>=0 and a2{i}>=0)|(a1{i}>=0 and a3{i}>=0)|(a2{i}>=0 and a3{i}>=0) then a4{i}=min(a1{i},a2{i},a3{i});
 end;
 drop i mccdate1_bf1-mccdate1_bf24 mccdate1_at1-mccdate1_at12 mccdate1_af1-mccdate1_af12 
        mccdate2_bf1-mccdate2_bf24 mccdate2_at1-mccdate2_at12 mccdate2_af1-mccdate2_af12
        mccdate3_bf1-mccdate3_bf24 mccdate3_at1-mccdate3_at12 mccdate3_af1-mccdate3_af12;
run;
data mcc3e;set mcc3d;
 array a1{12} mccdate_at1-mccdate_at12;
 array a2{48} mccdate_bf1-mccdate_bf24 mccdate_at1-mccdate_at12 mccdate_af1-mccdate_af12;
 array a3{48} mccdate_bf25-mccdate_bf48 mccdate_at21-mccdate_at32 mccdate_af21-mccdate_af32;

if dxy=. and dxybf1>. then dxy=dxybf1;if dxy=. and dxybf2>. then dxy=dxybf2;if dxy=. and dxyaf>. then dxy=dxyaf;
if dxm=. and dxmbf1>. then dxm=dxmbf1;if dxm=. and dxmbf2>. then dxm=dxmbf2;if dxm=. and dxmaf>. then dxm=dxmaf;

do n=1 to 48; if a2{n}>0 then a3{n}=1;if a2{n} in (0,.) then a3{n}=0;end;
cc=0;dd=0;ee=0;ff=0;
do i=1 to 12;
 if i=dxm then do;if a1{i}>0 then meddx=1;if a1{i} in (0,.) then meddx=0;  /*medicaid enrollment at diagnosis*/
  do l=(dxm+23) to (dxm+12) by -1;   /*count months of Medicaid enrollment in 12 months before cancer diagnosis*/
   cc=sum(cc,a3{l});
  end;
  do m=(dxm+25) to (dxm+36);    /*count months of Medicaid enrollment in 12 months after diagnosis*/
   dd=sum(dd,a3{m});
  end;
  do h=(dxm+25) to (dxm+30);    /*count months of Medicaid enrollment in 6 months after diagnosis*/
   ee=sum(ee,a3{h});
  end;
  do k=(dxm+25) to (dxm+26);    /*count months of Medicaid enrollment in 2 months after diagnosis*/
   ff=sum(ff,a3{k});
  end;
 end;
end;
if cc>10 then pattern1=1;if 0<cc<=10 then pattern1=0;if cc=0 then pattern1=2;label pattern1='0=discontinuous enroll within 1 year before diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if dd>10 then pattern2=1;if 0<dd<=10 then pattern2=0;if dd=0 then pattern2=2;label pattern2='0=discontinuous enroll withn 1 year after diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if ee=6 then pattern3=1;if 0<ee<6 then pattern3=0;if ee=0 then pattern3=2;label pattern3='0=discontinuous enroll within 6 months after diagnosis,1=continuous enroll (=6 mons),2=no enroll';
if ff=2 then pattern4=1;if 0<ff<2 then pattern4=0;if ff=0 then pattern4=2;label pattern4='0=discontinuous enroll within 2 months after diagnosis,1=continuous enroll (=2 mons),2=no enroll';
if meddx=1 then do;
 if pattern1=1 and pattern2=1 then pattern5=1;if pattern1=1 and pattern2=0 then pattern5=2;if pattern1=0 and pattern2=1 then pattern5=3; 
 if pattern1=0 and pattern2=0 then pattern5=4;if pattern1=2 and pattern2=1 then pattern5=5;if pattern1=2 and pattern2=0 then pattern5=6;
 if pattern1=1 and pattern2=2 then pattern5=7;if pattern1=0 and pattern2=2 then pattern5=8;if pattern1=2 and pattern2=2 then pattern5=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern2=1 then pattern5=9;if pattern1=1 and pattern2=0 then pattern5=10;if pattern1=0 and pattern2=1 then pattern5=11; 
 if pattern1=0 and pattern2=0 then pattern5=12;if pattern1=2 and pattern2=1 then pattern5=13;if pattern1=2 and pattern2=0 then pattern5=14;
 if pattern1=1 and pattern2=2 then pattern5=15;if pattern1=0 and pattern2=2 then pattern5=16;if pattern1=2 and pattern2=2 then pattern5=17;
end;
label pattern5='1=continuous enroll 1yr before,at and 1yr after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern3=1 then pattern6=1;if pattern1=1 and pattern3=0 then pattern6=2;if pattern1=0 and pattern3=1 then pattern6=3; 
 if pattern1=0 and pattern3=0 then pattern6=4;if pattern1=2 and pattern3=1 then pattern6=5;if pattern1=2 and pattern3=0 then pattern6=6;
 if pattern1=1 and pattern3=2 then pattern6=7;if pattern1=0 and pattern3=2 then pattern6=8;if pattern1=2 and pattern3=2 then pattern6=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern3=1 then pattern6=9;if pattern1=1 and pattern3=0 then pattern6=10;if pattern1=0 and pattern3=1 then pattern6=11; 
 if pattern1=0 and pattern3=0 then pattern6=12;if pattern1=2 and pattern3=1 then pattern6=13;if pattern1=2 and pattern3=0 then pattern6=14;
 if pattern1=1 and pattern3=2 then pattern6=15;if pattern1=0 and pattern3=2 then pattern6=16;if pattern1=2 and pattern3=2 then pattern6=17;
end;
label pattern6='1=continuous enroll 1yr before,at and 6months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern4=1 then pattern7=1;if pattern1=1 and pattern4=0 then pattern7=2;if pattern1=0 and pattern4=1 then pattern7=3; 
 if pattern1=0 and pattern4=0 then pattern7=4;if pattern1=2 and pattern4=1 then pattern7=5;if pattern1=2 and pattern4=0 then pattern7=6;
 if pattern1=1 and pattern4=2 then pattern7=7;if pattern1=0 and pattern4=2 then pattern7=8;if pattern1=2 and pattern4=2 then pattern7=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern4=1 then pattern7=9;if pattern1=1 and pattern4=0 then pattern7=10;if pattern1=0 and pattern4=1 then pattern7=11; 
 if pattern1=0 and pattern4=0 then pattern7=12;if pattern1=2 and pattern4=1 then pattern7=13;if pattern1=2 and pattern4=0 then pattern7=14;
 if pattern1=1 and pattern4=2 then pattern7=15;if pattern1=0 and pattern4=2 then pattern7=16;if pattern1=2 and pattern4=2 then pattern7=17;
end;
label pattern7='1=continuous enroll 1yr before,at and 2months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 

gg=0;
do o=(dxm+11) to dxm by -1;   /*count months of Medicaid enrollment in 13-24 months before cancer diagnosis*/
 gg=sum(gg,a3{o});
end;

bb=0;
if pattern5 in (1,2,7,9,10,15) then do;
 do j=1 to 12;
  if j=dxm then do;
   if gg>10 and a3{j}>0 then medenr_mo=a2{j};
   if gg>10 and a3{j}=0 then medenr_mo=a2{j+1};
   if 0<gg<=10 then do;
    if a3{j+11}=0 and a3{j+12}>0 then medenr_mo=a2{j+12};
	if a3{j+11}=0 and a3{j+12}=0 and a3{j+13}>0 then medenr_mo=a2{j+13};
	if a3{j+11}>0 then do;
     do v=(dxm+11) to dxm by -1 until(a3{v}=0);
	  if a2{v}>0 then bb=a2{v};
     end;
	medenr_mo=bb;
    end;
   end;
   if gg=0 and a3{j+13}>0 then medenr_mo=a2{j+13};
   if gg=0 and a3{j+13}=0 then medenr_mo=a2{j+14};
  end;
 end;
end;
if pattern5 in (3,4,8) then do;
 do p=1 to 12;
  if p=dxm then do;
   if a3{p+23}=0 then medenr_mo=a1{p};
   if a3{p+23}>0 then do;
    do q=(dxm+23) to (dxm+12) by -1 until(a3{q}=0);
     if a2{q}>0 then bb=a2{q};
    end;
    medenr_mo=bb;
   end;
  end;
 end;
end;
if pattern5 in (5,6) then do;
 do r=1 to 12;
  if r=dxm then medenr_mo=a1{r};
 end;
end;
if pattern5=11 then do;
 do c=1 to 12;
  if c=dxm then do;
   if a3{c+23}>0 then do;
    do d=(dxm+23) to (dxm+12) by -1 until(a3{d}=0);
     if a2{d}>0 then bb=a2{d};
    end;
    medenr_mo=bb;
   end;
   if a3{c+23}=0 and a3{c+25}>0 then medenr_mo=a2{c+25};
   if a3{c+23}=0 and a3{c+25}=0 then medenr_mo=a2{c+26};
  end;
 end; 
end;
if pattern5=12 then do;
 do e=1 to 12;
  if e=dxm then do;
   if a3{e+23}>0 then do;
    do f=(dxm+23) to (dxm+12) by -1 until(a3{f}=0);
     if a2{f}>0 then bb=a2{f};
    end;
   end;
   if a3{e+23}=0 then do;
    do g=(dxm+25) to (dxm+36) until(a3{g}>0);
     if a2{g}>0 then bb=a2{g};
    end;
   end;
   medenr_mo=bb;
  end;
 end;
end;
if pattern5 in (13,14) then do; 
 do s=1 to 12;
  if s=dxm then do;
   if a3{s+25}=0 then do;
    do t=(dxm+25) to (dxm+36) until(a3{t}>0);
     if a2{t}>0 then bb=a2{t};
    end;
    medenr_mo=bb;
   end;
   if a3{s+25}>0 then medenr_mo=a2{s+25};
  end;
 end;
end;
if pattern5=18 then do;
 do u=1 to 12;
  if u=dxm then medenr_mo=a1{u};
 end;
end;
drop c d e f g h i j k l m n o p q r s t u v bb cc dd ee ff gg;
run;
data mcc3f;set mcc3e;mccenr_mo=medenr_mo;mccpattern5=pattern5;mccpattern6=pattern6;mccpattern7=pattern7;
keep patient_id mccenr_mo mccpattern5-mccpattern7;run;proc sort;by patient_id;run;

data mcc4a;set taf4a;
 array a1{12} mccdate13-mccdate24;
 array a2{12} at1-at12;
 array a3{12} bf1-bf12;
 array a4{12} bf13-bf24;
 array a5{12} af1-af12;
 array a6{48} bf1-bf24 at1-at12 af1-af12;
 array a7{48} ff1-ff48;
cc=0;dd=0;ee=0;ff=0;
do i=1 to 12;
 if diff=0 then a2{i}=a1{i};
 if diff=1 then a4{i}=a1{i};
 if diff=2 then a3{i}=a1{i};
 if diff=-1 then a5{i}=a1{i};
end;
do n=1 to 48;
 if a6{n}>0 then a7{n}=1;if a6{n} in (0,.) then a7{n}=0;
end;
do h=1 to 12;
 if h=dxm and a2{h}>0 then meddx=1;
end;
do k=(dxm+23) to (dxm+12) by -1;  /*count months of Medicaid enrollment in 12 months before cancer diagnosis*/
 cc=sum(cc,a7{k});
end;
do l=(dxm+25) to (dxm+36);  /*count months of Medicaid enrollment in 12 months after diagnosis*/
 dd=sum(dd,a7{l});
end;
do m=(dxm+25) to (dxm+30); /*count months of Medicaid enrollment in 6 months after diagnosis*/
 ee=sum(ee,a7{m});
end;
do o=(dxm+25) to (dxm+26); /*count months of Medicaid enrollment in 2 months after diagnosis*/
 ff=sum(ff,a7{o});
end;

if cc>10 then pattern1=1;if 0<cc<=10 then pattern1=0;if cc=0 then pattern1=2;label pattern1='0=discontinuous enroll within 1 year before diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if dd>10 then pattern2=1;if 0<dd<=10 then pattern2=0;if dd=0 then pattern2=2;label pattern2='0=discontinuous enroll withn 1 year after diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if ee=6 then pattern3=1;if 0<ee<6 then pattern3=0;if ee=0 then pattern3=2;label pattern3='0=discontinuous enroll within 6 months after diagnosis,1=continuous enroll (=6 mons),2=no enroll';
if ff=2 then pattern4=1;if 0<ff<2 then pattern4=0;if ff=0 then pattern4=2;label pattern4='0=discontinuous enroll within 2 months after diagnosis,1=continuous enroll (=2 mons),2=no enroll';
if meddx=1 then do;
 if pattern1=1 and pattern2=1 then pattern5=1;if pattern1=1 and pattern2=0 then pattern5=2;if pattern1=0 and pattern2=1 then pattern5=3; 
 if pattern1=0 and pattern2=0 then pattern5=4;if pattern1=2 and pattern2=1 then pattern5=5;if pattern1=2 and pattern2=0 then pattern5=6;
 if pattern1=1 and pattern2=2 then pattern5=7;if pattern1=0 and pattern2=2 then pattern5=8;if pattern1=2 and pattern2=2 then pattern5=18;
end;
if meddx in (0,.) then do;
 if pattern1=1 and pattern2=1 then pattern5=9;if pattern1=1 and pattern2=0 then pattern5=10;if pattern1=0 and pattern2=1 then pattern5=11; 
 if pattern1=0 and pattern2=0 then pattern5=12;if pattern1=2 and pattern2=1 then pattern5=13;if pattern1=2 and pattern2=0 then pattern5=14;
 if pattern1=1 and pattern2=2 then pattern5=15;if pattern1=0 and pattern2=2 then pattern5=16;if pattern1=2 and pattern2=2 then pattern5=17;
end;
label pattern5='1=continuous enroll 1yr before,at and 1yr after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern3=1 then pattern6=1;if pattern1=1 and pattern3=0 then pattern6=2;if pattern1=0 and pattern3=1 then pattern6=3; 
 if pattern1=0 and pattern3=0 then pattern6=4;if pattern1=2 and pattern3=1 then pattern6=5;if pattern1=2 and pattern3=0 then pattern6=6;
 if pattern1=1 and pattern3=2 then pattern6=7;if pattern1=0 and pattern3=2 then pattern6=8;if pattern1=2 and pattern3=2 then pattern6=18;
end;
if meddx in (0,.) then do;
 if pattern1=1 and pattern3=1 then pattern6=9;if pattern1=1 and pattern3=0 then pattern6=10;if pattern1=0 and pattern3=1 then pattern6=11; 
 if pattern1=0 and pattern3=0 then pattern6=12;if pattern1=2 and pattern3=1 then pattern6=13;if pattern1=2 and pattern3=0 then pattern6=14;
 if pattern1=1 and pattern3=2 then pattern6=15;if pattern1=0 and pattern3=2 then pattern6=16;if pattern1=2 and pattern3=2 then pattern6=17;
end;
label pattern6='1=continuous enroll 1yr before,at and 6months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern4=1 then pattern7=1;if pattern1=1 and pattern4=0 then pattern7=2;if pattern1=0 and pattern4=1 then pattern7=3; 
 if pattern1=0 and pattern4=0 then pattern7=4;if pattern1=2 and pattern4=1 then pattern7=5;if pattern1=2 and pattern4=0 then pattern7=6;
 if pattern1=1 and pattern4=2 then pattern7=7;if pattern1=0 and pattern4=2 then pattern7=8;if pattern1=2 and pattern4=2 then pattern7=18;
end;
if meddx in (0,.) then do;
 if pattern1=1 and pattern4=1 then pattern7=9;if pattern1=1 and pattern4=0 then pattern7=10;if pattern1=0 and pattern4=1 then pattern7=11; 
 if pattern1=0 and pattern4=0 then pattern7=12;if pattern1=2 and pattern4=1 then pattern7=13;if pattern1=2 and pattern4=0 then pattern7=14;
 if pattern1=1 and pattern4=2 then pattern7=15;if pattern1=0 and pattern4=2 then pattern7=16;if pattern1=2 and pattern4=2 then pattern7=17;
end;
label pattern7='1=continuous enroll 1yr before,at and 2months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 

gg=0;
do o=(dxm+11) to dxm by -1;   /*count months of Medicaid enrollment in 13-24 months before cancer diagnosis*/
 gg=sum(gg,a7{o});
end;

bb=0;
if pattern5 in (1,2,7,9,10,15) then do;
 do j=1 to 12;
  if j=dxm then do;
   if gg>10 and a7{j}>0 then medenr_mo=a6{j};
   if gg>10 and a7{j}=0 then medenr_mo=a6{j+1};
   if 0<gg<=10 then do;
    if a7{j+11}=0 and a7{j+12}>0 then medenr_mo=a6{j+12};
	if a7{j+11}=0 and a7{j+12}=0 and a7{j+13}>0 then medenr_mo=a6{j+13};
	if a7{j+11}>0 then do;
     do v=(dxm+11) to dxm by -1 until(a7{v}=0);
	  if a6{v}>0 then bb=a6{v};
     end;
	medenr_mo=bb;
    end;
   end;
   if gg=0 and a7{j+13}>0 then medenr_mo=a6{j+13};
   if gg=0 and a7{j+13}=0 then medenr_mo=a6{j+14};
  end;
 end;
end;
if pattern5 in (3,4,8) then do;
 do p=1 to 12;
  if p=dxm then do;
   if a7{p+23}=0 then medenr_mo=a2{p};
   if a7{p+23}>0 then do;
    do q=(dxm+23) to (dxm+12) by -1 until(a7{q}=0);
     if a6{q}>0 then bb=a6{q};
    end;
    medenr_mo=bb;
   end;
  end;
 end;
end;
if pattern5 in (5,6) then do;
 do r=1 to 12;
  if r=dxm then medenr_mo=a2{r};
 end;
end;
if pattern5=11 then do;
 do c=1 to 12;
  if c=dxm then do;
   if a7{c+23}>0 then do;
    do d=(dxm+23) to (dxm+12) by -1 until(a7{d}=0);
     if a6{d}>0 then bb=a6{d};
    end;
    medenr_mo=bb;
   end;
   if a7{c+23}=0 and a7{c+25}>0 then medenr_mo=a6{c+25};
   if a7{c+23}=0 and a7{c+25}=0 then medenr_mo=a6{c+26};
  end;
 end; 
end;
if pattern5=12 then do;
 do e=1 to 12;
  if e=dxm then do;
   if a7{e+23}>0 then do;
    do f=(dxm+23) to (dxm+12) by -1 until(a7{f}=0);
     if a6{f}>0 then bb=a6{f};
    end;
   end;
   if a7{e+23}=0 then do;
    do g=(dxm+25) to (dxm+36) until(a7{g}>0);
     if a6{g}>0 then bb=a6{g};
    end;
   end;
   medenr_mo=bb;
  end;
 end;
end;
if pattern5 in (13,14) then do; 
 do s=1 to 12;
  if s=dxm then do;
   if a7{s+25}=0 then do;
    do t=(dxm+25) to (dxm+36) until(a7{t}>0);
     if a6{t}>0 then bb=a6{t};
    end;
    medenr_mo=bb;
   end;
   if a7{s+25}>0 then medenr_mo=a6{s+25};
  end;
 end;
end;
if pattern5=18 then do;
 do u=1 to 12;
  if u=dxm then medenr_mo=a2{u};
 end;
end;
drop c d e f g h i j k l m n o p q r s t u v bb cc dd ee ff gg;
run;
data mcc4b;set mcc4a;mccenr_mo=medenr_mo;mccpattern5=pattern5;mccpattern6=pattern6;mccpattern7=pattern7;
keep patient_id mccenr_mo mccpattern5-mccpattern7;run;proc sort;by patient_id;run;

data mcc5;set mcc1a mcc2e mcc3f mcc4b;run;proc sort;by patient_id;run;  /*medicaid managed care enrollment timeline*/

/*FFS enrollment timeline*/
data ffs1;set taf10;
 array a1{12} ffsdate_at1-ffsdate_at12;
 array a2{48} ffsdate_bf1-ffsdate_bf24 ffsdate_at1-ffsdate_at12 ffsdate_af1-ffsdate_af12;
 array a3{48} ffsdate_bf25-ffsdate_bf48 ffsdate_at13-ffsdate_at24 ffsdate_af13-ffsdate_af24;

if dxy=. and dxybf1>. then dxy=dxybf1;if dxy=. and dxybf2>. then dxy=dxybf2;if dxy=. and dxyaf>. then dxy=dxyaf;
if dxm=. and dxmbf1>. then dxm=dxmbf1;if dxm=. and dxmbf2>. then dxm=dxmbf2;if dxm=. and dxmaf>. then dxm=dxmaf;

do n=1 to 48; if a2{n}>0 then a3{n}=1;if a2{n} in (0,.) then a3{n}=0;end;
cc=0;dd=0;ee=0;ff=0;
do i=1 to 12;
 if i=dxm then do;if a1{i}>0 then meddx=1;if a1{i} in (0,.) then meddx=0;  /*medicaid enrollment at diagnosis*/
  do l=(dxm+23) to (dxm+12) by -1;   /*count months of Medicaid enrollment in 12 months before cancer diagnosis*/
   cc=sum(cc,a3{l});
  end;
  do m=(dxm+25) to (dxm+36);    /*count months of Medicaid enrollment in 12 months after diagnosis*/
   dd=sum(dd,a3{m});
  end;
  do h=(dxm+25) to (dxm+30);    /*count months of Medicaid enrollment in 6 months after diagnosis*/
   ee=sum(ee,a3{h});
  end;
  do k=(dxm+25) to (dxm+26);    /*count months of Medicaid enrollment in 2 months after diagnosis*/
   ff=sum(ff,a3{k});
  end;
 end;
end;
if cc>10 then pattern1=1;if 0<cc<=10 then pattern1=0;if cc=0 then pattern1=2;label pattern1='0=discontinuous enroll within 1 year before diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if dd>10 then pattern2=1;if 0<dd<=10 then pattern2=0;if dd=0 then pattern2=2;label pattern2='0=discontinuous enroll withn 1 year after diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if ee=6 then pattern3=1;if 0<ee<6 then pattern3=0;if ee=0 then pattern3=2;label pattern3='0=discontinuous enroll within 6 months after diagnosis,1=continuous enroll (=6 mons),2=no enroll';
if ff=2 then pattern4=1;if 0<ff<2 then pattern4=0;if ff=0 then pattern4=2;label pattern4='0=discontinuous enroll within 2 months after diagnosis,1=continuous enroll (=2 mons),2=no enroll';
if meddx=1 then do;
 if pattern1=1 and pattern2=1 then pattern5=1;if pattern1=1 and pattern2=0 then pattern5=2;if pattern1=0 and pattern2=1 then pattern5=3; 
 if pattern1=0 and pattern2=0 then pattern5=4;if pattern1=2 and pattern2=1 then pattern5=5;if pattern1=2 and pattern2=0 then pattern5=6;
 if pattern1=1 and pattern2=2 then pattern5=7;if pattern1=0 and pattern2=2 then pattern5=8;if pattern1=2 and pattern2=2 then pattern5=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern2=1 then pattern5=9;if pattern1=1 and pattern2=0 then pattern5=10;if pattern1=0 and pattern2=1 then pattern5=11; 
 if pattern1=0 and pattern2=0 then pattern5=12;if pattern1=2 and pattern2=1 then pattern5=13;if pattern1=2 and pattern2=0 then pattern5=14;
 if pattern1=1 and pattern2=2 then pattern5=15;if pattern1=0 and pattern2=2 then pattern5=16;if pattern1=2 and pattern2=2 then pattern5=17;
end;
label pattern5='1=continuous enroll 1yr before,at and 1yr after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern3=1 then pattern6=1;if pattern1=1 and pattern3=0 then pattern6=2;if pattern1=0 and pattern3=1 then pattern6=3; 
 if pattern1=0 and pattern3=0 then pattern6=4;if pattern1=2 and pattern3=1 then pattern6=5;if pattern1=2 and pattern3=0 then pattern6=6;
 if pattern1=1 and pattern3=2 then pattern6=7;if pattern1=0 and pattern3=2 then pattern6=8;if pattern1=2 and pattern3=2 then pattern6=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern3=1 then pattern6=9;if pattern1=1 and pattern3=0 then pattern6=10;if pattern1=0 and pattern3=1 then pattern6=11; 
 if pattern1=0 and pattern3=0 then pattern6=12;if pattern1=2 and pattern3=1 then pattern6=13;if pattern1=2 and pattern3=0 then pattern6=14;
 if pattern1=1 and pattern3=2 then pattern6=15;if pattern1=0 and pattern3=2 then pattern6=16;if pattern1=2 and pattern3=2 then pattern6=17;
end;
label pattern6='1=continuous enroll 1yr before,at and 6months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern4=1 then pattern7=1;if pattern1=1 and pattern4=0 then pattern7=2;if pattern1=0 and pattern4=1 then pattern7=3; 
 if pattern1=0 and pattern4=0 then pattern7=4;if pattern1=2 and pattern4=1 then pattern7=5;if pattern1=2 and pattern4=0 then pattern7=6;
 if pattern1=1 and pattern4=2 then pattern7=7;if pattern1=0 and pattern4=2 then pattern7=8;if pattern1=2 and pattern4=2 then pattern7=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern4=1 then pattern7=9;if pattern1=1 and pattern4=0 then pattern7=10;if pattern1=0 and pattern4=1 then pattern7=11; 
 if pattern1=0 and pattern4=0 then pattern7=12;if pattern1=2 and pattern4=1 then pattern7=13;if pattern1=2 and pattern4=0 then pattern7=14;
 if pattern1=1 and pattern4=2 then pattern7=15;if pattern1=0 and pattern4=2 then pattern7=16;if pattern1=2 and pattern4=2 then pattern7=17;
end;
label pattern7='1=continuous enroll 1yr before,at and 2months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 

gg=0;
do o=(dxm+11) to dxm by -1;   /*count months of Medicaid enrollment in 13-24 months before cancer diagnosis*/
 gg=sum(gg,a3{o});
end;

bb=0;
if pattern5 in (1,2,7,9,10,15) then do;
 do j=1 to 12;
  if j=dxm then do;
   if gg>10 and a3{j}>0 then medenr_mo=a2{j};
   if gg>10 and a3{j}=0 then medenr_mo=a2{j+1};
   if 0<gg<=10 then do;
    if a3{j+11}=0 and a3{j+12}>0 then medenr_mo=a2{j+12};
	if a3{j+11}=0 and a3{j+12}=0 and a3{j+13}>0 then medenr_mo=a2{j+13};
	if a3{j+11}>0 then do;
     do v=(dxm+11) to dxm by -1 until(a3{v}=0);
	  if a2{v}>0 then bb=a2{v};
     end;
	medenr_mo=bb;
    end;
   end;
   if gg=0 and a3{j+13}>0 then medenr_mo=a2{j+13};
   if gg=0 and a3{j+13}=0 then medenr_mo=a2{j+14};
  end;
 end;
end;
if pattern5 in (3,4,8) then do;
 do p=1 to 12;
  if p=dxm then do;
   if a3{p+23}=0 then medenr_mo=a1{p};
   if a3{p+23}>0 then do;
    do q=(dxm+23) to (dxm+12) by -1 until(a3{q}=0);
     if a2{q}>0 then bb=a2{q};
    end;
    medenr_mo=bb;
   end;
  end;
 end;
end;
if pattern5 in (5,6) then do;
 do r=1 to 12;
  if r=dxm then medenr_mo=a1{r};
 end;
end;
if pattern5=11 then do;
 do c=1 to 12;
  if c=dxm then do;
   if a3{c+23}>0 then do;
    do d=(dxm+23) to (dxm+12) by -1 until(a3{d}=0);
     if a2{d}>0 then bb=a2{d};
    end;
    medenr_mo=bb;
   end;
   if a3{c+23}=0 and a3{c+25}>0 then medenr_mo=a2{c+25};
   if a3{c+23}=0 and a3{c+25}=0 then medenr_mo=a2{c+26};
  end;
 end; 
end;
if pattern5=12 then do;
 do e=1 to 12;
  if e=dxm then do;
   if a3{e+23}>0 then do;
    do f=(dxm+23) to (dxm+12) by -1 until(a3{f}=0);
     if a2{f}>0 then bb=a2{f};
    end;
   end;
   if a3{e+23}=0 then do;
    do g=(dxm+25) to (dxm+36) until(a3{g}>0);
     if a2{g}>0 then bb=a2{g};
    end;
   end;
   medenr_mo=bb;
  end;
 end;
end;
if pattern5 in (13,14) then do; 
 do s=1 to 12;
  if s=dxm then do;
   if a3{s+25}=0 then do;
    do t=(dxm+25) to (dxm+36) until(a3{t}>0);
     if a2{t}>0 then bb=a2{t};
    end;
    medenr_mo=bb;
   end;
   if a3{s+25}>0 then medenr_mo=a2{s+25};
  end;
 end;
end;
if pattern5=18 then do;
 do u=1 to 12;
  if u=dxm then medenr_mo=a1{u};
 end;
end;
drop c d e f g h i j k l m n o p q r s t u v bb cc dd ee ff gg;
run;
data ffs1a;set ffs1;ffsenr_mo=medenr_mo;ffspattern5=pattern5;ffspattern6=pattern6;ffspattern7=pattern7;
keep patient_id ffsenr_mo ffspattern5-ffspattern7;run;proc sort;by patient_id;run;

data ffs2a;set taf13;if mod(_n_,2)=1;
 array a1{48} ffsdate_bf1-ffsdate_bf24 ffsdate_at1-ffsdate_at12 ffsdate_af1-ffsdate_af12;
 array a2{48} ffsdate1_bf1-ffsdate1_bf24 ffsdate1_at1-ffsdate1_at12 ffsdate1_af1-ffsdate1_af12;
 do i=1 to 48;a2{i}=a1{i};end;
 drop ffsdate_bf1-ffsdate_bf24 ffsdate_at1-ffsdate_at12 ffsdate_af1-ffsdate_af12 i;
run;proc sort;by patient_id;run;
data ffs2b;set taf13;if mod(_n_,2)=0;
 array a1{48} ffsdate_bf1-ffsdate_bf24 ffsdate_at1-ffsdate_at12 ffsdate_af1-ffsdate_af12;
 array a2{48} ffsdate2_bf1-ffsdate2_bf24 ffsdate2_at1-ffsdate2_at12 ffsdate2_af1-ffsdate2_af12;
 do i=1 to 48;a2{i}=a1{i};end;
 keep patient_id ffsdate2_bf1-ffsdate2_bf24 ffsdate2_at1-ffsdate2_at12 ffsdate2_af1-ffsdate2_af12;
run;proc sort;by patient_id;run;
data ffs2c;merge ffs2a ffs2b;by patient_id;
 array a1{48} ffsdate1_bf1-ffsdate1_bf24 ffsdate1_at1-ffsdate1_at12 ffsdate1_af1-ffsdate1_af12;
 array a2{48} ffsdate2_bf1-ffsdate2_bf24 ffsdate2_at1-ffsdate2_at12 ffsdate2_af1-ffsdate2_af12;
 array a3{48} ffsdate_bf1-ffsdate_bf24 ffsdate_at1-ffsdate_at12 ffsdate_af1-ffsdate_af12;
 do i=1 to 48;
  if a1{i}=.|a2{i}=. then a3{i}=sum(a1{i},a2{i});
  if a1{i}>=0 and a2{i}>=0 then a3{i}=min(a1{i},a2{i});
 end;
 drop i ffsdate1_bf1-ffsdate1_bf24 ffsdate1_at1-ffsdate1_at12 ffsdate1_af1-ffsdate1_af12 ffsdate2_bf1-ffsdate2_bf24 ffsdate2_at1-ffsdate2_at12 ffsdate2_af1-ffsdate2_af12;
run;
data ffs2d;set ffs2c;
 array a1{12} ffsdate_at1-ffsdate_at12;
 array a2{48} ffsdate_bf1-ffsdate_bf24 ffsdate_at1-ffsdate_at12 ffsdate_af1-ffsdate_af12;
 array a3{48} ffsdate_bf25-ffsdate_bf48 ffsdate_at21-ffsdate_at32 ffsdate_af21-ffsdate_af32;

if dxy=. and dxybf1>. then dxy=dxybf1;if dxy=. and dxybf2>. then dxy=dxybf2;if dxy=. and dxyaf>. then dxy=dxyaf;
if dxm=. and dxmbf1>. then dxm=dxmbf1;if dxm=. and dxmbf2>. then dxm=dxmbf2;if dxm=. and dxmaf>. then dxm=dxmaf;

do n=1 to 48; if a2{n}>0 then a3{n}=1;if a2{n} in (0,.) then a3{n}=0;end;
cc=0;dd=0;ee=0;ff=0;
do i=1 to 12;
 if i=dxm then do;if a1{i}>0 then meddx=1;if a1{i} in (0,.) then meddx=0;  /*medicaid enrollment at diagnosis*/
  do l=(dxm+23) to (dxm+12) by -1;   /*count months of Medicaid enrollment in 12 months before cancer diagnosis*/
   cc=sum(cc,a3{l});
  end;
  do m=(dxm+25) to (dxm+36);    /*count months of Medicaid enrollment in 12 months after diagnosis*/
   dd=sum(dd,a3{m});
  end;
  do h=(dxm+25) to (dxm+30);    /*count months of Medicaid enrollment in 6 months after diagnosis*/
   ee=sum(ee,a3{h});
  end;
  do k=(dxm+25) to (dxm+26);    /*count months of Medicaid enrollment in 2 months after diagnosis*/
   ff=sum(ff,a3{k});
  end;
 end;
end;
if cc>10 then pattern1=1;if 0<cc<=10 then pattern1=0;if cc=0 then pattern1=2;label pattern1='0=discontinuous enroll within 1 year before diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if dd>10 then pattern2=1;if 0<dd<=10 then pattern2=0;if dd=0 then pattern2=2;label pattern2='0=discontinuous enroll withn 1 year after diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if ee=6 then pattern3=1;if 0<ee<6 then pattern3=0;if ee=0 then pattern3=2;label pattern3='0=discontinuous enroll within 6 months after diagnosis,1=continuous enroll (=6 mons),2=no enroll';
if ff=2 then pattern4=1;if 0<ff<2 then pattern4=0;if ff=0 then pattern4=2;label pattern4='0=discontinuous enroll within 2 months after diagnosis,1=continuous enroll (=2 mons),2=no enroll';
if meddx=1 then do;
 if pattern1=1 and pattern2=1 then pattern5=1;if pattern1=1 and pattern2=0 then pattern5=2;if pattern1=0 and pattern2=1 then pattern5=3; 
 if pattern1=0 and pattern2=0 then pattern5=4;if pattern1=2 and pattern2=1 then pattern5=5;if pattern1=2 and pattern2=0 then pattern5=6;
 if pattern1=1 and pattern2=2 then pattern5=7;if pattern1=0 and pattern2=2 then pattern5=8;if pattern1=2 and pattern2=2 then pattern5=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern2=1 then pattern5=9;if pattern1=1 and pattern2=0 then pattern5=10;if pattern1=0 and pattern2=1 then pattern5=11; 
 if pattern1=0 and pattern2=0 then pattern5=12;if pattern1=2 and pattern2=1 then pattern5=13;if pattern1=2 and pattern2=0 then pattern5=14;
 if pattern1=1 and pattern2=2 then pattern5=15;if pattern1=0 and pattern2=2 then pattern5=16;if pattern1=2 and pattern2=2 then pattern5=17;
end;
label pattern5='1=continuous enroll 1yr before,at and 1yr after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern3=1 then pattern6=1;if pattern1=1 and pattern3=0 then pattern6=2;if pattern1=0 and pattern3=1 then pattern6=3; 
 if pattern1=0 and pattern3=0 then pattern6=4;if pattern1=2 and pattern3=1 then pattern6=5;if pattern1=2 and pattern3=0 then pattern6=6;
 if pattern1=1 and pattern3=2 then pattern6=7;if pattern1=0 and pattern3=2 then pattern6=8;if pattern1=2 and pattern3=2 then pattern6=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern3=1 then pattern6=9;if pattern1=1 and pattern3=0 then pattern6=10;if pattern1=0 and pattern3=1 then pattern6=11; 
 if pattern1=0 and pattern3=0 then pattern6=12;if pattern1=2 and pattern3=1 then pattern6=13;if pattern1=2 and pattern3=0 then pattern6=14;
 if pattern1=1 and pattern3=2 then pattern6=15;if pattern1=0 and pattern3=2 then pattern6=16;if pattern1=2 and pattern3=2 then pattern6=17;
end;
label pattern6='1=continuous enroll 1yr before,at and 6months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern4=1 then pattern7=1;if pattern1=1 and pattern4=0 then pattern7=2;if pattern1=0 and pattern4=1 then pattern7=3; 
 if pattern1=0 and pattern4=0 then pattern7=4;if pattern1=2 and pattern4=1 then pattern7=5;if pattern1=2 and pattern4=0 then pattern7=6;
 if pattern1=1 and pattern4=2 then pattern7=7;if pattern1=0 and pattern4=2 then pattern7=8;if pattern1=2 and pattern4=2 then pattern7=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern4=1 then pattern7=9;if pattern1=1 and pattern4=0 then pattern7=10;if pattern1=0 and pattern4=1 then pattern7=11; 
 if pattern1=0 and pattern4=0 then pattern7=12;if pattern1=2 and pattern4=1 then pattern7=13;if pattern1=2 and pattern4=0 then pattern7=14;
 if pattern1=1 and pattern4=2 then pattern7=15;if pattern1=0 and pattern4=2 then pattern7=16;if pattern1=2 and pattern4=2 then pattern7=17;
end;
label pattern7='1=continuous enroll 1yr before,at and 2months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 

gg=0;
do o=(dxm+11) to dxm by -1;   /*count months of Medicaid enrollment in 13-24 months before cancer diagnosis*/
 gg=sum(gg,a3{o});
end;

bb=0;
if pattern5 in (1,2,7,9,10,15) then do;
 do j=1 to 12;
  if j=dxm then do;
   if gg>10 and a3{j}>0 then medenr_mo=a2{j};
   if gg>10 and a3{j}=0 then medenr_mo=a2{j+1};
   if 0<gg<=10 then do;
    if a3{j+11}=0 and a3{j+12}>0 then medenr_mo=a2{j+12};
	if a3{j+11}=0 and a3{j+12}=0 and a3{j+13}>0 then medenr_mo=a2{j+13};
	if a3{j+11}>0 then do;
     do v=(dxm+11) to dxm by -1 until(a3{v}=0);
	  if a2{v}>0 then bb=a2{v};
     end;
	medenr_mo=bb;
    end;
   end;
   if gg=0 and a3{j+13}>0 then medenr_mo=a2{j+13};
   if gg=0 and a3{j+13}=0 then medenr_mo=a2{j+14};
  end;
 end;
end;
if pattern5 in (3,4,8) then do;
 do p=1 to 12;
  if p=dxm then do;
   if a3{p+23}=0 then medenr_mo=a1{p};
   if a3{p+23}>0 then do;
    do q=(dxm+23) to (dxm+12) by -1 until(a3{q}=0);
     if a2{q}>0 then bb=a2{q};
    end;
    medenr_mo=bb;
   end;
  end;
 end;
end;
if pattern5 in (5,6) then do;
 do r=1 to 12;
  if r=dxm then medenr_mo=a1{r};
 end;
end;
if pattern5=11 then do;
 do c=1 to 12;
  if c=dxm then do;
   if a3{c+23}>0 then do;
    do d=(dxm+23) to (dxm+12) by -1 until(a3{d}=0);
     if a2{d}>0 then bb=a2{d};
    end;
    medenr_mo=bb;
   end;
   if a3{c+23}=0 and a3{c+25}>0 then medenr_mo=a2{c+25};
   if a3{c+23}=0 and a3{c+25}=0 then medenr_mo=a2{c+26};
  end;
 end; 
end;
if pattern5=12 then do;
 do e=1 to 12;
  if e=dxm then do;
   if a3{e+23}>0 then do;
    do f=(dxm+23) to (dxm+12) by -1 until(a3{f}=0);
     if a2{f}>0 then bb=a2{f};
    end;
   end;
   if a3{e+23}=0 then do;
    do g=(dxm+25) to (dxm+36) until(a3{g}>0);
     if a2{g}>0 then bb=a2{g};
    end;
   end;
   medenr_mo=bb;
  end;
 end;
end;
if pattern5 in (13,14) then do; 
 do s=1 to 12;
  if s=dxm then do;
   if a3{s+25}=0 then do;
    do t=(dxm+25) to (dxm+36) until(a3{t}>0);
     if a2{t}>0 then bb=a2{t};
    end;
    medenr_mo=bb;
   end;
   if a3{s+25}>0 then medenr_mo=a2{s+25};
  end;
 end;
end;
if pattern5=18 then do;
 do u=1 to 12;
  if u=dxm then medenr_mo=a1{u};
 end;
end;
drop c d e f g h i j k l m n o p q r s t u v bb cc dd ee ff gg;
run;
data ffs2e;set ffs2d;ffsenr_mo=medenr_mo;ffspattern5=pattern5;ffspattern6=pattern6;ffspattern7=pattern7;
keep patient_id ffsenr_mo ffspattern5-ffspattern7;run;proc sort;by patient_id;run;

data ffs3a;set taf14;if mod(_n_,3)=1;
 array a1{48} ffsdate_bf1-ffsdate_bf24 ffsdate_at1-ffsdate_at12 ffsdate_af1-ffsdate_af12;
 array a2{48} ffsdate1_bf1-ffsdate1_bf24 ffsdate1_at1-ffsdate1_at12 ffsdate1_af1-ffsdate1_af12;
 do i=1 to 48;a2{i}=a1{i};end;
 drop ffsdate_bf1-ffsdate_bf24 ffsdate_at1-ffsdate_at12 ffsdate_af1-ffsdate_af12 i;
run;proc sort;by patient_id;run;
data ffs3b;set taf14;if mod(_n_,3)=2;
 array a1{48} ffsdate_bf1-ffsdate_bf24 ffsdate_at1-ffsdate_at12 ffsdate_af1-ffsdate_af12;
 array a2{48} ffsdate2_bf1-ffsdate2_bf24 ffsdate2_at1-ffsdate2_at12 ffsdate2_af1-ffsdate2_af12;
 do i=1 to 48;a2{i}=a1{i};end;
 keep patient_id ffsdate2_bf1-ffsdate2_bf24 ffsdate2_at1-ffsdate2_at12 ffsdate2_af1-ffsdate2_af12;
run;proc sort;by patient_id;run;
data ffs3c;set taf14;if mod(_n_,3)=0;
 array a1{48} ffsdate_bf1-ffsdate_bf24 ffsdate_at1-ffsdate_at12 ffsdate_af1-ffsdate_af12;
 array a2{48} ffsdate3_bf1-ffsdate3_bf24 ffsdate3_at1-ffsdate3_at12 ffsdate3_af1-ffsdate3_af12;
 do i=1 to 48;a2{i}=a1{i};end;
 keep patient_id ffsdate3_bf1-ffsdate3_bf24 ffsdate3_at1-ffsdate3_at12 ffsdate3_af1-ffsdate3_af12;
run;proc sort;by patient_id;run;
data ffs3d;merge ffs3a ffs3b ffs3c;by patient_id;
 array a1{48} ffsdate1_bf1-ffsdate1_bf24 ffsdate1_at1-ffsdate1_at12 ffsdate1_af1-ffsdate1_af12;
 array a2{48} ffsdate2_bf1-ffsdate2_bf24 ffsdate2_at1-ffsdate2_at12 ffsdate2_af1-ffsdate2_af12;
 array a3{48} ffsdate3_bf1-ffsdate3_bf24 ffsdate3_at1-ffsdate3_at12 ffsdate3_af1-ffsdate3_af12;
 array a4{48} ffsdate_bf1-ffsdate_bf24 ffsdate_at1-ffsdate_at12 ffsdate_af1-ffsdate_af12;
 do i=1 to 48;
  if (a1{i}=. and a2{i}=. and a3{i}>=0)|(a1{i}=. and a2{i}>=0 and a3{i}=.)|(a1{i}>=0 and a2{i}=. and a3{i}=.) then a4{i}=sum(a1{i},a2{i},a3{i});
  if (a1{i}>=0 and a2{i}>=0)|(a1{i}>=0 and a3{i}>=0)|(a2{i}>=0 and a3{i}>=0) then a4{i}=min(a1{i},a2{i},a3{i});
 end;
 drop i ffsdate1_bf1-ffsdate1_bf24 ffsdate1_at1-ffsdate1_at12 ffsdate1_af1-ffsdate1_af12 
        ffsdate2_bf1-ffsdate2_bf24 ffsdate2_at1-ffsdate2_at12 ffsdate2_af1-ffsdate2_af12
        ffsdate3_bf1-ffsdate3_bf24 ffsdate3_at1-ffsdate3_at12 ffsdate3_af1-ffsdate3_af12;
run;
data ffs3e;set ffs3d;
 array a1{12} ffsdate_at1-ffsdate_at12;
 array a2{48} ffsdate_bf1-ffsdate_bf24 ffsdate_at1-ffsdate_at12 ffsdate_af1-ffsdate_af12;
 array a3{48} ffsdate_bf25-ffsdate_bf48 ffsdate_at21-ffsdate_at32 ffsdate_af21-ffsdate_af32;

if dxy=. and dxybf1>. then dxy=dxybf1;if dxy=. and dxybf2>. then dxy=dxybf2;if dxy=. and dxyaf>. then dxy=dxyaf;
if dxm=. and dxmbf1>. then dxm=dxmbf1;if dxm=. and dxmbf2>. then dxm=dxmbf2;if dxm=. and dxmaf>. then dxm=dxmaf;

do n=1 to 48; if a2{n}>0 then a3{n}=1;if a2{n} in (0,.) then a3{n}=0;end;
cc=0;dd=0;ee=0;ff=0;
do i=1 to 12;
 if i=dxm then do;if a1{i}>0 then meddx=1;if a1{i} in (0,.) then meddx=0;  /*medicaid enrollment at diagnosis*/
  do l=(dxm+23) to (dxm+12) by -1;   /*count months of Medicaid enrollment in 12 months before cancer diagnosis*/
   cc=sum(cc,a3{l});
  end;
  do m=(dxm+25) to (dxm+36);    /*count months of Medicaid enrollment in 12 months after diagnosis*/
   dd=sum(dd,a3{m});
  end;
  do h=(dxm+25) to (dxm+30);    /*count months of Medicaid enrollment in 6 months after diagnosis*/
   ee=sum(ee,a3{h});
  end;
  do k=(dxm+25) to (dxm+26);    /*count months of Medicaid enrollment in 2 months after diagnosis*/
   ff=sum(ff,a3{k});
  end;
 end;
end;
if cc>10 then pattern1=1;if 0<cc<=10 then pattern1=0;if cc=0 then pattern1=2;label pattern1='0=discontinuous enroll within 1 year before diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if dd>10 then pattern2=1;if 0<dd<=10 then pattern2=0;if dd=0 then pattern2=2;label pattern2='0=discontinuous enroll withn 1 year after diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if ee=6 then pattern3=1;if 0<ee<6 then pattern3=0;if ee=0 then pattern3=2;label pattern3='0=discontinuous enroll within 6 months after diagnosis,1=continuous enroll (=6 mons),2=no enroll';
if ff=2 then pattern4=1;if 0<ff<2 then pattern4=0;if ff=0 then pattern4=2;label pattern4='0=discontinuous enroll within 2 months after diagnosis,1=continuous enroll (=2 mons),2=no enroll';
if meddx=1 then do;
 if pattern1=1 and pattern2=1 then pattern5=1;if pattern1=1 and pattern2=0 then pattern5=2;if pattern1=0 and pattern2=1 then pattern5=3; 
 if pattern1=0 and pattern2=0 then pattern5=4;if pattern1=2 and pattern2=1 then pattern5=5;if pattern1=2 and pattern2=0 then pattern5=6;
 if pattern1=1 and pattern2=2 then pattern5=7;if pattern1=0 and pattern2=2 then pattern5=8;if pattern1=2 and pattern2=2 then pattern5=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern2=1 then pattern5=9;if pattern1=1 and pattern2=0 then pattern5=10;if pattern1=0 and pattern2=1 then pattern5=11; 
 if pattern1=0 and pattern2=0 then pattern5=12;if pattern1=2 and pattern2=1 then pattern5=13;if pattern1=2 and pattern2=0 then pattern5=14;
 if pattern1=1 and pattern2=2 then pattern5=15;if pattern1=0 and pattern2=2 then pattern5=16;if pattern1=2 and pattern2=2 then pattern5=17;
end;
label pattern5='1=continuous enroll 1yr before,at and 1yr after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern3=1 then pattern6=1;if pattern1=1 and pattern3=0 then pattern6=2;if pattern1=0 and pattern3=1 then pattern6=3; 
 if pattern1=0 and pattern3=0 then pattern6=4;if pattern1=2 and pattern3=1 then pattern6=5;if pattern1=2 and pattern3=0 then pattern6=6;
 if pattern1=1 and pattern3=2 then pattern6=7;if pattern1=0 and pattern3=2 then pattern6=8;if pattern1=2 and pattern3=2 then pattern6=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern3=1 then pattern6=9;if pattern1=1 and pattern3=0 then pattern6=10;if pattern1=0 and pattern3=1 then pattern6=11; 
 if pattern1=0 and pattern3=0 then pattern6=12;if pattern1=2 and pattern3=1 then pattern6=13;if pattern1=2 and pattern3=0 then pattern6=14;
 if pattern1=1 and pattern3=2 then pattern6=15;if pattern1=0 and pattern3=2 then pattern6=16;if pattern1=2 and pattern3=2 then pattern6=17;
end;
label pattern6='1=continuous enroll 1yr before,at and 6months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern4=1 then pattern7=1;if pattern1=1 and pattern4=0 then pattern7=2;if pattern1=0 and pattern4=1 then pattern7=3; 
 if pattern1=0 and pattern4=0 then pattern7=4;if pattern1=2 and pattern4=1 then pattern7=5;if pattern1=2 and pattern4=0 then pattern7=6;
 if pattern1=1 and pattern4=2 then pattern7=7;if pattern1=0 and pattern4=2 then pattern7=8;if pattern1=2 and pattern4=2 then pattern7=18;
end;
if meddx=0 then do;
 if pattern1=1 and pattern4=1 then pattern7=9;if pattern1=1 and pattern4=0 then pattern7=10;if pattern1=0 and pattern4=1 then pattern7=11; 
 if pattern1=0 and pattern4=0 then pattern7=12;if pattern1=2 and pattern4=1 then pattern7=13;if pattern1=2 and pattern4=0 then pattern7=14;
 if pattern1=1 and pattern4=2 then pattern7=15;if pattern1=0 and pattern4=2 then pattern7=16;if pattern1=2 and pattern4=2 then pattern7=17;
end;
label pattern7='1=continuous enroll 1yr before,at and 2months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 

gg=0;
do o=(dxm+11) to dxm by -1;   /*count months of Medicaid enrollment in 13-24 months before cancer diagnosis*/
 gg=sum(gg,a3{o});
end;

bb=0;
if pattern5 in (1,2,7,9,10,15) then do;
 do j=1 to 12;
  if j=dxm then do;
   if gg>10 and a3{j}>0 then medenr_mo=a2{j};
   if gg>10 and a3{j}=0 then medenr_mo=a2{j+1};
   if 0<gg<=10 then do;
    if a3{j+11}=0 and a3{j+12}>0 then medenr_mo=a2{j+12};
	if a3{j+11}=0 and a3{j+12}=0 and a3{j+13}>0 then medenr_mo=a2{j+13};
	if a3{j+11}>0 then do;
     do v=(dxm+11) to dxm by -1 until(a3{v}=0);
	  if a2{v}>0 then bb=a2{v};
     end;
	medenr_mo=bb;
    end;
   end;
   if gg=0 and a3{j+13}>0 then medenr_mo=a2{j+13};
   if gg=0 and a3{j+13}=0 then medenr_mo=a2{j+14};
  end;
 end;
end;
if pattern5 in (3,4,8) then do;
 do p=1 to 12;
  if p=dxm then do;
   if a3{p+23}=0 then medenr_mo=a1{p};
   if a3{p+23}>0 then do;
    do q=(dxm+23) to (dxm+12) by -1 until(a3{q}=0);
     if a2{q}>0 then bb=a2{q};
    end;
    medenr_mo=bb;
   end;
  end;
 end;
end;
if pattern5 in (5,6) then do;
 do r=1 to 12;
  if r=dxm then medenr_mo=a1{r};
 end;
end;
if pattern5=11 then do;
 do c=1 to 12;
  if c=dxm then do;
   if a3{c+23}>0 then do;
    do d=(dxm+23) to (dxm+12) by -1 until(a3{d}=0);
     if a2{d}>0 then bb=a2{d};
    end;
    medenr_mo=bb;
   end;
   if a3{c+23}=0 and a3{c+25}>0 then medenr_mo=a2{c+25};
   if a3{c+23}=0 and a3{c+25}=0 then medenr_mo=a2{c+26};
  end;
 end; 
end;
if pattern5=12 then do;
 do e=1 to 12;
  if e=dxm then do;
   if a3{e+23}>0 then do;
    do f=(dxm+23) to (dxm+12) by -1 until(a3{f}=0);
     if a2{f}>0 then bb=a2{f};
    end;
   end;
   if a3{e+23}=0 then do;
    do g=(dxm+25) to (dxm+36) until(a3{g}>0);
     if a2{g}>0 then bb=a2{g};
    end;
   end;
   medenr_mo=bb;
  end;
 end;
end;
if pattern5 in (13,14) then do; 
 do s=1 to 12;
  if s=dxm then do;
   if a3{s+25}=0 then do;
    do t=(dxm+25) to (dxm+36) until(a3{t}>0);
     if a2{t}>0 then bb=a2{t};
    end;
    medenr_mo=bb;
   end;
   if a3{s+25}>0 then medenr_mo=a2{s+25};
  end;
 end;
end;
if pattern5=18 then do;
 do u=1 to 12;
  if u=dxm then medenr_mo=a1{u};
 end;
end;
drop c d e f g h i j k l m n o p q r s t u v bb cc dd ee ff gg;
run;
data ffs3f;set ffs3e;ffsenr_mo=medenr_mo;ffspattern5=pattern5;ffspattern6=pattern6;ffspattern7=pattern7;
keep patient_id ffsenr_mo ffspattern5-ffspattern7;run;proc sort;by patient_id;run;

data ffs4a;set taf4a;
 array a1{12} ffsdate13-ffsdate24;
 array a2{12} at1-at12;
 array a3{12} bf1-bf12;
 array a4{12} bf13-bf24;
 array a5{12} af1-af12;
 array a6{48} bf1-bf24 at1-at12 af1-af12;
 array a7{48} ff1-ff48;
cc=0;dd=0;ee=0;ff=0;
do i=1 to 12;
 if diff=0 then a2{i}=a1{i};
 if diff=1 then a4{i}=a1{i};
 if diff=2 then a3{i}=a1{i};
 if diff=-1 then a5{i}=a1{i};
end;
do n=1 to 48;
 if a6{n}>0 then a7{n}=1;if a6{n} in (0,.) then a7{n}=0;
end;
do h=1 to 12;
 if h=dxm and a2{h}>0 then meddx=1;
end;
do k=(dxm+23) to (dxm+12) by -1;  /*count months of Medicaid enrollment in 12 months before cancer diagnosis*/
 cc=sum(cc,a7{k});
end;
do l=(dxm+25) to (dxm+36);  /*count months of Medicaid enrollment in 12 months after diagnosis*/
 dd=sum(dd,a7{l});
end;
do m=(dxm+25) to (dxm+30); /*count months of Medicaid enrollment in 6 months after diagnosis*/
 ee=sum(ee,a7{m});
end;
do o=(dxm+25) to (dxm+26); /*count months of Medicaid enrollment in 2 months after diagnosis*/
 ff=sum(ff,a7{o});
end;

if cc>10 then pattern1=1;if 0<cc<=10 then pattern1=0;if cc=0 then pattern1=2;label pattern1='0=discontinuous enroll within 1 year before diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if dd>10 then pattern2=1;if 0<dd<=10 then pattern2=0;if dd=0 then pattern2=2;label pattern2='0=discontinuous enroll withn 1 year after diagnosis,1=continuous enroll (>=11 mons),2=no enroll';
if ee=6 then pattern3=1;if 0<ee<6 then pattern3=0;if ee=0 then pattern3=2;label pattern3='0=discontinuous enroll within 6 months after diagnosis,1=continuous enroll (=6 mons),2=no enroll';
if ff=2 then pattern4=1;if 0<ff<2 then pattern4=0;if ff=0 then pattern4=2;label pattern4='0=discontinuous enroll within 2 months after diagnosis,1=continuous enroll (=2 mons),2=no enroll';
if meddx=1 then do;
 if pattern1=1 and pattern2=1 then pattern5=1;if pattern1=1 and pattern2=0 then pattern5=2;if pattern1=0 and pattern2=1 then pattern5=3; 
 if pattern1=0 and pattern2=0 then pattern5=4;if pattern1=2 and pattern2=1 then pattern5=5;if pattern1=2 and pattern2=0 then pattern5=6;
 if pattern1=1 and pattern2=2 then pattern5=7;if pattern1=0 and pattern2=2 then pattern5=8;if pattern1=2 and pattern2=2 then pattern5=18;
end;
if meddx in (0,.) then do;
 if pattern1=1 and pattern2=1 then pattern5=9;if pattern1=1 and pattern2=0 then pattern5=10;if pattern1=0 and pattern2=1 then pattern5=11; 
 if pattern1=0 and pattern2=0 then pattern5=12;if pattern1=2 and pattern2=1 then pattern5=13;if pattern1=2 and pattern2=0 then pattern5=14;
 if pattern1=1 and pattern2=2 then pattern5=15;if pattern1=0 and pattern2=2 then pattern5=16;if pattern1=2 and pattern2=2 then pattern5=17;
end;
label pattern5='1=continuous enroll 1yr before,at and 1yr after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern3=1 then pattern6=1;if pattern1=1 and pattern3=0 then pattern6=2;if pattern1=0 and pattern3=1 then pattern6=3; 
 if pattern1=0 and pattern3=0 then pattern6=4;if pattern1=2 and pattern3=1 then pattern6=5;if pattern1=2 and pattern3=0 then pattern6=6;
 if pattern1=1 and pattern3=2 then pattern6=7;if pattern1=0 and pattern3=2 then pattern6=8;if pattern1=2 and pattern3=2 then pattern6=18;
end;
if meddx in (0,.) then do;
 if pattern1=1 and pattern3=1 then pattern6=9;if pattern1=1 and pattern3=0 then pattern6=10;if pattern1=0 and pattern3=1 then pattern6=11; 
 if pattern1=0 and pattern3=0 then pattern6=12;if pattern1=2 and pattern3=1 then pattern6=13;if pattern1=2 and pattern3=0 then pattern6=14;
 if pattern1=1 and pattern3=2 then pattern6=15;if pattern1=0 and pattern3=2 then pattern6=16;if pattern1=2 and pattern3=2 then pattern6=17;
end;
label pattern6='1=continuous enroll 1yr before,at and 6months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO'; 
if meddx=1 then do;
 if pattern1=1 and pattern4=1 then pattern7=1;if pattern1=1 and pattern4=0 then pattern7=2;if pattern1=0 and pattern4=1 then pattern7=3; 
 if pattern1=0 and pattern4=0 then pattern7=4;if pattern1=2 and pattern4=1 then pattern7=5;if pattern1=2 and pattern4=0 then pattern7=6;
 if pattern1=1 and pattern4=2 then pattern7=7;if pattern1=0 and pattern4=2 then pattern7=8;if pattern1=2 and pattern4=2 then pattern7=18;
end;
if meddx in (0,.) then do;
 if pattern1=1 and pattern4=1 then pattern7=9;if pattern1=1 and pattern4=0 then pattern7=10;if pattern1=0 and pattern4=1 then pattern7=11; 
 if pattern1=0 and pattern4=0 then pattern7=12;if pattern1=2 and pattern4=1 then pattern7=13;if pattern1=2 and pattern4=0 then pattern7=14;
 if pattern1=1 and pattern4=2 then pattern7=15;if pattern1=0 and pattern4=2 then pattern7=16;if pattern1=2 and pattern4=2 then pattern7=17;
end;
label pattern7='1=continuous enroll 1yr before,at and 2months after diagnosis,2=CE E DE,3=DE E CE,4=DE E DE,5=NO E CE,6=NO E DE,7=CE E NO,8=DE E NO,
9=CE NO CE,10=CE NO DE,11=DE NO CE,12=DE NO DE,13=NO NO CE,14=NO NO DE,15=CE NO NO,16=DE NO NO,17=NO NO NO,18=NO E NO'; 

gg=0;
do o=(dxm+11) to dxm by -1;   /*count months of Medicaid enrollment in 13-24 months before cancer diagnosis*/
 gg=sum(gg,a7{o});
end;

bb=0;
if pattern5 in (1,2,7,9,10,15) then do;
 do j=1 to 12;
  if j=dxm then do;
   if gg>10 and a7{j}>0 then medenr_mo=a6{j};
   if gg>10 and a7{j}=0 then medenr_mo=a6{j+1};
   if 0<gg<=10 then do;
    if a7{j+11}=0 and a7{j+12}>0 then medenr_mo=a6{j+12};
	if a7{j+11}=0 and a7{j+12}=0 and a7{j+13}>0 then medenr_mo=a6{j+13};
	if a7{j+11}>0 then do;
     do v=(dxm+11) to dxm by -1 until(a7{v}=0);
	  if a6{v}>0 then bb=a6{v};
     end;
	medenr_mo=bb;
    end;
   end;
   if gg=0 and a7{j+13}>0 then medenr_mo=a6{j+13};
   if gg=0 and a7{j+13}=0 then medenr_mo=a6{j+14};
  end;
 end;
end;
if pattern5 in (3,4,8) then do;
 do p=1 to 12;
  if p=dxm then do;
   if a7{p+23}=0 then medenr_mo=a2{p};
   if a7{p+23}>0 then do;
    do q=(dxm+23) to (dxm+12) by -1 until(a7{q}=0);
     if a6{q}>0 then bb=a6{q};
    end;
    medenr_mo=bb;
   end;
  end;
 end;
end;
if pattern5 in (5,6) then do;
 do r=1 to 12;
  if r=dxm then medenr_mo=a2{r};
 end;
end;
if pattern5=11 then do;
 do c=1 to 12;
  if c=dxm then do;
   if a7{c+23}>0 then do;
    do d=(dxm+23) to (dxm+12) by -1 until(a7{d}=0);
     if a6{d}>0 then bb=a6{d};
    end;
    medenr_mo=bb;
   end;
   if a7{c+23}=0 and a7{c+25}>0 then medenr_mo=a6{c+25};
   if a7{c+23}=0 and a7{c+25}=0 then medenr_mo=a6{c+26};
  end;
 end; 
end;
if pattern5=12 then do;
 do e=1 to 12;
  if e=dxm then do;
   if a7{e+23}>0 then do;
    do f=(dxm+23) to (dxm+12) by -1 until(a7{f}=0);
     if a6{f}>0 then bb=a6{f};
    end;
   end;
   if a7{e+23}=0 then do;
    do g=(dxm+25) to (dxm+36) until(a7{g}>0);
     if a6{g}>0 then bb=a6{g};
    end;
   end;
   medenr_mo=bb;
  end;
 end;
end;
if pattern5 in (13,14) then do; 
 do s=1 to 12;
  if s=dxm then do;
   if a7{s+25}=0 then do;
    do t=(dxm+25) to (dxm+36) until(a7{t}>0);
     if a6{t}>0 then bb=a6{t};
    end;
    medenr_mo=bb;
   end;
   if a7{s+25}>0 then medenr_mo=a6{s+25};
  end;
 end;
end;
if pattern5=18 then do;
 do u=1 to 12;
  if u=dxm then medenr_mo=a2{u};
 end;
end;
drop c d e f g h i j k l m n o p q r s t u v bb cc dd ee ff gg;
run;
data ffs4b;set ffs4a;ffsenr_mo=medenr_mo;ffspattern5=pattern5;ffspattern6=pattern6;ffspattern7=pattern7;
keep patient_id ffsenr_mo ffspattern5-ffspattern7;run;proc sort;by patient_id;run;

data ffs5;set ffs1a ffs2e ffs3f ffs4b;run;proc sort;by patient_id;run;  /*FFS enrollment timeline*/

data taf2a;set taf2;if year=dxy;keep patient_id dual13 medicaid1-medicaid5;run;proc sort;by patient_id;run;  
data taf2b taf2c;set taf2a;by patient_id;if first.patient_id and last.patient_id then output taf2b;else output taf2c;run;  /*taf2c includes cases with 2+ observations*/
proc sort data=taf2c;by patient_id;run;
proc freq data=taf2c noprint;tables patient_id/out=taf2d;run;
data taf2e;set taf2d;if count=2;b=1;run;proc sort;by patient_id;run;
data taf2f;merge taf2c taf2e;by patient_id;if b=1;run;
data taf2g;set taf2f;if mod(_n_,2)=1;
 array a1{6} dual13 medicaid1-medicaid5;
 array a2{6} b1-b6;
 do i=1 to 6;a2{i}=a1{i};end;
 keep patient_id b1-b6;
 run;proc sort;by patient_id;run;
data taf2h;set taf2f;if mod(_n_,2)=0;
 array a1{6} dual13 medicaid1-medicaid5;
 array a2{6} c1-c6;
 do i=1 to 6;a2{i}=a1{i};end;
 keep patient_id c1-c6;
 run;proc sort;by patient_id;run;
data taf2i;merge taf2g taf2h;by patient_id;
array a1{6} b1-b6;
array a2{6} c1-c6;
array a3{6} dual13 medicaid1-medicaid5;
do i=1 to 6;
 if a1{i}=.|a2{i}=. then a3{i}=sum(a1{i},a2{i});
 if a1{i}>=0 and a2{i}>=0 then a3{i}=min(a1{i},a2{i});
end;
keep patient_id dual13 medicaid1-medicaid5;
run;proc sort;by patient_id;run;
data taf2j;set taf2d;if count=3;b=1;run;proc sort;by patient_id;run;
data taf2k;merge taf2c taf2j;by patient_id;if b=1;run;
data taf2l;set taf2k;if mod(_n_,3)=1;
 array a1{6} dual13 medicaid1-medicaid5;
 array a2{6} b1-b6;
 do i=1 to 6;a2{i}=a1{i};end;
 keep patient_id b1-b6;
 run;proc sort;by patient_id;run;
data taf2m;set taf2k;if mod(_n_,3)=2;
 array a1{6} dual13 medicaid1-medicaid5;
 array a2{6} c1-c6;
 do i=1 to 6;a2{i}=a1{i};end;
 keep patient_id c1-c6;
 run;proc sort;by patient_id;run;
 data taf2n;set taf2k;if mod(_n_,3)=0;
 array a1{6} dual13 medicaid1-medicaid5;
 array a2{6} d1-d6;
 do i=1 to 6;a2{i}=a1{i};end;
 keep patient_id d1-d6;
 run;proc sort;by patient_id;run;
data taf2o;merge taf2l taf2m taf2n;by patient_id;
array a1{6} b1-b6;
array a2{6} c1-c6;
array a3{6} d1-d6;
array a4{6} dual13 medicaid1-medicaid5;
do i=1 to 6;
  if (a1{i}=. and a2{i}=. and a3{i}>=0)|(a1{i}=. and a2{i}>=0 and a3{i}=.)|(a1{i}>=0 and a2{i}=. and a3{i}=.) then a4{i}=sum(a1{i},a2{i},a3{i});
  if (a1{i}>=0 and a2{i}>=0)|(a1{i}>=0 and a3{i}>=0)|(a2{i}>=0 and a3{i}>=0) then a4{i}=min(a1{i},a2{i},a3{i});
end;
keep patient_id dual13 medicaid1-medicaid5;
run;proc sort;by patient_id;run;

/******link Medicaid enrollment patterns to SEER data******/
proc sort data=taf2b;by patient_id;run;
proc sort data=seer2;by patient_id;run;
data seer3;merge seer2 taf15 mcc5 ffs5 taf2b taf2i taf2o;by patient_id;if aa=1;
/*insurance at diagnosis*/
if medicaid1=. and pattern5 in (9,10,11,12,13,14,15,16,17) then medicaid1=0;label medicaid1='0=no enrollment at diagnosis 1=enrolled at diagnosis';
if medicaid4=. and mccpattern5 in (9,10,11,12,13,14,15,16,17) then medicaid4=0;label medicaid4='0=no managed care enrollment at diagnosis 1=enrolled at diagnosis';
if medicaid5=. and ffspattern5 in (9,10,11,12,13,14,15,16,17) then medicaid5=0;label medicaid5='0=no FFS enrollment at diagnosis 1=enrolled at diagnosis';

if 1<=pattern5<17|pattern5=18 then medicaid22=1;else medicaid22=0;label medicaid22='0=no enrollment 1 yr before to 1yr after diagnosis,1=enrolled';
if pattern5 in (1,2,7,9,10,15) then conenroll_bf=2;if pattern5 in (3,4,8,11,12,16) then conenroll_bf=1;if pattern5 in (5,6,13,14,17,18) then conenroll_bf=0;
label conenroll_bf='continuous enrollment before diagnosis 0=no E 1=discontinuous E 2=continuous E';
if pattern6 in (1,3,5,9,11,13) then conenroll_af6=2;if pattern6 in (2,4,6,10,12,14,18) then conenroll_af6=1;if pattern6 in (7,8,15,16,17) then conenroll_af6=0;
if pattern7 in (1,3,5,9,11,13) then conenroll_af2=2;if pattern7 in (2,4,6,10,12,14,18) then conenroll_af2=1;if pattern7 in (7,8,15,16,17) then conenroll_af2=0;
label conenroll_af6='continuous enrollment 6mon after diagnosis 0=no E 1=discontinuous E 2=continuous E';
label conenroll_af2='continuous enrollment 2mon after diagnosis 0=no E 1=discontinuous E 2=continuous E';
label dual13='1=fuall dual eligible 2=partial dual 4=medicare not medicaid 5=not medicare';

if insurance2 in (1,2) and medicaid22=0 then insurdx=9;
if insurance2 in (1,2,10) and medicaid22=1 then insurdx=2;
if insurance2 in (20,21) and medicaid22=0 then insurdx=8;
if insurance2 in (20,21) and medicaid22=1 then insurdx=10;
if insurance2 in (31,35) and medicaid22=1 then insurdx=2;
if insurance2 in (31,35) and medicaid22=0 then insurdx=9;
if insurance2 in (60,61,64) and medicaid22=1 then insurdx=3;
if insurance2 in (60,61,64) and medicaid22=0 then insurdx=1;
if insurance2=62 and medicaid22=1 then insurdx=7;
if insurance2=62 and medicaid22=0 then insurdx=6;
if insurance2=63 and medicaid22=0 then insurdx=4;
if insurance2=63 and medicaid22=1 then insurdx=5;
if insurance2=99 and medicaid22=1 then insurdx=2;
if insurance2 in (65,66,67,68) then insurdx=12;
if insurdx=. then insurdx=13;
label insurdx='insurance at diagnosis (medicaid status based on 1yr to 1yr after diagnosis:
               1=FFS,2=Medicaid only,3=dual,4=FFS+Private,5=FFS+Private+Medicaid,6=MA,7=MA+Medicaid,8=private only,
               9=uninsured,10=private+medicaid,12=other,13=missing';

if 1<=mccpattern5<17|mccpattern5=18 then mccmedicaid22=1;else mccmedicaid22=0;label mccmedicaid22='0=no managed care enrollment 1 yr before to 1yr after diagnosis,1=enrolled';
if mccpattern5 in (1,2,7,9,10,15) then mccconenroll_bf=2;if mccpattern5 in (3,4,8,11,12,16) then mccconenroll_bf=1;
if mccpattern5 in (5,6,13,14,17,18) then mccconenroll_bf=0;
if mccpattern6 in (1,3,5,9,11,13) then mccconenroll_af6=2;if mccpattern6 in (2,4,6,10,12,14,18) then mccconenroll_af6=1;
if mccpattern6 in (7,8,15,16,17) then mccconenroll_af6=0;
if mccpattern7 in (1,3,5,9,11,13) then mccconenroll_af2=2;if mccpattern7 in (2,4,6,10,12,14,18) then mccconenroll_af2=1;
if mccpattern7 in (7,8,15,16,17) then mccconenroll_af2=0;
label mccconenroll_af6='continuous managed care enrollment 6mon after diagnosis 0=no E 1=discontinuous E 2=continuous E';
label mccconenroll_af2='continuous managed care enrollment 2mon after diagnosis 0=no E 1=discontinuous E 2=continuous E';
label mccconenroll_bf='continuous managed care enrollment in a year before diagnosis 0=no E 1=discontinuous E 2=continuous E';

if 1<=ffspattern5<17|ffspattern5=18 then ffsmedicaid22=1;else ffsmedicaid22=0;label ffsmedicaid22='0=no FFS enrollment 1 yr before to 1yr after diagnosis,1=enrolled';
if ffspattern5 in (1,2,7,9,10,15) then ffsconenroll_bf=2;if ffspattern5 in (3,4,8,11,12,16) then ffsconenroll_bf=1;
if ffspattern5 in (5,6,13,14,17,18) then ffsconenroll_bf=0;
if ffspattern6 in (1,3,5,9,11,13) then ffsconenroll_af6=2;if ffspattern6 in (2,4,6,10,12,14,18) then ffsconenroll_af6=1;
if ffspattern6 in (7,8,15,16,17) then ffsconenroll_af6=0;
if ffspattern7 in (1,3,5,9,11,13) then ffsconenroll_af2=2;if ffspattern7 in (2,4,6,10,12,14,18) then ffsconenroll_af2=1;
if ffspattern7 in (7,8,15,16,17) then ffsconenroll_af2=0;
label ffsconenroll_af6='continuous FFS enrollment 6mon after diagnosis 0=no E 1=discontinuous E 2=continuous E';
label ffsconenroll_af2='continuous FFS enrollment 2mon after diagnosis 0=no E 1=discontinuous E 2=continuous E';
label ffsconenroll_bf='continuous FFS enrollment in a year before diagnosis 0=no E 1=discontinuous E 2=continuous E';

/*continuous enrollment pattern based on different cutoffs*/
 dx_months=(dxy-2005)*12+dxm;
 if medicaid22=1 then do;conenroll=dx_months-medenr_mo; 
  if conenroll>1 then enrollstatus1=2;if 0<conenroll<=1 then enrollstatus1=1;if .<conenroll<=0 then enrollstatus1=0;if enrollstatus1=. and pattern5=16 then enrollstatus1=0;
  if conenroll>3 then enrollstatus3=2;if 0<conenroll<=3 then enrollstatus3=1;if .<conenroll<=0 then enrollstatus3=0;if enrollstatus3=. and pattern5=16 then enrollstatus3=0;
  if conenroll>6 then enrollstatus6=2;if 0<conenroll<=6 then enrollstatus6=1;if .<conenroll<=0 then enrollstatus6=0;if enrollstatus6=. and pattern5=16 then enrollstatus6=0;
  if conenroll>9 then enrollstatus9=2;if 0<conenroll<=9 then enrollstatus9=1;if .<conenroll<=0 then enrollstatus9=0;if enrollstatus9=. and pattern5=16 then enrollstatus9=0;
  if conenroll>12 then enrollstatus12=2;if 0<conenroll<=12 then enrollstatus12=1;if .<conenroll<=0 then enrollstatus12=0;if enrollstatus12=. and pattern5=16 then enrollstatus12=0;
 end;

 if medicaid22=1 then do;mccconenroll=dx_months-mccenr_mo; 
  if mccconenroll>1 then mccenrollstatus1=2;if 0<mccconenroll<=1 then mccenrollstatus1=1;if .<mccconenroll<=0 then mccenrollstatus1=0;if mccenrollstatus1=. and mccpattern5=16 then mccenrollstatus1=0;
  if mccconenroll>3 then mccenrollstatus3=2;if 0<mccconenroll<=3 then mccenrollstatus3=1;if .<mccconenroll<=0 then mccenrollstatus3=0;if mccenrollstatus3=. and mccpattern5=16 then mccenrollstatus3=0;
  if mccconenroll>6 then mccenrollstatus6=2;if 0<mccconenroll<=6 then mccenrollstatus6=1;if .<mccconenroll<=0 then mccenrollstatus6=0;if mccenrollstatus6=. and mccpattern5=16 then mccenrollstatus6=0;
  if mccconenroll>9 then mccenrollstatus9=2;if 0<mccconenroll<=9 then mccenrollstatus9=1;if .<mccconenroll<=0 then mccenrollstatus9=0;if mccenrollstatus9=. and mccpattern5=16 then mccenrollstatus9=0;
  if mccconenroll>12 then mccenrollstatus12=2;if 0<mccconenroll<=12 then mccenrollstatus12=1;if .<mccconenroll<=0 then mccenrollstatus12=0;if mccenrollstatus12=. and mccpattern5=16 then mccenrollstatus12=0;
 end;

 if medicaid22=1 then do;ffsconenroll=dx_months-ffsenr_mo; 
  if ffsconenroll>1 then ffsenrollstatus1=2;if 0<ffsconenroll<=1 then ffsenrollstatus1=1;if .<ffsconenroll<=0 then ffsenrollstatus1=0;if ffsenrollstatus1=. and ffspattern5=16 then ffsenrollstatus1=0;
  if ffsconenroll>3 then ffsenrollstatus3=2;if 0<ffsconenroll<=3 then ffsenrollstatus3=1;if .<ffsconenroll<=0 then ffsenrollstatus3=0;if ffsenrollstatus3=. and ffspattern5=16 then ffsenrollstatus3=0;
  if ffsconenroll>6 then ffsenrollstatus6=2;if 0<ffsconenroll<=6 then ffsenrollstatus6=1;if .<ffsconenroll<=0 then ffsenrollstatus6=0;if ffsenrollstatus6=. and ffspattern5=16 then ffsenrollstatus6=0;
  if ffsconenroll>9 then ffsenrollstatus9=2;if 0<ffsconenroll<=9 then ffsenrollstatus9=1;if .<ffsconenroll<=0 then ffsenrollstatus9=0;if ffsenrollstatus9=. and ffspattern5=16 then ffsenrollstatus9=0;
  if ffsconenroll>12 then ffsenrollstatus12=2;if 0<ffsconenroll<=12 then ffsenrollstatus12=1;if .<ffsconenroll<=0 then ffsenrollstatus12=0;if ffsenrollstatus12=. and ffspattern5=16 then ffsenrollstatus12=0;
 end;

/*histology classification based on Blom EF. Ann Am Thorac Soc 2020 with modifications + Potter AL BMJ 2022 + Howlader N N Engl J Med 2021*/
if histo3 in (8015,8050,8140,8141,8143,8144,8145,8147,8200,8201,8230,8250,8251,8252,8253,8254,8255,8256,8257,8260,8263,8265,8290,8310,8323,8333,8440,
              8470,8471,8480,8481,8490,8550,8551,8570,8571,8572,8573,8574,8575,8576) then his1=0;
if histo3 in (8023,8051,8052,8070,8071,8072,8073,8074,8075,8076,8078,8082,8083,8084,8094,8123) then his1=1;
if histo3 in (8012,8013,8014,8021,8034) then his1=2;
if histo3 in (8003,8004,8005,8011,8022,8030,8031,8032,8033,8035,8046,8430,8560,8562,8980) then his1=3;
if histo3 in (8002,8041,8042,8043,8044,8045) then his1=4;
IF HISTO3 in (8240,8241,8242,8243,8244,8245,8246) then his1=5;
if his1=. and histo3>. then his1=6;
if his1=. then his1=7;
label his1='histology 0=adenocarcinoma 1=squamous cell 2=large cell 3=NSC others 4=small cell 5=typical carcinoid 6=other 7=missing';
if his1 in (0,1,2,3,5) then his2=1;if his1=4 then his2=2;
if histo3 in (8250,8251,8252,8255,8256,8257) then his3=1;
label his3='NSCLC formerly classified as Bronchiolo-alveolar Adenocarcinoma';
 
/*age*/
IF 18<=AGE<50 then agegp=1;
IF 50<=AGE<55 then agegp=2;
IF 55<=AGE<60 then agegp=3;
IF 60<=AGE<65 then agegp=4;
IF 65<=AGE<70 then agegp=5;
IF 70<=AGE<75 then agegp=6;
IF 75<=AGE<80 then agegp=7;
IF 80<=AGE<85 then agegp=8;
IF 85<=AGE then agegp=9;

if agegp=1 then agegp2=1;
if agegp in (2,3) then agegp2=2;
if agegp in (4,5) then agegp2=3;
if agegp in (6,7) then agegp2=4;
if agegp in (8,9) then agegp2=5;
label agegp='age at dx: 1=18-49 2=50-54 3=55-59 4=60-64 5=65-69 6=70-74 7=75-79 8=80-84 9=85+';

if 18<=age<65 then agegp3=0;if 65<=age<=120 then agegp3=1;

if 18<=age<20 then agegp4=4;
if 20<=age<25 then agegp4=5;
if 25<=age<30 then agegp4=6;
if 30<=age<35 then agegp4=7;
if 35<=age<40 then agegp4=8;
if 40<=age<45 then agegp4=9;
if 45<=age<50 then agegp4=10;
IF 50<=AGE<55 then agegp4=11;
IF 55<=AGE<60 then agegp4=12;
IF 60<=AGE<65 then agegp4=13;
IF 65<=AGE<70 then agegp4=14;
IF 70<=AGE<75 then agegp4=15;
IF 75<=AGE<80 then agegp4=16;
IF 80<=AGE<85 then agegp4=17;
IF 85<=AGE then agegp4=18;

/*cancer stage: stage5/stage6=regional was defined as stage III; MA NY ID TX have no detailed staging info */
if stage7 in ('0','0A','0IS') then stage=0;
if stage7 in ('1','1A','1B','1B1') then stage=1;
if stage7 in ('2','2A','2A2','2B','2C') then stage=2;
if stage7 in ('3','3A','3B','3C','3C1') then stage=3;
if stage7 in ('4','4A','4B','4C') then stage=4;
if stage=. then do;
 if stage8 in ('0','0a','0is') then stage=0;
 if stage8 in ('1','1A','1A1','1A2','1A3','1B','1B1','1E') then stage=1;
 if stage8 in ('2','2A','2A2','2B','2C','2E') then stage=2;
 if stage8 in ('3','3A','3B','3C') then stage=3;
 if stage8 in ('4','4A','4B','4C') then stage=4;
end;
if stage=. then do;
 if 0<=stage1<=2 then stage=0;
 if 10<=stage1<30 then stage=1;
 if 30<=stage1<50 then stage=2;
 if 50<=stage1<70 then stage=3;
 if 70<=stage1<=74 then stage=4;
end;
if stage=. then do;
 if 0<=stage2<=20 then stage=0;
 if 100<=stage2<300 then stage=1;
 if 300<=stage2<500 then stage=2;
 if 500<=stage2<700 then stage=3;
 if 700<=stage2<=740 then stage=4;
end;
if stage=. then do;
 if stage3=0 then stage=0;
 if stage3=7 then stage=4;
 if stage3=1 then stage=1;
 if stage3=2 then stage=2;
 if stage3 in (3,4) then stage=3;
end;
if stage=. then do;
 if stage6='0' then stage=0;
 if stage6='1' then stage=1;
 if stage6='7' then stage=4;
 if stage6='2' then stage=3;
end;
if stage=. then stage=5;

/*convert AJCC 6th edition and EOD18 to AJCC 7th edition*/
 if stage2=0 then stage11=0;
 if stage2 in (120,150) then stage11=1;
 if stage2 in (300,320,330) then stage11=2;
 if stage2 in (500,520,530) then stage11=3;
 if stage2 in (700,720,730) then stage11=4;
 if stage11=. then do;
  if stage7 in ('0') then stage11=0;
  if stage7 in ('1A','1B') then stage11=1;
  if stage7 in ('2','2A','2B') then stage11=2;
  if stage7 in ('3','3A','3B') then stage11=3;
  if stage7 in ('4') then stage11=4;
 end;
 if stage11=. then do;
  if stage1=0 then stage11=0;
  if stage1=12 then stage11=1;
  if stage1=15 and (tnmother3<=50|tnmother3 in (991,992,993,994,995)) then stage11=1;
  if stage1=15 and 50<tnmother3<=989 then stage11=2; 
  if stage1=15 and stage11=. then stage11=1;
  if stage1 in (30,32,33) then stage11=2;
  if stage1=30 and tnmt1='40' and tnmn1 in ('00','10') and tnmm1='00' then stage11=3;
  if stage1 in (52,53) then stage11=3;
  if 70<=stage1<=74 then stage11=4;  
  if stage1=70 and tnmt1='40' and tnmn1 in ('00','10') and tnmm1='00' then stage11=3;
 end;
 if stage11=. then do;
  if stage8 in ('0','0a','0is') then stage11=0;
  if stage8 in ('1','1A','1A1','1A2','1A3','1B','1B1','1E') then stage11=1;
  if stage8 in ('2','2A','2A2','2B','2C','2E') then stage11=2;
  if stage8 in ('3','3A','3B','3C') then stage11=3;
  if stage8 in ('4','4A','4B','4C') then stage11=4;
 end;
 if stage11=. then stage11=5;

/*radiation*/
if rad1 in (0,7) then rad=0;
if rad1 in (1,2,3,4,5,6) then rad=1;
if rad=. and treatseq1>0 then rad=1;
if rad=. and treatseq1=0 then rad=0;
if rad=. then rad=2;
label rad='radiation 0=No 1=yes 2=missing';

/*chemo*/
if chemo=0 and treatseq2>0 then chemo=1;
if chemo=. and treatseq2>0 then chemo=1;
if chemo=. and treatseq2=0 then chemo=0;
if chemo=. and stage in (0,1) then chemo=0;

/*surgery*/
if surg1=0 then surg=0;  /*no surgery*/
if 10<surg1<30 then surg=1;  /*<1 lobe*/
if 30<=surg1<55 then surg=2; /*1 lobe<= but < l lung*/
if 55<=surg1<80 then surg=3; /*pneumonectomy*/
if surg1 in (80,90) then surg=4; /*NOS*/
if surg=. then do;
 if treatseq1='0' and treatseq2='0' then surg=0;
 if treatseq1 in ('2','3','4','5','6','7','9')|treatseq2 in ('2','3','4','5','6','7','9') then surg=4;
end;
label surg='surgery 0=no surg 1=<1 lobe 2=1 lobe-<1 lung 3=pneumonectomy 4=NOS 5=missing';
if surg=0 then surgery=0;if surg in (1,2,3,4) then surgery=1;

if 10<=nodeex<95 then surgq=1;if .<nodeex<10|nodeex=95 then surgq=2;if surgq=. then surgq=3;
label surgq='1=10+ LN examined for stage IA, IB, IIA, IIB NSCLC 2=<10 3=missing';

/*survival*/
if event1=1 then outcome1=1;if event2=1 then outcome1=2;if outcome1=. and event3=0 then outcome1=2;if outcome1=. then outcome1=0; 
if event3=0 then outcome2=1;if outcome2=. then outcome2=0; /*overall survival*/
label outcome1='outcome 0=alive 1=dead due to LC 2=dead due to other causes';

/*guideline-concordant treatment, Blom EF. Ann Am Thorac Soc 2020*/
if his2=1 then do;   /*treat1=0 nonconcordant treat1=1 less intensive than recommended*/
 if stage in (1,2) then do;
  if surg in (1,2,3,4)|rad1 in (1,4) then treat1=2;
  if treat1=. and rad=1 then treat1=1;
  if treat1=. and surg=0 and rad=0 then treat1=0;;
 end;
 if stage=3 then do;
  if (surg in (1,2,3,4) and chemo=1)|(rad=1 and chemo=1) then treat1=2;
  if treat1=. and (surg in (1,2,3,4)|chemo=1|rad=1) then treat1=1;
  if treat1=. and surg=0 and rad=0 and chemo=0 then treat1=0;
 end;
 if stage=4 then do;
  if chemo=1 then treat1=2;
  if treat1=. and (surg in (1,2,3,4)|rad=1) then treat1=1;
  if treat1=. and chemo=0 and rad=0 and surg=0 then treat1=0; 
 end; 
 if stage=0 then do;if surgery=1 then treat1=2;if surgery=0 then treat1=0;end;
end;
if his2=2 then do;
 if stage in (1,2,3) then do;
  if (surg in (1,2,3,4) and chemo=1)|(rad=1 and chemo=1) then treat1=2;
  if treat1=. and (surg in (1,2,3,4)|chemo=1|rad=1) then treat1=1;
  if treat1=. and surg=0 and rad=0 and chemo=0 then treat1=0;
 end;
 if stage=4 then do;
  if chemo=1 then treat1=2;
  if treat1=. and (surg in (1,2,3,4)|rad=1) then treat1=1;
  if treat1=. and chemo=0 and rad=0 and surg=0 then treat1=0; 
 end;
end;
label treat1='treatment adherence 0=no treat 1=less intensive 2=adherent';

/*follow up months: study cutoff is December 2020; registry 61,62,63,66 had no outcome info*/
dxmonth=(dxy-2006)*12+dxm;
followmonth=(followy-2006)*12+followm;
deathmonth=(deathy-2006)*12+deathm;
if event1 in (1,8)|event2 in (1,8)|event3=0 then event4=1;  /*death*/
if time1=. or time1=9999 then do;
 if event4=1 then time1=deathmonth-dxmonth;
 if event4 ne 1 then time1=followmonth-dxmonth;
 if .<time1<0 then time1=0;
end;
treatdelay1=txdate-dxdate;
if treatdelay1>2 then delay1=1;if .<treatdelay1<=2 then delay1=0;
if treatdelay1>3 then delay2=1;if .<treatdelay1<=3 then delay2=0;
label delay1='treatment initiation 1=>2mons 0=within 2mons';

/*urban/rural*/
if rural13 in (1,2,3) then rural=0;if rural13 in (4,5,6,7,8,9) then rural=1;if rural=. then rural=0;

/*marital status*/
if mari in ('2','6') then marital=0;if mari in ('1','3','4','5') then marital=1;if marital=. then marital=2;

/*Medicaid FFS coverage of lung cancer screening in 2018-2019
https://medquest.hawaii.gov/en/members-applicants/fee-for-service.html
*/
if reg1 in ('06','09','13','16','19','21','25','26','36','53') then coverscr=1;
if reg1 in ('22','48','49') then coverscr=0;
label coverscr='FFS Medicaid covering LCS: 0=No 1=Yes';
if reg1 in ('06','26','36','53') then screlig=1;
if reg1 in ('09','13','19','21','25') then screlig=2;
if reg1 in ('16') then screlig=3;
label screlig='eligibility for LC screening: 1=USPSTF,2=Other,3=Medicare';
if reg1 in ('06','13','16','21','25','26') then author=0;
if reg1 in ('09','19','36','53') then author=1;
label author='prior authorization: 0=No 1=Yes';
if coverscr=0 then coverauth=0;
if coverscr=1 and author=1 then coverauth=1;
if coverscr=1 and author=0 then coverauth=2;
label coverauth='FFS Medicaid covering LCS and prior authorization: 0=No coverage,1=coverage and prior author,2=coverage no prior author';
run;proc sort;by cty;run;

/******link neighborhood-level measures to SEER data******/
data ses;set seer.sesrural;cty=FIPs;keep cty ses10 ses00 ses00cat4 ses00cat5 ses10cat4 ses10cat5;run;proc sort;by cty;run;
data access1;set acce.acc2605_cty;run;proc sort;by cty;run;
data seer4;merge seer3 ses access1;by cty;if aa=1;
/*socioeconomic deprivation index*/
 if .<dxy<2010 then do;ses=ses00;ses4=ses00cat4;ses5=ses00cat5;end;
 if dxy>=2010 then do;ses=ses10;ses4=ses10cat4;ses5=ses10cat5;end;
 if ses5=. then ses5=5;  /*5 cases in New Mexico*/
 if ses4=. then ses4=4;
 if ses10cat5=. then ses10cat5=5;
if dxy<2013 then do;
 if rural03 in ('01','02','03') then lcs1=lcs10q;
 if rural03 in ('04','05','06','07') then lcs1=lcs30q;
 if rural03 in ('08','09') then lcs1=lcs60s;
end;
if dxy>=2013 then do;
 if rural13 in ('01','02','03') then lcs1=lcs10q;
 if rural13 in ('04','05','06','07') then lcs1=lcs30q;
 if rural13 in ('08','09') then lcs1=lcs60s;
end;
run;

/******data used for the analysis of geographic access to LCS centers and early-stage diagnosis******/
data LCS1;set seer4;
if 55<=age<=80;
if dxy in (2015,2016,2017,2018,2019);
if report in (6,7) then delete;
if lcs1>.;

/*define early-stage diagnosis. If SEER Combined Summary Stage is missing, use variables stage to determine early-stage.*/
if stage6 in (0,1) then stage102=1;if stage6 in (2,7) then stage102=0;
if stage102=. then do; if stage in (0,1,2) then stage102=1;if stage in (3,4) then stage102=0;end;

if race4 in ('3','9') then race4='9';
if ses4 in (1,2) then ses2=0;if ses4 in (3,4) then ses2=1;
if pattern5 in (1,2,3,4,5,6,7,8,18) then medicaid23=1; 
if pattern5 in (9,10,11,12,15,16) and conenroll>=0 then medicaid23=1;
if medicaid23=. then medicaid23=0;/*medicaid enrollment at diagnosis*/
if insurance2 in (1,2) and medicaid23=0 then insurdxn=9;
if insurance2 in (1,2,10) and medicaid23=1 then insurdxn=2;
if insurance2 in (20,21) and medicaid23=0 then insurdxn=8;
if insurance2 in (20,21) and medicaid23=1 then insurdxn=10;
if insurance2 in (31,35) and medicaid23=1 then insurdxn=2;
if insurance2 in (31,35) and medicaid23=0 then insurdxn=9;
if insurance2 in (60,61,64) and medicaid23=1 then insurdxn=3;
if insurance2 in (60,61,64) and medicaid23=0 then insurdxn=1;
if insurance2=62 and medicaid23=1 then insurdxn=7;
if insurance2=62 and medicaid23=0 then insurdxn=6;
if insurance2=63 and medicaid23=0 then insurdxn=4;
if insurance2=63 and medicaid23=1 then insurdxn=5;
if insurance2=99 and medicaid23=1 then insurdxn=2;
if insurance2 in (65,66,67) then insurdxn=11;
if insurance2=68 then insurdxn=12;
if insurdxn=. then insurdxn=13;
label insurdxn='insurance at diagnosis (medicaid status based on at and 1yr before diagnosis:
               1=FFS,2=Medicaid only,3=dual,4=FFS+Private,5=FFS+Private+Medicaid,6=MA,7=MA+Medicaid,8=private only,
               9=uninsured,10=private+medicaid,11=military related insurance 12=Indian,13=missing'; 
if insurdxn in (1,6) then insurdxn1=1;if insurdxn in (3,7) then insurdxn1=3;if insurdxn in (4,5) then insurdxn1=4;if insurdxn in (8,10) then insurdxn1=8;
if insurdxn in (2,9,11,12) then insurdxn1=insurdxn;if insurdxn1=. then insurdxn1=13;
label insurdxn1='insurance at diagnosis: 1=FFS/MA 2=Medicaid only,3=dual,4=FFS+private,8=private/private+Medicaid,9=uninsured,11=military related insurance,
                 12=indian 13=missing';
if insurdxn in (1,6,3,7,4,5) then insurdxn2=1;if insurdxn in (8,10) then insurdxn2=8;if insurdxn in (2,9,11) then insurdxn2=insurdxn;
if insurdxn2=. then insurdxn2=12;
label insurdxn2='insurance at diagnosis: 1=Medicare 2=Medicaid only,8=private/private+Medicaid,9=uninsured,11=military related insurance,12=missing';
if 55<=age<70 then agegp5=0;if 70<=age<=80 then agegp5=1;
if his1=6 then his2=3;if his1=7 then his2=4;
if reg in (1,21,25,26,31,35,41,23,61) then region=1;  /*West*/
 if reg in (2,44,62,63) then region=2;  /*Northeast*/
 if reg in (20,22) then region=3;  /*Midwest*/
 if reg in (27,37,42,43,47,66) then region=5;  /*South*/

if mccconenroll_af6=2 and ffsconenroll_af6 in (1,0) then medicaidma=1;if mccconenroll_af6=2 and ffsconenroll_af6=2 then medicaidma=3;
if mccconenroll_af6=1 and ffsconenroll_af6=0 then medicaidma=1;if mccconenroll_af6=1 and ffsconenroll_af6=1 then medicaidma=3;
if mccconenroll_af6 in (0,1) and ffsconenroll_af6=2 then medicaidma=2;
if mccconenroll_af6=0 and ffsconenroll_af6=1 then medicaidma=2;
if mccconenroll_af6=0 and ffsconenroll_af6=0 then do;
  if medicaid4=1 then medicaidma=1;if medicaid5=1 then medicaidma=2;  /*managed care at diagnosis*/
  if medicaidma=. then do;
   if mccconenroll_bf in (1,2) and ffsconenroll_bf=0 then medicaidma=1;
   if mccconenroll_bf=0 and ffsconenroll_bf in (1,2) then medicaidma=2;
   if mccconenroll_bf=2 and ffsconenroll_bf=1 then medicaidma=1;
   if mccconenroll_bf=1 and ffsconenroll_bf=2 then medicaidma=2;
  end;
end;
if medicaidma=. then medicaidma=4;
label medicaidma='1=managed care 2=FFS 3=both 4=unknown';
agegp6=agegp;if agegp in (7,8) then agegp6=7;

run;
data lcs2;set lcs1;if stage102 in (0,1);run;
proc rank data=lcs2 out=lcs2a groups=4 ties=low;
 var carownership lcs1;
 ranks cargp4 lcs1gp4;
 run;
proc rank data=lcs2a out=lcs2b groups=2 ties=low;
 var carownership lcs1;
 ranks cargp2 lcs1gp2;
run;

/******Table 1******/
proc freq data=lcs2b;tables agegp6 sex race4 insurdxn2 dxy his2 rural stage102;run;
proc means data=lcs2b median mean std;var age time1;run;
proc means data=lcs2b mean std;class lcs1gp4;var age;run;
proc freq data=lcs2b;tables agegp6 sex race4 insurdxn2 dxy his2 ses5 rural cargp2 region;run;
proc freq data=lcs2b;tables lcs1gp4*(agegp6 sex race4 insurdxn2 dxy his2 ses5 rural cargp2 region)/chisq;run;
data lcs2b1;set lcs2b;if race4=9;run;proc freq;tables race1-race3;run;
data lcs2b2;set lcs2b;if insurdxn2=12;run;proc freq;tables insurdxn;run;
data lcs2b3;set lcs2b;if his2=3;run;proc freq;tables his1;run;

/******Figure 1******/
proc freq data=lcs2b;tables lcs1gp4*stage102 agegp5*lcs1gp2*stage102 
insurdxn2*lcs1gp2*stage102 race4*lcs1gp2*stage102 rural*lcs1gp2*stage102 cargp2*lcs1gp2*stage102 ses2*lcs1gp2*stage102 sex*lcs1gp2*stage102
region*lcs1gp2*stage102;
run;
/*Modified Poisson regression (Zou log-linear method) with a robust variance estimator to estimate prevalence ratio*/
proc sort data=lcs2b;by cty;run;
proc genmod data=lcs2b;   
class agegp6 sex rural ses5 marital lcs1gp4 medicaidma race4 stage cty cargp2 insurdxn2 reg1 his2 dxy region patient_id/ref=last;
model stage102=lcs1gp4 agegp6 sex race4 insurdxn2 dxy his2 ses5 rural/dist=poisson link=log type3;
repeated subject=cty/type=ind; 
/*applying a robust (sandwich) error variance estimator. This converts standard Poisson into a modified Poisson model. 
It ensures p value and 95% CI are statistically valid for binary data lacking a time component.
When run modified Poisson, type=ind should be used*/
lsmeans lcs1gp4/pdiff cl exp;
run;  
proc genmod data=lcs2b;
class agegp6 sex rural ses5 marital lcs1gp4 medicaidma race4 stage cty patient_id cargp2 insurdxn2 reg1 his2 region/ref=last;
model stage102=lcs1 agegp6 sex race4 insurdxn2 dxy his2 ses5 rural/dist=poisson link=log type3;
repeated subject=cty/type=ind;
run;  
proc genmod data=lcs2b descending;
class lcs1gp2 agegp5 sex rural ses5 marital medicaidma race4 stage cty patient_id cargp2 insurdxn2 reg1 dxy his2 region/ref=last;
model stage102=lcs1gp2|agegp5 sex race4 insurdxn2 ses5 rural dxy his2/dist=poisson link=log type3;
repeated subject=cty/type=ind;
lsmeans lcs1gp2*agegp5/pdiff cl exp;
run;  
proc genmod data=lcs2b descending;
class lcs1gp2 agegp6 sex rural ses5 marital medicaidma race4 stage cty patient_id cargp2 insurdxn2 reg1 dxy his2 region/ref=last;
model stage102=lcs1gp2|sex agegp6 race4 insurdxn2 ses5 rural dxy his2/dist=poisson link=log type3;
repeated subject=cty/type=ind;
lsmeans lcs1gp2*sex/pdiff cl exp;
run;  
proc genmod data=lcs2b descending;
class lcs1gp2 agegp6 sex ses5 marital medicaidma race4 stage cty patient_id cargp2 insurdxn2 reg1 his2 dxy region rural/ref=last;
model stage102=lcs1gp2|rural his2 insurdxn2 sex agegp6 race4 ses5 dxy/dist=poisson link=log type3;
repeated subject=cty/type=ind;
lsmeans lcs1gp2*rural/pdiff cl exp;
run;  
proc genmod data=lcs2b descending;
class lcs1gp2 agegp6 sex rural ses5 marital medicaidma race4 stage cty patient_id cargp2 insurdxn2 reg1 his2 ses2 dxy region/ref=last;
model stage102=lcs1gp2|ses2 his2 insurdxn2 sex agegp6 race4 rural dxy /dist=poisson link=log type3;
repeated subject=cty/type=ind;
lsmeans lcs1gp2*ses2/pdiff cl exp;
run;  
data lcs2c;set lcs2b;if race4 in (1,2,4,5);run;proc sort;by cty;run;
proc genmod data=lcs2c descending;
class agegp6 sex rural ses5 marital lcs1gp2 medicaidma race4 stage cty patient_id cargp2 insurdxn2 reg1 dxy his2 region/ref=last;
model stage102=lcs1gp2|race4 agegp6 sex insurdxn2 ses5 rural dxy his2/dist=poisson link=log type3;
repeated subject=cty/type=ind;
lsmeans lcs1gp2*race4/pdiff cl exp;
lsmestimate lcs1gp2*race4 'DID in Asian vs NHW'  -1 0 1 0 1 0 -1 0;
lsmestimate lcs1gp2*race4 'DID in Hispanic vs NHW'  -1 0 0 1 1 0 0 -1;
lsmestimate lcs1gp2*race4 'DID in NHB vs NHW'  -1 1 0 0 1 -1 0 0;
run;  
proc genmod data=lcs2c descending;
class agegp6 sex rural ses5 marital lcs1gp2 medicaidma race4 stage cty patient_id cargp2 insurdxn2 reg1 dxy his2 region/ref=last;
model stage102=lcs1gp2|race4 agegp6 sex insurdxn2 ses5 rural dxy his2 reg1/dist=poisson link=log type3;
repeated subject=cty/type=ind;
lsmeans lcs1gp2*race4/pdiff cl exp;
lsmestimate lcs1gp2*race4 'DID in Asian vs NHW'  -1 0 1 0 1 0 -1 0;
lsmestimate lcs1gp2*race4 'DID in Hispanic vs NHW'  -1 0 0 1 1 0 0 -1;
lsmestimate lcs1gp2*race4 'DID in NHB vs NHW'  -1 1 0 0 1 -1 0 0;
run;  
proc genmod data=lcs2c descending;
class agegp6 sex rural ses5 marital lcs1gp2 medicaidma race4 stage cty patient_id cargp2 insurdxn2 reg1 dxy his2 region/ref=last;
model stage102=lcs1gp2|race4 agegp6 sex insurdxn2 dxy his2 ses5 rural/dist=poisson link=log type3;
repeated subject=patient_id(cty)/type=ind;
lsmeans lcs1gp2*race4/pdiff cl exp;
lsmestimate lcs1gp2*race4 'DID in Asian vs NHW'  -1 0 1 0 1 0 -1 0;
lsmestimate lcs1gp2*race4 'DID in Hispanic vs NHW'  -1 0 0 1 1 0 0 -1;
lsmestimate lcs1gp2*race4 'DID in NHB vs NHW'  -1 1 0 0 1 -1 0 0;
run;  
proc genmod data=lcs2c descending;
class agegp6 sex rural ses5 marital lcs1gp2 medicaidma race4 stage cty patient_id cargp2 insurdxn2 reg1 dxy his2 region/ref=last;
model stage102=lcs1gp2|race4 agegp6 sex insurdxn2 dxy his2 ses5 rural/dist=binomial link=logit type3;
repeated subject=cty/type=ind;
lsmeans lcs1gp2*race4/pdiff cl exp;
lsmestimate lcs1gp2*race4 'DID in Asian vs NHW'  -1 0 1 0 1 0 -1 0;
lsmestimate lcs1gp2*race4 'DID in Hispanic vs NHW'  -1 0 0 1 1 0 0 -1;
lsmestimate lcs1gp2*race4 'DID in NHB vs NHW'  -1 1 0 0 1 -1 0 0;
run;  
proc genmod data=lcs2b descending;   
class agegp6 sex rural ses5 marital lcs1gp4 medicaidma race4 stage cty cargp2 insurdxn2 reg1 his2 dxy region patient_id/ref=last;
model stage102=lcs1gp4 agegp6 sex race4 insurdxn2 dxy his2 ses5 rural/dist=binomial link=logit type3;
repeated subject=cty/type=ind; 
lsmeans lcs1gp4/pdiff cl exp;
run;  
data lcs2d;set lcs2b;if insurdxn2 in (1,2,8,9,11);run;proc sort;by cty;run;
proc genmod data=lcs2d descending;
class agegp6 sex rural ses5 marital lcs1gp2 medicaidma race4 stage cty patient_id cargp2 insurdxn2 reg1 dxy his2 region/ref=last;
model stage101=lcs1gp2|insurdxn2 sex agegp6 race4 ses5 rural his2 dxy/dist=poisson link=log type3;
repeated subject=cty/type=ind;
lsmeans lcs1gp2*insurdxn2/pdiff cl exp;
lsmestimate lcs1gp2*insurdxn2 'DID in Medicare vs private'  1 0 -1 0 0 -1 0 1 0 0;
lsmestimate lcs1gp2*insurdxn2 'DID in no insurance vs private' 0 0 -1 1 0 0 0 1 -1 0;
lsmestimate lcs1gp2*insurdxn2 'DID in Medicaid vs private' 0 1 -1 0 0 0 -1 1 0 0;
lsmestimate lcs1gp2*insurdxn2 'DID in military vs private' 0 0 -1 0 1 0 0 1 0 -1;
run;  
proc genmod data=lcs2b descending;
class agegp6 sex rural ses5 marital lcs1gp2 medicaidma race4 stage cty patient_id cargp2 insurdxn2 reg1 dxy his2 region/ref=last;
model stage102=lcs1gp2|region agegp6 sex race4 insurdxn2 ses5 rural dxy his2/dist=poisson link=log type3;
repeated subject=cty/type=ind;
lsmeans lcs1gp2*region/pdiff cl exp;
lsmestimate lcs1gp2*region 'DID in northeast vs south'  0 1 0 -1 0 -1 0 1;
lsmestimate lcs1gp2*region 'DID in south vs west'  -1 0 0 1 1 0 0 -1;
lsmestimate lcs1gp2*region 'DID in midwest vs south'  0 0 1 -1 0 0 -1 1;
run;  

/******Figure 2******/
proc freq data=lcs2b;tables cargp2*lcs1gp4*stage102;run;
proc genmod data=lcs2b descending;
class lcs1gp4 agegp6 sex rural ses5 marital medicaidma race4 stage cty patient_id cargp2 insurdxn2 reg1 dxy his2 region/ref=last;
model stage102=lcs1gp4|cargp2 agegp6 sex race4 insurdxn2 ses5 rural dxy his2/dist=poisson link=log type3;
repeated subject=cty/type=ind;
lsmeans lcs1gp4*cargp2/pdiff cl exp;
run;  
proc genmod data=lcs2b descending;
class lcs1gp2 agegp6 sex rural ses5 marital medicaidma race4 stage cty patient_id cargp2 insurdxn2 reg1 dxy his2 region/ref=last;
model stage102=lcs1gp2|cargp2 agegp6 sex race4 insurdxn2 ses5 rural dxy his2/dist=poisson link=log type3;
repeated subject=cty/type=ind;
lsmeans lcs1gp2*cargp2/pdiff cl exp;
run;  
proc sort data=lcs2b;by cargp2;run;
proc genmod data=lcs2b descending;by cargp2;
class lcs1gp4 agegp6 sex rural ses5 marital medicaidma race4 stage cty patient_id cargp2 insurdxn2 reg1 dxy his2 region/ref=last;
model stage102=lcs1 agegp6 sex race4 insurdxn2 ses5 rural dxy his2/dist=poisson link=log type3;
repeated subject=cty/type=ind;
run;  

/******Table S1******/
data lcs2e;set lcs2b;if rural=0;run;proc sort;by cty;run;
/******read in the dataset containing density of public transit stops******/
data stop1;set da38605p1;tract1=TRACT_FIPS*1;cty=int(tract1/1000000);run;proc sort;by cty;run;
proc means data=stop1 sum;by cty;var CENSUS_TRACT_AREA COUNT_NTM_STOPS TOT_POP_2010;
output out=stop2 sum(CENSUS_TRACT_AREA)=area sum(COUNT_NTM_STOPS)=count sum(TOT_POP_2010)=pop;
run;
data stop3;set stop2;countbypop=count/pop;countbyarea=count/area;keep cty countbypop countbyarea;run;proc sort;by cty;run;
data lcs2f;merge lcs2e stop3;by cty;if aa=1;run;
proc rank data=lcs2f out=lcs2f1 groups=4 ties=low;
 var countbyarea countbypop;
 ranks countbyarea4 countbypop4;
 run;
proc rank data=lcs2f1 out=lcs2f2 groups=2 ties=low;
 var countbyarea countbypop;
 ranks countbyarea2 countbypop2;
 run;
proc rank data=lcs2f2 out=lcs2f3 groups=3 ties=low;
 var countbyarea countbypop;
 ranks countbyarea3 countbypop3;
 run;
data lcs2f4;set lcs2f3;
if countbyarea2=0 and cargp2=0 then transport4=0;if countbyarea2=1|cargp2=1 then transport4=1;
run;
proc freq data=lcs2f4;tables transport4*lcs1gp4*stage102;run; 
proc genmod data=lcs2f4;
class lcs1gp4 agegp6 sex rural ses5 marital medicaidma race4 stage cty patient_id cargp2 insurdxn2 reg1 dxy his2 region transport4/ref=last;
model stage102=lcs1gp4|transport4 agegp6 sex race4 ses5 dxy his2 insurdxn2/dist=poisson link=log type3;
repeated subject=cty/type=ind;
lsmeans lcs1gp4*transport4/pdiff cl exp;
run;  
proc sort data=lcs2f4;by transport4;run;
proc genmod data=lcs2f4;by transport4;
class lcs1gp4 agegp6 sex rural ses5 marital medicaidma race4 stage cty patient_id cargp2 insurdxn2 reg1 dxy his2 region transport4/ref=last;
model stage102=lcs1 agegp6 sex race4 ses5 dxy his2 insurdxn2/dist=poisson link=log type3;
repeated subject=cty/type=ind;
run;  

/******Figure 3******/
data lcs2f6;set lcs2f4;if insurdxn2 in (2);run;
proc freq data=lcs2f6;tables transport4*lcs1gp4*stage102;run;
proc genmod data=lcs2f6;
class lcs1gp4 agegp6 sex rural ses5 marital medicaidma race4 stage cty patient_id cargp2 insurdxn2 reg1 dxy his2 region transport4/ref=last;
model stage102=lcs1gp4|transport4 agegp6 sex race4 ses5 dxy his2/dist=poisson link=log type3;
repeated subject=cty/type=ind;
lsmeans lcs1gp4*transport4/pdiff cl exp;
run;  
proc sort data=lcs2f6;by transport4;run;
proc genmod data=lcs2f6;by transport4;
class lcs1gp4 agegp6 sex rural ses5 marital medicaidma race4 stage cty patient_id cargp2 insurdxn2 reg1 dxy his2 region transport4/ref=last;
model stage102=lcs1 agegp6 sex race4 ses5 dxy his2/dist=poisson link=log type3;
repeated subject=cty/type=ind;
run;  

/******Table S2******/
data lcs2f7;set lcs2f6;
if reg1 in (48,49) then delete;  /*delete 2 states where Medicaid FFS did not cover LCS and no Medicaid expansion*/
/*medicaid expansion start dates:
california/connecticut/Hawaii/Iowa/Kentucky/Massachussett/NJ/NM/NY/WA 1/1/2014   
georgia/TX no 
Idaho/Utah 1/1/2020
Louisiana 7/1/2016
Michiga 4/1/2014
*/
run;
proc freq data=lcs2f7;tables transport4*lcs1gp4*stage102;run;
proc genmod data=lcs2f7;
class lcs1gp4 agegp6 sex rural ses5 marital medicaidma race4 stage cty patient_id cargp2 insurdxn2 reg1 dxy his2 region transport4/ref=last;
model stage102=lcs1gp4|transport4 age sex race4 dxy his2 ses5/dist=poisson link=log type3;
repeated subject=cty/type=ind;
lsmeans lcs1gp4*transport4/pdiff cl exp;
run;  
proc sort data=lcs2f7;by transport4;run;
proc genmod data=lcs2f7;by transport4;
class lcs1gp2 agegp6 sex rural ses5 marital medicaidma race4 stage cty patient_id cargp2 insurdxn2 reg1 dxy his2 region countbyarea2/ref=last;
model stage102=lcs1 agegp6 sex race4 ses5 dxy his2/dist=poisson link=log type3;
repeated subject=cty/type=ind;
run;  



