physiological.tab <- fread("~/project/gtps/gtps/doc/magnet_redcap.txt")
#physiological.tab <- physiological.tab[sample_name %in% pedtrait.info$sample_name ]
physiological.tab <- physiological.tab[1:1809]

#write.table(file="~/project/gtps/doc/physiological.tab",x = physiological.tab , row.names = F, 
	#col.names =T,  sep="\t", quote=F )
physio.subset <- physiological.tab
physio.subset <- physio.subset[,list(sample_name, site, age, gender, race, ethnicity, patient_weight_kg, height_cm, etiology, heart_weight_grams, prior_cabg, prior_valve_surgery, prior_device,
		     is_the_patient_on_dialysis, history_of_afib_aflutter,type_of_afib, prior_cardioversion, prior_ablation, maze_procedure, history_of_vt_vf,
		    prior_icd_shock, prior_sustained_vt, prior_vt_ablation, history_of_diabetes, insulin, oral_hypoglycemics, history_of_hypertension,
		    ace_inhibitor, aldosterone_antagonist, amiodarone, angiotensin_ii_antagonist, aspirin, beta_blocker, calcium_channel_blocker, aniti_coagulnats,
		    digoxin, diurectics, hydralazine, nitrate, other_antiarrythmic, anti_platelet, pde5_inhibitors, statin, other_lipid_lowering_drug, thyroid_hormone,
		    is_patient_dobutamine, is_patient_dopamine, is_patient_on_milrinone, lvedd, lvesd,lvef , mitral_regurgitation, tricuspid_regurgitation, hemos_avail, hemos_rap,
		    hemos_systolic, hemos_diastolic, hemos_pwcp, hemos_co, cardiac_index, sbp, dbp, creatinine_level, bnp) ]

physio.subset[,site.p:=sapply(site, function(xx) switch(xx, "Penn"=0, "Cleveland"=1, NA))]
physio.subset[, age.p:=as.numeric(age)]
physio.subset[,gender.p:=sapply(gender, function(xx) switch(xx, "Male"=0, "Female"=1, NA))]
physio.subset[,race.p:=sapply(race, function(xx) switch(xx, "White/Caucasian"=0, "Other"=1, NA))]
physio.subset[,ethnicity.p:=sapply(ethnicity, function(xx) switch(xx, "Not Hispanic or Latino"=0, "Hispanic or Latino"=1, NA))]
physio.subset[,race.p:=sapply(race, function(xx) switch(xx, "White/Caucasian"=0, "Other"=1, NA))]
physio.subset[,patient_weight_kg.p:=as.numeric(patient_weight_kg)]
physio.subset[,patient_weight_kg.p:=ifelse((patient_weight_kg.p==0),NA, patient_weight_kg) ]

physio.subset[,height_cm.p:=as.numeric(height_cm)]
physio.subset[,height_cm.p:=ifelse((height_cm.p==0),NA, height_cm.p)]

physio.subset[,etiology.Idiopathic.p:=sapply(etiology, function(xx) switch(xx, "Idiopathic Dilated CMP"=1, 0))]

physio.subset[,etiology.Ischemic.p:=sapply(etiology, function(xx) switch(xx, "Ischemic"=1, 0))]

physio.subset[,etiology.donor.p:=sapply(etiology, function(xx) switch(xx, "donor"=1, 0))]

physio.subset[,etiology.donor.p:=sapply(etiology, function(xx) switch(xx, "other"=1, 0))]

physio.subset[,etiology.Valvular.p:=sapply(etiology, function(xx) switch(xx, "Valvular"=1, 0))]

physio.subset[,heart_weight_grams.p:=as.numeric(heart_weight_grams)]
physio.subset[,heart_weight_grams.p:=ifelse(heart_weight_grams.p==0,NA, heart_weight_grams.p)]


physio.subset[,prior_cabg.p:=sapply(prior_cabg, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]

physio.subset[,prior_valve_surgery.p:=sapply(prior_valve_surgery, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]

physio.subset[,prior_device.ICD.p:=sapply(prior_device, function(xx) switch(xx, "ICD"=1, 0))]

physio.subset[,prior_device.BIVICD.p:=sapply(prior_device, function(xx) switch(xx, "BIVICD"=1, 0))]

physio.subset[,prior_device.pacemaker.p:=sapply(prior_device, function(xx) switch(xx, "Pacemaker"=1, 0))]

physio.subset[,prior_device.pacemaker.p:=sapply(prior_device, function(xx) switch(xx, "Pacemaker"=1, 0))]

physio.subset[,history_of_afib_aflutter.p:=sapply(history_of_afib_aflutter, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]

#physio.subset[,type_of_afib.p:=sapply(type_of_afib, function(xx) switch(xx, "Paroxsymal"=1, "Permanent"=0, NA))]
#
#
#physio.subset[,prior_cardioversion.p:=sapply(prior_cardioversion, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
#
physio.subset[,prior_ablation.p:=sapply(prior_ablation, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]

physio.subset[,prior_icd_shock.p:=sapply(prior_icd_shock, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]

#physio.subset[,prior_sustained_vt.p:=sapply(prior_sustained_vt, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
#
physio.subset[,prior_vt_ablation.p:=sapply(prior_vt_ablation, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,history_of_diabetes.p:=sapply(history_of_diabetes, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,insulin.p:=sapply(insulin, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,oral_hypoglycemics.p:=sapply(oral_hypoglycemics, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,history_of_hypertension.p:=sapply(history_of_hypertension, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,ace_inhibitor.p:=sapply(ace_inhibitor, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,aldosterone_antagonist.p:=sapply(aldosterone_antagonist, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]

physio.subset[,amiodarone.p:=sapply(amiodarone, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,angiotensin_ii_antagonist.p:=sapply(angiotensin_ii_antagonist, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,beta_blocker.p:=sapply(beta_blocker, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,calcium_channel_blocker.p:=sapply(calcium_channel_blocker, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,aniti_coagulnats.p:=sapply(aniti_coagulnats, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,digoxin.p:=sapply(digoxin, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,diurectics.p:=sapply(diurectics, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,hydralazine.p:=sapply(hydralazine, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,nitrate.p:=sapply(nitrate, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,other_antiarrythmic.p:=sapply(other_antiarrythmic, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]

physio.subset[,anti_platelet.p:=sapply(anti_platelet, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,pde5_inhibitors.p:=sapply(pde5_inhibitors, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,statin.p:=sapply(statin, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,other_lipid_lowering_drug.p:=sapply(other_lipid_lowering_drug, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,thyroid_hormone.p:=sapply(thyroid_hormone, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,is_patient_dobutamine.p:=sapply(is_patient_dobutamine, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,is_patient_dopamine.p:=sapply(is_patient_dopamine, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,is_patient_on_milrinone.p:=sapply(is_patient_on_milrinone, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
#physio.subset[,lvedd.p:=as.numeric(lvedd)]
physio.subset[, lvedd.p:=ifelse(grepl(lvedd, pattern="normal", ignore.case=T),4.8 , ifelse(grepl(lvedd, pattern="-"),
					sapply(lvedd, function(tt) mean(as.numeric(strsplit(unique(tt)[1],split="\\-")[[1]]))), 
					ifelse( grepl(lvedd, pattern="<25%"), 25, ifelse(grepl(lvedd, pattern="<"), as.numeric(gsub(lvedd, pattern="<", replacement="")),
									     as.numeric(lvedd)))))]

physio.subset[,lvedd.p:=ifelse(lvedd.p==0,NA,lvedd.p)]
physio.subset[,lvedd.p:=ifelse(lvedd.p> 20, lvedd.p/10,lvedd.p)]
physio.subset[,lvesd.p:=as.numeric(lvesd)]
physio.subset[,lvesd.p:=ifelse(lvesd.p==0,NA,lvesd.p)]
physio.subset[,lvesd.p:=ifelse(lvesd.p> 20, lvesd.p/10,lvesd.p)]

physio.subset[, lvef:=ifelse(grepl(lvef, pattern="<"), gsub(lvef, pattern="<", replacement=""), lvef)  ]
physio.subset[, lvef:=ifelse(grepl(lvef, pattern=">"), gsub(lvef, pattern=">", replacement=""), lvef)  ]
physio.subset[, lvef:=ifelse(grepl(lvef, pattern="%"), gsub(lvef, pattern="%", replacement=""), lvef)  ]
physio.subset[, lvef.p:=ifelse(grepl(lvef, pattern="normal", ignore.case=T), 58, ifelse(grepl(lvef, pattern="-"),
					sapply(lvef, function(tt) mean(as.numeric(strsplit(unique(tt)[1],split="\\-")[[1]]))), 
					ifelse( grepl(lvef, pattern="<25%"), 25, ifelse(grepl(lvef, pattern="<"), as.numeric(gsub(lvef, pattern="<", replacement="")),
									     as.numeric(lvef)))))]
physio.subset[,lvef.p:=ifelse(lvef.p==0,NA,lvef.p)]
physio.subset[,lvef.p:=ifelse(lvef.p<1, lvef.p * 100,lvef.p)]


#sapply( as.numeric(strsplit(unique(physio.subset$lvef)[1],split="\\-")[[1]])
physio.subset[,mitral_regurgitation.p:=sapply(mitral_regurgitation, function(xx) switch(xx, "none" = 0, "trace/mild"=1, "mild to moderate"=2, "moderate"=3, "moderate to severe"=4, "severe"=5, NA))]
physio.subset[,tricuspid_regurgitation.p:=sapply(tricuspid_regurgitation, function(xx) switch(xx, "none" = 0, "trace/mild"=1, "mild to moderate"=2, "moderate"=3, "moderate to severe"=4, "severe"=5, NA))]

#physio.subset[,hemos_avail.p:=sapply(hemos_avail, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
physio.subset[,hemos_rap.p:=as.numeric(hemos_rap)]
physio.subset[,hemos_rap.p:=ifelse(hemos_rap.p==0,NA,hemos_rap.p)]
physio.subset[,hemos_systolic.p:=as.numeric(hemos_systolic)]
physio.subset[,hemos_systolic.p:=ifelse(hemos_systolic.p==0,NA,hemos_systolic.p)]
physio.subset[,hemos_diastolic.p:=as.numeric(hemos_diastolic) ] 
physio.subset[,hemos_diastolic.p:=ifelse(hemos_diastolic.p==0,NA,hemos_diastolic.p)]
physio.subset[,hemos_pwcp.p:=as.numeric(hemos_pwcp)]
physio.subset[,hemos_pwcp.p:=ifelse(hemos_pwcp.p==0,NA,hemos_pwcp.p)]
physio.subset[,hemos_co.p:=as.numeric(hemos_co)]
physio.subset[,hemos_co.p:=ifelse(hemos_co.p==0,NA,hemos_co.p)]
physio.subset[,cardiac_index.p:=as.numeric(cardiac_index)]
physio.subset[,cardiac_index.p:=ifelse(cardiac_index.p==0,NA,cardiac_index.p)]
physio.subset[,cardiac_index.p:=ifelse(cardiac_index.p> 10, cardiac_index.p/10,cardiac_index.p)]
physio.subset[,sbp.p:=as.numeric(sbp)]
physio.subset[,sbp.p:=ifelse(sbp.p==0,NA,sbp.p)]
physio.subset[,dbp.p:=as.numeric(dbp)]
physio.subset[,dbp.p:=ifelse(dbp.p==0,NA,dbp.p)]
physio.subset[,creatinine_level.p:=as.numeric(creatinine_level)]
physio.subset[,creatinine_level.p:=ifelse(creatinine_level.p==0,NA,creatinine_level.p)]
physio.subset[,creatinine_level.p:=ifelse(creatinine_level.p> 20, creatinine_level.p/10,creatinine_level.p)]
physio.subset[,bnp.p:=as.numeric(bnp)]
physio.subset[,bnp.p:=ifelse(bnp.p==0,NA,bnp.p)]


 physio.subset[,disease.p:=ifelse(etiology=="donor",1,2)]
physio.reduced <- physio.subset[,list(site.p, age.p, gender.p, race.p, ethnicity.p, patient_weight_kg.p, height_cm.p, heart_weight_grams.p, prior_cabg.p,
				   history_of_afib_aflutter.p, history_of_diabetes.p, history_of_hypertension.p, ace_inhibitor.p, lvedd.p, lvesd.p, lvef.p,
				   mitral_regurgitation.p, tricuspid_regurgitation.p, hemos_rap.p, hemos_systolic.p,hemos_pwcp.p,cardiac_index.p, creatinine_level.p,
				   disease.p, etiology)]
physio.reduced <- as.data.frame(physio.reduced)
