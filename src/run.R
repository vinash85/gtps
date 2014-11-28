load("/cbcbhomes/vinash85/shortcuts/lasso/data/Data/MikeData/magnet_PedTraits.RData")
library(data.table)
require(mi) 
pedtrait.info <- pedtrait[, 1:5,with=F]
expression <- pedtrait[, grep(x=colnames(pedtrait), pattern="^X"), with=F]
expression.z <- scale(expression, scale=apply(expression, 2, sd), center=T)
rownames(expression.z)  = pedtrait.info$sample_name


# analysis with limma

library(limma)
exp.mat  <- t(as.matrix(expression))
design <- cbind(grp1=1, disease=pedtrait$disease) 
fit  <-  lmFit(exp.mat, design)
# Fold-change thresholding
fit2 <- treat(fit,lfc=0.1)
topTreat(fit2,coef=2)
fit2 <- eBayes(fit2) 
date=Sys.Date()
outDir <- paste0("result/", date)
dir.create(outDir)

topGene <- topTreat(fit2,coef=2,200)
topGene$affyid <- gsub(rownames(topGene), pattern="^X",replacement="")
if(F){
  library(biomaRt)
  ensembl <- useMart("ensembl", dataset= "hsapiens_gene_ensembl")
  affy2ensembl = getBM(attributes = c("chromosome_name", "start_position", "end_position","strand","ensembl_gene_id", "affy_hugene_1_0_st_v1"), mart=ensembl)
  affy2ensembl = affy2ensembl[!is.na(affy2ensembl$affy_hugene_1_0_st_v1) ,]
  topGene.ensembl  <- merge(x=topGene, y=affy2ensembl, by.x = "affyid", by.y="affy_hugene_1_0_st_v1", sort = F)


  write.fit(fit2,file=paste0(outDir, "/limMaFit.MagNet.Genes.txt"))
  pdf(paste0(outDir, "/limMaFit.MagNet.Genes.pdf"))
  volcanoplot(fit2,coef=2,highlight=100, names=rownames(exp.mat))
  dev.off()

  write.table(file="limMaFit.diff.200.Genes.txt" , x =topTreat(fit2,coef=2,200),  row.names = T, 
	      col.names =F,  sep="\t", quote=F )

  write.table(file="limMaFit.diff.200.ensemblGenes.txt" , x =topGene.ensembl,  row.names = F, 
	      col.names =T,  sep="\t", quote=F )
  # linear model of phenotype variance by gene expression

  #gene.df  <- copy(expression)
  gene.df  <- copy(expression)
  gene.df$disease  <- diseaseid$disease -1 
  annova.gene  <-  aov(disease ~ ., data=gene.df)
  temp3 <- cbind(gene.df[,c(1:311, 313), with=F], gene.df[,disease])
  lm.gene <- lm(V2 ~., data=temp3)
  #explained variance
  xx <- as.matrix(expression.z[, top.diff.genes]) ; yy = gene.df$disease
  yy <- scale(yy, center=T,scale=sd(yy))
  library(MASS)
  #to big to execute
  bb  <- ginv(t(xx) %*% xx)  %*% xx %*% yet 
  ev.gene <- 1- var( yy - xx %*% bb)/var(yy)
  # 
  #manova.gene  <-  manova(disease ~ ., data=gene.df)
}




# pca of physiological factors.
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

impute=F
library(mi)
if(impute){ 
  reduced.info = mi.info( physio.reduced) 
reduced.info <- update(reduced.info, "type", list(mitral_regurgitation.p="count", tricuspid_regurgitation.p = "count"))
#reduced.info <- update(reduced.info, "type", list())
physio.transformed <- mi.preprocess(physio.reduced, reduced.info)
IMP.reduced <- mi(physio.transformed, n.iter=1000,  max.minutes=500, check.coef.convergence=T,  add.noise=noise.control(post.run.iter=20))
save(file="IMP.reduced.RData", IMP)
}

load("IMP.reduced.RData")

imputed.matrix = mi.completed(IMP.reduced)
imputed1.mat <-imputed.matrix[[1]]
imputed1.mat <- imputed1.mat[!duplicated(physio.subset$sample_name[1:1809]), ]
#rownames(imputed1.mat)  <- unique(physio.subset$sample_name[1:1809])
physio.imputed <- as.data.table(imputed1.mat)
physio.imputed$sample_name  <- physio.subset$sample_name
setkey(physio.subset, sample_name)
setkey(pedtrait, sample_name)
setkey(physio.imputed, sample_name)
gene.names <- grep(x=colnames(pedtrait), pattern="^X", value=T)
physio.imputed.pedtrait <- merge(x=pedtrait[,c("sample_name", gene.names),with=F], y=physio.imputed, by = "sample_name", sort = F)
physio.imputed.pedtrait$etiology  <- (physio.subset[physio.imputed.pedtrait$sample_name, ])$etiology
#summary((physio.imputed.pedtrait[etiology=="donor"])$tricuspid_regurgitation.p)



if(F){
library(doMC)
registerDoMC(cores=8)
library(glmnet)
xx <- as.matrix(physio.imputed.pedtrait[,gene.names,with=F])
xx.healthy <- xx[physio.imputed.pedtrait$disease.p==1,]
xx.failure <- xx[physio.imputed.pedtrait$disease.p==2,]
healthy.inx <- which(physio.imputed.pedtrait$disease.p==1)
failure.inx <- which(physio.imputed.pedtrait$disease.p==2)
inx.curr <- physio.imputed.pedtrait$disease.p==2


#xx.curr <- xx.failure

inx.curr <- seq(nrow(xx)) 
  xx.curr <- xx
yy <- physio.imputed.pedtrait$lvef.p[inx.curr] 
inx1 <- which(!is.na(yy)  )
inx <- sample(inx1, .8*length(inx1))
#inx <- inx1
cvobj = cv.glmnet(xx.curr[inx,], yy[inx], parallel=T, family="gaussian")
aa  <-  coef(cvobj) 
gene.selected  <- aa[which(aa !=0),,drop=F]
gene.sel.names <- gsub(x=rownames(gene.selected)[-1], pattern="^X", replace="")
#gene.sel.hf  = gene.selected

test  <-  inx1[!( inx1 %in% inx) ]
#test  <-  inx1
test <- seq(nrow(xx))
inx.curr <- seq(nrow(xx))
xx.curr <- xx
mydf <- physio.imputed.pedtrait[inx.curr,list(sample_name, lvef.p)]
mydf$etiology=(physio.subset[mydf$sample_name, ])$etiology
mydf <- mydf[test]
mydf$predicted.lvef<-  predict(cvobj,newx=xx.curr[test,  ], s="lambda.min")
cor.test(mydf$lvef.p, mydf$predicted.lvef)
cor.test(mydf$lvef.p[healthy.inx], mydf$predicted.lvef[healthy.inx])
cor.test(mydf$lvef.p[failure.inx], mydf$predicted.lvef[failure.inx])
cat("R-square : ", 1 - min(cvobj$cvm)/var( (mydf$lvef)), "\n") 
p  <-  ggplot(data=mydf, aes( lvef.p, predicted.lvef))
p + geom_point(aes(color=etiology)) + ylab("predicted lvef by training on all samples") + xlab("lvef")

ggsave("prediction_cardiac_index_all_add_hf.pdf")




#### fig 1

yy.all <- physio.imputed.pedtrait$tricuspid_regurgitation.p
mydf <- physio.imputed.pedtrait[,list(sample_name, tricuspid_regurgitation.p)]
mydf$etiology=(physio.subset[mydf$sample_name, ])$etiology

out1 <-  model.gene(xx, yy.all, seq(nrow(xx)))
mydf$predicted <- out1$predicted 
p  <-  ggplot(data=mydf, aes( tricuspid_regurgitation.p, predicted))
p + geom_point(aes(color=etiology)) + ylab("predicted tricuspid_regurgitation") + xlab("tricuspid_regurgitation")
ggsave("tricuspid_regurgitation_predicted_imputed.pdf")

out2 <-  model.gene(xx, yy.all, healthy.inx)
mydf$predicted <- out2$predicted 
p  <-  ggplot(data=mydf, aes( tricuspid_regurgitation.p, predicted))
p + geom_point(aes(color=etiology)) + ylab("predicted tricuspid_regurgitation (trained on Donors) ") + xlab("tricuspid_regurgitation") 
ggsave("tricuspid_regurgitation_predicted_imputed_healthy.pdf")

out3 <-  model.gene(xx, yy.all, failure.inx)
mydf$predicted <- out3$predicted 
p  <-  ggplot(data=mydf, aes( tricuspid_regurgitation.p, predicted))
p + geom_point(aes(color=etiology)) + ylab("predicted tricuspid_regurgitation (trained on Heart failure cases) ") + xlab("tricuspid_regurgitation") 
ggsave("tricuspid_regurgitation_predicted_imputed_failure.pdf")
out1$r2
out2$r2
out3$r2
out1$correlation
out2$correlation
out3$correlation


 # binary feature 
yy.all <- physio.imputed.pedtrait$disease.p.p
mydf <- physio.imputed.pedtrait[,list(sample_name, disease.p.p)]
mydf$etiology=(physio.subset[mydf$sample_name, ])$etiology

out1 <-  model.gene(xx, factor(yy.all), seq(nrow(xx)), family="binomial")
mydf$predicted <- out1$predicted 
p  <-  ggplot(data=mydf, aes( disease.p.p, predicted))
p + geom_point(aes(color=etiology)) + ylab("predicted disease.p") + xlab("disease.p")
ggsave("disease.p_predicted_imputed.pdf")

out2 <-  model.gene(xx, yy.all, healthy.inx)
mydf$predicted <- out2$predicted 
p  <-  ggplot(data=mydf, aes( disease.p.p, predicted))
p + geom_point(aes(color=etiology)) + ylab("predicted disease.p (trained on Donors) ") + xlab("disease.p") 
ggsave("disease.p_predicted_imputed_healthy.pdf")

out3 <-  model.gene(xx, yy.all, failure.inx)
mydf$predicted <- out3$predicted 
p  <-  ggplot(data=mydf, aes( disease.p.p, predicted))
p + geom_point(aes(color=etiology)) + ylab("predicted disease.p (trained on Heart failure cases) ") + xlab("disease.p") 
ggsave("disease.p_predicted_imputed_failure.pdf")
out1$r2
out2$r2
out3$r2
out1$correlation
out2$correlation
out3$correlation

}



