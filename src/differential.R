###### Transcriptomic space of donors #######
library(data.table)
load("~/shortcuts/Data/MikeData/MAGnet_combat_eset.Rdata") 
eset.combat
expression <- exprs(eset.combat)
source("~/project/gtps/gtps/src/physioProcess.R")
physio.reduced <- physio.subset[,list(sample_name, site.p, age.p, gender.p, race.p, ethnicity.p, patient_weight_kg.p, height_cm.p, heart_weight_grams.p, prior_cabg.p,
				   history_of_afib_aflutter.p, history_of_diabetes.p, history_of_hypertension.p, ace_inhibitor.p, lvedd.p, lvesd.p, lvef.p,
				   mitral_regurgitation.p, tricuspid_regurgitation.p, hemos_rap.p, hemos_systolic.p,hemos_pwcp.p,cardiac_index.p, creatinine_level.p,
				   disease.p, etiology)] 
setkey(physio.reduced, sample_name)
expression = expression[, colnames(expression) %in% physio.reduced$sample_name]
physio = physio.reduced[colnames(expression)]
library(limma)
exp.mat  <- (as.matrix(expression))
design <- cbind(grp1=1, disease=physio$disease.p) 
fit  <-  lmFit(exp.mat, design)
# Fold-change thresholding
fit2 <- treat(fit,lfc=0.1)
topTreat(fit2,coef=2)
fit2 <- eBayes(fit2) 
# date=Sys.Date()
# outDir <- paste0("result/", date)
# dir.create(outDir)

topGene <- topTreat(fit2,coef=2,2000)

# topGene$affyid <- gsub(rownames(topGene), pattern="^X",replacement="")
topGene.dt = data.table(topGene)
topGene.dt$affyid = rownames(topGene)
write.table(file="topDifferentialGene.txt",x = topGene.dt , row.names = F,col.names =T,  sep="\t", quote=F )