## cluster and see if there is pattern due to etiology 
library(gplots)
require("avinash")
expression.disease <- expression.z[physio.pedtrait$disease.ped ==2,] 
etiology  <-  as.factor(physio.pedtrait$etiology[physio.pedtrait$disease.ped ==2])
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R")
mydatascale <- t(expression.disease)
hr <- hclust(as.dist(1-cor(t(mydatascale), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(mydatascale, method="spearman")), method="complete")
mycl <- cutree(hc, k=3); myrowhc <- sample(rainbow(256)); myrowhc <- myrowhc[as.vector(mycl)]; 
etiology.col <- as.factor(physio.pedtrait$etiology[physio.pedtrait$disease.ped ==2])
levels(etiology.col) = sample(rainbow(256))[1:4]
colcolor <- cbind( matrix(myrowhc, ncol=1),as.character(etiology.col) )
colnames(colcolor) <- c("cluster","etiology")
pdf(file="etiology_heatmap.pdf", height=12, width=16)
heatmap.3(mydatascale, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=my.colorFct(), scale="row", trace="none",ColSideColors= colcolor, dendrogram="col", key=F)
dev.off()



###index####
healthy.inx <- which(physio.pedtrait$disease.ped==1) 
failure.inx <- which(physio.pedtrait$disease.ped==2) 
physio.healthy <- physio.pedtrait[healthy.inx]
physio.failure <- physio.pedtrait[failure.inx]

physio.etiology <- physio.pedtrait[etiology %in% c("donor", "Idiopathic Dilated CMP", "Ischemic")]
#### differential expression between etiology ###
# analysis with limma
inx <- which(etiology %in% c("Idiopathic Dilated CMP", "Ischemic")) 
library(limma)
exp.mat  <- t(as.matrix(expression[inx]))
design <- ifelse( etiology[inx]=="Ischemic",1,-1)
design <- cbind(etiology=design) 
fit2 <- treat(fit,lfc=0.1)
topTreat(fit2,coef=1)
fit2 <- eBayes(fit2) 
#outDir <- paste0("result/", date)
#write.fit(fit2,file=paste0(outDir, "/limMaFit.MagNet.Genes.txt"))
pdf("limMaFit.eitology.Genes.pdf")
volcanoplot(fit2,coef=1,highlight=100, names=rownames(exp.mat))
dev.off()
top.genes.tab <- topTable(fit2, coef=1, number=10)
top.diff.genes <- rownames(top.genes.tab)



##########weight

library(ggplot2)

p <- ggplot(physio.pedtrait, aes(factor(etiology), patient_weight_kg)) 
p + geom_boxplot() + xlab("Etiology") + ylab("Weight (in kg)")
ggsave("etiology_weight.pdf")
mydf <- cbind(patient_weight_kg=physio.pedtrait$patient_weight_kg, expression)
mydf <- mydf[!is.na(mydf$patient_weight_kg), ]
weight.h.fit  <- lm(patient_weight_kg ~ ., data=mydf[healthy.inx, 1:100, with=F] ) 
weight.f.fit  <- lm(patient_weight_kg ~ ., data=mydf[failure.inx, ] ) 
library(genlasso)
p <-  ncol(expression)
D = diag(1,p)
mydf.curr <- mydf[healthy.inx]
out = genlasso(mydf$patient_weight_kg , X=as.matrix(mydf[,2:ncol(mydf), with=F]), D=D)

p <- ggplot(physio.pedtrait, aes(factor(etiology), age.ped)) 
p + geom_boxplot() + xlab("Etiology") + ylab("Age")
ggsave("etiology_age.pdf")

p <- ggplot(physio.pedtrait, aes(factor(etiology), height_cm)) 
p + geom_boxplot() + xlab("Etiology") + ylab("Height(cm)")
ggsave("etiology_height.pdf")

p <- ggplot(physio.pedtrait, aes(factor(etiology), patient_weight_kg/((height_cm/100)^2) )) 
p + geom_boxplot() + xlab("Etiology") + ylab("BMI")
ggsave("etiology_bmi.pdf")


p <- ggplot(physio.pedtrait, aes(factor(etiology),heart_weight_grams.p )) 
p + geom_boxplot() + xlab("Etiology") + ylab("Heart weight (grams)")
ggsave("etiology_heart_weight.pdf")

p <- ggplot(physio.pedtrait, aes(factor(etiology),heart_weight_grams.p/patient_weight_kg.p )) 
p + geom_boxplot() + xlab("Etiology") + ylab("Heart weight / patient_weight ")
ggsave("etiology_heartweight_BY_weight.pdf")


p <- ggplot(physio.pedtrait, aes(factor(etiology),heart_weight_grams.p/(patient_weight_kg/((height_cm/100)^2)) )) 
p + geom_boxplot() + xlab("Etiology") + ylab("Heart weight / BMI  ")
ggsave("etiology_heartweight_BY_BMI.pdf")


p <- ggplot(physio.pedtrait, aes(factor(etiology), cardiac_index.p)) 
p + geom_boxplot() + xlab("Etiology") + ylab("Cardiac index  ")
ggsave("etiology_cardiac.pdf")
#### binary prediction of demographics #####

p <- ggplot(physio.pedtrait, aes(factor(etiology), history_of_hypertension.p )) 
p + geom_violin() + xlab("Etiology") + ylab("History of hypertension  ")
ggsave("etiology_hypertension.pdf")

p <- ggplot(physio.pedtrait, aes(factor(etiology), history_of_diabetes.p )) 
p + geom_violin() + xlab("Etiology") + ylab("History of diabetes  ")
ggsave("etiology_diabetes.pdf")


p <- ggplot(physio.etiology[history_of_diabetes.p %in% 0:1], aes(factor(history_of_diabetes), patient_weight_kg.p)) 
p + geom_boxplot() + facet_grid(.~ etiology ) + xlab("Diabetes") + ylab("Weight (in kg) ")
ggsave("diabetes_weight.pdf")


p <- ggplot(physio.pedtrait[failure.inx], aes(factor(ace_inhibitor), patient_weight_kg )) 
p + geom_boxplot() + xlab("Etiology") + ylab("Ace inhibitor ")
ggsave("etiology_diabetes.pdf")


p <- ggplot(physio.etiology[history_of_hypertension %in% c("Yes", "No")], aes(factor(history_of_hypertension), patient_weight_kg )) 
p + geom_boxplot() + facet_grid(.~ etiology ) + xlab("History of hypertension") + ylab("Weight (in kg) ")
ggsave("history_hypertension_weight_etiology.pdf")

mf_labeller <- function(var, value){
    value <- as.character(value)
    if (var=="history_of_hypertension") { 
        value[value=="Yes"] <- "Hypertension"
        value[value=="No"]   <- "No hypertension"
    }
    if (var=="history_of_diabetes") { 
        value[value=="Yes"] <- "Diabetes"
        value[value=="No"]   <- "No diabetes"
    }
    return(value)
}

p <- ggplot(physio.etiology[(history_of_hypertension %in% c("Yes", "No")) & (history_of_diabetes.p %in% 0:1) ], aes(factor(etiology), patient_weight_kg )) 
p + geom_boxplot() + facet_grid(history_of_diabetes ~history_of_hypertension  , labeller=mf_labeller ) + xlab("etiology") + ylab("Weight (in kg) ")
ggsave("etiology_diabetes_weight_hypertension.pdf")


p <- ggplot(physio.etiology[(history_of_hypertension %in% c("Yes", "No")) & (history_of_diabetes.p %in% 0:1) ], aes(factor(history_of_diabetes), patient_weight_kg )) 
p + geom_boxplot() + facet_grid(history_of_hypertension ~ etiology, labeller=mf_labeller ) + xlab("History of Diabetes") + ylab("Weight (in kg) ")
ggsave("hypertension_weight_diabetes_etiology.pdf")

p <- ggplot(physio.etiology[(history_of_hypertension %in% c("Yes", "No")) & (history_of_diabetes.p %in% 0:1) ], aes(factor(history_of_hypertension), patient_weight_kg )) 
p + geom_boxplot() + facet_grid(history_of_diabetes ~ etiology, labeller=mf_labeller ) + xlab("History of hypertension ") + ylab("Weight (in kg) ")
ggsave("diabetes_weight_hypertension_etiology.pdf")



p <- ggplot(physio.etiology[(history_of_hypertension %in% c("Yes", "No")) & (history_of_diabetes.p %in% 0:1) ], aes(factor(history_of_diabetes), heart_weight_grams.p )) 
 p + geom_boxplot() + facet_grid(history_of_hypertension  ~ etiology, , labeller=mf_labeller ) + xlab("History of Diabetes") + ylab("heart Weight (in gm) ")
ggsave("etiology_diabetes_heart_weight_hypertension.pdf")

p <- ggplot(physio.etiology[(history_of_hypertension %in% c("Yes", "No")) & (history_of_diabetes.p %in% 0:1) ], aes(factor(history_of_diabetes), sbp.p )) 
 p + geom_boxplot() + facet_grid(history_of_hypertension  ~ etiology, , labeller=mf_labeller ) + xlab("History of Diabetes") + ylab("sbp ")
ggsave("etiology_diabetes_sbp_hypertension.pdf")

p <- ggplot(physio.etiology[(history_of_hypertension %in% c("Yes", "No")) & (history_of_diabetes.p %in% 0:1) ], aes(factor(history_of_diabetes), dbp.p )) 
 p + geom_boxplot() + facet_grid(history_of_hypertension  ~ etiology, , labeller=mf_labeller ) + xlab("History of Diabetes") + ylab("dbp ")
ggsave("etiology_diabetes_dbp_hypertension.pdf")


p <- ggplot(physio.etiology, aes(factor(etiology), lvedd.p)) 
p + geom_boxplot() + xlab("Etiology") + ylab("lvedd")
ggsave("etiology_lvedd.pdf")

p <- ggplot(physio.etiology, aes(factor(etiology), lvesd.p)) 
p + geom_boxplot() + xlab("Etiology") + ylab("lvesd")
ggsave("etiology_lvesd.pdf")

p <- ggplot(physio.etiology, aes(factor(etiology), lvef.p)) 
p + geom_boxplot() + xlab("Etiology") + ylab("lvef")
ggsave("etiology_lvef.pdf")
p <- ggplot(physio.etiology[(history_of_hypertension %in% c("Yes", "No")) & (history_of_diabetes.p %in% 0:1) ], aes(factor(history_of_diabetes), sbp.p )) 
 p + geom_boxplot() + facet_grid(history_of_hypertension  ~ etiology, , labeller=mf_labeller ) + xlab("History of Diabetes") + ylab("sbp ")
ggsave("etiology_diabetes_sbp_hypertension.pdf")

p <- ggplot(physio.etiology[(history_of_hypertension %in% c("Yes", "No")) & (history_of_diabetes.p %in% 0:1) ], aes(factor(history_of_diabetes), dbp.p )) 
 p + geom_boxplot() + facet_grid(history_of_hypertension  ~ etiology, , labeller=mf_labeller ) + xlab("History of Diabetes") + ylab("dbp ")
ggsave("etiology_diabetes_dbp_hypertension.pdf")



p <- ggplot(physio.etiology, aes(factor(etiology),mitral_regurgitation.p)) 
p + geom_boxplot() + xlab("Etiology") + ylab("mitral_regurgitation")
ggsave("etiology_mitral_regurgitation.pdf")

p <- ggplot(physio.etiology, aes(factor(etiology),tricuspid_regurgitation.p)) 
p + geom_boxplot() + xlab("Etiology") + ylab("tricuspid_regurgitation")
ggsave("etiology_tricuspid_regurgitation.pdf")


p <- ggplot(physio.etiology, aes(factor(etiology), hemos_rap.p)) 
p + geom_boxplot() + xlab("Etiology") + ylab("hemos_rap")
ggsave("etiology_hemos_rap.pdf")

p <- ggplot(physio.etiology, aes(factor(etiology), hemos_systolic.p)) 
p + geom_boxplot() + xlab("Etiology") + ylab("hemos_systolic")
ggsave("etiology_hemos_systolic.pdf")
p <- ggplot(physio.etiology, aes(factor(etiology), hemos_pwcp.p)) 
p + geom_boxplot() + xlab("Etiology") + ylab("hemos_pwcp")
ggsave("etiology_hemos_pwcp.pdf")

p <- ggplot(physio.etiology, aes(factor(etiology), sbp.p)) 
p + geom_boxplot() + xlab("Etiology") + ylab("sbp")
ggsave("etiology_sbp.pdf")

p <- ggplot(physio.etiology, aes(factor(etiology), dbp.p)) 
p + geom_boxplot() + xlab("Etiology") + ylab("dbp")
ggsave("etiology_dbp.pdf")

p <- ggplot(physio.etiology, aes(factor(etiology), creatinine_level.p)) 
p + geom_boxplot() + xlab("Etiology") + ylab("creatinine_level")
ggsave("etiology_creatinine_level.pdf")

p <- ggplot(physio.etiology, aes(factor(etiology), bnp.p)) 
p + geom_boxplot() + xlab("Etiology") + ylab("bnp")
ggsave("etiology_bnp.pdf")

#aa <- physio.pedtrait[, list(lvedd.p, lvesd.p,lvef.p,  hemos_rap.p,
		    #hemos_systolic.p, hemos_diastolic.p, hemos_pwcp.p, hemos_co.p, cardiac_index.p, sbp.p, dbp.p, creatinine_level.p)]

### explaining variables

library(doMC)
registerDoMC(cores=8)
library(glmnet)
setkey(pedtrait, sample_name)
xx <- as.matrix(pedtrait[physio.pedtrait$sample_name, grep(x=colnames(pedtrait), pattern="^X"), with=F ])
xx.healthy <- xx[healthy.inx,]
xx.failure <- xx[failure.inx,]
inx1 <- which(!is.na(physio.pedtrait$disease.ped))
inx <- sample(inx1, 250)
#inx <- inx1
cvobj = cv.glmnet(xx[inx,], factor(physio.pedtrait$disease.ped[inx]), parallel=T, family="binomial", type.measure="class")
aa  <-  coef(cvobj, s="lambda.min") 
gene.selected  <- aa[which(aa !=0),,drop=F]
gene.sel.names <- gsub(x=rownames(gene.selected)[-1], pattern="^X", replace="")
gene.sel.hf  = gene.selected

test  <-  inx1[!( inx1 %in% inx) ]
mydf <- physio.pedtrait[test,list(sample_name, etiology=factor(etiology), disease.ped)]
mydf$predicted.hf<-  predict(cvobj,newx=xx[test,  ], s="lambda.min")
p  <-  ggplot(data=mydf, aes( etiology, predicted.hf))
p + geom_point(aes(color=etiology)) + ylab("predicted cardiac_index by training on all samples") + xlab("cardiac_index")
ggsave("prediction_cardiac_index_all_add_hf.pdf")




predicted.hf <-  predict(cvobj,newx=xx, s="lambda.min")
#xx1 <-  cbind(hypertenstion=physio.pedtrait$history_of_hypertension.p , diabetes= physio.pedtrait$history_of_diabetes.p,  hf=physio.pedtrait$disease.ped,  xx)
xx1 <-  cbind( hf=physio.pedtrait$disease.ped,  xx)
#xx2 <- xx1[1:200,]
inx1 <- which(!is.na(physio.pedtrait$cardiac_index.p))
inx <- sample(inx1, .7 * length(inx1))
cvobj1 = cv.glmnet(xx1[inx,], physio.pedtrait$cardiac_index.p[inx], parallel=T, family="gaussian", type.measure="auc" )
aa  <-  coef(cvobj1) 
gene.selected1  <- aa[which(aa !=0),,drop=F]
gene.sel.names1 <- gsub(x=rownames(gene.selected1)[-1], pattern="^X", replace="")

test  <-  inx1[!(inx1 %in% inx) ]
mydf <- physio.pedtrait[test,list(sample_name, etiology=factor(etiology), cardiac_index.p)]
mydf$predicted.hf<-  predict(cvobj1,newx=xx1[test,  ], s="lambda.min")
p  <-  ggplot(data=mydf, aes( cardiac_index.p, predicted.hf))
p + geom_point(aes(color=etiology)) + ylab("predicted cardiac_index by training on all samples") + xlab("cardiac_index")


inx <- which(!is.na(physio.healthy$cardiac_index.p))
cvobj2 = cv.glmnet(xx.healthy[inx,], scale(physio.healthy$cardiac_index.p[inx]), parallel=T, family="gaussian", nfolds=length(inx))
aa  <-  coef(cvobj2) 
gene.selected2  <- aa[which(aa !=0),,drop=F]
gene.sel.names2 <- gsub(x=rownames(gene.selected2)[-1], pattern="^X", replace="")

#glmobj <- glmnet(xx.healthy[inx,], scale(physio.healthy$cardiac_index.p[inx]), family="gaussian")

inx <- which(!is.na(physio.failure$cardiac_index.p))
cvobj3 = cv.glmnet(xx.failure[inx,], scale(physio.failure$cardiac_index.p[inx]), parallel=T, family="gaussian")
aa  <-  coef(cvobj3) 
gene.selected3  <- aa[which(aa !=0),,drop=F]
gene.sel.names3 <- gsub(x=rownames(gene.selected3[-1]), pattern="^X", replace="")

gene.sel.mat1 <- data.table(gene=rownames(gene.selected)[-1], coef=gene.selected1[-1])
setkey(gene.sel.mat1, gene)
go.tabl <-  read.table("GO.hf.regression.all.tab", header=T, sep="\t", comment.char = "", stringsAsFactors =F,  strip.white =T)
go.tabl$regression.coef <-  gene.sel.mat1[paste0("X", go.tabl$ID), ]$coef
write.table(file="GO.cardiac_index.regression.all.txt",x = go.tabl[,c(1,ncol(go.tabl), (2:ncol(go.tabl) -1)), with=F] , row.names = F, 
	col.names =T,  sep="\t", quote=F )


inx <- which(!is.na(physio.pedtrait$disease.ped))
cvobj = cv.glmnet(xx[inx,], as.factor(physio.pedtrait$disease.ped[inx]), parallel=T, family="binomial" )
aa  <-  coef(cvobj) 
gene.selected  <- aa[which(aa !=0),,drop=F]
gene.sel.names <- gsub(x=rownames(gene.selected)[-1], pattern="^X", replace="")

mydf <- physio.pedtrait[,list(sample_name, etiology=factor(etiology), cardiac_index.p)]

mydf$predicted.cardiac_index.healthy <-  predict(cvobj2,newx=xx, s="lambda.min")
p  <-  ggplot(data=mydf, aes(cardiac_index.p, predicted.cardiac_index.healthy))
p + geom_point(aes(color=etiology)) + ylab("predicted cardiac_index by training on donor") + xlab("cardiac_index")
ggsave("prediction_cardiac_index_healthy.pdf")


mydf$predicted.cardiac_index.all <-  predict(cvobj1,newx=xx1, s="lambda.min")
p  <-  ggplot(data=mydf, aes(cardiac_index.p, predicted.cardiac_index.all))
p + geom_point(aes(color=etiology)) + ylab("predicted cardiac_index by training on all samples") + xlab("cardiac_index")
ggsave("prediction_cardiac_index_all_add_hf.pdf")

#temp.gene <- c(8160168,8123446,8160670,7899187,8031550,8123104,7914075,8074991,7968242,7992987,7984813,7980970,7923347,8165277,7989146,8164177,7918558,7976496,8156278,7930139)
#0.03486717
#-0.09089275
#-0.17719294
#-0.09089275
#0.20541076
#-0.03146813
#-0.40317097
#-0.21424678
#-0.21424678
#-0.17719294
#-0.40317097
#-0.09089275
#-0.17719294
#-0.21424678
#0.03486717
#0.20541076
#0.03486717
#-0.03146813
#-0.40317097
#0.20541076
#temp  <-  predict(cvobj,newx=xx.failure, s="lambda.min") 
#inx1 <- which((!is.na(physio.failure$cardiac_index.p)) & (!(physio.failure$etiology == "Ischemic")))
inx1 <- which((!is.na(physio.failure$cardiac_index.p)) )
inx <- inx1[1:100]
 cc = apply(xx.failure[inx,], 2 , function(tt) (cor.test(tt, physio.failure$cardiac_index.p[inx] ))$p.value)
bb <- order(cc)
f.inx <- which(cc < 0.01)
#glmobj3 = cv.glmnet(xx.failure[inx,bb[1:200]], physio.failure$cardiac_index.p[inx],  family="gaussian", parallel=T, nfolds=10, nlambda= 1000,  maxit=1e6, thresh = 1e-08)
cvobj1 = cv.glmnet(xx.failure[inx,f.inx], physio.failure$cardiac_index.p[inx],  family="gaussian", parallel=T, nfolds=length(inx), alpha=0.5, lambda=la)
#gene.selected3  <- which(glmobj3$beta[,40] !=0) 
#gene.sel.names3 <- gsub(x=names(gene.selected3), pattern="^X", replace="")
aa <- coef(cvobj1)
gene.selected3  <- aa[which(aa !=0),,drop=F]
gene.sel.names3 <- gsub(x=rownames(gene.selected3)[-1], pattern="^X", replace="")
print(gene.sel.names3)
plot(glmobj3)
#temp <-  predict(glmobj3, xx.failure[inx1[71:90],], s=0.001,type="response")

#aa = apply(temp,2, function(tt) (cor.test(tt, physio.failure$cardiac_index.p[inx1[51:76]]))$p.value)


mydf <- physio.pedtrait[,list(sample_name, etiology=factor(etiology), cardiac_index.p)]

mydf$predicted.cardiac_index<-  predict(cvobj1,newx=xx[,f.inx], s="lambda.min")
p  <-  ggplot(data=mydf, aes(cardiac_index.p, predicted.cardiac_index))
p + geom_point(aes(color=etiology)) + ylab("predicted cardiac_index by heart failure model") + xlab("cardiac_index")
ggsave("prediction_cardiac_index_heart_failure_model.pdf")

####epi_eQTL
inx1 <- which((!is.na(physio.failure$disease.ped)) )
inx <- inx1[1:150]
nFeat <- ncol(xx.failure)
yy <- as.matrix( physio.failure$disease.ped[inx], ncol=1)
out <- epi_eQTL(x=xx.failure[inx, ], y=yy, feature=matrix(0, ncol=1,nrow=nFeat), pairFeature= matrix(0, ncol=1,nrow=nFeat), mask = cbind(rep(1, nFeat), 1:nFeat,1:nFeat), 
		alpHa= c(logit(0.5), rep(0,2)), ratio =1,  
		gamMa =rep(0, nFeat) ,  estimate_alpha=F, estimate_beta=T, B_inv_alpHa= matrix(1, ncol=1,nrow=1) , itermax = 100, thin= 5, burnIn=40, threads=8, gamMa_thres = .9,
		beTa_thres=0, balance=T, use_raoblackwell =F, logistic_variable_selection=F, num_logit_train=10000, regulator_prior=0.5, accIter=10, prior=1, rho=(nFeat/5)^2, verbose=T)


inx1 <- which((!is.na(physio.pedtrait$disease.ped)) )
inx <- sample(inx1, length=100)
#inx <- inx1[]
nFeat <- ncol(xx)
yy <- as.matrix( physio.pedtrait$disease.ped[inx], ncol=1)
out <- epi_eQTL(x=xx[inx, ], y=yy, feature=matrix(0, ncol=1,nrow=nFeat), pairFeature= matrix(0, ncol=1,nrow=nFeat), mask = cbind(rep(1, nFeat), 1:nFeat,1:nFeat), 
		alpHa= c(logit(0.5), rep(0,2)), ratio =1,  
		gamMa =rep(0, nFeat) ,  estimate_alpha=F, estimate_beta=T, B_inv_alpHa= matrix(1, ncol=1,nrow=1) , itermax = 100, thin= 5, burnIn=40, threads=8, gamMa_thres = .9,
		beTa_thres=0, balance=T, use_raoblackwell =F, logistic_variable_selection=F, num_logit_train=10000, regulator_prior=0.5, accIter=10, prior=1, rho=(nFeat/25)^2, verbose=T)
gamMa = apply(out$gamMa.mat,2,mean)
beTa = apply(out$beTa.mat,2, mean)
eeSNP = which(gamMa > 0.5)


#svm
#cardiac index
library(avinash)
xx.norm <- apply(xx,2, function(tt) ((tt -  min(tt))/(max(tt) - min(tt))))
neg.th <- quantile(physio.pedtrait$cardiac_index.p, na.rm=T, probs=.25)
neg.inx <- which(physio.pedtrait$cardiac_index.p <= neg.th)
xx.neg <- xx.norm[neg.inx, ]
yy.neg <- physio.pedtrait$cardiac_index.p[neg.inx]

pos.th <- quantile(physio.pedtrait$cardiac_index.p, na.rm=T, probs=.75)
pos.inx <- which(physio.pedtrait$cardiac_index.p >=  pos.th)
xx.pos <- xx.norm[pos.inx, ]
yy.pos <- physio.pedtrait$cardiac_index.p[pos.inx]

dataDir  <- "svm" 
writeSVMFeature(paste(dataDir,"trainset", sep="/"), neg=xx.neg, pos=xx.pos) 

#heart-failure
xx.norm <- apply(expression.z,2, function(tt) ((tt -  min(tt))/(max(tt) - min(tt))))
xx.neg <- xx.norm[as.character(pedtrait$sample_name[pedtrait$disease ==1]), ]
xx.pos <- xx.norm[as.character(pedtrait$sample_name[pedtrait$disease ==2]), ]
writeSVMFeature(paste(dataDir,"trainset", sep="/"), neg=xx.neg, pos=xx.pos)



#imputation 

physiological.tab <- fread("~/project/gtps/gtps/doc/magnet_redcap.txt")
physiological.tab <- physiological.tab[1:1809]
#physio.subset.temp  <- copy(physio.subset)
physio.subset  <- physiological.tab
####use run.R
 physio.subset[,disease.p:=ifelse(etiology=="donor",1,2)]
physiological.parse <- data.frame(physio.subset[,grep(colnames(physio.subset), pattern="\\.p$"),with=F ])
physiological.parse <-  physiological.parse[,-(8:11)]
physiological.parse$etiology  <- as.factor(physio.subset$etiology)
info = mi.info( physiological.parse) 
info <- update(info, "type", list(mitral_regurgitation.p="count"))
info <- update(info, "type", list(  tricuspid_regurgitation.p = "count"))
# run the imputation with data transformation
dat.transformed <- mi.preprocess(physiological.parse, info)
IMP2 <- mi(dat.transformed, n.iter=1000,max.minutes=500, add.noise=F)
IMP1 <- mi(dat.transformed, n.iter=1000,  max.minutes=200, check.coef.convergence=T,  add.noise=noise.control(post.run.iter=10))
IMP1 <- mi(IMP1, n.iter=1000,  max.minutes=360 )
#MP <- mi(MP, n.iter=10,  max.minutes=20, add.noise=noise.control(post.run.iter=5))
#MP <- mi(at.transformed, n.iter=6, check.coef.convergence=TRUE, add.noise=noise.control(post.run.iter=6))

# removing variables
#info <- update(info, "type", list(  tricuspid_regurgitation.p = "count"))

