load("/cbcbhomes/vinash85/shortcuts/lasso/data/Data/MikeData/magnet_PedTraits.RData")
library(data.table)
require(mi) 
setkey(pedtrait, sample_name)
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

soruce("physioProcess.R")

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



