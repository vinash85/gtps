library(R.matlab)

# lincs = readMat('~/project/project2-scratch/data/lincs/l1000/level2/gex_deltap_n53094x978.matrix.mat') # lincMat, cid and rid
lincs.dim = readMat('~/project/project2-scratch/data/lincs/l1000/level2/gex_epsilon_n1429794x978.dim.mat')
lincs = fread("gex_epsilon_n1429794x978.matrix.csv")
save(file="gex_epsilon_n1429794x978.matrix.RData", lincs, lincs.dim)

load("gex_epsilon_n1429794x978.matrix.RData")

library(avinash)
tab = fread("../metadata/inst.info")
##### pertubation category can be obtained from  http://support.lincscloud.org/hc/en-us/articles/202216073####
cid = unlist(lincs.dim$cid)
rid = unlist(lincs.dim$rid)
# sum(lincs$cid %in% tab$distil_id)
 # skl is skeletal myocytes.
# tab.sub = tab[distil_id %in%  lincs$cid]
tab.mcf7 = tab[cell_id=="MCF7"] 
tab.mcf7 = tab.mcf7[pert_type %in% c("trt_cp", "ctl_untrt")]
tab.mcf7 = tab.mcf7[pert_time == 24]			

 library(biomaRt)
 ensembl <- useMart("ensembl", dataset= "hsapiens_gene_ensembl")
 affy2ensembl = getBM(attributes = c("affy_hg_u133a", "start_position", "end_position", "affy_hugene_1_0_st_v1"), mart=ensembl)
affy2ensembl = affy2ensembl[affy2ensembl$affy_hg_u133a %in%  rid, ]
affy2ensembl = affy2ensembl[!is.na(affy2ensembl$affy_hugene_1_0_st_v1),]
affy2ensembl.dt = as.data.table(affy2ensembl)
load("~/shortcuts/Data/MikeData/MAGnet_combat_eset.Rdata") 
eset.combat
expression <- exprs(eset.combat)
physiological.tab <- fread("~/project/gtps/gtps/doc/magnet_redcap.txt")
setkey(physiological.tab, sample_name)
expression.annot = physiological.tab[colnames(expression)]
expression.annot[,disease.p:=ifelse(etiology=="donor",1,2)]
expression.annot[,history_of_afib_aflutter.p:=sapply(history_of_afib_aflutter, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
expression.annot[,history_of_diabetes.p:=sapply(history_of_diabetes, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
expression.annot[,history_of_hypertension.p:=sapply(history_of_hypertension, function(xx) switch(xx, "Yes"=1, "No"=0, NA))]
expression = t(expression)
affy2ensembl.dt = affy2ensembl.dt[affy_hugene_1_0_st_v1 %in% colnames(expression)]
affy2ensembl.dt = affy2ensembl.dt[,c(1,4),with=F]
affy2ensembl.dt = affy2ensembl.dt[!duplicated(affy2ensembl.dt$affy_hugene_1_0_st_v1)]
expression = expression[,as.character(affy2ensembl.dt$affy_hugene_1_0_st_v1)]

inx = cid %in% tab.mcf7$distil_id
setkey(tab.mcf7,distil_id)
tab.mcf7 = tab.mcf7[cid[inx] ]
lincs1 = as.matrix(lincs[,which(inx),with=F])
rm(lincs, lincs.dim)
rownames(lincs1) = rid
lincs1 = lincs1[affy2ensembl.dt$affy_hg_u133a,]
exp.trt = lincs1[  , tab.mcf7$pert_type == "trt_cp"]
exp.untrt= lincs1[, tab.mcf7$pert_type == "ctl_untrt"]
 tab.untrt = tab.mcf7[pert_type== "trt_cp"]
qnorm.col <- function(gene.exp.resid)
{
  gene.exp.norm <- gene.exp.resid
  for( sl in seq(ncol(gene.exp.resid ))) {
    mat = gene.exp.resid[,sl];
    mat = rank(mat,  rank, ties.method = "average");
    mat = qnorm(mat / (length(mat)+1));
    gene.exp.norm[,sl] = mat;
  }
   gene.exp.norm
}

exp.trt = qnorm.col(exp.trt)
exp.untrt = qnorm.col(exp.untrt)
magnet.exp = qnorm.col(t(expression))
exp.donor = apply(magnet.exp[,expression.annot$disease.p ==1], 1, mean)
exp.hf = apply(magnet.exp[,expression.annot$disease.p ==2], 1, mean)
exp.avg.untrt = rowMeans(exp.untrt[,1:6])
zero.exp = min(c(exp.donor, exp.hf))
exp.donor.trt= exp.trt
exp.hf.trt = exp.trt
# pertubations = unique(tab.mcf7$pert_desc)
pertubations = seq(ncol(exp.trt))
for (pert in pertubations) {
	exp.del = exp.trt[,pert] - exp.avg.untrt
	exp.donor.trt[,pert] = ifelse( exp.del + exp.donor > zero.exp, exp.del + exp.donor, zero.exp)
	exp.hf.trt[,pert] = ifelse( exp.del + exp.hf > zero.exp, exp.del + exp.hf, zero.exp)
}
# exp.donor.df = do.call(exp.donor.trt,rbind)
# exp.hf.df = do.call(exp.hf.trt_cp,rbind)


## transcriptomic space.

library(limma)
exp.mat  <- magnet.exp
rownames(exp.mat) <- affy2ensembl.dt$affy_hugene_1_0_st_v1
design <- cbind(grp1=1, disease=expression.annot$disease.p) 
fit  <-  lmFit(exp.mat, design)
# Fold-change thresholding
fit2 <- treat(fit,lfc=0.1)
topTreat(fit2,coef=2)
fit2 <- eBayes(fit2) 
topGene <- topTreat(fit2,coef=2,2000)
expression.diff.hf  <- t(exp.mat) # take al genes 
# expression.diff.hf= scale(expression.diff.hf,center=T, scale=apply(expression.diff.hf,2,sd)) 
source("~/project/gtps/gtps/src/fisher.score.R")
library(FactoMineR)
gene.weight <- apply(expression.diff.hf ,  2,function(tt) fisher.score(tt, expression.annot$disease) )
# see the google group discussion for more information https://groups.google.com/forum/#!topic/factominer-users/N8xzRCxCfQM
expression.join = rbind(expression.diff.hf, t(exp.donor.trt), t(exp.hf.trt))
rownames(expression.join) = c(rownames(expression.diff.hf), 
	paste0("donor", colnames(exp.donor.trt)), 
	paste0("hf", colnames(exp.hf.trt)) )
ind.sup = (nrow(expression.diff.hf) + 1): nrow(expression.join)
pca.gene <- PCA(expression.join, scale.unit=T, ind.sup=ind.sup , col.w=gene.weight, ncp=ncol(expression.diff.hf), graph=F) 
gene.proj.df <- as.data.table(pca.gene$ind.sup$coord)
eigen.gene <- colnames(pca.gene$ind.sup$coord)

magnet.gene.df = as.data.table(pca.gene$ind$coord)
donor.eigen.gene <- colMeans( magnet.gene.df[expression.annot$disease.p==1, eigen.gene, with=F])
expression.annot$tdi = apply(magnet.gene.df[,eigen.gene, with=F], 1 , function(tt) sum((tt - donor.eigen.gene)^2)) 
tdi <- apply(gene.proj.df[,eigen.gene, with=F], 1 , function(tt) sum((tt - donor.eigen.gene)^2)) 

