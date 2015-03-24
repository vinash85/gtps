############ pca using missing values ##########
#library(pcaMethods)
library(FactoMineR)

## use only cardiac measurement as physiology ###
setkey(physio.subset, sample_name)
setkey(pedtrait.info, sample_name)

physio.pedtrait <- merge(x=physio.subset, y=pedtrait.info, by = "sample_name", sort = F)
#setnames(physio.pedtrait, 1:5, paste0(colnames(pedtrait.info), ".ped"))

# princomp of gene expression  
top.genes.tab <- topTable(fit2, coef=2, number=2000)
top.diff.genes <- rownames(top.genes.tab)
if(F){
pca.gene <- prcomp(expression.z[, top.diff.genes], retx = TRUE, scale = FALSE, center=TRUE)
proj.gene  <- expression.z[, top.diff.genes] %*%  pca.gene$rotation
eigen.gene <- colnames(proj.gene)

expression.df <- as.data.table(proj.gene)
expression.df$sample_name  <- pedtrait$sample_name
expression.df$disease <- pedtrait$disease
setkey(expression.df, sample_name)
physio.pedtrait1 <- copy(physio.pedtrait)
setkey(physio.pedtrait1, sample_name)
expression.merge <- merge(x=expression.df, y=physio.pedtrait, by = "sample_name", sort = F)
donor.eigen.gene <- colMeans( expression.merge[disease.ped==1, eigen.gene, with=F])
expression.merge$tdi <- apply(expression.merge[,eigen.gene, with=F], 1 , function(tt) sum((tt - donor.eigen.gene))) 

library(ggplot2)
p <- ggplot(data=expression.merge, aes(x=PC1, y=PC2))
p+geom_point(aes(col=factor(etiology)))
ggsave("tdi_hf.pdf")

p <- ggplot(data=expression.merge, aes(x=tdi, y=lvef.p))
p+geom_point(aes(col=factor(etiology)))
ggsave("tdi_lvef.pdf")

p <- ggplot(data=expression.merge, aes(x=tdi, y=lvedd.p))
p+geom_point(aes(col=factor(etiology)))
ggsave("tdi_lvedd.pdf")

p <- ggplot(data=expression.merge, aes(x=tdi, y=lvesd.p))
p+geom_point(aes(col=factor(etiology)))
ggsave("tdi_lvesd.pdf")
}

source("~/project/gtps/gtps/src/fisher.score.R")
library(limma)
topGene <- topTreat(fit2,coef=2,2000)
expression.diff.hf  <-  expression.z[, rownames(topGene)]
gene.weight <- apply(expression.diff.hf ,  2,function(tt) fisher.score(tt, pedtrait.info$disease) )
pca.gene <- PCA(expression.diff.hf, col.w=gene.weight, ncp=ncol(expression.diff.hf), graph=F)
gene.proj.df <- as.data.table(pca.gene$ind$coord)
eigen.gene <- colnames(pca.gene$ind$coord)
gene.proj.df$sample_name  <-  rownames(pca.gene$ind$coord)



###without imputation #### 
sample.name.pedtrait <-  pedtrait.info$sample_name[ pedtrait.info$sample_name %in%  physio.subset$sample_name]
physio.pedtrait <- physio.subset[as.character(sample.name.pedtrait)]

physio.reduced <- physio.pedtrait[, list(site.p, age.p, gender.p, race.p, ethnicity.p, patient_weight_kg.p, height_cm.p, heart_weight_grams.p, prior_cabg.p,
				   history_of_afib_aflutter.p, history_of_diabetes.p, history_of_hypertension.p, ace_inhibitor.p, lvedd.p, lvesd.p, lvef.p,
				   mitral_regurgitation.p, tricuspid_regurgitation.p, hemos_rap.p, 
				   hemos_systolic.p,hemos_pwcp.p,cardiac_index.p, creatinine_level.p,disease.p, etiology)] 
aa =sapply(physio.reduced, function(tt) sum(is.na(tt)))
physio.reduced <-  physio.pedtrait[, list(sample_name, age.p, gender.p, race.p,  patient_weight_kg.p, height_cm.p, heart_weight_grams.p, prior_cabg.p,
				   history_of_afib_aflutter.p, history_of_diabetes.p, history_of_hypertension.p, ace_inhibitor.p, lvedd.p, lvesd.p, lvef.p,
				   mitral_regurgitation.p, tricuspid_regurgitation.p, cardiac_index.p, creatinine_level.p,disease.p, etiology)]



setkey(pedtrait.info, sample_name)
expression.merge <- merge(x=physio.reduced, y=gene.proj.df, by = "sample_name", sort = F)
donor.eigen.gene <- colMeans( expression.merge[disease.p==1, eigen.gene, with=F])
expression.merge$tdi <- apply(expression.merge[,eigen.gene, with=F], 1 , function(tt) sum((tt - donor.eigen.gene)^2)) 

p <- ggplot(data=expression.merge, aes(x=Dim.1, y=Dim.2))
p+geom_point(aes(col=factor(history_of_diabetes.p)))
ggsave("pca_weighted_hf.pdf")

params = c( "lvef.p","lvedd.p","lvesd.p","mitral_regurgitation.p","tricuspid_regurgitation.p","cardiac_index.p","creatinine_level.p") 

for(curr.param in params){
  curr.name <- toupper(gsub(curr.param, pattern="\\.p$", replace=""))
#donor.di <- colMeans( expression.merge[disease.p==1, curr.param, with=F])
expression.merge$pdi <- eval(parse(text=paste0("expression.merge$",curr.param)))
aa <-  cor.test(expression.merge$pdi, expression.merge$tdi, ,method="spearman")
p <- ggplot(data=expression.merge, aes(x=rank(tdi), y=rank(pdi)))
p+geom_point(aes(col=factor(etiology)) ) + geom_smooth(method=lm) +
labs(y= curr.name, x = "Transcriptomic Deviation index (Rank)",   title=paste0(curr.name , "\n correlation =",  format( aa$estimate, digits=3) , " \n p =", format(aa$p.value, digits=3, scientific=T)))   
ggsave(paste0( "tdi_weighted_",curr.name, ".pdf"))
}


if(F) {
  p <- ggplot(data=expression.merge[!etiology%in%c("Valvular", "Other")], aes(y=tdi, x=factor(history_of_hypertension.p)))
  p+geom_boxplot() + facet_grid(.~ etiology )
}


### with imputation###
library(missMDA)
physio.reduced.df <- as.data.frame(physio.subset[,list( age.p, gender.p, race.p,  patient_weight_kg.p, height_cm.p, heart_weight_grams.p, prior_cabg.p,
				   history_of_afib_aflutter.p, history_of_diabetes.p, history_of_hypertension.p, ace_inhibitor.p, lvedd.p, lvesd.p, lvef.p,
				   mitral_regurgitation.p, tricuspid_regurgitation.p, cardiac_index.p, creatinine_level.p,disease.p, etiology.Ischemic.p, etiology.Idiopathic.p) ])
imp <-  imputePCA(physio.reduced.df, ncp=10)

#physio.reduced.df <- as.data.frame(physio.reduced[,list( age.p, gender.p, race.p,  patient_weight_kg.p, height_cm.p, heart_weight_grams.p, prior_cabg.p,
				   #history_of_afib_aflutter.p, history_of_diabetes.p, history_of_hypertension.p, ace_inhibitor.p, lvedd.p, lvesd.p, lvef.p,
				   #mitral_regurgitation.p, tricuspid_regurgitation.p, cardiac_index.p, creatinine_level.p,disease.p, etiology) ])
#imp <-  imputePCA(physio.reduced.df, ncp=10)
physio.imputed  <-  as.data.table(imp$completeObs) 
physio.imputed$sample_name = physio.subset$sample_name
physio.imputed$etiology =   physio.subset$etiology
expression.merge <- merge(x=physio.imputed, y=gene.proj.df, by = "sample_name", sort = F, suffixes= c("", ".gene"))
donor.eigen.gene <- colMeans( expression.merge[disease.p==1, eigen.gene, with=F])
expression.merge$tdi <- apply(expression.merge[,eigen.gene, with=F], 1 , function(tt) sum((tt - donor.eigen.gene)^2)) 

for(curr.param in params){
  curr.name <- toupper(gsub(curr.param, pattern="\\.p$", replace=""))
#donor.di <- colMeans( expression.merge[disease.p==1, curr.param, with=F])
expression.merge$pdi <- eval(parse(text=paste0("expression.merge$",curr.param)))
aa <-  cor.test(expression.merge$pdi, expression.merge$tdi, ,method="spearman")
p <- ggplot(data=expression.merge, aes(x=rank(tdi), y=rank(pdi)))
p+geom_point(aes(col=factor(etiology)) ) + geom_smooth(method=lm) +
labs(y= curr.name, x = "Transcriptomic Deviation index (Rank)",   title=paste0(curr.name , "\n correlation =",  format( aa$estimate, digits=3) , " \n p =", format(aa$p.value, digits=3, scientific=T)))   
ggsave(paste0( "tdi_weighted_imputed_",curr.name, ".pdf"))
}

####local pdi####
for(curr.param in params){
  curr.name <- toupper(gsub(curr.param, pattern="\\.p$", replace=""))
donor.di <- colMeans( expression.merge[disease.p==1, curr.param, with=F])
expression.merge$pdi <- apply(expression.merge[,curr.param, with=F], 1 , function(tt) sum((tt - donor.di)^2))
aa <-  cor.test(expression.merge$pdi, expression.merge$tdi, ,method="spearman")
p <- ggplot(data=expression.merge, aes(x=rank(tdi), y=rank(pdi)))
p+geom_point(aes(col=factor(etiology)) ) + geom_smooth(method=lm) +
labs(y="Phenotypic deviation index (Rank) ", x = "Transcriptomic Deviation index (Rank)",   title=paste0(curr.name , "\n correlation =",  format( aa$estimate, digits=3) , " \n p =", format(aa$p.value, digits=3, scientific=T)))   
ggsave(paste0( "tdi_weighted_imputed_",curr.name, "pdi.pdf"))
}

for(curr.param in params){
  curr.name <- toupper(gsub(curr.param, pattern="\\.p$", replace=""))
donor.di <- colMeans( expression.merge[disease.p==1, curr.param, with=F])
expression.merge$pdi <- apply(expression.merge[,curr.param, with=F], 1 , function(tt) sum((tt - donor.di)^2))
aa <-  cor.test(expression.merge$pdi, expression.merge$tdi)
p <- ggplot(data=expression.merge, aes(x=tdi, y=pdi))
p+geom_point(aes(col=factor(etiology)) ) + geom_smooth(method=lm) +
labs(y="Phenotypic deviation index ", x = "Transcriptomic Deviation index",   title=paste0(curr.name , "\n correlation =",  format( aa$estimate, digits=3) , " \n p =", format(aa$p.value, digits=3, scientific=T)))   
ggsave(paste0( "tdi_weighted_imputed_",curr.name, ".value.pdi.pdf"))
}

#### GPDI

physio.imputed1 <- physio.imputed[!duplicated(sample_name)]
physio.imputed.params <-  as.matrix( physio.imputed1[,list( lvedd.p, lvesd.p, lvef.p,
				   mitral_regurgitation.p, tricuspid_regurgitation.p, cardiac_index.p, creatinine_level.p)])

rownames(physio.imputed.params)  <-  physio.imputed1$sample_name
physio.weight <- apply(physio.imputed.params,  2,function(tt) fisher.score(tt, physio.imputed1$disease.p) )
pca.physio <- PCA(physio.imputed.params, col.w=physio.weight, ncp=ncol(physio.imputed.params), graph=F)
physio.proj.df <- as.data.table(pca.physio$ind$coord)
eigen.physio <- colnames(pca.physio$ind$coord)
eigen.physio <-  paste0(eigen.physio, ".p")
setnames(physio.proj.df, 1:ncol(physio.proj.df), eigen.physio)
physio.proj.df$sample_name  <-  rownames(pca.physio$ind$coord)
physio.proj.df$etiology  <-  physio.imputed1$etiology
physio.proj.df$disease.p <-  physio.imputed1$disease.p
setkey(physio.proj.df, sample_name)
setkey(gene.proj.df, sample_name)
expression.merge <- merge(x=gene.proj.df, y=physio.proj.df, by = "sample_name", sort = F, suffixes== c(".gene", ""))
#eigen.physio <- paste0(eigen.physio, ".x") 
curr.name <- "GPDI"
curr.param <- "gpdi"
donor.eigen.physio <- colMeans( expression.merge[disease.p==1, eigen.physio, with=F])
expression.merge$gpdi <- apply(expression.merge[,eigen.physio, with=F], 1 , function(tt) sum((tt - donor.eigen.physio)^2))
donor.eigen.gene <- colMeans( expression.merge[disease.p==1, eigen.gene, with=F])
expression.merge$tdi <- apply(expression.merge[,eigen.gene, with=F], 1 , function(tt) sum((tt - donor.eigen.gene)^2)) 

aa <-  cor.test(expression.merge$gpdi, expression.merge$tdi, method="spearman")
p <- ggplot(data=expression.merge, aes(x=rank(tdi), y=rank(gpdi)))
p <- p+geom_point(aes(col=factor(etiology)) ) + geom_smooth(method=lm) +
labs(y="Global phenotypic deviation index (Rank) ", x = "Transcriptomic Deviation index (Rank)",   title=paste0(curr.name , "\n correlation =",  format( aa$estimate, digits=3) , " \n p =", format(aa$p.value, digits=3, scientific=T)))   
ggsave(paste0( "tdi_weighted_imputed_",curr.name, "pdi.pdf"), p)

aa <-  cor.test(expression.merge$gpdi, expression.merge$tdi)
p <- ggplot(data=expression.merge, aes(x=tdi, y=gpdi))
p <- p+geom_point(aes(col=factor(etiology)) ) + geom_smooth(method=lm) +
labs(y="Global phenotypic deviation index ", x = "Transcriptomic Deviation index ",   title=paste0(curr.name , "\n correlation =",  format( aa$estimate, digits=3) , " \n p =", format(aa$p.value, digits=3, scientific=T)))   
ggsave(paste0( "tdi_weighted_imputed_",curr.name, ".value.pdi.pdf"), p)




### genetic co-ordination####
## find associated snp with gwas ... 
## run lasso 
## run correlation 

snpid <- grep(colnames(pedtrait) ,pattern="SNP_A\\.", value=T)  
snp.df <- as.matrix(pedtrait[,snpid, with=F])
rownames(snp.df)  <-  pedtrait$sample_name
iter <- length(snpid)
iters <- 1:iter
disease.vec <-  pedtrait$disease - 1
gwas <- mclapply(iters ,  function(tt) {
		       aa = cor.test(snp.df[,tt], disease.vec);
		       return(c(aa$p.value, aa$estimate))
}, mc.cores=8)
snp.corr <- t(do.call(cbind, gwas))
aa  <-  order(snp.corr[,1], decreasing=F)
snp.select <- aa[1:10000] 
snp.sel.df  <- snp.df[,snp.select]


#### impute snp
imp <-  imputePCA(snp.sel.df , ncp=300)
snp.imputed <-  as.data.table(imp$completeObs) 
snp.weight <- apply(snp.imputed ,  2,function(tt) fisher.score(tt, pedtrait.info$disease) )
pca.snp <- PCA(snp.imputed, col.w=snp.weight, ncp=ncol(snp.imputed), graph=F)
snp.proj.df <- as.data.table(pca.snp$ind$coord)
eigen.snp <- colnames(pca.snp$ind$coord)
snp.proj.df$sample_name  <-  rownames(snp.sel.df)
snp.proj.df$disease.p <- pedtrait.info$disease
donor.eigen.snp <- colMeans( snp.proj.df[disease.p==1, eigen.snp, with=F])
snp.proj.df$gdi <- apply(snp.proj.df[,eigen.snp, with=F], 1 , function(tt) sum((tt - donor.eigen.snp)^2)) 

gene.proj.df$disease.p  <- pedtrait.info$disease
donor.eigen.gene <- colMeans( gene.proj.df[disease.p==1, eigen.gene, with=F])
gene.proj.df$tdi <- apply(gene.proj.df[,eigen.gene, with=F], 1 , function(tt) sum((tt - donor.eigen.gene)^2))
gene.proj.df$gdi <- snp.proj.df$gdi
expression.merge <- merge(x=physio.imputed, y=gene.proj.df, by = "sample_name", sort = F, suffixes = c("", ".gene") )

donor.eigen.physio <- colMeans( physio.proj.df[disease.p==1, eigen.physio, with=F])
physio.proj.df$gpdi <- apply(physio.proj.df[,eigen.physio, with=F], 1 , function(tt) sum((tt - donor.eigen.physio)^2))
expression.merge$gpdi <- (physio.proj.df[expression.merge$sample_name])$gpdi


curr.name <- "GDI"
curr.param <- "gdi"
aa <-  cor.test(expression.merge$gdi, expression.merge$tdi, method="spearman")
p <- ggplot(data=expression.merge, aes(x=rank(tdi), y=rank(gdi)))
p+geom_point(aes(col=factor(etiology)) ) + geom_smooth(method=lm) +
labs(y="Genotypic deviation index (Rank) ", x = "Transcriptomic Deviation index (Rank)",   title=paste0(curr.name , "\n correlation =",  format( aa$estimate, digits=3) , " \n p =", format(aa$p.value, digits=3, scientific=T)))   
ggsave(paste0( "tdi_weighted_imputed_",curr.name, "pdi.pdf"))


curr.name <- "GDI"
curr.param <- "gdi"
aa <-  cor.test(expression.merge$gdi, expression.merge$gpdi, ,method="spearman")
p <- ggplot(data=expression.merge, aes(x=rank(gdi), y=rank(gpdi)))
p+geom_point(aes(col=factor(etiology)) ) + geom_smooth(method=lm) +
labs(x="Genotypic deviation index (Rank) ", y = "Global phenotype Deviation index (Rank)",   title=paste0(curr.name , "\n correlation =",  format( aa$estimate, digits=3) , " \n p =", format(aa$p.value, digits=3, scientific=T)))   
ggsave(paste0( "gpdi_weighted_imputed_",curr.name, "pdi.pdf"))



curr.name <- "GDI"
curr.param <- "gdi"
expression.merge.disease <- expression.merge[disease.p==2]
aa <-  cor.test(expression.merge.disease$gdi, expression.merge.disease$gpdi )
p <- ggplot(data=expression.merge.disease, aes(x=gdi, y=gpdi))
p+geom_point(aes(col=factor(etiology)) ) + geom_smooth(method=lm) +
labs(x="Genotypic deviation index (Rank) ", y = "Global phenotype Deviation index (Rank)",   title=paste0(curr.name , "\n correlation =",  format( aa$estimate, digits=3) , " \n p =", format(aa$p.value, digits=3, scientific=T)))   
####local pdi with genetic####

for(curr.param in params){
  curr.name <- toupper(gsub(curr.param, pattern="\\.p$", replace=""))
donor.di <- colMeans( expression.merge[disease.p==1, curr.param, with=F])
expression.merge$pdi <- apply(expression.merge[,curr.param, with=F], 1 , function(tt) sum((tt - donor.di)^2))
aa <-  cor.test(expression.merge$pdi, expression.merge$tdi, ,method="spearman")
p <- ggplot(data=expression.merge, aes(x=rank(gdi), y=rank(pdi)))
p+geom_point(aes(col=factor(etiology)) ) + geom_smooth(method=lm) +
labs(x="Genotypic deviation index (Rank) ", y = "Phenotypic deviation index (Rank)",   title=paste0(curr.name , "\n correlation =",  format( aa$estimate, digits=3) , " \n p =", format(aa$p.value, digits=3, scientific=T)))   
ggsave(paste0( "gdi_weighted_imputed_",curr.name, "pdi.pdf"))
}


for(curr.param in params){
  curr.name <- toupper(gsub(curr.param, pattern="\\.p$", replace=""))
donor.di <- colMeans( expression.merge[disease.p==1, curr.param, with=F])
expression.merge$pdi <- apply(expression.merge[,curr.param, with=F], 1 , function(tt) sum((tt - donor.di)^2))
aa <-  cor.test(expression.merge$pdi, expression.merge$tdi)
p <- ggplot(data=expression.merge, aes(x=gdi, y=pdi))
p+geom_point(aes(col=factor(etiology)) ) + geom_smooth(method=lm) +
labs(x="Genotypic deviation index ", y = "Phenotypic deviation index ",   title=paste0(curr.name , "\n Pearson correlation =",  format( aa$estimate, digits=3) , " \n p =", format(aa$p.value, digits=3, scientific=T)))   
ggsave(paste0( "gdi_weighted_imputed_",curr.name, "value_pdi.pdf"))
}


## independent analysis to check for over-training###
source("pca_independent.R")

### combining genetic and transcriptomic state-space ###
exp.snp.merge <- merge(x=gene.proj.df, y=snp.proj.df[,c("sample_name", eigen.snp), with=F], by="sample_name", sort = F, suffixes = c(".gene", ".snp"))
eigen.gene.m <- paste0(eigen.gene,".gene")
eigen.snp.m <- paste0(eigen.gene,".snp")



donor.eigen.snp.gene <- colMeans( exp.snp.merge[disease.p==1, c(eigen.snp.m, eigen.gene.m, eigen.snp.m) , with=F])
exp.snp.merge$gtdi <- apply(exp.snp.merge[, c(eigen.snp.m, eigen.gene.m, eigen.snp.m), with=F], 1 , function(tt) sum((tt - donor.eigen.snp.gene)^2)) 
exp.snp.merge <- merge(x=exp.snp.merge, y=physio.proj.df, by="sample_name", sort = F, suffixes = c(".snp.gene", ".p"))
#exp.snp.merge$gpdi <- (physio.proj.df[exp.snp.merge$sample_name])$gpdi
#exp.snp.merge$etiology <- (physio.proj.df[exp.snp.merge$sample_name])$etiology

curr.name <- "GTDI"
curr.param <- "gtdi"
aa <-  cor.test(exp.snp.merge$gtdi, exp.snp.merge$gpdi)
p <- ggplot(data=exp.snp.merge, aes(x=rank(gtdi), y=rank(gpdi)))
p+geom_point(aes(col=factor(etiology)) ) + geom_smooth(method=lm) +
labs(x="Genotypic-transcriptomic deviation index (Rank) ", y = "Global phenotype Deviation index (Rank)",   title=paste0(curr.name , "\n correlation =",  format( aa$estimate, digits=3) , " \n p =", format(aa$p.value, digits=3, scientific=T)))   
ggsave(paste0( "gtdi_weighted_imputed_",curr.name, "pdi.pdf"))



####  combined state - space ######  

snp.gene <- cbind(expression.diff.hf,snp.imputed)
snp.gene.weight <- apply(snp.gene ,  2,function(tt) fisher.score(tt, pedtrait.info$disease) )
pca.snp.gene <- PCA(snp.gene, col.w=snp.gene.weight, ncp=ncol(snp.gene), graph=F)
snp.gene.proj.df <- as.data.table(pca.snp.gene$ind$coord)
eigen.snp.gene <- colnames(pca.snp.gene$ind$coord)
snp.gene.proj.df$sample_name  <-  rownames(snp.sel.df)
snp.gene.proj.df$disease.p <- pedtrait.info$disease
donor.eigen.snp.gene <- colMeans( snp.gene.proj.df[disease.p==1, eigen.snp.gene, with=F])
snp.gene.proj.df$gtdi <- apply(snp.gene.proj.df[,eigen.snp.gene, with=F], 1 , function(tt) sum((tt - donor.eigen.snp.gene)^2)) 
exp.snp.merge <- merge(x=snp.gene.proj.df, y=physio.proj.df, by="sample_name", sort = F, suffixes = c(".snp.gene", ".p"))

aa <-  cor.test(exp.snp.merge$gtdi, exp.snp.merge$gpdi )


#### hierarchical combining

exp.snp.mat <- merge(x=gene.proj.df, y=snp.proj.df[,c("sample_name", eigen.snp), with=F], by="sample_name", sort = F, suffixes = c(".gene", ".snp"))
eigen.gene.m <- paste0(eigen.gene,".gene")
eigen.snp.m <- paste0(eigen.gene,".snp")
exp.snp.gene.pca1 <- as.matrix(exp.snp.mat[,c(eigen.gene.m , eigen.snp.m ), with=F])
snp.gene.weight <- apply(exp.snp.gene.pca1 ,  2,function(tt) fisher.score(tt, pedtrait.info$disease) )
pca.snp.gene <- PCA(exp.snp.gene.pca1 , col.w=snp.gene.weight, ncp=ncol(snp.gene), graph=F)

snp.gene.proj.df <- as.data.table(pca.snp.gene$ind$coord)
eigen.snp.gene <- colnames(pca.snp.gene$ind$coord)
snp.gene.proj.df$sample_name  <-  rownames(snp.sel.df)
snp.gene.proj.df$disease.p <- pedtrait.info$disease
donor.eigen.snp.gene <- colMeans( snp.gene.proj.df[disease.p==1, eigen.snp.gene, with=F])
snp.gene.proj.df$gtdi <- apply(snp.gene.proj.df[,eigen.snp.gene, with=F], 1 , function(tt) sum((tt - donor.eigen.snp.gene)^2)) 
exp.snp.merge <- merge(x=snp.gene.proj.df, y=physio.proj.df, by="sample_name", sort = F, suffixes = c(".snp.gene", ".p"))

aa <-  cor.test(exp.snp.merge$gtdi, exp.snp.merge$gpdi )


curr.name <- "GTDI"
curr.param <- "gtdi"
aa <-  cor.test(exp.snp.merge$gtdi, exp.snp.merge$gpdi)
p <- ggplot(data=exp.snp.merge, aes(x=rank(gtdi), y=rank(gpdi)))
p+geom_point(aes(col=factor(etiology)) ) + geom_smooth(method=lm) +
labs(x="Genotypic-transcriptomic deviation index (Rank) ", y = "Global phenotype Deviation index (Rank)",   title=paste0(curr.name , "\n Pearson correlation =",  format( aa$estimate, digits=3) , " \n p =", format(aa$p.value, digits=3, scientific=T)))   
ggsave(paste0( "gtdi_weighted_imputed_",curr.name, "pdi.pdf"))

######################################################
### limiting to ischemic 
############################################

setkey(physio.imputed, sample_name)
physio.imputed.sub <- physio.imputed[as.character(sample.name.pedtrait)]
pedtrait.merge <- merge(x=pedtrait, y=physio.imputed.sub, by = "sample_name", sort = F)
ishemic.inx <- which(pedtrait.merge$etiology %in% c("donor", "Ischemic"))
pedtrait.ishemic <- pedtrait.merge[ishemic.inx]
expression <- pedtrait.ishemic[, grep(x=colnames(pedtrait), pattern="^X"), with=F]
expression.z <- scale(expression, scale=apply(expression, 2, sd), center=T)
rownames(expression.z)  = pedtrait.ishemic$sample_name
# analysis with limma
library(limma)
exp.mat  <- t(as.matrix(expression))
design <- cbind(grp1=1, disease=pedtrait.ishemic$disease) 
fit  <-  lmFit(exp.mat, design)
# Fold-change thresholding
fit2 <- treat(fit,lfc=0.1)
fit2 <- eBayes(fit2) 
topGene <- topTreat(fit2,coef=2,2000)
expression.diff.hf  <-  expression.z[, rownames(topGene)]
gene.weight <- apply(expression.diff.hf ,  2,function(tt) fisher.score(tt, pedtrait.ishemic$disease) )
pca.gene <- PCA(expression.diff.hf, col.w=gene.weight, ncp=ncol(expression.diff.hf), graph=F)
gene.proj.df <- as.data.table(pca.gene$ind$coord)
eigen.gene <- colnames(pca.gene$ind$coord)
gene.proj.df$sample_name  <-  rownames(pca.gene$ind$coord)

setkey(pedtrait.ishemic, sample_name)
expression.merge <- merge(x=pedtrait.ishemic, y=gene.proj.df, by = "sample_name", sort = F)
donor.eigen.gene <- colMeans( expression.merge[disease.p==1, eigen.gene, with=F])
expression.merge$tdi <- apply(expression.merge[,eigen.gene, with=F], 1 , function(tt) sum((tt - donor.eigen.gene)^2)) 

#####tdi 
expression.merge.copy <- copy(expression.merge)
expression.merge  <- expression.merge.copy[, c(params, "tdi", "disease", "sample_name", eigen.gene, "etiology", "disease.p"), with=F]
donor.eigen.physio <- colMeans( physio.proj.df[disease.p==1, eigen.physio, with=F])
physio.proj.df$gpdi <- apply(physio.proj.df[,eigen.physio, with=F], 1 , function(tt) sum((tt - donor.eigen.physio)^2))
setkey(physio.proj.df, sample_name)
setkey(snp.proj.df, sample_name)
expression.merge$gpdi <- (physio.proj.df[expression.merge$sample_name])$gpdi
expression.merge$gdi <- snp.proj.df[expression.merge$sample_name]$gdi


