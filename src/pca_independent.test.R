setwd("independent_test")
sample.name <- rownames(expression.z)  
train <- sample(sample.name, size=200, replace=F)  
test <- sample.name[!( sample.name %in%  train)]
expression.train  <- expression.diff.hf[train, ]; 
expression.test  <- expression.diff.hf[test, ];
pedtrait.info.train <- pedtrait.info[train]
pedtrait.info.test <- pedtrait.info[test]
gene.weight <- apply(expression.train ,  2,function(tt) fisher.score(tt, pedtrait.info.train$disease) )
pca.gene <- PCA(expression.train, col.w=gene.weight, ncp=ncol(expression.train), graph=F)
gene.proj.df <- as.data.table(expression.test %*% pca.gene$var$coord)
eigen.gene <- colnames(pca.gene$ind$coord)
#proj.gene  <- expression.diff.hf[test,] %*%  pca.gene$rotation
gene.proj.df$sample_name  <-  test

expression.merge <- merge(x=physio.reduced, y=gene.proj.df, by = "sample_name", sort = F)
donor.eigen.gene <- colMeans( expression.merge[disease.p==1, eigen.gene, with=F])
expression.merge$tdi <- apply(expression.merge[,eigen.gene, with=F], 1 , function(tt) sum((tt - donor.eigen.gene)^2)) 

params = c( "lvef.p","lvedd.p","lvesd.p","mitral_regurgitation.p","tricuspid_regurgitation.p","cardiac_index.p","creatinine_level.p") 

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


# snp
snp.diff.hf <- as.matrix(snp.imputed)
rownames(snp.diff.hf)  <-  pedtrait$sample_name
snp.train  <- snp.diff.hf[train, ]; 
snp.test  <- snp.diff.hf[test, ];
snp.weight <- apply(snp.train ,  2,function(tt) fisher.score(tt, pedtrait.info.train$disease) )
pca.snp <- PCA(snp.train, col.w=snp.weight, ncp=ncol(snp.train), graph=F)
snp.proj.df <- as.data.table(snp.test %*% pca.snp$var$coord)
eigen.snp <- colnames(pca.snp$ind$coord)
#proj.snp  <- snp.diff.hf[test,] %*%  pca.snp$rotation
snp.proj.df$sample_name  <-  test

snp.proj.df$disease.p <- pedtrait.info.test$disease
donor.eigen.snp <- colMeans( snp.proj.df[disease.p==1, eigen.snp, with=F])
snp.proj.df$gdi <- apply(snp.proj.df[,eigen.snp, with=F], 1 , function(tt) sum((tt - donor.eigen.snp)^2)) 

#gene.proj.df$disease.p  <- pedtrait.info$disease
donor.eigen.gene <- colMeans( gene.proj.df[disease.p==1, eigen.gene, with=F])
gene.proj.df$tdi <- apply(gene.proj.df[,eigen.gene, with=F], 1 , function(tt) sum((tt - donor.eigen.gene)^2))
gene.proj.df$gdi <- snp.proj.df$gdi
expression.merge <- merge(x=physio.imputed, y=gene.proj.df, by = "sample_name", sort = F)

donor.eigen.physio <- colMeans( physio.proj.df[disease.p==1, eigen.physio, with=F])
physio.proj.df$gpdi <- apply(physio.proj.df[,eigen.physio, with=F], 1 , function(tt) sum((tt - donor.eigen.physio)^2))
expression.merge$gpdi <- (physio.proj.df[expression.merge$sample_name])$gpdi


curr.name <- "GDI"
curr.param <- "gdi"
aa <-  cor.test(expression.merge$gdi, expression.merge$tdi, ,method="spearman")
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

setwd("../")
